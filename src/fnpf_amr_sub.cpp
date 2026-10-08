/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Architect: Hans Bihs
--------------------------------------------------------------------*/


#include"fnpf_amr.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"slice4.h"
#include"fnpf_fsfbc.h"
#include"fnpf_sigma.h"
#include"fnpf_fsf_update.h"
#include"fnpf_bed_update.h"
#include"fnpf_laplace_cds2.h"
#include"reefmg_core.h"
#include"reefamr_krylov.h"
#include"fnpf_amr_fill.h"
#include"fnpf_body.h"
#include"fnpf_6DOF.h"
#include<cmath>
#include<mpi.h>
#include<algorithm>

//  Subcycling in time for REEF3D::FNPF on the patch hierarchy (G 7 1): level l takes 2^l steps of
//  dt/2^l per level-0 step.
//
//  The free-surface conditions of FNPF (kfsfbc, dfsfbc) carry no fluxes, and Fi is solved from
//  Fifsf in every stage: there are no flux registers and no synchronisation projection; after the
//  two steps of a level the coarser grids take the restricted eta, Fifsf and Fi in their covered
//  columns and rebuild their free-surface derivatives, sigma grid and surface velocities there.
//
//  - level 0 steps alone (fnpf_RK3, the hooks off), its Laplace solve over all its columns, the
//    covered ones included (its own solver)
//  - then every level its two steps, recursively: per stage the tendencies of the patches, the
//    stage values, the cells around the patches from the parent linear in time between the start
//    and the end of the parent step (eta, Fifsf, and the Fi columns of the edge: the parent
//    columns the fills read are kept at the start of the parent step and swapped in), the sigma
//    grid, the Laplace solve over the level-l patches with the parent columns fixed at the time
//    of the stage output (the composite solver on a level window; its preconditioner restricts the
//    residual through the levels below to a level-0 V-cycle), the vertical velocity at the surface
//  - resolved bodies (fnpf_6DOF): the finest level advances the body (fnpf_6DOF_sub.cpp)
//  - time step: the patch cells limit the level-0 step as cells 2^l times larger

namespace
{
typedef reefamr_comms_off comms_off;

inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

enum { NR=REEFAMR_NR, NRH=REEFAMR_NRH, NPV=REEFAMR_NPV, NVV=REEFAMR_NVV, NS=REEFAMR_NS, NT=REEFAMR_NT,
       NPH=REEFAMR_NPH, NSH=REEFAMR_NSH, NRP=REEFAMR_NTMP, NPRE=REEFAMR_NVEC };

const int PSWEEP = 1;
}

void fnpf_amr::attach_level0(fnpf_fsf *pf, fnpf_sigma *ps, fnpf_fsf_update *pu, fnpf_bed_update *pb)
{
    pf0 = pf;
    psig0 = ps;
    pfu0 = pu;
    pbu0 = pb;
}

// --------------------------------------------------------------------- parent columns
// values per column: eta, Fifsf, Fi per node
int fnpf_amr::tc_nv(int g)
{
    return 2 + glex(g)->knoz + 1;
}

// the columns of grid g that the fills of the next level read: 5x5 around every source cell
void fnpf_amr::tc_build(int g)
{
    tcol &T = tc[g+1];
    if(T.layout==layout_id)
    return;

    T.layout = layout_id;
    T.ci.clear();
    T.cj.clear();

    const int lc = (g<0) ? 0 : P[g]->lev;
    if(lc>=maxlev)
    return;

    lexer *q = glex(g);
    vector<char> mark((size_t)q->imax*q->jmax,0);

    auto add = [&](const reefamr_fill &f)
    {
        if(f.kind!=1 || f.g!=g)
        return;
        for(int a=-2; a<=2; ++a)
        for(int b=-2; b<=2; ++b)
        {
            const int ii = f.si+a, jj = f.sj+b;
            if(ii<q->imin || ii>=q->imin+q->imax || jj<q->jmin || jj>=q->jmin+q->jmax)
            continue;
            char &mk = mark[lij(q,ii,jj)];
            if(mk)
            continue;
            mk = 1;
            T.ci.push_back(ii);
            T.cj.push_back(jj);
        }
    };

    for(int id : lev[lc+1])
    for(const reefamr_fill &f : P[id]->fill)
    add(f);

    for(const reefamr_fill &f : gserve[lc+1])
    add(f);
}

void fnpf_amr::tc_pack(int g, int n, double *v)
{
    lexer *q = glex(g);
    fdm_fnpf *c = gfd(g);
    const tcol &T = tc[g+1];
    const int ii = T.ci[n], jj = T.cj[n];
    v[0] = c->eta(ii,jj);
    v[1] = c->Fifsf(ii,jj);
    for(int k=0; k<=q->knoz; ++k)
    v[2+k] = c->Fi[fidx(q,ii,jj,k)];
}

void fnpf_amr::tc_unpack(int g, int n, const double *v)
{
    lexer *q = glex(g);
    fdm_fnpf *c = gfd(g);
    const tcol &T = tc[g+1];
    const int ii = T.ci[n], jj = T.cj[n];
    c->eta(ii,jj) = v[0];
    c->Fifsf(ii,jj) = v[1];
    for(int k=0; k<=q->knoz; ++k)
    c->Fi[fidx(q,ii,jj,k)] = v[2+k];
}

// state of grid g at the start of its step (the columns the fills of the next level read)
void fnpf_amr::sub_snapshot(int g)
{
    if((int)tc.size()<(int)P.size()+1)
    tc.resize(P.size()+1);
    tc_build(g);
    tcol &T = tc[g+1];
    const int nv = tc_nv(g);
    T.old.resize(T.ci.size()*(size_t)nv);
    for(size_t n=0; n<T.ci.size(); ++n)
    tc_pack(g,(int)n,&T.old[n*nv]);
}

// the parent columns of level lp at time theta of their step (0: start, 1: end = current)
void fnpf_amr::sub_swap_in(int lp, double theta)
{
    vector<int> grids;
    if(lp==0)
    grids.push_back(-1);
    else
    grids = lev[lp];

    for(int g : grids)
    {
        tcol &T = tc[g+1];
        const int nv = tc_nv(g);
        T.live.resize(T.ci.size()*(size_t)nv);
        vector<double> v(nv);

        for(size_t n=0; n<T.ci.size(); ++n)
        {
            double *L = &T.live[n*nv];
            const double *O = &T.old[n*nv];
            tc_pack(g,(int)n,L);
            for(int m=0; m<nv; ++m)
            v[m] = (1.0-theta)*O[m] + theta*L[m];
            tc_unpack(g,(int)n,&v[0]);
        }
    }
}

void fnpf_amr::sub_swap_out(int lp)
{
    vector<int> grids;
    if(lp==0)
    grids.push_back(-1);
    else
    grids = lev[lp];

    for(int g : grids)
    {
        tcol &T = tc[g+1];
        const int nv = tc_nv(g);
        for(size_t n=0; n<T.ci.size(); ++n)
        tc_unpack(g,(int)n,&T.live[n*nv]);
    }
}

// the columns around the level-l patches of psi_0 (body loads, f[g+1] of grid g): phi_t of the
// parent step (at fixed sigma: linear in time, as the fills of Fi), the siblings' current values
void fnpf_amr::sub_psi_edges(int l, double **f)
{
    vector<int> grids;
    if(l==1)
    grids.push_back(-1);
    else
    grids = lev[l-1];

    for(int g : grids)
    {
        tcol &T = tc[g+1];
        lexer *q = glex(g);
        fdm_fnpf *c = gfd(g);
        const size_t n7 = (size_t)q->imax*q->jmax*(q->kmax+2);
        if((size_t)T.ptn!=n7)
        {
            delete [] T.pt;
            T.pt = new double[n7]();
            T.ptn = (int)n7;
        }
        const int nv = tc_nv(g);
        const double r = 1.0/dtlev[l-1];
        for(size_t n=0; n<T.ci.size(); ++n)
        {
            const double *O = &T.old[n*nv];
            for(int k=0; k<=q->knoz; ++k)
            {
                const int m = fidx(q,T.ci[n],T.cj[n],k);
                T.pt[m] = (c->Fi[m] - O[2+k])*r;
            }
        }
    }

    fill_col(l,7900+l,[&](int g) -> double*
    {
        if(g>=0 && P[g]->lev==l)
        return f[g+1];
        return tc[g+1].pt;
    });
}

// --------------------------------------------------------------------- the step
// end of the level-0 step (fnpf_RK3, its stages done): the finer levels, two steps each, then the
// restriction into level 0
void fnpf_amr::sub_step(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    const double t = p->simtime, dt = p->dt;

    fnpf_6DOF *b = dynamic_cast<fnpf_6DOF*>(body);
    if(b!=nullptr)
    b->amr_restore(p,pgc);

    dtlev.assign(maxlev+1,0.0);
    dtlev[0] = dt;
    bfin = 0;

    sub_level(p,pgc,1,0,t,0.5*dt);
    sub_level(p,pgc,1,1,t+0.5*dt,0.5*dt);

    sub_sync(p,pgc,0,1.0);

    slev = 0;
    tsub = t;
}

// step k (0, 1) of the level-l patches from time t with dt, within the step of level l-1; then
// the steps of the finer levels and the restriction into level l
void fnpf_amr::sub_level(lexer *p, ghostcell *pgc, int l, int k, double t, double dt)
{
    dtlev[l] = dt;
    fnpf_6DOF *b = dynamic_cast<fnpf_6DOF*>(body);

    // the state at the start of the step: the parent columns of the fills of level l+1
    if(l<maxlev)
    for(int id : lev[l])
    sub_snapshot(id);

    if(b!=nullptr && l<maxlev)
    b->amr_save();

    for(int s=0; s<3; ++s)
    {
        slev = l;
        tsub = t;
        sub_stage(p,pgc,l,k,s,dt);
    }

    if(b!=nullptr && l<maxlev)
    b->amr_restore(p,pgc);
    if(l==maxlev)
    ++bfin;

    if(l>=maxlev)
    return;

    sub_level(p,pgc,l+1,0,t,0.5*dt);
    sub_level(p,pgc,l+1,1,t+0.5*dt,0.5*dt);

    sub_sync(p,pgc,l,0.5*double(k+1));
}

// stage s of step k of the level-l patches (as fnpf_RK3 with the patch tendencies)
void fnpf_amr::sub_stage(lexer *p, ghostcell *pgc, int l, int k, int s, double dt)
{
    double t0 = MPI_Wtime();
    fnpf_6DOF *b = dynamic_cast<fnpf_6DOF*>(body);

    // tendencies
    {
    comms_off guard(pgc);
    for(int id : lev[l])
    {
        fnpf_amr_patch *pc = FP(id);
        lexer *pp = pc->pp;
        fdm_fnpf *cc = pc->c;
        pp->dt = dt;
        pp->dt_old = dt;
        pp->simtime = tsub;
        pp->count = p->count;

        slice4 &ek = *pc->ek, &fk = *pc->fk;
        slice4 &erk1 = *pc->erk1, &erk2 = *pc->erk2;
        slice &Ein = (s==0) ? static_cast<slice&>(cc->eta) : (s==1) ? static_cast<slice&>(erk1) : static_cast<slice&>(erk2);

        lexer *p = pp;
        fdm_fnpf *c = cc;

        pc->pf->kfsfbc(p,c,pgc);
        SLICELOOP4
        ek(i,j) = c->K(i,j);

        pc->pf->dfsfbc(p,c,pgc,Ein);
        SLICELOOP4
        fk(i,j) = c->K(i,j);
    }
    }
    tm[0] += MPI_Wtime()-t0;

    // the body: the finest level advances it (loads of the stage, unless those of level 0 at the
    // start of the level-0 step), the levels in between a predicted copy
    if(b!=nullptr)
    {
        if(l==maxlev)
        b->amr_sub_finest(p,c0,pgc,l,s,tsub,dt,!(bfin==0 && s==0),bfin==(1<<maxlev)-1);
        else
        b->amr_sub_predict(p,pgc,s,dt);
    }

    t0 = MPI_Wtime();

    // stage values, the footprint of the body
    for(int id : lev[l])
    {
        fnpf_amr_patch *pc = FP(id);
        {
        comms_off guard(pgc);
        fdm_fnpf *cc = pc->c;
        slice4 &ek = *pc->ek, &fk = *pc->fk;
        slice4 &erk1 = *pc->erk1, &erk2 = *pc->erk2, &frk1 = *pc->frk1, &frk2 = *pc->frk2;
        lexer *p = pc->pp;
        fdm_fnpf *c = cc;

        if(s==0)
        {
            SLICELOOP4
            erk1(i,j) = c->eta(i,j) + p->dt*ek(i,j);
            SLICELOOP4
            frk1(i,j) = c->Fifsf(i,j) + p->dt*fk(i,j);
        }
        if(s==1)
        {
            SLICELOOP4
            erk2(i,j) = 0.75*c->eta(i,j) + 0.25*erk1(i,j) + 0.25*p->dt*ek(i,j);
            SLICELOOP4
            frk2(i,j) = 0.75*c->Fifsf(i,j) + 0.25*frk1(i,j) + 0.25*p->dt*fk(i,j);
        }
        if(s==2)
        {
            SLICELOOP4
            c->eta(i,j) = (1.0/3.0)*c->eta(i,j) + (2.0/3.0)*erk2(i,j) + (2.0/3.0)*p->dt*ek(i,j);
            SLICELOOP4
            c->Fifsf(i,j) = (1.0/3.0)*c->Fifsf(i,j) + (2.0/3.0)*frk2(i,j) + (2.0/3.0)*p->dt*fk(i,j);
        }
        }

        if(body!=nullptr)
        body->amr_surface(pc->pp,pgc,id,sval(id,s,0),sval(id,s,1));
    }

    // the cells around the patches at the time of the stage output, the parent interpolated in
    // time; the Fi columns of the edge as well (fixed in the Laplace solve)
    const double theta = 0.5*(double(k)+rko(s));
    const bool swap = (theta<1.0);
    if(swap)
    sub_swap_in(l-1,theta);

    auto sel = [&](int g, int m) -> slice&
    {
        if(g>=0 && P[g]->lev==l)
        return sval(g,s,m);
        return (m==0) ? static_cast<slice&>(gfd(g)->eta) : static_cast<slice&>(gfd(g)->Fifsf);
    };
    fill_sl(l,2,7950+l,sel);

    for(int id : lev[l])
    {
        walls_sl(*FP(id),sel(id,0),gcval_eta);
        walls_sl(*FP(id),sel(id,1),gcval_fifsf);
    }

    fill_col(l,7960+l,[&](int g) -> double* { return gfd(g)->Fi; });

    if(swap)
    sub_swap_out(l-1);

    tm[0] += MPI_Wtime()-t0;

    // the Laplace equation of the stage on the level-l patches
    t0 = MPI_Wtime();
    {
    comms_off guard(pgc);
    for(int id : lev[l])
    {
        fnpf_amr_patch *pc = FP(id);
        lexer *pp = pc->pp;
        fdm_fnpf *cc = pc->c;
        slice &Se = sval(id,s,0);
        slice &Sf = sval(id,s,1);

        pc->pf->fsfdisc(pp,cc,pgc,Se,Sf);
        pc->psig->sigma_update(pp,cc,pgc,pc->pf,Se);
        pc->pfu->fsfbc_sig(pp,cc,pgc,Sf,cc->Fi);
        pc->pbu->bedbc_sig(pp,cc,pgc,cc->Fi,pc->pf);
        walls_fi(*pc,cc->Fi);
    }
    }

    if(b!=nullptr)
    b->amr_geometry_level(p,pgc,l);

    {
    comms_off guard(pgc);
    for(int id : lev[l])
    {
        fnpf_amr_patch *pc = FP(id);
        pc->plap->start(pc->pp,pc->c,pgc,nullptr,pc->pf,pc->c->Fi,sval(id,s,1));
    }
    }

    ltgt.assign(P.size()+1,nullptr);
    ltgt[0] = c0->Fi;
    for(int id=0; id<(int)P.size(); ++id)
    ltgt[id+1] = FP(id)->c->Fi;

    wlo = whi = l;
    lap_edge = true;
    lap_kind = 0;
    lap_core(p,pgc);
    lap_edge = false;
    wlo = 0;
    whi = -1;

    if(b!=nullptr)
    b->amr_post_solve_level(p,pgc,l);

    {
    comms_off guard(pgc);
    for(int id : lev[l])
    {
        fnpf_amr_patch *pc = FP(id);
        walls_fi(*pc,pc->c->Fi);
        pc->pf->fsfwvel(pc->pp,pc->c,pgc,sval(id,s,0),sval(id,s,1));
    }
    }

    tm[1] += MPI_Wtime()-t0;
}

// a psi solve of the body loads on the level-l patches
void fnpf_amr::lap_solve_psi_level(lexer *p, ghostcell *pgc, int l, double **f, slice **D)
{
    const double t0 = MPI_Wtime();

    {
    comms_off guard(pgc);
    for(int id : lev[l])
    {
        fnpf_amr_patch *pc = FP(id);
        pc->plap->start(pc->pp,pc->c,pgc,nullptr,pc->pf,f[id+1],*D[id+1]);
    }
    }

    ltgt.assign(f,f+P.size()+1);

    wlo = whi = l;
    lap_edge = true;
    lap_kind = 1;
    lap_core(p,pgc);
    lap_kind = 0;
    lap_edge = false;
    wlo = 0;
    whi = -1;

    tm[1] += MPI_Wtime()-t0;
}

// after the two steps of level l+1: eta, Fifsf and Fi of the finer levels into their covered
// columns, then the free-surface derivatives, sigma grid and surface velocities of level l (its
// parent at time theta of the parent step)
void fnpf_amr::sub_sync(lexer *p, ghostcell *pgc, int l, double theta)
{
    const double t0 = MPI_Wtime();

    wlo = l;
    whi = -1;
    restrict_sl(2,[&](int g, int m) -> slice& { return (m==0) ? static_cast<slice&>(gfd(g)->eta) : static_cast<slice&>(gfd(g)->Fifsf); });
    restrict_col([&](int g) -> double* { return gfd(g)->Fi; });
    wlo = 0;

    sub_derived(p,pgc,l,theta);

    tm[0] += MPI_Wtime()-t0;
}

// the derived fields of the level-l grids after the restriction
void fnpf_amr::sub_derived(lexer *p, ghostcell *pgc, int l, double theta)
{
    if(l==0)
    {
        fdm_fnpf *c = c0;
        pgc->gcsl_start4(p,c->eta,gcval_eta);
        pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
        pf0->fsfdisc(p,c,pgc,c->eta,c->Fifsf);
        psig0->sigma_update(p,c,pgc,pf0,c->eta);
        pfu0->fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);
        pbu0->bedbc_sig(p,c,pgc,c->Fi,pf0);
        pgc->start7V(p,c->Fi,c->bc,250);
        pf0->fsfwvel(p,c,pgc,c->eta,c->Fifsf);
        pfu0->velcalc_sig(p,c,pgc,c->Fi);
        return;
    }

    // the cells around the level-l patches, the parent at this time
    const bool swap = (theta<1.0);
    if(swap)
    sub_swap_in(l-1,theta);
    fill_sl(l,2,7970+l,[&](int g, int m) -> slice& { return (m==0) ? static_cast<slice&>(gfd(g)->eta) : static_cast<slice&>(gfd(g)->Fifsf); });
    fill_col(l,7975+l,[&](int g) -> double* { return gfd(g)->Fi; });
    if(swap)
    sub_swap_out(l-1);

    comms_off guard(pgc);
    for(int id : lev[l])
    {
        fnpf_amr_patch *pc = FP(id);
        lexer *pp = pc->pp;
        fdm_fnpf *cc = pc->c;
        walls_sl(*pc,cc->eta,gcval_eta);
        walls_sl(*pc,cc->Fifsf,gcval_fifsf);
        pc->pf->fsfdisc(pp,cc,pgc,cc->eta,cc->Fifsf);
        pc->psig->sigma_update(pp,cc,pgc,pc->pf,cc->eta);
        pc->pfu->fsfbc_sig(pp,cc,pgc,cc->Fifsf,cc->Fi);
        pc->pbu->bedbc_sig(pp,cc,pgc,cc->Fi,pc->pf);
        walls_fi(*pc,cc->Fi);
        pc->pf->fsfwvel(pp,cc,pgc,cc->eta,cc->Fifsf);
    }
}

// --------------------------------------------------------------------- preconditioner
// the level solve of the window wlo..wtop (wlo > 0, its parent columns fixed): the FAC sweep of
// lap_prec over the window, with the coarse correction of a level-0 V-cycle on the window residual
// restricted through the levels below (in the vectors of the grids below the window, which have
// no rows in this solve: scratch, 0 outside the covered columns)
void fnpf_amr::lap_prec_win(int kr, int kz)
{
    const int wl = wlo;
    auto selz = [&](int g) -> double* { return lvec(g,kz); };
    auto inwin = [&](int g) { const int lg = (g<0) ? 0 : P[g]->lev; return lg>=wl && lg<=wtop(); };

    // 1. pre-smoothing on the window patches
    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        const size_t n = (size_t)q->imax*q->jmax*(q->kmax+2);
        std::fill(lvec(g,kz),lvec(g,kz)+n,0.0);
        if(!inwin(g))
        std::fill(lvec(g,NRP),lvec(g,NRP)+n,0.0);
    }

    for(int id=0; id<(int)P.size(); ++id)
    if(inwin(id))
    lap_local(id,kr,kz);

    // 2. window residual, z after the pre-smoothing kept
    lap_apply(kz,NRP);
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *r = lvec(g,kr);
        double *t = lvec(g,NRP);
        for(int qq : lg[g+1].lq)
        t[qq] = r[qq] - t[qq];
        for(int qq : lg[g+1].cq)
        t[qq] = 0.0;
    }

    for(int id=0; id<(int)P.size(); ++id)
    if(inwin(id))
    {
        lexer *q = glex(id);
        std::copy(lvec(id,kz),lvec(id,kz)+q->imax*q->jmax*(q->kmax+2),lvec(id,NPRE));
    }

    // 3. restricted down to level 0 through all levels, one level-0 V-cycle
    wlo = 0;
    restrict_col([&](int g) -> double* { return lvec(g,NRP); });

    {
        sc_level &L = mg0->fine();
        const double *r = lvec(-1,NRP);
        double *z = lvec(-1,kz);

        std::fill(L.u.begin(),L.u.end(),0.0);
        std::fill(L.f.begin(),L.f.end(),0.0);
        for(int ii=0; ii<p0->knox; ++ii)
        for(int jj=0; jj<p0->knoy; ++jj)
        for(int kk=0; kk<p0->knoz; ++kk)
        {
            const long lq = L.idx(ii,jj,kk);
            if(L.act[lq])
            L.f[lq] = r[fidx(p0,ii,jj,kk)];
        }

        mg0->vcycle(0,1,1);

        for(int ii=0; ii<p0->knox; ++ii)
        for(int jj=0; jj<p0->knoy; ++jj)
        for(int kk=0; kk<p0->knoz; ++kk)
        {
            const long lq = L.idx(ii,jj,kk);
            z[fidx(p0,ii,jj,kk)] = L.act[lq] ? L.u[lq] : 0.0;
        }

        if(p0->mpi_size>1)
        {
            pgc0->gcparax7(p0,z,7);
            pgc0->gcparax7co(p0,z,7);
        }
    }

    // 4. coarse to fine: below the window the correction interpolated (scratch), from the window
    //    on as lap_prec (the increment since the pre-smoothing, post-smoothing)
    for(int l=1; l<=wtop(); ++l)
    {
        if(l<wl)
        {
            prolong_interior_col(l,[](reefamr_patch*) { return true; },
                                 [&](int g) -> double* { return lvec(g,kz); },
                                 [&](int id) -> double* { return lvec(id,kz); });
            fill_col(l,7980+l,selz);
            continue;
        }

        if(l>wl)
        for(int id : lev[l-1])
        {
            lexer *q = glex(id);
            const double *z = lvec(id,kz);
            double *d = lvec(id,NPRE);
            const int n = q->imax*q->jmax*(q->kmax+2);
            for(int m=0; m<n; ++m)
            d[m] = z[m] - d[m];
        }

        prolong_interior_col(l,[](reefamr_patch*) { return true; },
                             [&](int g) -> double* { return (g<0 || P[g]->lev<wl) ? lvec(g,kz) : lvec(g,NPRE); },
                             [&](int id) -> double* { return lvec(id,NRP); });

        for(int id : lev[l])
        {
            fnpf_amr_patch *c = FP(id);
            lexer *pp = c->pp;
            double *z = lvec(id,kz);
            const double *t = lvec(id,NRP);
            for(int ii=EXT; ii<EXT+c->nx; ++ii)
            for(int jj=EXT; jj<EXT+c->ny; ++jj)
            for(int kk=0; kk<=pp->knoz; ++kk)
            {
                const int qq = fidx(pp,ii,jj,kk);
                z[qq] += t[qq];
            }
        }

        // the window: its lowest level with the fixed parent columns (0 in z)
        wlo = wl;
        lap_dir = true;
        fill_col(l,7700+l,selz);

        for(int id : lev[l])
        lap_local(id,kr,kz);

        if(l<wtop())
        fill_col(l,7700+l,selz);
        wlo = 0;
    }

    wlo = wl;
}
