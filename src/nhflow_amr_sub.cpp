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

#include"nhflow_amr.h"
#include"nhflow_amr_fill.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_pressure.h"
#include"6DOF_nhflow.h"
#include<cmath>
#include<mpi.h>
#include<algorithm>

//  Subcycling in time for NHFLOW (G 7 1): Berger-Oliger with refluxing, level projections and
//  synchronisation projections (decision A of PLAN_step12c).
//
//  Level 0 takes its RK step alone with its own pressure solve (its covered columns included).
//  Its faces next to a patch keep their own flux, layer by layer; flux_hook sums it with the RK
//  weight of the stage times dt into the coarse register of the face (the face depth is the fine
//  one, as in the synchronous coupling).
//
//  Then every level l >= 1 takes two steps of half the step of level l-1 (sub_level, recursive).
//  The columns around a level-l patch are filled at the stage input times from the parent state
//  linear in time between the start of the parent step and its end: the parent columns the fills
//  read (5x5 around the source cell) are kept at the start of the parent step (tcol) and for the
//  fill replaced by the interpolated values (sub_swap_in/out), so fill_stage itself is unchanged.
//  Wet where either end was wet and the interpolated level is above A 544.  The patches sum their
//  box face fluxes into fine registers.  The pressure of a level-l stage is one solve over all
//  level-l patches (the composite solver restricted to the level, the parent columns fixed at the
//  values of the fill: sub_press_level).
//
//  After the two steps of level l+1 (sub_sync): the coarse cells of level l next to the patches
//  take the difference of the summed fine and coarse fluxes (continuity into WL, layer by layer
//  with the sigma thickness, momentum UH, VH, WH per layer), the levels above l are restricted
//  into level l, then the synchronisation projection: one composite correction solve over the
//  levels l..L for the divergence the coupled velocities have now (the parent columns of level l
//  fixed at 0), the velocity correction on all those grids, the correction added to the pressure.
//
//  The time step of level 0 is limited by the cells of level l as if they were 2^l times larger
//  (dt_cell_size).  The regrid runs at the end of a level-0 step.  Not with floating bodies; the
//  implicit diffusion per grid (G 31 1 needs all grids at the same time).

namespace
{
typedef reefamr_comms_off comms_off;
using nhflow_amr_detail::lij;
}

double nhflow_amr::rkw(int s) const
{
    if(mom0->stages()==2)
    return 0.5;
    return (s==2) ? 2.0/3.0 : 1.0/6.0;
}

double nhflow_amr::rkc(int s) const
{
    if(mom0->stages()==2)
    return (s==0) ? 0.0 : 1.0;
    return (s==0) ? 0.0 : ((s==1) ? 1.0 : 0.5);
}

double nhflow_amr::rko(int s) const
{
    if(mom0->stages()==2)
    return 1.0;
    return (s==1) ? 0.5 : 1.0;
}

// --------------------------------------------------------------------- parent columns
// values per column: eta, wet, deep, vb, U, V, W per layer, P per node
int nhflow_amr::tc_nv(int g)
{
    const int K = glex(g)->knoz;
    return 4 + 3*K + K+1;
}

// the columns of grid g that the fills of the next level read: 5x5 around every source cell
void nhflow_amr::tc_build(int g)
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

void nhflow_amr::tc_pack(int g, int n, double *v)
{
    lexer *q = glex(g);
    fdm_nhf *d = gfd(g);
    const tcol &T = tc[g+1];
    const int ii = T.ci[n], jj = T.cj[n];
    const int K = q->knoz;

    v[0] = d->eta(ii,jj);
    v[1] = q->wet[lij(q,ii,jj)];
    v[2] = q->deep[lij(q,ii,jj)];
    v[3] = (p0->A550==1) ? d->vb(ii,jj) : 0.0;
    for(int k=0; k<K; ++k)
    {
        const int c = cidx(q,ii,jj,k);
        v[4+k] = d->U[c];
        v[4+K+k] = d->V[c];
        v[4+2*K+k] = d->W[c];
    }
    for(int k=0; k<=K; ++k)
    v[4+3*K+k] = d->P[fidx(q,ii,jj,k)];
}

void nhflow_amr::tc_unpack(int g, int n, const double *v)
{
    lexer *q = glex(g);
    fdm_nhf *d = gfd(g);
    const tcol &T = tc[g+1];
    const int ii = T.ci[n], jj = T.cj[n];
    const int K = q->knoz;

    d->eta(ii,jj) = v[0];
    q->wet[lij(q,ii,jj)] = (int)v[1];
    q->deep[lij(q,ii,jj)] = (int)v[2];
    if(p0->A550==1)
    d->vb(ii,jj) = v[3];
    for(int k=0; k<K; ++k)
    {
        const int c = cidx(q,ii,jj,k);
        d->U[c] = v[4+k];
        d->V[c] = v[4+K+k];
        d->W[c] = v[4+2*K+k];
    }
    for(int k=0; k<=K; ++k)
    d->P[fidx(q,ii,jj,k)] = v[4+3*K+k];
}

// state of grid g at the start of its step (the columns the fills of the next level read)
void nhflow_amr::sub_snapshot(int g)
{
    tc_build(g);
    tcol &T = tc[g+1];
    const int nv = tc_nv(g);
    T.old.resize(T.ci.size()*(size_t)nv);
    for(size_t n=0; n<T.ci.size(); ++n)
    tc_pack(g,(int)n,&T.old[n*nv]);
}

// the parent columns of level lp at time theta of their step (0: start, 1: end = current)
// the grids of level l: level 0 the native grid and the blocks 1.. of a multi-block grid
vector<int> nhflow_amr::level_grids(int l) const
{
    vector<int> grids;
    if(l==0)
    grids.push_back(-1);
    for(int id : lev[l])
    grids.push_back(id);
    return grids;
}

void nhflow_amr::sub_swap_in(int lp, double theta)
{
    const vector<int> grids = level_grids(lp);

    for(int g : grids)
    {
        tc_build(g);
        tcol &T = tc[g+1];
        lexer *q = glex(g);
        fdm_nhf *d = gfd(g);
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

            // flags: either end wet and the interpolated level above the wet-dry depth; deep
            // where both ends are, otherwise the nearer end
            const double wl = v[0] + d->depth(T.ci[n],T.cj[n]);
            const int wo = (int)O[1], wn = (int)L[1];
            const int wt = (theta==0.0) ? wo : (((wo==1 || wn==1) && wl>q->A544+1.0e-6) ? 1 : 0);
            const int dp = (theta==0.0) ? (int)O[2] : ((O[2]==1.0 && L[2]==1.0) ? 1 : (int)(theta<0.5 ? O[2] : L[2]));
            v[1] = wt;
            v[2] = (wt==1) ? dp : 0;

            tc_unpack(g,(int)n,&v[0]);
        }
    }
}

void nhflow_amr::sub_swap_out(int lp)
{
    const vector<int> grids = level_grids(lp);

    for(int g : grids)
    {
        tcol &T = tc[g+1];
        const int nv = tc_nv(g);
        for(size_t n=0; n<T.ci.size(); ++n)
        tc_unpack(g,(int)n,&T.live[n*nv]);
    }
}

// the columns around the level-l patches for stage s (s<0: end of the step), the parent at time
// theta of its step
void nhflow_amr::sub_fill(ghostcell *pgc, int l, int s, double theta)
{
    const double t0 = MPI_Wtime();

    const bool swap = (theta<1.0);
    if(swap)
    sub_swap_in(l-1,theta);

    fill_stage(pgc,l,s);

    if(swap)
    sub_swap_out(l-1);

    tm[0] += MPI_Wtime()-t0;
}

// --------------------------------------------------------------------- registers
void nhflow_amr::sub_creg_reset(int g)
{
    const int t = g+1;
    const int K = glex(g)->knoz;
    cregL[t].assign(match[t].size()*4*(size_t)K,0.0);
    cregR[t].assign(rmatch[t].size()*4*(size_t)K,0.0);
}

// fine face depths of the patch boxes for the coarse faces next to them (the coarse grids run
// before the patches in a step): the patch's own face depth of its last stage, a fresh patch the
// mean depth (as the reconstruction)
void nhflow_amr::sub_dfx()
{
    for(auto q : P)
    {
        nhflow_amr_patch &c = *NP(q);
        fdm_nhf *d = c.d;
        const int il = EXT-1, ih = EXT+c.nx-1;
        const int jl = EY-1, jh = EY+c.ny-1;
        const bool fresh = c.rec[0][0].empty();

        for(int side=0; side<4; ++side)
        c.rec[0][side].resize(side<2 ? c.ny : c.nx);

        auto fx = [&](int a, int b) { return fresh ? 0.5*(d->depth(a+1,b)+d->depth(a,b)) : d->dfx(a,b); };
        auto fy = [&](int a, int b) { return fresh ? 0.5*(d->depth(a,b+1)+d->depth(a,b)) : d->dfy(a,b); };

        for(int r=0; r<c.ny; ++r)
        {
            c.rec[0][0][r] = fx(il,EY+r);
            c.rec[0][1][r] = fx(ih,EY+r);
        }
        for(int r=0; r<c.nx; ++r)
        {
            c.rec[0][2][r] = fy(EXT+r,jl);
            c.rec[0][3][r] = fy(EXT+r,jh);
        }
    }

    for(int l=1; l<=maxlev; ++l)
    exchange_fluxes(l);
}

// --------------------------------------------------------------------- the step
void nhflow_amr::sub_step(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_momentum_func *mom, nhflow_stage_obj &S0)
{
    const int ns = mom->stages();
    const double t = p->simtime, dt = p->dt;

    dtlev.assign(maxlev+1,0.0);
    dtlev[0] = dt;

    tc.resize(P.size()+1);
    cregL.resize(P.size()+1);
    cregR.resize(P.size()+1);
    rflux.resize(P.size()+1);

    double t0 = MPI_Wtime();
    sub_dfx();
    tm[2] += MPI_Wtime()-t0;

    S0p = &S0;
    mom->step_begin(p,d,pgc,S0);

    for(int g : level_grids(0))
    {
        sub_snapshot(g);
        sub_creg_reset(g);
    }

    // floating bodies: level 0 with a predicted copy (the finest level advances the body)
    const bool bpred = (b6!=nullptr && maxlev>0);
    if(bpred)
    {
        for(int nb=0; nb<b6->objects(); ++nb)
        b6->object(nb)->amr_save();
        b6->amr_predict = true;
    }
    bfin = 0;

    // level 0 alone, with its own pressure (multi-block: the native grid and the blocks, sub_level0_blocks)
    if(nblk>0)
    sub_level0_blocks(p,d,pgc,mom,S0);
    else
    for(int s=0; s<ns; ++s)
    {
        cur_stage = s;
        hstage = s;
        mom->phase_F(p,d,pgc,S0,s);
        mom->phase_M(p,d,pgc,S0,s);
        mom->phase_P(p,d,pgc,S0,s);
        mom->phase_E(p,d,pgc,S0,s);
    }
    cur_stage = -1;

    if(bpred)
    {
        b6->amr_predict = false;
        for(int nb=0; nb<b6->objects(); ++nb)
        b6->object(nb)->amr_restore(p,pgc);
    }

    // the finer levels, two steps each, then refluxing, restriction and the synchronisation
    for(int id : lev[1])
    {
        nhflow_amr_patch *c = NP(id);
        for(int ip=1; ip<=4; ++ip)
        for(int side=0; side<4; ++side)
        c->freg[ip][side].assign((size_t)(side<2 ? c->ny : c->nx)*c->pp->knoz,0.0);
    }

    sub_level(p,pgc,1,0,t,0.5*dt);
    sub_level(p,pgc,1,1,t+0.5*dt,0.5*dt);

    sub_sync(p,pgc,0);

    if(bpred)
    sub_body_end(p,pgc);
}

// step k (0, 1) of the level-l patches from time t with dt, within the step of level l-1; then
// the steps of the finer levels and the synchronisation of level l with them
// G 7 1 on a multi-block grid: the step of level 0 (the native grid and the blocks 1..) alone, as the stages of
// nhflow_amr::step without the patches: before every stage the cells around the blocks and the ring of the native
// grid from the cells they stand for (fill_level0), the continuity and momentum of the blocks and of the native
// grid (with breaking, A 550 1: the continuity of all, the breaking of level 0 over the blocks, then the momentum),
// the pressure as one solve over the grids of level 0 (the composite solve with the level window 0: the columns
// covered by the patches are unknowns of level 0 as on a uniform grid), the end of the stage.  The implicit
// diffusion per grid (G 31 1 needs all grids at the same time, as for the patches with G 7 1).  The turbulence of
// the blocks as in turb_step (level 0 has made its step in the driver)
void nhflow_amr::sub_level0_blocks(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_momentum_func *mom, nhflow_stage_obj &S0)
{
    const int ns = mom->stages();

    {
    comms_off guard(pgc);
    for(int id : lev[0])
    {
        nhflow_amr_patch *c = NP(id);
        lexer *pp = c->pp;
        pp->dt = p->dt;
        pp->dt_old = p->dt_old;
        pp->simtime = p->simtime;
        pp->count = p->count;
        pp->amr_follow(*p);
        pscope ps(pgc,c->d,d0);
        c->pmom->step_begin(pp,c->d,pgc,c->S);
    }
    }

    // A 560: the blocks step their models with the cells around them from level 0 and from each other
    if(turb)
    {
        const double t0 = MPI_Wtime();
        turb_fill0(pgc);
        for(int id : lev[0])
        turb_run(pgc,*NP(id));
        turb_fill0(pgc);
        for(int id : lev[0])
        turb_save(pgc,*NP(id));
        tturb += MPI_Wtime()-t0;
    }

    for(int s=0; s<ns; ++s)
    {
        cur_stage = s;
        hstage = s;

        double t0 = MPI_Wtime();
        fill_level0(pgc,s);
        tm[0] += MPI_Wtime()-t0;

        if(brk0)
        {
            blocks_phase(pgc,s,0);
            mom->phase_F(p,d,pgc,S0,s);
            t0 = MPI_Wtime();
            breaking_level0(pgc,-1);
            tm[1] += MPI_Wtime()-t0;
            blocks_phase(pgc,s,1);
            mom->phase_M(p,d,pgc,S0,s);
        }
        else
        {
            blocks_phase(pgc,s,0);
            blocks_phase(pgc,s,1);
            mom->phase_F(p,d,pgc,S0,s);
            mom->phase_M(p,d,pgc,S0,s);
        }

        halo0(pgc,s);

        mom->phase_P1(p,d,pgc,S0,s);
        {
        comms_off guard(pgc);
        for(int id : lev[0])
        {
            nhflow_amr_patch *c = NP(id);
            pscope ps(pgc,c->d,d0);
            c->pmom->phase_P1(c->pp,c->d,pgc,c->S,s);
        }
        }

        // the pressure of the grids of level 0 in one solve
        wlo = whi = 0;
        press_solve(p,pgc,s);
        wlo = 0;
        whi = -1;

        {
        comms_off guard(pgc);
        for(int id : lev[0])
        {
            nhflow_amr_patch *c = NP(id);
            pscope ps(pgc,c->d,d0);
            c->pmom->phase_P2(c->pp,c->d,pgc,c->S,s);
        }
        }
        mom->phase_P2(p,d,pgc,S0,s);

        {
        comms_off guard(pgc);
        for(int id : lev[0])
        {
            nhflow_amr_patch *c = NP(id);
            pscope ps(pgc,c->d,d0);
            c->pmom->phase_E(c->pp,c->d,pgc,c->S,s);
        }
        }
        mom->phase_E(p,d,pgc,S0,s);
    }
    cur_stage = -1;

    // the cells around the blocks and the ring at the end of the step (the parents of level 1)
    fill_level0(pgc,-1);
}

void nhflow_amr::sub_level(lexer *p, ghostcell *pgc, int l, int k, double t, double dt)
{
    const int ns = mom0->stages();
    double t0;
    dtlev[l] = dt;

    {
    comms_off guard(pgc);
    for(int n : lev[l])
    {
        nhflow_amr_patch *c = NP(n);
        lexer *pp = c->pp;
        pp->dt = dt;
        pp->dt_old = dt;
        pp->simtime = t;
        pp->count = p->count;
        pp->amr_follow(*p);
        pscope ps(pgc,c->d,d0);
        c->pmom->step_begin(pp,c->d,pgc,c->S);
    }
    }

    // floating bodies: a level coarser than the finest steps with a predicted copy
    const bool bpred = (b6!=nullptr && l<maxlev);
    if(bpred)
    for(int nb=0; nb<b6->objects(); ++nb)
    b6->object(nb)->amr_save();

    for(int s=0; s<ns; ++s)
    {
        // columns around the patches at the stage input time, the parents interpolated in time
        sub_fill(pgc,l,s,0.5*(double(k)+rkc(s)));

        // start state of the step (with the filled columns) for the fills of the finer level
        if(s==0 && l<maxlev)
        for(int n : lev[l])
        {
            sub_snapshot(n);
            sub_creg_reset(n);
        }

        cur_stage = s;
        hstage = s;

        // the bodies at the stage output of this level
        if(b6!=nullptr)
        sub_body_stage(p,pgc,l,s,t,dt);

        t0 = MPI_Wtime();
        {
        comms_off guard(pgc);
        for(int n : lev[l])
        {
            nhflow_amr_patch *c = NP(n);
            pscope ps(pgc,c->d,d0);
            c->pmom->phase_F(c->pp,c->d,pgc,c->S,s);
            c->pmom->phase_M(c->pp,c->d,pgc,c->S,s);
            c->pmom->phase_P1(c->pp,c->d,pgc,c->S,s);
        }
        }
        tm[1] += MPI_Wtime()-t0;

        sub_press_level(p,pgc,l,k,s);

        t0 = MPI_Wtime();
        {
        comms_off guard(pgc);
        for(int n : lev[l])
        {
            nhflow_amr_patch *c = NP(n);
            pscope ps(pgc,c->d,d0);
            c->pmom->phase_P2(c->pp,c->d,pgc,c->S,s);
            c->pmom->phase_E(c->pp,c->d,pgc,c->S,s);
        }
        }
        tm[1] += MPI_Wtime()-t0;

        // the finest level: the loads after its stage, every hull triangle on the finest grid at
        // the stage output (amr_grid_nhflow with cur_stage); the force output at the end of the
        // level-0 step
        if(b6!=nullptr && l==maxlev)
        {
            // the force file labels the loads with the start of the finest step, as with one
            // step for all levels
            const bool fin = (s==ns-1 && bfin==(1<<maxlev)-1);
            const double ts = p->simtime;
            p->simtime = t;
            for(int nb=0; nb<b6->objects(); ++nb)
            b6->object(nb)->hydrodynamic_forces_nhflow(p,d0,pgc,d0->WL,fin);
            p->simtime = ts;
            
            // external solver (REEF3D::CORAL): the finest step is its fluid step
            const double dts = p->dt;
            p->simtime = t;
            p->dt = dt;
            b6->amr_body_loads(p,pgc,s,ns);
            p->simtime = ts;
            p->dt = dts;
        }
    }
    cur_stage = -1;

    if(b6!=nullptr && l==maxlev)
    {
        const double dts = p->dt;
        p->dt = dt;
        for(int nb=0; nb<b6->objects(); ++nb)
        b6->object(nb)->saveTimeStep(p,ns-1);
        p->dt = dts;
        ++bfin;
    }

    if(bpred)
    for(int nb=0; nb<b6->objects(); ++nb)
    b6->object(nb)->amr_restore(p,pgc);

    if(l>=maxlev)
    return;

    // columns around the patches at the end of the step: the end state of the parents of l+1
    sub_fill(pgc,l,-1,0.5*double(k+1));

    for(int id : lev[l+1])
    {
        nhflow_amr_patch *c = NP(id);
        for(int ip=1; ip<=4; ++ip)
        for(int side=0; side<4; ++side)
        c->freg[ip][side].assign((size_t)(side<2 ? c->ny : c->nx)*c->pp->knoz,0.0);
    }

    sub_level(p,pgc,l+1,0,t,0.5*dt);
    sub_level(p,pgc,l+1,1,t+0.5*dt,0.5*dt);

    sub_sync(p,pgc,l);
}

// --------------------------------------------------------------------- floating bodies
// the synchronisation leaves the cells next to the body forcing alone: the divergence the
// reforcing of the last stage leaves there is not a mismatch of the levels, the next stage
// projects it as on one grid.  The forcing acts where the hull level set FB < psi (the half width
// of the smoothed Heaviside, A 526); a cell with a forced neighbour has FB < psi + one cell (FB is
// a distance), from its own value only: the same on any partition
void nhflow_amr::sub_body_mask(lexer *p, fdm_nhf *d)
{
    int i,j,k;
    size_t n=0;
    vector<double> &V = d->rhsvec.V;
    LOOP
    {
        const double psi = p->A526*0.5*(p->DXN[IP] + p->DYN[JP]);
        if(d->FB[IJK] < psi + MAX(p->DXN[IP],p->DYN[JP]))
        if(n<V.size())
        V[n] = 0.0;
        ++n;
    }
}

// stage s of the bodies with the step of level l (from t with dt): the finest level advances the
// body (external forces, the loads of its last stage), a coarser one its predicted copy (the
// loads frozen); the patches cast the hull and take their forcing in phase_P1
void nhflow_amr::sub_body_stage(lexer *p, ghostcell *pgc, int l, int s, double t, double dt)
{
    const double dts = p->dt, ts = p->simtime;
    p->dt = dt;
    p->simtime = t;

    b6->amr_body_stage(p,d0,pgc,s,l<maxlev);

    p->dt = dts;
    p->simtime = ts;
}

// end of the level-0 step: the hull on level 0 at the final pose (FB, the position line) and the
// output of the bodies, once per level-0 step: the state at the end of the level-0 step, labelled
// with that time (as the fluid output after the driver advances simtime); the vtp/stl print
// interval logic still sees the start of the last finest step
void nhflow_amr::sub_body_end(lexer *p, ghostcell *pgc)
{
    const double ts = p->simtime;
    p->simtime = ts + dtlev[0] - dtlev[maxlev];

    for(int nb=0; nb<b6->objects(); ++nb)
    {
        sixdof_obj_nhflow *o = b6->object(nb);
        o->update_position_nhflow(p,d0,pgc,d0->fs,true);

        if(p->X50==1)
        o->print_vtp(p,pgc);

        if(p->X50==2)
        o->print_stl(p,pgc);

        o->print_parameter(p,pgc,ts + dtlev[0]);
    }

    p->simtime = ts;
}

// the pressure of stage s of step k of level l: one solve over the level-l patches, the parent
// columns fixed at the parent pressure at the time of the stage output (interpolated in the parent
// step; pjm_corr: PCORR there is that pressure minus the pressure of the stage fill)
void nhflow_amr::sub_press_level(lexer *p, ghostcell *pgc, int l, int k, int s)
{
    if(p->A520==0)
    return;

    const double t0 = MPI_Wtime();
    const double alpha = mom0->stage_alpha(s);

    {
    comms_off guard(pgc);
    for(int n : lev[l])
    {
        nhflow_amr_patch *c = NP(n);
        pscope ps(pgc,c->d,d0);
        c->ptgt = c->ppress->amr_prepare(c->pp,c->d,pgc,alpha);
    }
    }

    // the parent columns at the time of the stage output (siblings: overwritten by the solve)
    {
        const double th = 0.5*(double(k)+rko(s));
        const bool swap = (th<1.0);
        if(swap)
        sub_swap_in(l-1,th);

        fill_col(l,7800+l,[&](int g) -> double* { return (g>=0 && P[g]->lev==l) ? NP(g)->ptgt : gfd(g)->P; });

        if(swap)
        sub_swap_out(l-1);

        if(p->A520==2)
        for(int n : lev[l])
        {
            nhflow_amr_patch *c = NP(n);
            lexer *pp = c->pp;
            for(const reefamr_fill &f : c->fill)
            for(int kk=0; kk<=pp->knoz; ++kk)
            {
                const int q = fidx(pp,f.di,f.dj,kk);
                c->ptgt[q] -= c->d->P[q];
            }
        }
    }

    wlo = whi = l;
    pr_edge = true;
    pr_core(p,pgc);
    pr_edge = false;
    wlo = 0;
    whi = -1;
    sub_lv_it += pr_it_last;
    ++sub_lv_n;

    {
    comms_off guard(pgc);
    for(int n : lev[l])
    {
        nhflow_amr_patch *c = NP(n);
        stg O = stage_out(n,s);
        pscope ps(pgc,c->d,d0);
        c->ppress->amr_finish(c->pp,c->d,pgc,*O.WL,O.UH,O.VH,O.WH,alpha);
    }
    }

    tm[4] += MPI_Wtime()-t0;
}

// --------------------------------------------------------------------- synchronisation
// level l after the two steps of level l+1: refluxing, restriction, synchronisation projection
void nhflow_amr::sub_sync(lexer *p, ghostcell *pgc, int l)
{
    const int last = mom0->stages()-1;
    double t0 = MPI_Wtime();

    sub_reflux(l);

    wlo = l;
    whi = maxlev;
    {
    comms_off guard(pgc);
    restrict_surface(last);
    restrict_momentum(last,false);
    }
    restrict_flags(pgc,last);
    wlo = 0;
    whi = -1;

    if(l==0)
    halo0(pgc,last);

    tm[3] += MPI_Wtime()-t0;

    if(p->A520>0)
    sub_project(p,pgc,l);
}

// coarse cells of level l next to the level-(l+1) patches: the summed fine face fluxes replace the
// summed coarse ones, layer by layer
void nhflow_amr::sub_reflux(int l)
{
    const int K = klev(l);
    const int Kf = klev(l+1);
    const int fz = Kf/K;
    const int NR = 4*K;

    for(size_t t=0; t<rmatch.size() && t<rflux.size(); ++t)
    rflux[t].assign(rmatch[t].size()*(size_t)NR,0.0);

    // the mean of the fine faces r, r+1 in the coarse layer k
    auto fmean = [&](const vector<double> &R, int side, int r, int k)
    {
        if(fz==1)
        return 0.5*(R[r*Kf+k]+R[(r+fr1(side))*Kf+k]);
        return 0.25*(R[r*Kf+2*k]+R[r*Kf+2*k+1]+R[(r+fr1(side))*Kf+2*k]+R[(r+fr1(side))*Kf+2*k+1]);
    };

    face_run_at(l+1,NR,7280+l+1,
                [&](reefamr_patch *q, int side, int r, double *sb)
                {
                    nhflow_amr_patch *c = NP(q);
                    for(int ip=1; ip<=4; ++ip)
                    for(int k=0; k<K; ++k)
                    sb[(ip-1)*K+k] = fmean(c->freg[ip][side],side,r,k);
                },
                [&](int tgt, int idx, const double *v)
                {
                    for(int m=0; m<NR; ++m)
                    rflux[tgt][(size_t)idx*NR+m] = v[m];
                });

    const vector<int> grids = level_grids(l);

    const double wd = p0->A544;
    vector<double> dlt(NR);

    for(int g : grids)
    {
        const int t = g+1;
        lexer *q = glex(g);
        fdm_nhf *d = gfd(g);

        auto fix = [&](const reefamr_match &m)
        {
            // the outside cell as seen on the target grid (multi-block: a turned neighbour block, side2)
            const int s2 = (m.side2>=0) ? m.side2 : m.side;
            const int ii = m.fi + (s2==1 ? 1 : 0);
            const int jj = m.fj + (s2==3 ? 1 : 0);
            const double sg = (s2==0 || s2==2) ? -1.0 : 1.0;
            // curvilinear grid: face length over cell area (the fluxes are per unit length)
            const double f = (q->geo_curv==1) ? sg*((m.dir==0) ? q->LEN1[lij(q,m.fi,m.fj)] : q->LEN2[lij(q,m.fi,m.fj)])/q->AREA2D[lij(q,ii,jj)]
                           : ((m.dir==0) ? sg/q->DXN[ii+marge] : sg*q->y_dir/q->DYN[jj+marge]);

            // continuity: the layer fluxes per unit sigma times the layer thickness
            double dw = 0.0;
            for(int k=0; k<K; ++k)
            dw += q->DZN[k+marge]*dlt[3*K+k];
            d->WL(ii,jj) += f*dw;

            for(int k=0; k<K; ++k)
            {
                const int c = cidx(q,ii,jj,k);
                d->UH[c] += f*dlt[k];
                d->VH[c] += f*dlt[K+k];
                d->WH[c] += f*dlt[2*K+k];
            }

            // derived values as after the restriction
            const double wl = d->WL(ii,jj);
            const double wlvl = fabs(wl)>wd ? wl : 1.0e20;
            d->eta(ii,jj) = wl - d->depth(ii,jj);
            for(int k=0; k<K; ++k)
            {
                const int c = cidx(q,ii,jj,k);
                d->U[c] = d->UH[c]/wlvl;
                d->V[c] = d->VH[c]/wlvl;
                d->W[c] = d->WH[c]/wlvl;
            }
        };

        for(size_t n=0; n<match[t].size(); ++n)
        {
            const reefamr_match &m = match[t][n];
            nhflow_amr_patch *c = NP(m.child);
            // (multi-block: the fine fluxes against the normal of the coarse face, sgn, as flux_hook)
            for(int ip=1; ip<=4; ++ip)
            for(int k=0; k<K; ++k)
            dlt[(ip-1)*K+k] = double(m.sgn)*fmean(c->freg[ip][m.side],m.side,m.r,k) - cregL[t][(4*n+ip-1)*K+k];
            fix(m);
        }

        for(size_t n=0; n<rmatch[t].size(); ++n)
        {
            const reefamr_match &m = rmatch[t][n];
            for(int m2=0; m2<NR; ++m2)
            dlt[m2] = double(m.sgn)*rflux[t][n*NR+m2] - cregR[t][n*NR+m2];
            fix(m);
        }
    }
}

// synchronisation projection of the levels l..L at the end of a level-l step: the divergence of
// the coupled velocities removed by one composite correction solve (the parent columns of level
// l > 0 fixed at 0), the correction added to the pressure
void nhflow_amr::sub_project(lexer *p, ghostcell *pgc, int l)
{
    const double t0 = MPI_Wtime();
    const int last = mom0->stages()-1;

    // the window grids (level 0: the native grid and the blocks of a multi-block grid)
    vector<int> grids = level_grids(l);
    for(int m=l+1; m<=maxlev; ++m)
    for(int id : lev[m])
    grids.push_back(id);

    // columns around the patches of the levels above l at this time (all synchronous); the
    // columns around level l keep the values of its last fill (the parent is fixed)
    for(int m=l+1; m<=maxlev; ++m)
    fill_stage(pgc,m,-1);

    // the step of level l on all window grids
    vector<double> dts(grids.size());
    for(size_t n=0; n<grids.size(); ++n)
    {
        lexer *q = glex(grids[n]);
        dts[n] = q->dt;
        q->dt = dtlev[l];
    }

    // pjm (A 520 1) solves for the pressure itself: here for the correction, the pressure kept
    const bool total = (p->A520==1);
    vector<vector<double>> Psave;
    if(total)
    {
        Psave.resize(grids.size());
        for(size_t n=0; n<grids.size(); ++n)
        {
            lexer *q = glex(grids[n]);
            fdm_nhf *dd = gfd(grids[n]);
            const size_t nf = (size_t)q->imax*q->jmax*q->kmaxF;
            Psave[n].assign(dd->P,dd->P+nf);
            std::fill(dd->P,dd->P+nf,0.0);
        }
    }

    for(int g : grids)
    {
        if(g<0)
        {
            ptgt0 = S0p->ppress->amr_prepare(p0,d0,pgc,1.0);
            continue;
        }
        nhflow_amr_patch *c = NP(g);
        comms_off guard(pgc);
        pscope ps(pgc,c->d,d0);
        c->ptgt = c->ppress->amr_prepare(c->pp,c->d,pgc,1.0);
    }

    // floating bodies: not the divergence the reforcing leaves next to the body
    if(b6!=nullptr)
    for(int g : grids)
    sub_body_mask(glex(g),gfd(g));

    wlo = l;
    whi = maxlev;
    pr_edge = (l>0);
    pr_core(p,pgc);
    pr_edge = false;
    sub_sy_it += pr_it_last;
    ++sub_sy_n;

    // the sharp hull (X 15 8): amr_finish stores the Rhie-Chow face corrections of the pressure it is given
    // (nhflow_membrane_rc_store, linear in the pressure), here those of the correction; the next continuity
    // flux needs the corrections of the last stage plus those of the synchronisation
    // (the stored ones set to 0 first, so that a cell amr_finish does not write, e.g. a halo cell of a patch, keeps
    // its value once; the correction is the difference to the pressure kept: A 520 1)
    vector<vector<double>> RCX(grids.size()), RCY(grids.size());
    if(hb && total)
    for(size_t n=0; n<grids.size(); ++n)
    {
        lexer *q = glex(grids[n]);
        fdm_nhf *dd = gfd(grids[n]);
        if(dd->MBETA==nullptr)
        continue;
        const size_t nc = (size_t)q->imax*q->jmax*q->kmax;
        RCX[n].assign(dd->MRCX,dd->MRCX+nc);
        std::fill(dd->MRCX,dd->MRCX+nc,0.0);
        if(q->j_dir==1)
        {
        RCY[n].assign(dd->MRCY,dd->MRCY+nc);
        std::fill(dd->MRCY,dd->MRCY+nc,0.0);
        }
    }

    for(int g : grids)
    {
        stg O = stage_out(g,last);
        if(g<0)
        {
            S0p->ppress->amr_finish(p0,d0,pgc,*O.WL,O.UH,O.VH,O.WH,1.0);
            continue;
        }
        nhflow_amr_patch *c = NP(g);
        comms_off guard(pgc);
        pscope ps(pgc,c->d,d0);
        c->ppress->amr_finish(c->pp,c->d,pgc,*O.WL,O.UH,O.VH,O.WH,1.0);
    }

    if(total)
    for(size_t n=0; n<grids.size(); ++n)
    {
        fdm_nhf *dd = gfd(grids[n]);
        for(size_t m=0; m<Psave[n].size(); ++m)
        dd->P[m] += Psave[n][m];
    }

    for(size_t n=0; n<grids.size(); ++n)
    {
        fdm_nhf *dd = gfd(grids[n]);
        for(size_t m=0; m<RCX[n].size(); ++m)
        dd->MRCX[m] += RCX[n][m];
        for(size_t m=0; m<RCY[n].size(); ++m)
        dd->MRCY[m] += RCY[n][m];
    }

    for(size_t n=0; n<grids.size(); ++n)
    glex(grids[n])->dt = dts[n];

    // velocities of the corrected momentum, then the covered cells and the halo
    for(int g : grids)
    {
        stg O = stage_out(g,last);
        if(g<0)
        {
            mom0->velcalc(p0,d0,pgc,O.UH,O.VH,O.WH,*O.WL,1.0);
            continue;
        }
        nhflow_amr_patch *c = NP(g);
        comms_off guard(pgc);
        pscope ps(pgc,c->d,d0);
        c->pmom->velcalc(c->pp,c->d,pgc,O.UH,O.VH,O.WH,*O.WL,1.0);
    }

    {
    comms_off guard(pgc);
    restrict_momentum(last,true);
    }
    wlo = 0;
    whi = -1;

    if(l==0)
    halo0(pgc,last);

    // multi-block: the ring and the cells around the blocks from the corrected state
    if(l==0 && nblk>0)
    fill_level0(pgc,-1);

    tm[4] += MPI_Wtime()-t0;
    tsync += MPI_Wtime()-t0;
}
