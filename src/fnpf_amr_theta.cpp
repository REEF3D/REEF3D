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
#include"fnpf_sz_weights.h"
#include"fnpf_amr_fill.h"
#include<cmath>
#include<mpi.h>

//  theta-Robin (A 310 5) on the patches of the mesh refinement (G 1).
//
//  The patches take the theta step of fnpf_theta with the formulas of level 0 (see fnpf_theta.h): the tendencies
//  at t^n, the nonlinear parts N = (kinematic - w, dynamic + g eta), the Adams-Bashforth 2 predictor of the sigma
//  grid (Euler in the first step of a patch), the Robin rows of the free surface eliminated into the top Laplace
//  rows, and the corrector passes with N at t^n+1/2 trapezoidal.  The patches have no relaxation zones (B 96),
//  coastline (fnpf_amr: the patches stay off it) or bodies (refused by fnpf_theta with G 1), so their Robin rows are
//  those of level 0 with r = 1, T = 0.  Breaking (A 350 1-3) as on level 0 (patch 0022): the damping with the
//  breaking viscosity at the start of the step, level by level (theta_damping), the detection and the filters at its
//  end (theta_breaking; an edit there: the consistency solve of fnpf_theta).  The rows of all grids are solved together (lap_core): the
//  in-column couplings of the sz stencils of 5 and 7 points (xcol) are in the operator (lap_apply) and lumped into
//  the k-1 coupling in the V-cycles of the preconditioner (lap_prepare), as in REEFMG.  The state of a step is
//  made consistent between the grids after each update (theta_surface: restriction into the covered columns, the
//  cells around the patches filled), as the RK3 stages do (stage_surface).

namespace
{
typedef reefamr_comms_off comms_off;
}

struct fnpf_amr_theta
{
    explicit fnpf_amr_theta(lexer *pp) : wn(pp),Ne(pp),Nf(pp),Neo(pp),Nfo(pp),eko(pp),fko(pp),nEm(pp),nFm(pp),R0(pp),Rs(pp),
                                         fnew(pp),wnew(pp),colA(pp),colB(pp),colD(pp),fsold(pp),szw(pp)
    {
        int np = 3;
        if(pp->A311==2 || pp->A311==4 || pp->A311==5) np = 5;
        if(pp->A311==6) np = 7;
        while(np>3 && pp->knoz<np-1) np -= 2;
        wsz = szw.weights(np);
        npt = np;
    }

    slice4 wn, Ne, Nf, Neo, Nfo, eko, fko, nEm, nFm, R0, Rs, fnew, wnew, colA, colB, colD, fsold;
    vector<double> mts, fiold;
    fnpf_sz_weights szw;            // the sz stencil of the patch's sigma grid (G 6 1: more layers)
    const double *wsz = nullptr;
    int npt = 0;
    int nstep = 0;                  // steps of the patch (0: Euler predictor)
    double dt_prev = 0.0;
};

void fnpf_amr::theta_delete(fnpf_amr_patch &c)
{
    delete c.th;
    c.th = nullptr;
}

// tendencies at t^n, N, the predictor (erk1, frk1) and R of pass 0 of every patch
void fnpf_amr::theta_predict(lexer *p0p, ghostcell *pgc, double theta, double beta)
{
    if(!active())
    return;

    const double t0 = MPI_Wtime();
    const double g = fabs(p0p->W22);

    comms_off guard(pgc);

    for(auto q : P)
    {
        fnpf_amr_patch *pc = FP(q);
        lexer *pp = pc->pp;
        fdm_fnpf *cc = pc->c;

        pp->dt = p0p->dt;
        pp->dt_old = p0p->dt_old;
        pp->simtime = p0p->simtime;
        pp->count = p0p->count;
        pp->amr_follow(*p0p);

        if(pc->th==nullptr)
        pc->th = new fnpf_amr_theta(pp);
        theta_mode = true;
        fnpf_amr_theta &T = *pc->th;

        if(T.wsz==nullptr)
        {
            cout<<"FNPF AMR, A 310 5: no sz stencil with "<<T.npt<<" points for "<<pp->knoz<<" sigma layers on a patch"<<endl;
            MPI_Abort(MPI_COMM_WORLD,-2770);
        }

        lexer *p = pp;
        fdm_fnpf *c = cc;
        const double dt = p->dt;
        slice4 &ek = *pc->ek, &fk = *pc->fk, &erk1 = *pc->erk1, &frk1 = *pc->frk1;

        // w^n (a new patch has no Fz of its own yet: from Fi as at the end of a step; Fi is that of the surface
        // before the damping, as w^n of level 0)
        pc->pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);

        // eta at the start of the step (breaking detection, after the damping)
        if(brk)
        {
            slice4 &en = *pc->en;
            SLICELOOP4
            en(i,j) = c->eta(i,j);
        }

        pc->pf->kfsfbc(p,c,pgc);
        SLICELOOP4
        ek(i,j) = c->K(i,j);

        pc->pf->dfsfbc(p,c,pgc,c->eta);
        SLICELOOP4
        fk(i,j) = c->K(i,j);

        double ab0=1.0, ab1=0.0;
        if(T.nstep>0 && T.dt_prev>0.0)
        {
            const double r = dt/T.dt_prev;
            ab0 = 1.0 + 0.5*r;
            ab1 = -0.5*r;
        }

        SLICELOOP4
        {
            T.wn(i,j) = c->Fz(i,j);
            T.Ne(i,j) = ek(i,j) - c->Fz(i,j);
            T.Nf(i,j) = fk(i,j) + g*c->eta(i,j);

            T.nEm(i,j) = ab0*T.Ne(i,j) + ab1*T.Neo(i,j);
            T.nFm(i,j) = ab0*T.Nf(i,j) + ab1*T.Nfo(i,j);

            T.R0(i,j) = c->Fifsf(i,j) - g*dt*c->eta(i,j) - g*(theta-beta)*dt*dt*T.wn(i,j);
            T.Rs(i,j) = T.R0(i,j) - g*theta*dt*dt*T.nEm(i,j) + dt*T.nFm(i,j);

            erk1(i,j) = c->eta(i,j) + dt*(ab0*ek(i,j) + ab1*T.eko(i,j));
            frk1(i,j) = c->Fifsf(i,j) + dt*(ab0*fk(i,j) + ab1*T.fko(i,j));

            T.eko(i,j) = ek(i,j);
            T.fko(i,j) = fk(i,j);
        }
    }

    tm[0] += MPI_Wtime()-t0;
}

// the state s of all grids consistent (Se, Sf: that of level 0): restriction into the covered columns, the cells
// around the patches from the coarser grid or a sibling
void fnpf_amr::theta_surface(lexer *p, ghostcell *pgc, slice &Se, slice &Sf, int s)
{
    Se0 = &Se;
    Sf0 = &Sf;
    stg = s;

    if(!active())
    return;

    const double t0 = MPI_Wtime();

    auto sel = [&](int g, int m) -> slice& { return sval(g,s,m); };

    restrict_sl(2,sel);
    restrict_flags(pgc);

    pgc->gcsl_start4(p,Se,gcval_eta);
    pgc->gcsl_start4(p,Sf,gcval_fifsf);

    for(int l=1; l<=maxlev; ++l)
    {
        fill_sl(l,2,7500+l,sel);

        for(int id : lev[l])
        {
            walls_sl(*FP(id),sel(id,0),gcval_eta);
            walls_sl(*FP(id),sel(id,1),gcval_fifsf);
        }
    }

    tm[0] += MPI_Wtime()-t0;
}

// the Robin rows of a patch (fnpf_theta::robin_build_mg with r = 1, T = 0, all columns wet): the coupling mt of the
// top row to the surface node, the surface node eliminated, the couplings beyond k-1 in xcol
void fnpf_amr::theta_rows(fnpf_amr_patch &pc, slice &fs, double beta)
{
    lexer *p = pc.pp;
    fdm_fnpf *c = pc.c;
    fnpf_amr_theta &T = *pc.th;
    const double g = fabs(p->W22);
    const int N = p->knoz;
    const double *w = T.wsz;

    double *Mp = c->M.p.data();
    double *Mb = c->M.b.data();
    double *rhs = c->rhsvec.V.data();

    c->M.xcol.clear();

    int n=0;
    ILOOP
    JLOOP
    {
        T.colD(i,j) = 0.0;

        KLOOP
        if(p->flag4[IJK]>0)
        {
            const int r = n;
            ++n;

            if(k!=N-1 || p->flag7[FIJK]<=0)
            continue;

            const double mt = T.mts[r];
            const double A = g*beta*p->dt*p->dt*p->sigz[IJ];
            const double D = 1.0 + A*w[0];
            const double B = T.Rs(i,j);

            T.colD(i,j) = D;
            T.colA(i,j) = A;
            T.colB(i,j) = B;

            rhs[r] += mt*fs(i,j) - mt*B/D;
            Mp[r] += -mt*A*w[1]/D;

            if(T.npt>2 && N>1)
            Mb[r] += -mt*A*w[2]/D;

            for(int m=3; m<T.npt; ++m)
            if(N-m>=0)
            {
                const double v = -mt*A*w[m]/D;
                if(v!=0.0)
                c->M.xcol.push_back({i,j,N-1,N-m,v});
            }
        }
    }
}

// the Robin problem of all grids: level 0 assembled by fnpf_theta (c->M, rhs, xcol), f its Fi; the patches on the
// sigma grid of their state s (0 pass 0, 1 corrector passes).  Fifsf^n+1 of a patch in fnew, w^n+1 in wnew; dmax,
// xmax: the largest change and value of Fi and Fs of the patch columns (local)
void fnpf_amr::theta_solve(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *f, int s, double beta, double &dmax, double &xmax)
{
    const double t0 = MPI_Wtime();

    dmax = xmax = 0.0;

    // the patches: free-surface derivatives, sigma grid, boundary conditions, the rows
    {
    comms_off guard(pgc);
    for(int id=0; id<(int)P.size(); ++id)
    {
        fnpf_amr_patch *pc = FP(id);
        lexer *pp = pc->pp;
        fdm_fnpf *cc = pc->c;
        fnpf_amr_theta &T = *pc->th;
        slice &Se = sval(id,s,0);
        slice &Sf = sval(id,s,1);

        pc->pf->fsfdisc(pp,cc,pgc,Se,Sf);
        pc->psig->sigma_update(pp,cc,pgc,pc->pf,Se);
        pc->pfu->fsfbc_sig(pp,cc,pgc,Sf,cc->Fi);
        pc->pbu->bedbc_sig(pp,cc,pgc,cc->Fi,pc->pf);
        walls_fi(*pc,cc->Fi);

        const double ts = MPI_Wtime();
        if(T.mts.size()!=cc->M.p.size())
        T.mts.assign(cc->M.p.size(),0.0);
        pc->plap->mtsurf = T.mts.data();
        pc->plap->start(pp,cc,pgc,nullptr,pc->pf,cc->Fi,Sf);
        pc->plap->mtsurf = nullptr;
        theta_rows(*pc,Sf,beta);
        lap_tsetup(MPI_Wtime()-ts);

        // Fi and Fs before the solve
        const size_t nf = (size_t)pp->imax*pp->jmax*(pp->kmax+2);
        if(T.fiold.size()!=nf)
        T.fiold.resize(nf);
        std::copy(cc->Fi,cc->Fi+nf,T.fiold.begin());
        lexer *p = pp;
        SLICELOOP4
        T.fsold(i,j) = Sf(i,j);
    }
    }

    ltgt.assign(P.size()+1,nullptr);
    ltgt[0] = f;
    for(int id=0; id<(int)P.size(); ++id)
    ltgt[id+1] = FP(id)->c->Fi;

    lap_kind = 0;
    lap_core(p,pgc);


    // the patches: Fs and w^n+1 of the Robin rows
    {
    comms_off guard(pgc);
    for(int id=0; id<(int)P.size(); ++id)
    {
        fnpf_amr_patch *pc = FP(id);
        lexer *pp = pc->pp;
        fdm_fnpf *cc = pc->c;
        fnpf_amr_theta &T = *pc->th;
        const int N = pp->knoz;
        walls_fi(*pc,cc->Fi);
        cc->M.xcol.clear();

        lexer *p = pp;
        fdm_fnpf *c = cc;
        SLICELOOP4
        {
            T.fnew(i,j) = T.fsold(i,j);
            T.wnew(i,j) = 0.0;
        }

        for(int ii=EXT; ii<EXT+pc->nx; ++ii)
        for(int jj=EY; jj<EY+pc->ny; ++jj)
        {
            i = ii;
            j = jj;
            if(T.colD(i,j)<=0.0)
            continue;

            double sm = 0.0;
            for(int m=1; m<T.npt; ++m)
            {
                k = N-m;
                sm += T.wsz[m]*c->Fi[FIJK];
            }

            const double Fs = (T.colB(i,j) - T.colA(i,j)*sm)/T.colD(i,j);
            T.fnew(i,j) = Fs;
            T.wnew(i,j) = p->sigz[IJ]*(T.wsz[0]*Fs + sm);

            dmax = MAX(dmax,fabs(Fs-T.fsold(i,j)));
            xmax = MAX(xmax,fabs(Fs));

            KLOOP
            if(p->flag4[IJK]>0)
            {
                dmax = MAX(dmax,fabs(c->Fi[FIJK]-T.fiold[FIJK]));
                xmax = MAX(xmax,fabs(c->Fi[FIJK]));
            }
        }
    }
    }

    tm[1] += MPI_Wtime()-t0;
}

// corrector state of the patches: eta^n+1 of the last solution (erk2) with Fifsf^n+1 (frk2 = fnew), all grids
// consistent (Se, Sf: level 0), then N^n+1, N^n+1/2 trapezoidal and R of the next pass
void fnpf_amr::theta_corrector(lexer *p0p, ghostcell *pgc, slice &Se, slice &Sf, double theta, double beta)
{
    if(!active())
    {
        Se0 = &Se;
        Sf0 = &Sf;
        stg = 1;
        return;
    }

    double t0 = MPI_Wtime();
    const double g = fabs(p0p->W22);

    {
    comms_off guard(pgc);
    for(auto q : P)
    {
        fnpf_amr_patch *pc = FP(q);
        fnpf_amr_theta &T = *pc->th;
        lexer *p = pc->pp;
        fdm_fnpf *c = pc->c;
        const double dt = p->dt;
        slice4 &erk2 = *pc->erk2, &frk2 = *pc->frk2;

        SLICELOOP4
        {
            erk2(i,j) = c->eta(i,j) + dt*(theta*T.wnew(i,j) + (1.0-theta)*T.wn(i,j) + T.nEm(i,j));
            frk2(i,j) = T.fnew(i,j);
            c->Fz(i,j) = T.wnew(i,j);
        }
    }
    }

    tm[0] += MPI_Wtime()-t0;

    theta_surface(p0p,pgc,Se,Sf,1);

    t0 = MPI_Wtime();

    comms_off guard(pgc);
    for(auto q : P)
    {
        fnpf_amr_patch *pc = FP(q);
        fnpf_amr_theta &T = *pc->th;
        lexer *p = pc->pp;
        fdm_fnpf *c = pc->c;
        const double dt = p->dt;
        slice4 &erk2 = *pc->erk2, &frk2 = *pc->frk2;

        pc->pf->fsfdisc(p,c,pgc,erk2,frk2);

        pc->pf->kfsfbc(p,c,pgc);
        SLICELOOP4
        T.nEm(i,j) = 0.5*(T.Ne(i,j) + c->K(i,j) - T.wnew(i,j));

        pc->pf->dfsfbc(p,c,pgc,erk2);
        SLICELOOP4
        {
            T.nFm(i,j) = 0.5*(T.Nf(i,j) + c->K(i,j) + g*erk2(i,j));
            T.Rs(i,j) = T.R0(i,j) - g*theta*dt*dt*T.nEm(i,j) + dt*T.nFm(i,j);
        }
    }

    tm[0] += MPI_Wtime()-t0;
}

// eta^n+1 and Fifsf^n+1 of the patches, the history of the predictor; then all grids consistent (theta_surface 2,
// called by fnpf_theta with level 0)
void fnpf_amr::theta_final(lexer *p0p, ghostcell *pgc, double theta)
{
    if(!active())
    return;

    const double t0 = MPI_Wtime();

    comms_off guard(pgc);
    for(auto q : P)
    {
        fnpf_amr_patch *pc = FP(q);
        fnpf_amr_theta &T = *pc->th;
        lexer *p = pc->pp;
        fdm_fnpf *c = pc->c;
        const double dt = p->dt;

        SLICELOOP4
        {
            c->eta(i,j) = c->eta(i,j) + dt*(theta*T.wnew(i,j) + (1.0-theta)*T.wn(i,j) + T.nEm(i,j));
            c->Fifsf(i,j) = T.fnew(i,j);
            T.Neo(i,j) = T.Ne(i,j);
            T.Nfo(i,j) = T.Nf(i,j);
        }

        T.dt_prev = dt;
        ++T.nstep;
    }

    tm[0] += MPI_Wtime()-t0;
}

// the derived fields of the patches at the end of the step: free-surface derivatives, the sigma grid of eta^n+1,
// the boundary values of Fi and w^n+1 (as fnpf_theta on level 0)
void fnpf_amr::theta_end(lexer *p0p, ghostcell *pgc)
{
    if(!active())
    return;

    const double t0 = MPI_Wtime();

    comms_off guard(pgc);
    for(auto q : P)
    {
        fnpf_amr_patch *pc = FP(q);
        lexer *pp = pc->pp;
        fdm_fnpf *cc = pc->c;

        pc->pf->fsfdisc(pp,cc,pgc,cc->eta,cc->Fifsf);
        pc->psig->sigma_update(pp,cc,pgc,pc->pf,cc->eta);
        pc->pfu->fsfbc_sig(pp,cc,pgc,cc->Fifsf,cc->Fi);
        pc->pbu->bedbc_sig(pp,cc,pgc,cc->Fi,pc->pf);
        walls_fi(*pc,cc->Fi);
        pc->pf->fsfwvel(pp,cc,pgc,cc->eta,cc->Fifsf);
    }

    tm[0] += MPI_Wtime()-t0;
}

// A 350 1-3: the damping of eta and Fifsf with the breaking viscosity at the start of the step, as fnpf_theta on level 0
// (patch 0022: a source of the step, no solve after it), level by level: the cells around the patches of a level hold
// the damped values of the coarser grid (level_damping), then all grids consistent and the derivatives of the patches
bool fnpf_amr::theta_damping(lexer *p, ghostcell *pgc)
{
    if(!active() || !brk)
    return false;

    double vbm = 0.0;
    for(auto q : P)
    {
        fnpf_amr_patch *pc = FP(q);
        for(int ii=EXT; ii<EXT+q->nx; ++ii)
        for(int jj=EY; jj<EY+q->ny; ++jj)
        vbm = MAX(vbm,fabs(pc->c->vb(ii,jj)));
    }
    vbm = pgc->globalmax(vbm);

    if(vbm<=0.0)
    return false;

    const double t0 = MPI_Wtime();

    auto sel = [&](int g, int m) -> slice& { return sval(g,2,m); };
    Se0 = &c0->eta;
    Sf0 = &c0->Fifsf;
    stg = 2;

    for(int l=1; l<=maxlev; ++l)
    {
        fill_sl(l,2,7640+l,sel);

        for(int id : lev[l])
        {
            walls_sl(*FP(id),sel(id,0),gcval_eta);
            walls_sl(*FP(id),sel(id,1),gcval_fifsf);
        }

        level_damping(l,pgc,2,0,gcval_eta,1.0);
        level_damping(l,pgc,2,1,gcval_fifsf,1.0);
    }

    tm[0] += MPI_Wtime()-t0;

    theta_surface(p,pgc,c0->eta,c0->Fifsf,2);

    comms_off guard(pgc);
    for(auto q : P)
    {
        fnpf_amr_patch *pc = FP(q);
        pc->pf->fsfdisc(pc->pp,pc->c,pgc,pc->c->eta,pc->c->Fifsf);
    }

    return true;
}

// A 350 1-3: the breaking detection and the filters (A 352) on the patches at the end of the step, as fnpf_theta on
// level 0 (the cells around the patches then take the viscosity of their source grid); true if eta or Fifsf of a patch
// changed: all grids consistent again, and fnpf_theta makes its consistency solve
bool fnpf_amr::theta_breaking(lexer *p, ghostcell *pgc)
{
    if(!active() || !brk)
    return false;

    {
    comms_off guard(pgc);
    for(auto q : P)
    {
        fnpf_amr_patch *pc = FP(q);
        lexer *p = pc->pp;
        fdm_fnpf *c = pc->c;
        slice4 &e1 = *pc->erk1, &f1 = *pc->frk1;
        SLICELOOP4
        {
            e1(i,j) = c->eta(i,j);
            f1(i,j) = c->Fifsf(i,j);
        }
    }
    }

    patch_breaking(-1,pgc,2,1.0);

    double emax = 0.0;
    for(auto q : P)
    {
        fnpf_amr_patch *pc = FP(q);
        fdm_fnpf *c = pc->c;
        slice4 &e1 = *pc->erk1, &f1 = *pc->frk1;
        for(int ii=EXT; ii<EXT+q->nx; ++ii)
        for(int jj=EY; jj<EY+q->ny; ++jj)
        emax = MAX(emax, MAX(fabs(c->eta(ii,jj)-e1(ii,jj)), fabs(c->Fifsf(ii,jj)-f1(ii,jj))));
    }
    emax = pgc->globalmax(emax);

    if(emax<=0.0)
    return false;

    theta_surface(p,pgc,c0->eta,c0->Fifsf,2);

    return true;
}
