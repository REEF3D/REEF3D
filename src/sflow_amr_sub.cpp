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

#include"sflow_amr.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"slice.h"
#include"sflow_eta.h"
#include"sflow_fsf.h"
#include"sflow_momentum_RK3.h"
#include<cmath>
#include<cstring>
#include<algorithm>
#include<mpi.h>

//  Subcycling in time (G 7 1, hydrostatic SFLOW A 220 0): Berger-Oliger with refluxing.
//
//  Level 0 takes its RK3 step alone (sflow_momentum_RK3::start, the patches are not touched in
//  the stage hooks).  Its faces next to a patch keep their own flux; hll_hook sums it over the
//  stages with the RK3 weights of the step (1/6, 1/6, 2/3) times dt into the coarse register of
//  the face (cregL local patches, cregR patches of other ranks).  The face depth is the fine one
//  as in the synchronous coupling, so the coarse cell stays well balanced.
//
//  Then every level l >= 1 takes two steps of half the step of level l-1 (sub_level, recursive:
//  after each of its steps the finer levels take their two).  The cells around a level-l patch
//  are filled at the stage times (RK3 stage inputs at t, t+dt, t+dt/2) from the parent state
//  linear in time between the start of the parent step (told, kept at its start) and its end
//  (the current parent state): tint = (k + c_s)/2 for step k = 0, 1 of the patch.  A patch sums
//  its box face fluxes into its fine register (freg) in the same way.
//
//  After the two steps of level l+1: the coarse cells of level l next to a level-(l+1) patch take
//  the difference of the summed fine flux (mean of the two fine faces) and the summed coarse
//  flux of the face (sub_reflux: continuity, both momentum components), then the level-(l+1)
//  patches are restricted into level l.  Mass and momentum of the face fluxes are conserved
//  across the coarse-fine interface (as on one grid, the wet-dry clipping and the sources are
//  not).
//
//  The time step of level 0 is limited by every level at its own CFL number times 2^l
//  (sflow_amr::timestep).  The regrid runs at the end of a level-0 step, all levels synchronous.

namespace
{
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

typedef reefamr_comms_off comms_off;

// times of the RK3 stage inputs in the step (SSP RK3: t, t+dt, t+dt/2)
const double rkc[3] = {0.0, 1.0, 0.5};
}

sflow_amr_told& sflow_amr::told(int g)
{
    if(g<0)
    return told0;
    return SP(g)->told;
}

void sflow_amr::told_free(sflow_amr_told &T)
{
    for(int k=0; k<4; ++k)
    {
        delete T.f[k];
        T.f[k] = nullptr;
    }
    T.wet.clear();
}

// state of grid g at the start of its step
void sflow_amr::sub_snapshot(int g)
{
    gh G = grid(g);
    sflow_amr_told &T = told(g);
    const size_t n = (size_t)G.q->imax*G.q->jmax;

    slice *src[4] = {&G.b->WL,&G.b->UH,&G.b->VH,&G.b->WH};
    for(int k=0; k<4; ++k)
    {
        if(T.f[k]==nullptr)
        T.f[k] = new slice(G.q);
        memcpy(T.f[k]->V,src[k]->V,n*sizeof(double));
    }
    T.wet.assign(G.q->wet,G.q->wet+n);
}

void sflow_amr::sub_creg_reset(int g)
{
    const int t = g+1;
    cregL[t].assign(5*match[t].size(),0.0);
    cregR[t].assign(5*rmatch[t].size(),0.0);
}

// fine face depths of the patch boxes (static bed) for the coarse faces next to them: the coarse
// grids run before the patches in a step
void sflow_amr::sub_dfx()
{
    for(auto q : P)
    {
        sflow_amr_patch &c = *SP(q);
        fdm2D *b = c.b;
        const int il = EXT-1, ih = EXT+c.nx-1;
        const int jl = EXT-1, jh = EXT+c.ny-1;

        for(int r=0; r<c.ny; ++r)
        {
            c.rec[0][0][r] = b->dfx(il,EXT+r);
            c.rec[0][1][r] = b->dfx(ih,EXT+r);
        }
        for(int r=0; r<c.nx; ++r)
        {
            c.rec[0][2][r] = b->dfy(EXT+r,jl);
            c.rec[0][3][r] = b->dfy(EXT+r,jh);
        }
    }

    for(int l=1; l<=maxlev; ++l)
    exchange_fluxes(l);
}

// start of a level-0 step (sflow_momentum_RK3::start, before the level-0 stages)
void sflow_amr::sub_begin(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double t0 = MPI_Wtime();

    cregL.resize(P.size()+1);
    cregR.resize(P.size()+1);
    rflux.resize(P.size()+1);
    for(size_t t=0; t<rflux.size(); ++t)
    rflux[t].assign(5*rmatch[t].size(),0.0);

    sub_dfx();

    sub_snapshot(-1);
    sub_creg_reset(-1);
    hstage = 0;
    ++sub_steps[0];

    tm[2] += MPI_Wtime()-t0;
}

// end of a level-0 step (step_end): the patches take their steps
void sflow_amr::sub_end(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const double t = p->simtime, dt = p->dt;

    for(int n : lev[1])
    for(int ip=0; ip<5; ++ip)
    for(int side=0; side<4; ++side)
    std::fill(SP(n)->freg[ip][side].begin(),SP(n)->freg[ip][side].end(),0.0);

    sub_level(p,pgc,1,0,t,0.5*dt);
    sub_level(p,pgc,1,1,t+0.5*dt,0.5*dt);

    sub_sync(p,pgc,0);
}

// step k (0, 1) of the level-l patches from time t with dt, within the step of level l-1; then
// the steps of the finer levels, refluxing and restriction into level l
void sflow_amr::sub_level(lexer *p, ghostcell *pgc, int l, int k, double t, double dt)
{
    double t0;

    {
    comms_off guard(pgc);
    for(int n : lev[l])
    {
        sflow_amr_patch *c = SP(n);
        c->pp->dt = dt;
        c->pp->dt_old = dt;
        c->pp->simtime = t;
        c->pp->count = p->count;
        c->pmom->inflow(c->pp,c->b,pgc,pflow_void);
    }
    }

    if(l<=7)
    ++sub_steps[l];

    for(int s=0; s<3; ++s)
    {
        // cells around the patches at the stage time, the parent interpolated in time
        t0 = MPI_Wtime();
        cache_stage(s);
        tint = 0.5*(double(k)+rkc[s]);
        fill_level(pgc,l,s);
        tint = -1.0;

        // start state of the step (with the filled cells) for the fills of the finer level
        if(s==0 && l<maxlev)
        for(int n : lev[l])
        {
            sub_snapshot(n);
            sub_creg_reset(n);
        }
        tm[0] += MPI_Wtime()-t0;

        t0 = MPI_Wtime();
        hstage = s;
        {
        comms_off guard(pgc);
        for(int n : lev[l])
        SP(n)->pmom->rk_stage(SP(n)->pp,SP(n)->b,pgc,s);
        }
        tm[1] += MPI_Wtime()-t0;
    }

    // end of the step as for level 0 in sflow_f::mainloop
    {
    comms_off guard(pgc);
    for(int n : lev[l])
    {
        sflow_amr_patch *c = SP(n);
        lexer *pp = c->pp;
        c->pmom->rk_finish(pp,c->b,pgc);

        if(p->A248==0)
        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        c->b->breaking(ii,jj)=0;

        c->pfsf->depth_update(pp,c->b,pgc,c->b->WL);
    }
    }

    if(l>=maxlev)
    return;

    // cells around the patches at the end of the step: the end state of the parents of level l+1
    t0 = MPI_Wtime();
    cache_stage(0);
    tint = 0.5*double(k+1);
    fill_level(pgc,l,0);
    tint = -1.0;
    tm[0] += MPI_Wtime()-t0;

    for(int n : lev[l+1])
    for(int ip=0; ip<5; ++ip)
    for(int side=0; side<4; ++side)
    std::fill(SP(n)->freg[ip][side].begin(),SP(n)->freg[ip][side].end(),0.0);

    sub_level(p,pgc,l+1,0,t,0.5*dt);
    sub_level(p,pgc,l+1,1,t+0.5*dt,0.5*dt);

    sub_sync(p,pgc,l);
}

// level l after the two steps of level l+1: refluxing, restriction (level 0: halo)
void sflow_amr::sub_sync(lexer *p, ghostcell *pgc, int l)
{
    double t0 = MPI_Wtime();

    sub_reflux(l);
    restrict_level(p0,l+1,2);

    if(l==0)
    exchange_level0(p0,b0,pgc,2);

    tm[3] += MPI_Wtime()-t0;
}

// coarse cells of level l next to the level-(l+1) patches: the summed fine face fluxes replace
// the summed coarse ones
void sflow_amr::sub_reflux(int l)
{
    // fine sums of patches on other ranks
    face_run_at(l+1,5,7250+l+1,
                [&](reefamr_patch *q, int side, int r, double *sb)
                {
                    sflow_amr_patch *c = SP(q);
                    for(int ip=0; ip<5; ++ip)
                    sb[ip] = 0.5*(c->freg[ip][side][r]+c->freg[ip][side][r+1]);
                },
                [&](int tgt, int idx, const double *v)
                {
                    for(int ip=0; ip<5; ++ip)
                    rflux[tgt][5*idx+ip] = v[ip];
                });

    vector<int> grids;
    if(l==0)
    grids.push_back(-1);
    else
    grids = lev[l];

    const double wd = p0->A244;

    for(int g : grids)
    {
        const int t = g+1;
        if(match[t].empty() && rmatch[t].empty())
        continue;

        gh G = grid(g);
        lexer *q = G.q;
        fdm2D *pb = G.b;

        // corrected cell: outside the patch, the face is its high side for side 0/2, its low side
        // for side 1/3; ipol 4 continuity (WL), 1 UH, 2 VH, 3 WH
        auto fix = [&](const reefamr_match &m, const double *d)
        {
            const int ii = m.fi + (m.side==1 ? 1 : 0);
            const int jj = m.fj + (m.side==3 ? 1 : 0);
            const double sg = (m.side==0 || m.side==2) ? -1.0 : 1.0;
            const double f = (m.dir==0) ? sg/q->DXN[ii+marge] : sg*p0->y_dir/q->DYN[jj+marge];

            pb->WL(ii,jj) += f*d[4];
            pb->UH(ii,jj) += f*d[1];
            pb->VH(ii,jj) += f*d[2];
            pb->WH(ii,jj) += f*d[3];

            // derived values as after the restriction
            const double wl = pb->WL(ii,jj);
            const int w = wl>wd+eps ? 1 : 0;
            q->wet[lij(q,ii,jj)] = w;
            const double wlvl = fabs(wl)>wd ? wl : 1.0e20;
            pb->eta(ii,jj) = wl - pb->depth(ii,jj);
            pb->U(ii,jj) = w==1 ? pb->UH(ii,jj)/wlvl : 0.0;
            pb->V(ii,jj) = w==1 ? pb->VH(ii,jj)/wlvl : 0.0;
            pb->W(ii,jj) = w==1 ? pb->WH(ii,jj)/wlvl : 0.0;
            pb->hp(ii,jj) = wl;
        };

        double d[5];
        for(size_t k=0; k<match[t].size(); ++k)
        {
            const reefamr_match &m = match[t][k];
            sflow_amr_patch &c = *SP(m.child);
            for(int ip=0; ip<5; ++ip)
            d[ip] = 0.5*(c.freg[ip][m.side][m.r]+c.freg[ip][m.side][m.r+1]) - cregL[t][5*k+ip];
            fix(m,d);
        }

        for(size_t k=0; k<rmatch[t].size(); ++k)
        {
            const reefamr_match &m = rmatch[t][k];
            for(int ip=0; ip<5; ++ip)
            d[ip] = rflux[t][5*k+ip] - cregR[t][5*k+ip];
            fix(m,d);
        }
    }
}
