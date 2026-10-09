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
#include"matrix2D.h"
#include"vec2D.h"
#include"sflow_pressure_nh.h"
#include"sflow_momentum_RK3.h"
#include<cmath>
#include<mpi.h>
#include<algorithm>

//  Subcycling in time (G 7 1) with the non-hydrostatic pressure (A 220 1-3), as nhflow_amr 0019:
//
//  Every level solves the pressure of its stages alone: level 0 in its own step (nh_solve0, called
//  by sflow_pjm_lin/quad), the level-l patches after each of their stages (nh_level), with the
//  composite solver of sflow_amr_nh.cpp on the window [l,l].  The cells of level l under level
//  l+1 are fixed at the restricted pressure of level l+1 (from the last synchronisation): as in the
//  composite solve, the cells next to the patches see the fine pressure, and the covered cells
//  are not solved with the restricted fine velocities, which do not satisfy the coarse constraint
//  (with them as rows the coarse pressure took up their residual in every level-0 step: the bar
//  wave at x = 11 m 0.094 mm rms from G 7 0 instead of 0.036 mm, and more with smaller steps).
//  The cells around the level-l patches that the parent fills are fixed at the parent pressure
//  linear in time between the start and the end of the parent step (at the time of the stage
//  output), siblings coupled.
//
//  After the two steps of level l+1, refluxing and restriction (sflow_amr_sub.cpp), the
//  synchronisation projection over the levels l..L (nh_project): one composite correction solve
//  for the change of the constraint residual (-(WL div u + 2 (w + u.grad d))/dt, the linear part of
//  the rows of sflow_pjm_lin/quad) since the end of each grid's own step: the residual its own
//  projection leaves is part of the scheme on one grid, the change comes from the coupling
//  (refluxed cells, restricted cells, the cells around the patches).  The parent cells of level
//  l > 0 are fixed at 0; the correction of UH, VH, WH as sflow_pjm_lin/quad without the bed
//  acceleration of the quadratic pressure; the pressure itself is kept.

namespace
{
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

typedef reefamr_comms_off comms_off;

// times of the RK3 stage outputs in the step (SSP RK3: t+dt, t+dt/2, t+dt)
const double rko[3] = {1.0, 0.5, 1.0};
}

// the pressure of the level-l stage s of step k (0, 1)
void sflow_amr::nh_level(lexer *p, ghostcell *pgc, int l, int k, int s)
{
    const double t0 = MPI_Wtime();

    // cells around the patches: the parent pressure at the time of the stage output (siblings:
    // the previous pressure, overwritten by the solve)
    const double th = 0.5*(double(k)+rko[s]);
    fill_run(l,1,7350+l,
             [&](const reefamr_fill &f, double *v) { v[0] = nh_eval_t(f,th,4); },
             [&](reefamr_patch *c, int n, const reefamr_fill &f, const double *w) { nh_vec(n,-1)(f.di,f.dj) = w[0]; });

    wlo = whi = l;
    nh_edge = true;
    nh_prepare(pgc,true);
    nh_core(p);
    nh_edge = false;
    wlo = 0;
    whi = -1;
    sub_lv_it += nh_it_last;
    ++sub_lv_n;

    // velocity correction and the rest of the stage
    {
    comms_off guard(pgc);
    for(int n : lev[l])
    {
        sflow_amr_patch *c = SP(n);
        sflow_momentum_RK3 *m = c->pmom;
        pgc->gcsl_start4(c->pp,c->b->press,c->pnh->gcval());
        c->pnh->correct(c->pp,c->b,*m->nhUH,*m->nhVH,*m->nhWH,*m->nhWL,m->nh_alpha);
        m->stage_finish(c->pp,c->b,pgc,*m->nhUH,*m->nhVH,*m->nhWH,*m->nhWL);
    }
    }

    tm[5] += MPI_Wtime()-t0;
}

// constraint residual WL div u + 2 (w + u.grad d) of the cells of grid g with a non-hydrostatic row
void sflow_amr::nh_resid(int g, vector<double> &R)
{
    gh G = grid(g);
    lexer *q = G.q;
    fdm2D *pb = G.b;
    sflow_pressure_nh *pn = (g<0) ? pnh0 : SP(g)->pnh;

    R.assign((size_t)q->imax*q->jmax,0.0);

    int i0,i1,j0,j1;
    if(g<0) { i0=0; i1=NX0-1; j0=0; j1=NY0-1; }
    else    { i0=EXT; i1=EXT+P[g]->nx-1; j0=EXT; j1=EXT+P[g]->ny-1; }

    for(int ii=i0; ii<=i1; ++ii)
    for(int jj=j0; jj<=j1; ++jj)
    {
        if(q->flagslice4[lij(q,ii,jj)]<=0 || pn->is_active(q,pb,ii,jj)==0)
        continue;

        const double dx = q->DXP[ii+marge] + q->DXP[ii-1+marge];
        const double dy = q->DYP[jj+marge] + q->DYP[jj-1+marge];
        const double dudx = (pb->U(ii+1,jj) - pb->U(ii-1,jj))/dx;
        const double dvdy = (pb->V(ii,jj+1) - pb->V(ii,jj-1))/dy*p0->y_dir;
        const double dddx = (pb->depth(ii+1,jj) - pb->depth(ii-1,jj))/dx;
        const double dddy = (pb->depth(ii,jj+1) - pb->depth(ii,jj-1))/dy*p0->y_dir;

        R[lij(q,ii,jj)] = pb->WL(ii,jj)*(dudx + dvdy) + 2.0*(pb->W(ii,jj) + pb->U(ii,jj)*dddx + pb->V(ii,jj)*dddy);
    }
}

void sflow_amr::nh_store_end(int g)
{
    nh_resid(g, (g<0) ? nh0_rend : SP(g)->rend);
}

// synchronisation projection of the levels l..L at the end of a level-l step (dt: its step)
void sflow_amr::nh_project(lexer *p, ghostcell *pgc, int l, double dt)
{
    const double t0 = MPI_Wtime();

    vector<int> grids;
    if(l==0)
    grids.push_back(-1);
    for(int m=MAX(l,1); m<=maxlev; ++m)
    for(int id : lev[m])
    grids.push_back(id);

    // cells around the patches of the levels above l at this time (all synchronous); the cells
    // around level l keep the values of its last fill (the parent is fixed)
    auto refill = [&]()
    {
        cache_stage(0);
        for(int m=l+1; m<=maxlev; ++m)
        fill_level(pgc,m,0);
    };
    refill();

    // rows of the window grids for the step of level l: the pressure of the grids is the unknown
    // (the correction, starting at 0; kept and restored afterwards)
    vector<vector<double>> psave(grids.size());
    vector<double> dts(grids.size());
    vector<double> Rnow;

    for(size_t n=0; n<grids.size(); ++n)
    {
        const int g = grids[n];
        gh G = grid(g);
        lexer *q = G.q;
        fdm2D *pb = G.b;
        sflow_pressure_nh *pn = (g<0) ? pnh0 : SP(g)->pnh;

        dts[n] = q->dt;
        q->dt = dt;

        psave[n].assign(pb->press.begin(),pb->press.end());
        std::ranges::fill(pb->press,0.0);

        if(g<0)
        pn->assemble(q,pb,pgc,pb->UH,pb->VH,pb->WL,pb->U,pb->V,1.0);
        else
        {
            comms_off guard(pgc);
            pn->assemble(q,pb,pgc,pb->UH,pb->VH,pb->WL,pb->U,pb->V,1.0);
        }

        // right-hand side: the change of the residual since the end of the grid's own step
        nh_resid(g,Rnow);
        const vector<double> &Re = (g<0) ? nh0_rend : SP(g)->rend;
        const bool have = (Re.size()==Rnow.size());

        int r=0;
        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        if(q->flagslice4[lij(q,ii,jj)]>0)
        {
            const int c = lij(q,ii,jj);
            pb->rhsvec.V[r] = (have && pn->is_active(q,pb,ii,jj)) ? -(Rnow[c]-Re[c])/dt : 0.0;
            ++r;
        }
    }

    wlo = l;
    whi = -1;
    nh_edge = (l>0);
    nh_prepare(pgc,false);
    nh_core(p);
    nh_edge = false;
    wlo = 0;
    whi = -1;
    sub_sy_it += nh_it_last;
    ++sub_sy_n;

    // velocity correction (as sflow_pjm_lin/quad, without the bed acceleration), velocities
    for(size_t n=0; n<grids.size(); ++n)
    {
        const int g = grids[n];
        gh G = grid(g);
        lexer *q = G.q;
        fdm2D *pb = G.b;
        sflow_pressure_nh *pn = (g<0) ? pnh0 : SP(g)->pnh;
        double cbq, cwq;
        pn->coef(cbq,cwq);

        auto correct = [&]()
        {
            pgc->gcsl_start4(q,pb->press,pn->gcval());

            int i0,i1,j0,j1;
            if(g<0) { i0=0; i1=NX0-1; j0=0; j1=NY0-1; }
            else    { i0=EXT; i1=EXT+P[g]->nx-1; j0=EXT; j1=EXT+P[g]->ny-1; }

            slice &Pc = pb->press;
            slice &WL = pb->WL;
            const double f = dt/p0->W1;

            for(int ii=i0; ii<=i1; ++ii)
            for(int jj=j0; jj<=j1; ++jj)
            {
                if(q->flagslice4[lij(q,ii,jj)]<=0 || pn->is_active(q,pb,ii,jj)==0)
                continue;

                const double dx = q->DXP[ii+marge] + q->DXP[ii-1+marge];
                pb->UH(ii,jj) -= f*((WL(ii+1,jj)*Pc(ii+1,jj) - WL(ii-1,jj)*Pc(ii-1,jj))/dx
                                    - cbq*Pc(ii,jj)*(pb->depth(ii+1,jj) - pb->depth(ii-1,jj))/dx);

                if(p0->j_dir==1)
                {
                const double dy = q->DYP[jj+marge] + q->DYP[jj-1+marge];
                pb->VH(ii,jj) -= f*((WL(ii,jj+1)*Pc(ii,jj+1) - WL(ii,jj-1)*Pc(ii,jj-1))/dy
                                    - cbq*Pc(ii,jj)*(pb->depth(ii,jj+1) - pb->depth(ii,jj-1))/dy);
                }

                pb->WH(ii,jj) += f*cwq*Pc(ii,jj);
            }

            G.m->velcalc(q,pb,pgc,pb->UH,pb->VH,pb->WH,pb->WL,2);
        };

        if(g<0)
        correct();
        else
        {
            comms_off guard(pgc);
            correct();
        }

        std::ranges::copy(psave[n],pb->press.begin());
        q->dt = dts[n];
    }

    // the levels above l into level l again, the level-0 halo, the cells around the patches
    for(int m=maxlev; m>=l+1; --m)
    restrict_level(p0,m,2);

    if(l==0)
    exchange_level0(p0,b0,pgc,2);

    refill();

    // the residual of the coupled state is the reference of the next synchronisation
    for(int g : grids)
    nh_store_end(g);

    tm[5] += MPI_Wtime()-t0;
}
