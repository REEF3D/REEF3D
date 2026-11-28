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

#include"sflow_curv.h"
#include"sflow_amr.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"slice.h"
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
//  (with them as rows the coarse pressure would take up their residual in every level-0 step).
//  The cells around the level-l patches that the parent fills are fixed at the parent pressure
//  linear in time between the start and the end of the parent step (at the time of the stage
//  output), siblings coupled.
//
//  The preconditioner of a level solve is the patch-local part of sflow_amr_nh.cpp, with or without the
//  level-0 part before it (the residual restricted down to level 0, one V-cycle there, the correction
//  interpolated back up).  The level-0 part saves iterations where the patches of the level form
//  a wide block and the pressure reaches far (deep water for the cell size), and costs a V-cycle of
//  the whole level-0 grid in every application; which of the two is cheaper depends on the case
//  and changes with the patches.  Every level measures it: of every NHCYC level solves the first
//  NHPRB run the other way, the next NHPRB the chosen way, and the choice changes where the other
//  way needs less work (iterations times the cells an iteration works on, nh_c0_choose).  The choice
//  and the counters are part of the checkpoint; the synchronisation projection of the levels l..L
//  (l > 0) has the level-0 part where one of its levels chose it.
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

// level solves per cycle of the choice of the level-0 part, and the solves measured each way (the
// three stages of a step)
const int NHCYC = 96, NHPRB = 3;
}

// the pressure of the level-l stage s of step k (0, 1)
void sflow_amr::nh_level(lexer *p, ghostcell *pgc, int l, int k, int s)
{
    const double t0 = MPI_Wtime();

    // cells around the patches: the parent pressure at the time of the stage output (siblings:
    // the previous pressure, overwritten by the solve)
    const double th = (double(k)+rko[s])/double(msub[l]);
    fill_run(l,1,7350+l,
             [&](const reefamr_fill &f, double *v) { v[0] = nh_eval_t(f,th,4); },
             [&](reefamr_patch *c, int n, const reefamr_fill &f, const double *w) { nh_vec(n,-1)(f.di,f.dj) = w[0]; });

    // the level-0 part of the preconditioner: the choice of the level, the other way in the first
    // solves of a cycle
    const int cyc = nh_cn[l];
    nh_coarse = ((cyc<NHPRB ? 1-nh_c0[l] : nh_c0[l])==1);

    wlo = whi = l;
    nh_edge = true;
    nh_prepare(pgc,true);
    nh_core(p);
    nh_edge = false;
    wlo = 0;
    whi = -1;
    nh_coarse = true;
    sub_lv_it += nh_it_last;
    ++sub_lv_n;

    if(cyc<2*NHPRB)
    nh_ci[l][cyc<NHPRB ? 0 : 1] += nh_it_last;
    ++nh_cn[l];
    if(nh_cn[l]==2*NHPRB)
    nh_c0_choose(pgc,l);
    if(nh_cn[l]>=NHCYC)
    nh_cn[l] = nh_ci[l][0] = nh_ci[l][1] = 0;

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

// the level-0 part in the preconditioner of the level-l solves: the work of a solve is its iterations
// times the cells an application of the preconditioner works on, the cells of the level-l patches
// in the patch-local passes and the rest of the iteration, and with the level-0 part the cells of
// level 0 and of the levels between (their V-cycle, restriction and interpolation: 1.5 per cell,
// measured against a patch-local pass)
void sflow_amr::nh_c0_choose(ghostcell *pgc, int l)
{
    vector<double> w(l+1,0.0);
    for(auto c : P)
    if(c->lev>=1 && c->lev<=l)
    w[c->lev] += double(c->nx)*double(c->ny);
    for(int k=1; k<=l; ++k)
    w[k] = pgc->globalsum(w[k]);
    w[0] = double(p0->cellnumtot2D);

    const double own = w[l]*double((nlevg[l]>1 ? 4 : 1) + 1);
    double below = 0.0;
    for(int k=0; k<l; ++k)
    below += 1.5*w[k];

    const double f[2] = {own, own+below};
    const int cur = nh_c0[l], oth = 1-cur;

    if(double(nh_ci[l][0])*f[oth] < 0.95*double(nh_ci[l][1])*f[cur])
    nh_c0[l] = oth;
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
    else    { i0=EXT; i1=EXT+P[g]->nx-1; j0=EY; j1=EY+P[g]->ny-1; }

    for(int ii=i0; ii<=i1; ++ii)
    for(int jj=j0; jj<=j1; ++jj)
    {
        if(q->flagslice4[lij(q,ii,jj)]<=0 || pn->is_active(q,pb,ii,jj)==0)
        continue;

        double dudx,dvdy,dddx,dddy;

        // curvilinear grid: chain-rule gradients as in sflow_pjm_lin/quad
        if(q->geo_curv==1)
        {
        double a1;
        sflow_curv::grad(q,pb->U,ii,jj,dudx,a1);
        sflow_curv::grad(q,pb->V,ii,jj,a1,dvdy);
        sflow_curv::grad(q,pb->depth,ii,jj,dddx,dddy);
        }
        else
        {
        const double dx = q->DXP[ii+marge] + q->DXP[ii-1+marge];
        const double dy = q->DYP[jj+marge] + q->DYP[jj-1+marge];
        dudx = (pb->U(ii+1,jj) - pb->U(ii-1,jj))/dx;
        dvdy = (pb->V(ii,jj+1) - pb->V(ii,jj-1))/dy*p0->y_dir;
        dddx = (pb->depth(ii+1,jj) - pb->depth(ii-1,jj))/dx;
        dddy = (pb->depth(ii,jj+1) - pb->depth(ii,jj-1))/dy*p0->y_dir;
        }

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

        const size_t nc = (size_t)q->imax*q->jmax;
        psave[n].assign(pb->press.V,pb->press.V+nc);
        std::fill(pb->press.V,pb->press.V+nc,0.0);

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
    nh_coarse = false;
    for(int m=MAX(l,1); l>0 && m<=maxlev && m<8; ++m)
    if(nh_c0[m]==1)
    nh_coarse = true;
    nh_prepare(pgc,false);
    nh_core(p);
    nh_coarse = true;
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
            else    { i0=EXT; i1=EXT+P[g]->nx-1; j0=EY; j1=EY+P[g]->ny-1; }

            slice &Pc = pb->press;
            slice &WL = pb->WL;
            const double f = dt/p0->W1;

            for(int ii=i0; ii<=i1; ++ii)
            for(int jj=j0; jj<=j1; ++jj)
            {
                if(q->flagslice4[lij(q,ii,jj)]<=0 || pn->is_active(q,pb,ii,jj)==0)
                continue;

                if(q->geo_curv==1)
                {
                double gx,gy,dx_,dy_;
                sflow_curv::grad5(q,ii,jj,WL(ii+1,jj)*Pc(ii+1,jj),WL(ii-1,jj)*Pc(ii-1,jj),
                                         WL(ii,jj+1)*Pc(ii,jj+1),WL(ii,jj-1)*Pc(ii,jj-1),gx,gy);
                sflow_curv::grad(q,pb->depth,ii,jj,dx_,dy_);
                pb->UH(ii,jj) -= f*(gx - cbq*Pc(ii,jj)*dx_);
                pb->VH(ii,jj) -= f*(gy - cbq*Pc(ii,jj)*dy_);
                pb->WH(ii,jj) += f*cwq*Pc(ii,jj);
                continue;
                }

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

        std::copy(psave[n].begin(),psave[n].end(),pb->press.V);
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
