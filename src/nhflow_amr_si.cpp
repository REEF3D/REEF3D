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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"nhflow_amr.h"
#include"reef3d_mpmd.h"
#include"nhflow_amr_fill.h"
#include"nhflow_momentum_SI.h"
#include"nhflow_tracer_SI.h"
#include"nhflow_longrun.h"
#include"si_eta2D.h"
#include"reefmg2D.h"
#include"vec2D.h"
#include"matrix2D.h"
#include"reefamr_krylov.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice4.h"
#include"nhflow_curv.h"
#include"sflow_curv.h"
#include<mpi.h>
#include<map>
#include<Eigen/Sparse>
#include<Eigen/SparseLU>
#include<iomanip>

// Hydrostatic semi-implicit free surface (A 534 1) on the multi-block grid (G 50 1, DIVEMesh R 20 1) and with mesh
// refinement (G 1).
//
// Every grid (the native grid -1, the blocks 1.. of lev[0], the patches of the levels) has its own
// nhflow_momentum_SI. The step runs its parts on all grids; between them the cells around the grids take the values
// of the cells they stand for (blocks, ring, siblings) or are interpolated from the coarser level, as the stages of
// the explicit integrators (fill_level0, fill_stage):
//   si_restore                       the water column of the step before
//   fill_level0, fill_stage, si_fillF   eta, WL, wet, U, V, W, UH, VH, WH, omega of t^n around the grids
//   si_begin                         t^n, face depths, bed friction, open boundaries
//   per predictor stage:  si_stage_vel, fill U, V, W, si_stage_rhs (finest level first: the coarse faces next to a
//                                    patch take the mean of its fine advective fluxes)
//   si_predictor_end, fill UHD, VHD, si_assemble   the face coefficients of every grid; the coarse faces next to a
//                                    patch take the fine q^E and q^n; si_rhs: the eta system of every grid
//   si_solve                         one eta system over the leaf cells of all grids (reefamr_si.cpp): BiCGStab,
//                                    preconditioned by one REEFMG2D V-cycle per level-0 grid on its own rows
//                                    (restricted additive Schwarz) and, with patches, the FAC sweep over the
//                                    levels; the wetting and drying iteration over all grids together
//   si_end_faces                     the layer face discharges; the coarse faces next to a patch take the fine
//                                    Q, T, q (the volume of the coarse cell changes by the fine fluxes)
//   si_end_cells                     WL, cells, omega
//   restriction                      the covered cells take the mean of their children (WL, UH, VH, WH, flags)
// A face between two blocks is computed by both from the same values (the filled cells hold the values of the cells
// they stand for), its discharge is the same on both sides.

using namespace nhflow_amr_detail;

namespace
{
typedef reefamr_comms_off comms_off;
}

nhflow_momentum_SI* nhflow_amr::gsi(int g)
{
    return static_cast<nhflow_momentum_SI*>(gmom(g));
}

// ---- the hooks of the eta system of all grids (reefamr_si.cpp)
bool nhflow_amr::rsi_rows_of(int g, rsi_rows &R)
{
    nhflow_momentum_SI *m = gsi(g);
    R.q = glex(g);
    R.e = m->psi;
    R.ax = &m->ax;
    R.ay = &m->ay;
    R.etaI = &m->etaI;
    R.rhsW = &m->rhsW;
    R.depth = &gfd(g)->depth;
    R.bdiag = &m->bdiag;
    return true;
}

// the water level boundaries of every grid: the boundary value in the ghost cell of eta^{n+1}
void nhflow_amr::rsi_bound()
{
    for(int g : rsi_grids())
    {
        nhflow_momentum_SI *m = gsi(g);

        for(const auto &f : m->bcf)
        if(f.type==nhflow_momentum_SI::bc_level)
        m->etaI(f.gi,f.gj) = m->etaB(f.gi,f.gj);
    }
}

void nhflow_amr::rsi_run(int g, const std::function<void()> &f)
{
    if(g<0)
    {
        f();
        return;
    }

    comms_off guard(sipgc);
    pscope ps(sipgc,NP(g)->d,d0);
    f();
}

// f(g, lexer, fdm, SI object, stage objects) on every level-0 grid: the native grid with the MPI exchanges of its
// partition, the blocks 1.. on their rank alone (pscope)
template<class F>
void nhflow_amr::si_each(ghostcell *pgc, F f)
{
    f(-1,p0,d0,gsi(-1),*S0p);

    comms_off guard(pgc);
    for(int id : lev[0])
    {
        nhflow_amr_patch *c = NP(id);
        pscope ps(pgc,c->d,d0);
        f(id,c->pp,c->d,gsi(id),c->S);
    }
}

// ... on the patches of level l (on their rank alone)
template<class F>
void nhflow_amr::si_level(ghostcell *pgc, int l, F f)
{
    if(l<1 || l>=(int)lev.size())
    return;

    comms_off guard(pgc);
    for(int id : lev[l])
    {
        nhflow_amr_patch *c = NP(id);
        pscope ps(pgc,c->d,d0);
        f(id,c->pp,c->d,gsi(id),c->S);
    }
}

template<class F>
void nhflow_amr::si_all(ghostcell *pgc, F f)
{
    si_each(pgc,f);

    if(si_lv)
    for(int l=1; l<=maxlev; ++l)
    si_level(pgc,l,f);
}

// nf layer arrays A(g,m) of every grid: the cells around the blocks and the ring of the native grid from the cells
// they stand for, then the partition halo of the native grid; with patches the cells around them, coarse to fine: a
// patch of the same level, or interpolated from the coarser level layer by layer (wet cells, as fill_stage)
template<class A>
void nhflow_amr::si_fill3(ghostcell *pgc, int nf, int tag, A arr)
{
    const int K = p0->knoz;
    const int nv = nf*K;

    if(nblk>0)
    fill_run0(nv,tag,
              [&](const reefamr_fill &f, double *v)
              {
                  if(f.kind!=0)
                  {
                      for(int m=0; m<nv; ++m)
                      v[m] = 0.0;
                      return;
                  }
                  lexer *q = glex(f.g);
                  for(int m=0; m<nf; ++m)
                  {
                      const double *src = arr(f.g,m);
                      for(int k=0; k<K; ++k)
                      v[m*K+k] = src[cidx(q,f.si,f.sj,k)];
                  }
              },
              [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
              {
                  lexer *pp = q->pp;
                  for(int m=0; m<nf; ++m)
                  {
                      double *dst = arr(id,m);
                      for(int k=0; k<K; ++k)
                      dst[cidx(pp,f.di,f.dj,k)] = w[m*K+k];
                  }
              },
              [&](const reefamr_fill &f, const double *w)
              {
                  for(int m=0; m<nf; ++m)
                  {
                      double *dst = arr(-1,m);
                      for(int k=0; k<K; ++k)
                      dst[cidx(p0,f.di,f.dj,k)] = w[m*K+k];
                  }
              });

    if(p0->mpi_size>1)
    for(int m=0; m<nf; ++m)
    {
        pgc->gcparaxV(p0,arr(-1,m),14);
        pgc->gcparacoxV(p0,arr(-1,m),14);
    }

    if(!si_lv)
    return;

    for(int l=1; l<=maxlev; ++l)
    fill_run(l,nv,tag+20+l,
             [&](const reefamr_fill &f, double *v)
             {
                 if(f.kind!=0 && f.kind!=1)
                 {
                     for(int m=0; m<nv; ++m)
                     v[m] = 0.0;
                     return;
                 }
                 lexer *q = glex(f.g);
                 if(f.kind==0)
                 {
                     for(int m=0; m<nf; ++m)
                     {
                         const double *src = arr(f.g,m);
                         for(int k=0; k<K; ++k)
                         v[m*K+k] = src[cidx(q,f.si,f.sj,k)];
                     }
                     return;
                 }
                 double w1[25];
                 pweights(q,f.si,f.sj,f.ox,f.oy,w1,1);
                 for(int m=0; m<nf; ++m)
                 for(int k=0; k<K; ++k)
                 v[m*K+k] = pq3w(arr(f.g,m),q,f.si,f.sj,k,true,w1);
             },
             [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
             {
                 lexer *pp = q->pp;
                 for(int m=0; m<nf; ++m)
                 {
                     double *dst = arr(id,m);
                     for(int k=0; k<K; ++k)
                     dst[cidx(pp,f.di,f.dj,k)] = w[m*K+k];
                 }
             });
}

// omega (sigma velocity, K+1 levels) of all grids: the cells around the blocks and the ring of the native grid from the
// cells they stand for (the face columns between two blocks take the omega of both cells, col_coef), then the
// partition halo of the native grid; the cells around the patches from a sibling or interpolated (pcol)
void nhflow_amr::si_fillF(ghostcell *pgc)
{
    const int K = p0->knoz;
    const int nv = K+1;

    if(nblk>0)
    {
    fill_run0(nv,7640,
              [&](const reefamr_fill &f, double *v)
              {
                  if(f.kind!=0)
                  {
                      for(int m=0; m<nv; ++m)
                      v[m] = 0.0;
                      return;
                  }
                  lexer *q = glex(f.g);
                  const double *src = gfd(f.g)->omegaF;
                  for(int k=0; k<=K; ++k)
                  v[k] = src[fidx(q,f.si,f.sj,k)];
              },
              [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
              {
                  lexer *pp = q->pp;
                  double *dst = NP(q)->d->omegaF;
                  for(int k=0; k<=K; ++k)
                  dst[fidx(pp,f.di,f.dj,k)] = w[k];
              },
              [&](const reefamr_fill &f, const double *w)
              {
                  for(int k=0; k<=K; ++k)
                  d0->omegaF[fidx(p0,f.di,f.dj,k)] = w[k];
              });

    if(p0->mpi_size>1)
    pgc->start7S(p0,d0->omegaF,17);
    }

    if(!si_lv)
    return;

    for(int l=1; l<=maxlev; ++l)
    fill_run(l,nv,7660+l,
             [&](const reefamr_fill &f, double *v)
             {
                 if(f.kind==0)
                 {
                     lexer *q = glex(f.g);
                     const double *src = gfd(f.g)->omegaF;
                     for(int k=0; k<=K; ++k)
                     v[k] = src[fidx(q,f.si,f.sj,k)];
                     return;
                 }
                 if(f.kind==1)
                 {
                     pcol(f.g,f.si,f.sj,f.ox,f.oy,gfd(f.g)->omegaF,K,v,1);
                     return;
                 }
                 for(int m=0; m<nv; ++m)
                 v[m] = 0.0;
             },
             [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
             {
                 lexer *pp = q->pp;
                 double *dst = NP(q)->d->omegaF;
                 for(int k=0; k<=K; ++k)
                 dst[fidx(pp,f.di,f.dj,k)] = w[k];
             });
}

void nhflow_amr::si_solve(lexer *p, ghostcell *pgc)
{
    rsi_result r;
    rsi_solve(p,pgc,p->A537==1 ? 20 : 1,r);

    for(int g : rsi_grids())
    {
        nhflow_momentum_SI *m = gsi(g);
        m->psi_iter = r.iter;
        m->psi_time = r.time;
        m->wd_iter = r.wd_iter;
        m->wd_cells = r.dry;
    }
}

// the advection hooks of the grids: a coarse face next to a finer patch takes the mean of the two fine fluxes of every
// layer (multi-block: with the sign of the coarse face), a patch records the fluxes of its box faces; the records of
// a level are sent after its stage (si_step: finest first)
void nhflow_amr::si_hooks(bool on)
{
    for(int g : rsi_grids())
    {
        nhflow_momentum_SI *m = gsi(g);

        if(!on)
        {
            m->fhook = nullptr;
            continue;
        }

        m->fhook = [this,g](lexer *q, int comp, double *Fx, double *Fy)
        {
            const int K = q->knoz;

            rsm_take(g,[&](int dir, int fi, int fj, const double *v, int sgn)
            {
                double *F = (dir==1) ? Fx : Fy;
                for(int k=0; k<K; ++k)
                F[cidx(q,fi,fj,k)] = sgn*v[comp*K+k];
            });

            if(g<0 || P[g]->lev<1)
            return;

            for(int side=0; side<4; ++side)
            for(int r=0; r<rsm_count(g,side); ++r)
            {
                int dir,fi,fj;
                rsm_face(g,side,r,dir,fi,fj);
                double *s = rsm_slot(g,side,r);
                const double *F = (dir==1) ? Fx : Fy;
                for(int k=0; k<K; ++k)
                s[comp*K+k] = F[cidx(q,fi,fj,k)];
            }
        };
    }
}

// q^E and q^n of the coarse faces next to the patches: the means of the fine faces (the right-hand side of the coarse
// cell takes the fine fluxes)
void nhflow_amr::si_match_faces()
{
    rsm_begin(2);

    for(int g : rsi_grids())
    {
        if(g<0 || P[g]->lev<1)
        continue;

        nhflow_momentum_SI *m = gsi(g);
        for(int side=0; side<4; ++side)
        for(int r=0; r<rsm_count(g,side); ++r)
        {
            int dir,fi,fj;
            rsm_face(g,side,r,dir,fi,fj);
            double *s = rsm_slot(g,side,r);
            s[0] = (dir==1) ? m->qxE(fi,fj) : m->qyE(fi,fj);
            s[1] = (dir==1) ? m->qxn(fi,fj) : m->qyn(fi,fj);
        }
    }

    for(int l=maxlev; l>=1; --l)
    rsm_send(l,8700);

    for(int g : rsi_grids())
    {
        nhflow_momentum_SI *m = gsi(g);
        rsm_take(g,[&](int dir, int fi, int fj, const double *v, int sgn)
        {
            if(dir==1)
            {
            m->qxE(fi,fj) = sgn*v[0];
            m->qxn(fi,fj) = sgn*v[1];
            }
            else
            {
            m->qyE(fi,fj) = sgn*v[0];
            m->qyn(fi,fj) = sgn*v[1];
            }
        });
    }
}

// Q, T and q of the coarse faces next to the patches after the face update: the means of the fine faces (the coarse
// cell gains the volume the fine faces carry; omega and the next step's Q^n follow)
void nhflow_amr::si_match_layers()
{
    const int K = p0->knoz;
    rsm_begin(2*K+1);

    for(int g : rsi_grids())
    {
        if(g<0 || P[g]->lev<1)
        continue;

        nhflow_momentum_SI *m = gsi(g);
        lexer *q = glex(g);
        for(int side=0; side<4; ++side)
        for(int r=0; r<rsm_count(g,side); ++r)
        {
            int dir,fi,fj;
            rsm_face(g,side,r,dir,fi,fj);
            double *s = rsm_slot(g,side,r);
            const double *Q = (dir==1) ? m->Qx : m->Qy;
            const double *T = (dir==1) ? m->Tx : m->Ty;
            for(int k=0; k<K; ++k)
            {
                s[k] = Q[cidx(q,fi,fj,k)];
                s[K+k] = T[cidx(q,fi,fj,k)];
            }
            s[2*K] = (dir==1) ? m->qx(fi,fj) : m->qy(fi,fj);
        }
    }

    for(int l=maxlev; l>=1; --l)
    rsm_send(l,8720);

    for(int g : rsi_grids())
    {
        nhflow_momentum_SI *m = gsi(g);
        lexer *q = glex(g);
        rsm_take(g,[&](int dir, int fi, int fj, const double *v, int sgn)
        {
            double *Q = (dir==1) ? m->Qx : m->Qy;
            double *T = (dir==1) ? m->Tx : m->Ty;
            for(int k=0; k<K; ++k)
            {
                Q[cidx(q,fi,fj,k)] = sgn*v[k];
                T[cidx(q,fi,fj,k)] = sgn*v[K+k];
            }
            if(dir==1)
            m->qx(fi,fj) = sgn*v[2*K];
            else
            m->qy(fi,fj) = sgn*v[2*K];
        });
    }
}

void nhflow_amr::si_step(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_momentum_func *mom, nhflow_stage_obj &S0)
{
    S0p = &S0;
    sipgc = pgc;
    rsipgc = pgc;

    // no blocks and no patches (the adaptive refinement before its first patches): the step of the native grid alone
    if(!active())
    {
        gsi(-1)->step_single(p,d,pgc,S0);
        if(maxlev>=1)
        regrid_step(p,pgc);
        return;
    }

    si_lv = (maxlev>=1 && patches_total>0);

    // the time of the blocks and the patches
    for(auto q : P)
    {
        lexer *pp = q->pp;
        pp->dt = p->dt;
        pp->dt_old = p->dt_old;
        pp->simtime = p->simtime;
        pp->count = p->count;
        pp->amr_follow(*p);
    }

    // A 560: the turbulence of the blocks and the patches (level 0 has made its step before the momentum)
    if(si_lv)
    turb_step(p,pgc);

    // the water column of the step before, then the state of t^n around the grids
    si_all(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        m->si_restore(q,dd,pgc);
    });

    if(nblk>0)
    fill_level0(pgc,-1);
    if(si_lv)
    for(int l=1; l<=maxlev; ++l)
    fill_stage(pgc,l,-1);
    si_fillF(pgc);

    si_all(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        m->si_begin(q,dd,pgc,S);
    });

    // the advection hooks (face matching with the patches)
    if(si_lv)
    si_hooks(true);

    // predictor: the velocities of a stage around the grids before its fluxes
    const int ns = gsi(-1)->si_stages();

    for(int s=0; s<ns; ++s)
    {
        si_all(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
        {
            m->si_stage_vel(q,dd,pgc,S,s);
        });

        si_fill3(pgc,3,7601+s,[&](int g, int c) -> double*
        {
            fdm_nhf *dd = gfd(g);
            return c==0 ? dd->U : (c==1 ? dd->V : dd->W);
        });

        // the explicit terms, finest level first: the coarse faces next to the patches take the fine fluxes
        if(si_lv)
        {
            rsm_begin(3*p0->knoz);
            for(int l=maxlev; l>=1; --l)
            {
                si_level(pgc,l,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
                {
                    m->si_stage_rhs(q,dd,pgc,S,s);
                });
                rsm_send(l,8740);
            }
        }

        si_each(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
        {
            m->si_stage_rhs(q,dd,pgc,S,s);
        });
    }

    if(si_lv)
    si_hooks(false);

    // the cell predictor around the grids, then the face coefficients and the eta system of every grid
    si_all(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        m->si_predictor_end(q,dd,pgc);
    });

    si_fill3(pgc,2,7610,[&](int g, int c) -> double*
    {
        return c==0 ? gsi(g)->UHD : gsi(g)->VHD;
    });

    si_all(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        m->si_assemble(q,dd,pgc);
    });

    if(si_lv)
    si_match_faces();

    si_all(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        m->si_rhs(q,dd,pgc);
    });

    // one eta system over all grids
    si_solve(p,pgc);

    // faces, water level, cells
    si_all(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        m->si_end_faces(q,dd,pgc);
    });

    if(si_lv)
    si_match_layers();

    si_all(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        m->si_end_cells(q,dd,pgc,S);
    });

    // the covered cells: the mean of their children
    if(si_lv)
    {
        {
        comms_off guard(pgc);
        restrict_surface(0);
        restrict_momentum(0,false);
        }
        restrict_flags(pgc,0);
        halo0(pgc,0);

        // the restricted water level is the start of the next step of the coarse cell
        for(int g : rsi_grids())
        {
            lexer *q = glex(g);
            nhflow_momentum_SI *m = gsi(g);
            fdm_nhf *dd = gfd(g);
            for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
            for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
            m->WLk(ii,jj) = dd->WL(ii,jj);
        }
    }

    // the state of t^{n+1} around the blocks and in the ring (output, time step)
    if(nblk>0)
    fill_level0(pgc,-1);
    si_fillF(pgc);

    // tracers (A 600) with the transport of the step
    if(p->A600>0)
    si_tracer(p,pgc);

    // the log of all grids: volume of the leaf cells, discharges over the open boundaries
    double loc[4] = {0.0,0.0,0.0,0.0}, all[4];

    for(int g : rsi_grids())
    {
        lexer *q = glex(g);
        fdm_nhf *dd = gfd(g);
        nhflow_momentum_SI *m = gsi(g);

        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        if(rsi_leaf(g,ii,jj))
        loc[0] += dd->WL(ii,jj)*m->carea(q,ii,jj);

        // the native grid has summed its discharges over the ranks already
        if(g>=0)
        {
        loc[1] += m->Qin;
        loc[2] += m->Qout;
        loc[3] += m->vrelax;
        }
    }

    MPI_Allreduce(loc,all,4,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    nhflow_momentum_SI *m0 = gsi(-1);
    m0->write_log(p,all[0],all[1]+m0->Qin,all[2]+m0->Qout,all[3]+m0->vrelax,m0->psi_iter,m0->wd_iter,m0->wd_cells);

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
    cout<<"SI eta ("<<(si_lv ? "all grids" : "multi-block")<<"): iter "<<m0->psi_iter<<"  time "<<setprecision(3)<<m0->psi_time;

    if(p->A537==1)
    cout<<"  wetdry iterations "<<m0->wd_iter<<"  dry cells "<<m0->wd_cells;

    cout<<endl;
    }

    // the refinement follows the flags; the cells around the patches at the end of the step (the output)
    if(maxlev>=1)
    {
        regrid_step(p,pgc);
        if(patches_total>0)
        fill_out(pgc);
        if(nblk>0)
        fill_level0(pgc,-1);
    }
}


// tracers (A 600) of all blocks: the substeps of every block in lockstep (the same number on all grids
// and ranks), the cells around the blocks and the ring filled before each; the log over the own cells,
// the gauges P 51 from the nearest own cell centre of all grids
void nhflow_amr::si_tracer(lexer *p, ghostcell *pgc)
{
    const int nt = gsi(-1)->ptrc->nt;

    auto arr = [&](int g, int c) -> double* { return nhflow_tracer_SI::field(gfd(g),c); };

    double r = 0.0;

    si_each(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        r = max(r,m->ptrc->prepare(q,dd));
    });

    double rg;
    MPI_Allreduce(&r,&rg,1,MPI_DOUBLE,MPI_MAX,pgc->mpi_comm);
    const int ns = nhflow_tracer_SI::substeps(p,rg);

    // the fluid mask around the blocks and in the halo (static, once): filled as a tracer
    if(!sitr_mask)
    {
        si_fill3(pgc,1,7674,[&](int g, int c) -> double* { return gsi(g)->ptrc->mask; });
        sitr_mask = true;
    }

    for(int s=0; s<ns; ++s)
    {
        si_fill3(pgc,nt,7670,arr);

        si_each(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
        {
            m->ptrc->flushed(q,dd,s,ns);
        });

        si_fill3(pgc,nt,7676,[&](int g, int c) -> double* { return gsi(g)->ptrc->Cx[c]; });

        si_each(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
        {
            m->ptrc->substep(q,dd,s,ns);
        });
    }

    si_fill3(pgc,nt,7672,arr);

    si_each(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        m->ptrc->ghosts(q,dd);
    });

    // the log: own cells of all grids
    nhflow_tracer_SI::stats st;
    nhflow_tracer_SI::stats_zero(st);

    for(int g : rsi_grids())
    gsi(g)->ptrc->sums(glex(g),gfd(g),[&](int ii, int jj){ return rsi_leaf(g,ii,jj); },st);

    nhflow_tracer_SI::reduce(p,st);
    gsi(-1)->ptrc->write_log(p,st,ns,rg);

    // gauges: the nearest own cell centre of the grids of this rank, then over the ranks
    if(p->P51>0)
    {
        const int ng = p->P51;

        if(sitr_dist.empty())
        {
            sitr_dist.assign(ng,1.0e30);
            sitr_cell.assign(ng,-1);
            sitr_grid.assign(ng,-2);

            for(int g : rsi_grids())
            {
                vector<double> dist;
                vector<int> cell;
                gsi(g)->ptrc->gauge_find(glex(g),[&](int ii, int jj){ return rsi_leaf(g,ii,jj); },dist,cell);

                for(int k=0; k<ng; ++k)
                if(cell[k]>=0 && dist[k]<sitr_dist[k])
                {
                sitr_dist[k] = dist[k];
                sitr_cell[k] = cell[k];
                sitr_grid[k] = g;
                }
            }
        }

        const size_t nv = 3*nt;
        vector<double> val((size_t)ng*nv,0.0), dist = sitr_dist;

        for(int g : rsi_grids())
        {
            vector<int> cell(ng,-1);
            for(int k=0; k<ng; ++k)
            if(sitr_grid[k]==g)
            cell[k] = sitr_cell[k];

            vector<double> v;
            gsi(g)->ptrc->gauge_values(glex(g),gfd(g),cell,v);

            for(int k=0; k<ng; ++k)
            if(sitr_grid[k]==g)
            for(size_t q=0; q<nv; ++q)
            val[k*nv+q] = v[k*nv+q];
        }

        nhflow_tracer_SI::gauge_combine(p,dist,val);
        nhflow_tracer_SI::gauge_write(p,val);
    }
}

// the time step of the semi-implicit step over the block cells (nhflow_timestep_SI::start): the
// advective limit per cell and the explicit horizontal diffusion of the predictor
void nhflow_amr::dt_si_cmin(double cfl, double &cmin, double &hmax)
{
    for(auto q : P)
    {
        nhflow_amr_patch *c = NP(q);
        lexer *pp = c->pp;
        fdm_nhf *dd = c->d;
        double *EVh = (pp->A561>0) ? dd->EVH : dd->EV;

        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EY; jj<EY+c->ny; ++jj)
        {
            const int n2 = lij(pp,ii,jj);

            if(pp->flagslice4[n2]<0 || pp->wet[n2]!=1 || dd->WL(ii,jj)<=pp->A544)
            continue;

            hmax = MAX(hmax,dd->WL(ii,jj));

            for(int kk=0; kk<pp->knoz; ++kk)
            {
                const int n = cidx(pp,ii,jj,kk);
                double sigma;

                if(pp->geo_curv==1)
                {
                const double sp = sqrt(dd->U[n]*dd->U[n] + dd->V[n]*dd->V[n]);
                sigma = sp/nhflow_curv::width1(pp,ii,jj) + sp/nhflow_curv::width2(pp,ii,jj);
                }

                else
                {
                sigma = fabs(dd->U[n])/pp->DXN[ii+marge];

                if(pp->j_dir==1)
                sigma += fabs(dd->V[n])/pp->DYN[jj+marge];
                }

                if(sigma>1.0e-20)
                cmin = MIN(cmin,1.0/sigma);

                if(pp->A512>0 || pp->A560>0 || pp->A547_c>0.0 || pp->A547_nu0>0.0)
                {
                const double visc = dd->VISC[n] + EVh[n];
                const double lam = pp->geo_curv==1 ? sflow_curv::inv2(pp,ii,jj)
                                 : 1.0/(pp->DXN[ii+marge]*pp->DXN[ii+marge]) + pp->y_dir/(pp->DYN[jj+marge]*pp->DYN[jj+marge]);

                cmin = MIN(cmin, 0.25/(visc*lam + 1.0e-20)/(2.0*cfl));
                }
            }
        }
    }
}

// long-run statistics (P 160): the blocks 1.. (their own cells)
void nhflow_amr::longrun_grids(nhflow_longrun *lr)
{
    for(int id : lev[0])
    {
        nhflow_amr_patch *c = NP(id);
        reefamr_patch *q = P[id];
        lr->add_grid(c->pp,c->d,EXT,EXT+q->nx,EY,EY+q->ny);
    }
}
