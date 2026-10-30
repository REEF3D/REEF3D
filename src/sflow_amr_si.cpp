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
#include"sflow_momentum_SI.h"
#include"sflow_momentum_RK3.h"
#include"sflow_curv.h"
#include"sflow_fsf.h"
#include"si_eta2D.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"slice4.h"
#include<cmath>
#include<mpi.h>
#include<iomanip>

//  Semi-implicit free surface (A 226 1) on a multi-block grid (G 50 1) and with mesh refinement (G 1).
//
//  Every grid (the native grid -1, the blocks 1.. of lev[0], the patches of the levels) has its own
//  sflow_momentum_SI; the level-0 object hands its step over (sflow_momentum_SI::start -> si_step).  The parts of
//  the step run on all grids; between them the cells around the grids take the values of the cells they stand for
//  (blocks, ring, siblings) or are interpolated from the coarser level (the fills of the explicit stages: well
//  balanced surface, wet-dry flags of the source cell):
//    si_restore                      the water column of the step before
//    fill_level0, fill_level         WL, UH, VH, eta, U, V, wet of t^n around the grids
//    si_begin                        t^n, face depths, bed friction, open boundaries
//    per predictor stage:  si_stage_vel, fill U, V, UH, VH of the stage, si_stage_rhs (finest level first: the
//                                    coarse faces next to a patch take the mean of its fine advective fluxes)
//    si_predictor_end, fill UH*, VH*, si_assemble    the face coefficients of every grid; the coarse faces next to
//                                    a patch take the fine q^E and q^n; si_rhs the eta system of every grid
//    rsi_solve                       one eta system over the leaf cells of all grids (reefamr_si.cpp)
//    si_end_faces                    the face discharges; the coarse faces next to a patch take the fine q (the
//                                    volume of the coarse cell changes by the fine fluxes: exact volume)
//    si_end_cells                    WL, cells, relaxation zones
//    restriction                     the covered cells take the mean of their children (WL, UH, VH, flags)
//  One global time step: the advective limit over the cells of all grids (si_timestep).  Not with subcycling (G 7 1).

namespace
{
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

inline double mmod(double a, double b)
{
    if(a*b<=0.0)
    return 0.0;
    return fabs(a)<fabs(b)?a:b;
}

typedef reefamr_comms_off comms_off;
}

void sflow_amr::si_attach(sflow_momentum_SI *m)
{
    si = 1;
    psim0 = m;
    m->prun = this;
    m->amrgrid = true;
}

sflow_momentum_SI* sflow_amr::gsi(int g)
{
    return (g<0) ? psim0 : SP(g)->psim;
}

bool sflow_amr::rsi_rows_of(int g, rsi_rows &R)
{
    sflow_momentum_SI *m = gsi(g);
    if(m==nullptr)
    return false;

    R.q = glex(g);
    R.e = m->psi;
    R.ax = &m->ax;
    R.ay = &m->ay;
    R.etaI = &m->etaI;
    R.rhsW = &m->rhsW;
    R.depth = (g<0) ? &b0->depth : &SP(g)->b->depth;
    R.bdiag = &m->bdiag;
    return true;
}

// the water level boundaries of every grid: the boundary value in the ghost cell of eta^{n+1}
void sflow_amr::rsi_bound()
{
    for(int g : rsi_grids())
    {
        sflow_momentum_SI *m = gsi(g);

        for(const auto &f : m->bcf)
        if(f.type==sflow_momentum_SI::bc_level)
        m->etaI(f.gi,f.gj) = m->etaB(f.gi,f.gj);
    }
}

// f(g, lexer, fdm, SI object) on the level-0 grids: the native grid with the MPI exchanges of its partition, the
// blocks 1.. on their rank alone
template<class F>
void sflow_amr::si_level0(ghostcell *pgc, F f)
{
    f(-1,p0,b0,psim0);

    comms_off guard(pgc);
    for(int id : lev[0])
    f(id,SP(id)->pp,SP(id)->b,SP(id)->psim);
}

template<class F>
void sflow_amr::si_level(ghostcell *pgc, int l, F f)
{
    if(l<1 || l>=(int)lev.size())
    return;

    comms_off guard(pgc);
    for(int id : lev[l])
    f(id,SP(id)->pp,SP(id)->b,SP(id)->psim);
}

template<class F>
void sflow_amr::si_all(ghostcell *pgc, F f)
{
    si_level0(pgc,f);

    if(si_lv)
    for(int l=1; l<=maxlev; ++l)
    si_level(pgc,l,f);
}

// limited linear interpolation of a cell field of the coarse cell (ic,jc) into its child (ox,oy): minmod slopes over
// the wet fluid neighbours (velocities of a dry cell are 0)
double sflow_amr::si_interp(lexer *q, fdm2D *b, slice &x, int ic, int jc, int ox, int oy)
{
    const double c = x(ic,jc);
    auto ok = [&](int a, int d) { return q->flagslice4[lij(q,a,d)]>0 && q->wet[lij(q,a,d)]==1; };

    if(!ok(ic,jc))
    return c;

    double sx=0.0, sy=0.0;
    if(ok(ic+1,jc) && ok(ic-1,jc))
    sx = mmod(x(ic+1,jc)-c,c-x(ic-1,jc));
    if(oy!=0 && ok(ic,jc+1) && ok(ic,jc-1))
    sy = mmod(x(ic,jc+1)-c,c-x(ic,jc-1));

    return c + 0.25*ox*sx + 0.25*oy*sy;
}

// nf fields F(g,k) of every grid: the cells around the blocks and the ring of the native grid from the cells they
// stand for, the partition halo of the native grid; with patches the cells around them, coarse to fine: a patch of
// the same level, or interpolated from the coarser level
void sflow_amr::si_fillS(ghostcell *pgc, int nf, int tag, const std::function<slice&(int,int)> &F)
{
    if(nblk>0)
    {
        fill_run0(nf,tag,
                  [&](const reefamr_fill &f, double *v)
                  {
                      for(int k=0; k<nf; ++k)
                      v[k] = (f.kind==0) ? F(f.g,k)(f.si,f.sj) : 0.0;
                  },
                  [&](reefamr_patch *c, int id, const reefamr_fill &f, const double *w)
                  {
                      for(int k=0; k<nf; ++k)
                      F(id,k)(f.di,f.dj) = w[k];
                  },
                  [&](const reefamr_fill &f, const double *w)
                  {
                      for(int k=0; k<nf; ++k)
                      F(-1,k)(f.di,f.dj) = w[k];
                  });

        if(p0->mpi_size>1)
        for(int k=0; k<nf; ++k)
        {
            pgc->gcslparax(p0,F(-1,k),4);
            pgc->gcslparacox(p0,F(-1,k),10);
        }
    }

    if(!si_lv)
    return;

    for(int l=1; l<=maxlev; ++l)
    fill_run(l,nf,tag+20+l,
             [&](const reefamr_fill &f, double *v)
             {
                 if(f.kind==0)
                 {
                     for(int k=0; k<nf; ++k)
                     v[k] = F(f.g,k)(f.si,f.sj);
                     return;
                 }
                 if(f.kind==1)
                 {
                     gh G = grid(f.g);
                     for(int k=0; k<nf; ++k)
                     v[k] = si_interp(G.q,G.b,F(f.g,k),f.si,f.sj,f.ox,f.oy);
                     return;
                 }
                 for(int k=0; k<nf; ++k)
                 v[k] = 0.0;
             },
             [&](reefamr_patch *c, int id, const reefamr_fill &f, const double *w)
             {
                 for(int k=0; k<nf; ++k)
                 F(id,k)(f.di,f.dj) = w[k];
             });
}

// the advection hooks: a coarse face next to a finer patch takes the mean of the two fine fluxes of UH and VH
// (multi-block: with the sign of the coarse face), a patch records the fluxes of its box faces
void sflow_amr::si_hooks(bool on)
{
    for(int g : rsi_grids())
    {
        sflow_momentum_SI *m = gsi(g);

        if(!on)
        {
            m->fhook = nullptr;
            continue;
        }

        m->fhook = [this,g](lexer *q, slice &fxU, slice &fxV, slice &fyU, slice &fyV)
        {
            rsm_take(g,[&](int dir, int fi, int fj, const double *v, int sgn)
            {
                if(dir==1)
                {
                fxU(fi,fj) = sgn*v[0];
                fxV(fi,fj) = sgn*v[1];
                }
                else
                {
                fyU(fi,fj) = sgn*v[0];
                fyV(fi,fj) = sgn*v[1];
                }
            });

            if(g<0 || P[g]->lev<1)
            return;

            for(int side=0; side<4; ++side)
            for(int r=0; r<rsm_count(g,side); ++r)
            {
                int dir,fi,fj;
                rsm_face(g,side,r,dir,fi,fj);
                double *s = rsm_slot(g,side,r);
                s[0] = (dir==1) ? fxU(fi,fj) : fyU(fi,fj);
                s[1] = (dir==1) ? fxV(fi,fj) : fyV(fi,fj);
            }
        };
    }
}

// the face values of the coarse faces next to the patches: the means of the fine faces
//  what 0: q^E and q^n before the right-hand sides of the eta systems
//  what 1: q after the face update (the volume of the coarse cell changes by the fine fluxes)
void sflow_amr::si_match(int what)
{
    const int nv = (what==0) ? 2 : 1;
    rsm_begin(nv);

    for(int g : rsi_grids())
    {
        if(g<0 || P[g]->lev<1)
        continue;

        sflow_momentum_SI *m = gsi(g);
        for(int side=0; side<4; ++side)
        for(int r=0; r<rsm_count(g,side); ++r)
        {
            int dir,fi,fj;
            rsm_face(g,side,r,dir,fi,fj);
            double *s = rsm_slot(g,side,r);
            if(what==0)
            {
            s[0] = (dir==1) ? m->qxE(fi,fj) : m->qyE(fi,fj);
            s[1] = (dir==1) ? m->qxn(fi,fj) : m->qyn(fi,fj);
            }
            else
            s[0] = (dir==1) ? m->qx(fi,fj) : m->qy(fi,fj);
        }
    }

    for(int l=maxlev; l>=1; --l)
    rsm_send(l,(what==0) ? 8760 : 8780);

    for(int g : rsi_grids())
    {
        sflow_momentum_SI *m = gsi(g);
        rsm_take(g,[&](int dir, int fi, int fj, const double *v, int sgn)
        {
            if(what==0)
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
                return;
            }
            if(dir==1)
            m->qx(fi,fj) = sgn*v[0];
            else
            m->qy(fi,fj) = sgn*v[0];
        });
    }
}

// --------------------------------------------------------------------- the step
void sflow_amr::si_step(lexer *p, fdm2D *b, ghostcell *pgc)
{
    sipgc = pgc;
    rsipgc = pgc;

    // no blocks and no patches (the adaptive refinement before its first patches): the step of the native grid alone
    if(!(nblk>0 || (maxlev>=1 && patches_total>0)))
    {
        psim0->step_single(p,b,pgc);

        if(maxlev>=1 && regrid_int>0 && p->count%regrid_int==0)
        regrid(p,pgc,false);
        load_step(p);
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

    // the water column of the step before, then the state of t^n around the grids
    si_all(pgc,[&](int g, lexer *q, fdm2D *bb, sflow_momentum_SI *m) { m->si_restore(q,bb,pgc); });

    cache_stage(0);
    if(nblk>0)
    fill_level0(pgc,0);
    if(si_lv)
    for(int l=1; l<=maxlev; ++l)
    fill_level(pgc,l,0);

    si_all(pgc,[&](int g, lexer *q, fdm2D *bb, sflow_momentum_SI *m) { m->si_begin(q,bb,pgc); });

    // predictor: the velocities and the input of a stage around the grids before its fluxes
    if(si_lv)
    si_hooks(true);

    const int ns = psim0->si_stages();

    for(int s=0; s<ns; ++s)
    {
        si_all(pgc,[&](int g, lexer *q, fdm2D *bb, sflow_momentum_SI *m) { m->si_stage_vel(q,bb,pgc,s); });

        si_fillS(pgc,4,8800+10*s,[&](int g, int k) -> slice&
        {
            sflow_momentum_SI *m = gsi(g);
            fdm2D *bb = (g<0) ? b0 : SP(g)->b;
            slice *UHs,*VHs,*UHo,*VHo;
            m->stage_io_si(s,UHs,VHs,UHo,VHo);
            return k==0 ? bb->U : (k==1 ? bb->V : (k==2 ? *UHs : *VHs));
        });

        // the explicit terms, finest level first: the coarse faces next to the patches take the fine fluxes
        if(si_lv)
        {
            rsm_begin(2);
            for(int l=maxlev; l>=1; --l)
            {
                si_level(pgc,l,[&](int g, lexer *q, fdm2D *bb, sflow_momentum_SI *m) { m->si_stage_rhs(q,bb,pgc,s); });
                rsm_send(l,8840);
            }
        }

        si_level0(pgc,[&](int g, lexer *q, fdm2D *bb, sflow_momentum_SI *m) { m->si_stage_rhs(q,bb,pgc,s); });
    }

    if(si_lv)
    si_hooks(false);

    // the cell predictor around the grids, then the face coefficients and the eta system of every grid
    si_all(pgc,[&](int g, lexer *q, fdm2D *bb, sflow_momentum_SI *m) { m->si_predictor_end(q,bb,pgc); });

    si_fillS(pgc,2,8860,[&](int g, int k) -> slice& { return k==0 ? gsi(g)->UHp : gsi(g)->VHp; });

    si_all(pgc,[&](int g, lexer *q, fdm2D *bb, sflow_momentum_SI *m) { m->si_assemble(q,bb,pgc); });

    if(si_lv)
    si_match(0);

    si_all(pgc,[&](int g, lexer *q, fdm2D *bb, sflow_momentum_SI *m) { m->si_rhs(q,bb,pgc); });

    // one eta system over the leaf cells of all grids
    rsi_result r;
    rsi_solve(p,pgc,p->A229==1 ? 20 : 1,r);

    for(int g : rsi_grids())
    {
        sflow_momentum_SI *m = gsi(g);
        m->psi->iter = r.iter;
        m->psi->time = r.time;
        m->wd_iter = r.wd_iter;
        m->wd_cells = r.dry;
    }

    // faces, water level, cells
    si_all(pgc,[&](int g, lexer *q, fdm2D *bb, sflow_momentum_SI *m) { m->si_end_faces(q,bb,pgc); });

    if(si_lv)
    si_match(1);

    si_all(pgc,[&](int g, lexer *q, fdm2D *bb, sflow_momentum_SI *m) { m->si_end_cells(q,bb,pgc); });

    // the covered cells: the mean of their children; the end of the step of the patches (as level 0 in the main loop)
    if(si_lv)
    {
        restrict_levels(p,2);
        exchange_level0(p,b,pgc,2);

        for(int g : rsi_grids())
        {
            lexer *q = glex(g);
            fdm2D *bb = (g<0) ? b0 : SP(g)->b;
            sflow_momentum_SI *m = gsi(g);
            for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
            for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
            m->WLk(ii,jj) = bb->WL(ii,jj);
        }

        comms_off guard(pgc);
        for(auto q : P)
        {
            if(q->lev<1)
            continue;
            sflow_amr_patch *c = SP(q);
            lexer *pp = c->pp;
            if(p->A248==0)
            for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
            for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
            c->b->breaking(ii,jj)=0;
            c->pfsf->depth_update(pp,c->b,pgc,c->b->WL);
        }
    }

    // the blocks 1..: the end of their step (as level 0 in the main loop)
    if(nblk>0)
    {
        comms_off guard(pgc);
        for(int id : lev[0])
        {
            sflow_amr_patch *c = SP(id);
            lexer *pp = c->pp;
            if(p->A248==0)
            for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
            for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
            c->b->breaking(ii,jj)=0;
            c->pfsf->depth_update(pp,c->b,pgc,c->b->WL);
        }
    }

    si_log(p,pgc);

    // the refinement follows the flags
    if(maxlev>=1 && ((regrid_int>0 && p->count%regrid_int==0) || replace_due(p)))
    {
    double t0 = MPI_Wtime();
    regrid(p,pgc,false);
    tm[4] += MPI_Wtime()-t0;
    }

    // multi-block: the ring of level 0 again from the blocks, for the time step, the output and the next step
    if(nblk>0)
    {
    cache_stage(0);
    fill_level0(pgc,0);
    }

    load_step(p);
}

// the log of all grids: the volume of the leaf cells, the discharges over the open boundaries, the eta solve. The
// volume from the water column of the step (WLk): sflow_eta::depth_update of the blocks and the patches above has
// given their dry cells the film of the wetting criterion again, which the next step takes back (si_restore)
void sflow_amr::si_log(lexer *p, ghostcell *pgc)
{
    double loc[4] = {0.0,0.0,0.0,0.0}, all[4];

    for(int g : rsi_grids())
    {
        lexer *q = glex(g);
        sflow_momentum_SI *m = gsi(g);

        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        if(rsi_leaf(g,ii,jj))
        loc[0] += m->WLk(ii,jj)*m->carea(q,ii,jj);

        // the native grid has summed its discharges over the ranks already
        if(g>=0)
        {
        loc[1] += m->Qin;
        loc[2] += m->Qout;
        loc[3] += m->vrelax;
        }
    }

    MPI_Allreduce(loc,all,4,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    sflow_momentum_SI *m0 = psim0;

    if(p->mpirank==0)
    m0->silog<<p->count<<" "<<std::setprecision(10)<<p->simtime+p->dt<<" "<<p->dt<<" "<<std::setprecision(16)<<all[0]<<" "
             <<all[1]+m0->Qin<<" "<<all[2]+m0->Qout<<" "<<all[3]+m0->vrelax<<" "<<m0->psi->iter<<" "<<m0->wd_iter<<" "<<m0->wd_cells<<std::endl;

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
    cout<<"SI eta ("<<(si_lv ? "all grids" : "multi-block")<<"): iter "<<m0->psi->iter<<"  time "<<std::setprecision(3)<<m0->psi->time;

    if(p->A229==1)
    cout<<"  wetdry iterations "<<m0->wd_iter<<"  dry cells "<<m0->wd_cells;

    cout<<endl;
    }
}

// the time step: the advective limit (and the explicit diffusion of A 212 1) over the cells of the blocks and the
// patches, as sflow_timestep_SI over level 0
void sflow_amr::si_timestep(lexer *p, ghostcell *pgc)
{
    const double cfl = p->N47;
    double cmin = 1.0e20;

    for(auto q : P)
    {
        sflow_amr_patch *c = SP(q);
        lexer *pp = c->pp;
        fdm2D *pb = c->b;
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EY; jj<EY+c->ny; ++jj)
        if(pp->wet[lij(pp,ii,jj)]==1 && pp->flagslice4[lij(pp,ii,jj)]>0)
        {
            double sigma;

            if(pp->geo_curv==1)
            sigma = sflow_curv::sigma(pp,pb,ii,jj,0.0);

            else
            {
            sigma = fabs(pb->U(ii,jj))/pp->DXN[ii+marge];

            if(p->j_dir==1)
            sigma += fabs(pb->V(ii,jj))/pp->DYN[jj+marge];
            }

            if(sigma>1.0e-20)
            cmin = MIN(cmin,1.0/sigma);

            if(p->A212==1)
            {
            double visc = p->W2 + pb->eddyv(ii,jj);

            if(p->A246==2 && pb->breaking(ii,jj)==1)
            visc += p->A250;

            const double lam = (pp->geo_curv==1) ? sflow_curv::inv2(pp,ii,jj)
                             : 1.0/(pp->DXN[ii+marge]*pp->DXN[ii+marge]) + p->y_dir/(pp->DYN[jj+marge]*pp->DYN[jj+marge]);

            cmin = MIN(cmin, 0.5/(visc*lam + 1.0e-20)/(2.0*cfl));
            }
        }
    }

    const double dtp = pgc->globalmin(cfl*2.0*cmin);

    if(p->N48==1)
    p->dt = MIN(p->dt,dtp);

    p->dt = pgc->timesync(p->dt);
}
