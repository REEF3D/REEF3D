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

#include"sflow_amr.h"
#include"sediment_f.h"
#include"sediment_fdm.h"
#include"sediment_exner.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"reefamr_sed.h"
#include"slice4.h"
#include"sflow_momentum_RK3.h"
#include<memory>
#include<cstdlib>
#include<cstring>
#include<sstream>
#include<fstream>

// Sediment transport (S 10) on the grids of the SFLOW mesh refinement, as nhflow_amr_sed.cpp: every grid has its
// own sediment_f (amrgrid 1 on the patches); a sediment step of level 0 runs all grids with the Exner face
// fluxes of the coarse faces next to finer patches taken from the patches (the bed volume is conserved over
// the leaf cells), one dtsed, the bed of the covered cells restricted (mean per area of the children). The
// filter (S 100) and the sand slide (S 90) run over all grids in steps with the exchanges between them
// (reefamr_sed.h): the interfaces of the grids of a level are open. The bed of the cells around the patches is
// the bed of the patch of the same level there, else linear from the coarser level; the still water depth of
// the filled cells (reefamr_fill::aux) follows it.

namespace
{
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

inline double mmod(double a, double b)
{
    if(a*b<=0.0)
    return 0.0;
    return fabs(a)<fabs(b)?a:b;
}

// weights of the 2x2 children (i0+a, j0+b), order (0,0) (1,0) (0,1) (1,1): their areas over the sum
inline void cw4(const lexer *q, int i0, int j0, double *w)
{
    if(q->geo_curv!=1)
    {
    w[0] = w[1] = w[2] = w[3] = 0.25;
    return;
    }
    const double a0 = q->AREA2D[lij(q,i0,j0)], a1 = q->AREA2D[lij(q,i0+1,j0)];
    const double a2 = q->AREA2D[lij(q,i0,j0+1)], a3 = q->AREA2D[lij(q,i0+1,j0+1)];
    const double s = a0 + a1 + a2 + a3;
    w[0] = a0/s; w[1] = a1/s; w[2] = a2/s; w[3] = a3/s;
}

// patch kernels: no exchange, reductions local to the patch, the domain edges those of the patch cells
// (the sediment code tests them with origin_i/j and gknox/y; a patch has no global position)
struct sedscope
{
    ghostcell *g;
    lexer *q;
    bool oc, ol;
    int oi, oj, gx, gy;
    sedscope(ghostcell *gg, lexer *pp) : g(gg), q(pp), oi(pp->origin_i), oj(pp->origin_j), gx(pp->gknox), gy(pp->gknoy)
    {
        oc = g->set_comms(false);
        ol = g->set_local(true);
        q->origin_i = q->origin_j = 0;
        q->gknox = q->gknoy = 1<<29;
    }
    ~sedscope()
    {
        q->origin_i = oi; q->origin_j = oj;
        q->gknox = gx; q->gknoy = gy;
        g->set_local(ol);
        g->set_comms(oc);
    }
};
}

void sflow_amr::sed_attach(sediment *ps)
{
    if(p0->S10==0)
    return;

    sed0 = dynamic_cast<sediment_f*>(ps);

    if(sed0!=nullptr && (maxlev>0 || nblk>0))
    sed0->grids = this;
}

// linear, mean preserving: the child in quadrant (ox,oy) of the coarse cell (ic,jc) of lexer q
double sflow_amr::sed_plin(slice &f, lexer *q, int ic, int jc, int ox, int oy)
{
    auto fl = [&](int a, int b) { return q->flagslice4[lij(q,a,b)]>0; };
    const double c = f(ic,jc);
    double sx=0.0, sy=0.0;
    if(fl(ic+1,jc) && fl(ic-1,jc))
    sx = mmod(f(ic+1,jc)-c,c-f(ic-1,jc));
    if(oy!=0 && fl(ic,jc+1) && fl(ic,jc-1))
    sy = mmod(f(ic,jc+1)-c,c-f(ic,jc-1));
    return c + 0.25*ox*sx + 0.25*oy*sy;
}

void sflow_amr::sed_patch_ini(ghostcell *pgc, sflow_amr_patch &c)
{
    if(c.sed==nullptr)
    return;

    lexer *pp = c.pp;
    fdm2D *b = c.b;
    sedscope ss(pgc,pp);

    // no solids on the patches and blocks (the bed is the erodible bed, sediment_f::update_sflow)
    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
    b->topobed(ii,jj) = b->bed(ii,jj);
    b->solidbed(ii,jj) = -1.0e10;
    }

    c.sed->ini_sflow(pp,b,pgc);

    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
    c.sed->s->bedzh(ii,jj) = b->topobed(ii,jj);
    c.sed->s->bedzh0(ii,jj) = b->topobed(ii,jj);
    }
}

// the bed of the cells around the level-l patches, the depths and the surface of the patches, the still water
// depth of the filled cells (local and served to other ranks)
void sflow_amr::sed_fill_bed(ghostcell *pgc, int l)
{
    auto bedv = [&](const reefamr_fill &f) -> double
    {
        if(f.kind==0)
        return grid(f.g).b->bed(f.si,f.sj);
        if(f.kind==1)
        return sed_plin(grid(f.g).b->bed,grid(f.g).q,f.si,f.sj,f.ox,f.oy);
        return p0->wd - f.aux;
    };

    fill_run(l,1,7900+l,
             [&](const reefamr_fill &f, double *v) { v[0] = bedv(f); },
             [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
             {
                 fdm2D *b = SP(q)->b;
                 b->bed(f.di,f.dj) = w[0];
                 b->topobed(f.di,f.dj) = w[0];
                 if(SP(q)->sed!=nullptr)
                 SP(q)->sed->s->bedzh(f.di,f.dj) = w[0];
             });

    // the served cells: the still water depth for their prolongation on this rank
    for(auto &f : gserve[l])
    f.aux = p0->wd - bedv(f);

    for(int id : lev[l])
    {
        sflow_amr_patch *c = SP(id);
        lexer *pp = c->pp;
        fdm2D *b = c->b;

        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            b->depth(ii,jj) = p0->wd - b->bed(ii,jj);
            b->eta(ii,jj) = b->WL(ii,jj) - b->depth(ii,jj);
        }

        {
        sedscope ss(pgc,pp);
        pgc->gcsl_start4(pp,b->depth,50);
        }

        for(int ii=pp->imin; ii<pp->imin+pp->imax-1; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        b->dfx(ii,jj) = 0.5*(b->depth(ii+1,jj)+b->depth(ii,jj));

        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax-1; ++jj)
        b->dfy(ii,jj) = 0.5*(b->depth(ii,jj+1)+b->depth(ii,jj));

        for(auto &f : c->fill)
        f.aux = b->depth(f.di,f.dj);
    }
}

// multi-block grid: the bed of the cells around the blocks and of the ring of the native grid, the depths and
// the surface of the blocks
void sflow_amr::sed_fill_bed0(ghostcell *pgc)
{
    fill_run0(1,7899,
              [&](const reefamr_fill &f, double *v)
              {
                  v[0] = (f.kind==0) ? grid(f.g).b->bed(f.si,f.sj) : p0->wd - f.aux;
              },
              [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
              {
                  fdm2D *pb = SP(q)->b;
                  pb->bed(f.di,f.dj) = pb->topobed(f.di,f.dj) = w[0];
                  SP(q)->sed->s->bedzh(f.di,f.dj) = w[0];
              },
              [&](const reefamr_fill &f, const double *w)
              {
                  b0->bed(f.di,f.dj) = b0->topobed(f.di,f.dj) = w[0];
                  sed0->s->bedzh(f.di,f.dj) = w[0];
                  b0->depth(f.di,f.dj) = p0->wd - w[0];
              });

    for(int id : lev[0])
    {
        sflow_amr_patch *c = SP(id);
        lexer *pp = c->pp;
        fdm2D *pb = c->b;

        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            pb->depth(ii,jj) = p0->wd - pb->bed(ii,jj);
            pb->eta(ii,jj) = pb->WL(ii,jj) - pb->depth(ii,jj);
        }

        {
        sedscope ss(pgc,pp);
        pgc->gcsl_start4(pp,pb->depth,50);
        }

        for(int ii=pp->imin; ii<pp->imin+pp->imax-1; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        pb->dfx(ii,jj) = 0.5*(pb->depth(ii+1,jj)+pb->depth(ii,jj));

        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax-1; ++jj)
        pb->dfy(ii,jj) = 0.5*(pb->depth(ii,jj+1)+pb->depth(ii,jj));

        for(auto &f : c->fill)
        f.aux = pb->depth(f.di,f.dj);
    }
}

// regrid: the bed of the fresh level-l patches from the current bed of their parents (linear, mean preserving;
// curvilinear grids per area), where an old patch of the level was its bed; the cells around them
void sflow_amr::sed_bed_regrid(ghostcell *pgc, int l, vector<reefamr_patch*> &oldP)
{
    block_down_if(l,5,7920+l,[](reefamr_patch *q) { return q->fresh; },
                  [&](const reefamr_block &B, int key, double *v)
                  {
                      gh G = grid(B.g);
                      for(int a=0; a<2; ++a)
                      for(int d=0; d<=sy; ++d)
                      v[2*a+d] = sed_plin(G.b->bed,G.q,B.ic,B.jc,a==0?-1:1,(sy==0) ? 0 : (d==0?-1:1));
                      v[4] = G.b->bed(B.ic,B.jc);
                  },
                  [&](reefamr_patch *q, int id, int k, const double *v)
                  {
                      if(!q->fresh)
                      return;
                      const int nby = q->ny>>sy;
                      const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                      double w[4];
                      cw4(q->pp,i0,j0,w);
                      const double sh = (q->pp->geo_curv==1) ? v[4] - (w[0]*v[0] + w[1]*v[2] + w[2]*v[1] + w[3]*v[3]) : 0.0;
                      fdm2D *b = SP(q)->b;
                      for(int a=0; a<2; ++a)
                      for(int d=0; d<=sy; ++d)
                      b->bed(i0+a,j0+d) = v[2*a+d] + sh;
                  });

    old_run(l,1,7940+l,oldP,
            [&](reefamr_patch *q, int io, int jo, double *v)
            {
                v[0] = SP(q)->b->bed(io,jo);
            },
            [&](reefamr_patch *q, int id, int ii, int jj, const double *v)
            {
                SP(q)->b->bed(ii,jj) = v[0];
            });

    for(int id : lev[l])
    if(P[id]->fresh)
    {
        sflow_amr_patch *c = SP(id);
        lexer *pp = c->pp;
        fdm2D *b = c->b;
        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            b->topobed(ii,jj) = b->bed0(ii,jj) = b->bed(ii,jj);
            b->solidbed(ii,jj) = -1.0e10;
        }
    }

    sed_fill_bed(pgc,l);
}


void sflow_amr::sediment_step_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow, sediment_f *s0)
{
    const int ng = (int)P.size();

    // 0. the cells around the patches with the current state (two patches that meet compute the same face
    //    flux), the face velocities of the patches with them
    cache_stage(0);
    if(nblk>0)
    fill_level0(pgc,0);
    for(int l=1; l<=maxlev; ++l)
    fill_level(pgc,l,0);

    for(int id=0; id<ng; ++id)
    {
        sflow_amr_patch *c = SP(id);
        lexer *pp = c->pp;
        fdm2D *pb = c->b;
        for(int ii=pp->imin; ii<pp->imin+pp->imax-1; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        pb->P(ii,jj) = (pp->wet[lij(pp,ii,jj)]==1 && pp->wet[lij(pp,ii+1,jj)]==1) ? 0.5*(pb->U(ii,jj)+pb->U(ii+1,jj)) : 0.0;
        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax-1; ++jj)
        pb->Q(ii,jj) = (pp->wet[lij(pp,ii,jj)]==1 && pp->wet[lij(pp,ii,jj+1)]==1) ? 0.5*(pb->V(ii,jj)+pb->V(ii,jj+1)) : 0.0;
    }

    // the lags of the bend morphology (S 95 2, S 64, sediment_exner_lag.cpp): the relaxation of every grid holds the
    // cells around it (filled from their sources) and its covered cells (the mean of their children)
    const bool lag = (p->S95==2 || p->S64>0.0);
    if(lag)
    {
        for(int g=-1; g<ng; ++g)
        {
        lexer *q = grid(g).q;
        slice4 &fx = gsed(g)->exner()->lagfix;
        for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
        for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
        fx(ii,jj) = 0.0;
        if(g>=0)
        for(auto &f : SP(g)->fill)
        fx(f.di,f.dj) = 1.0;
        if(g<0 && nblk>0)
        for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
        for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
        if(nring[lij(q,ii,jj)]==1)
        fx(ii,jj) = 1.0;
        }

        for(int l=1; l<=maxlev; ++l)
        {
        for(int id : lev[l])
        {
            reefamr_patch *c = P[id];
            for(size_t k=0; k<c->rgrid.size(); ++k)
            if(c->rgrid[k]>=-1)
            gsed(c->rgrid[k])->exner()->lagfix(c->ric[k],c->rjc[k]) = 1.0;
        }
        for(auto &B : bsrv[l])
        if(B.g>=-1)
        gsed(B.g)->exner()->lagfix(B.ic,B.jc) = 1.0;
        }
    }

    // 1. physics of every grid
    s0->amr_physics_sflow(p,b,pgc);

    for(int id=0; id<ng; ++id)
    {
        sflow_amr_patch *c = SP(id);
        sedscope ss(pgc,c->pp);
        c->sed->amr_physics_sflow(c->pp,c->b,pgc);
    }

    // the covered cells: qb of their children (mean per area): the cells around a patch take qb linear from
    // the coarser level, with the slopes of its neighbours, covered ones among them (step 2); again after every
    // pass of the lags
    auto sync = [&]()
    {
        for(int l=maxlev; l>=1; --l)
        block_up(l,sediment_exner::amr_nfield,8040+l,
                 [&](reefamr_patch *q, int id, int k, double *v)
                 {
                     sediment_f *sf = SP(q)->sed;
                     const int nby = q->ny>>sy;
                     const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                     double w[4];
                     cw4(q->pp,i0,j0,w);
                     for(int n=0; n<sediment_exner::amr_nfield; ++n)
                     {
                     auto f = [&](int a, int b) { return sf->exner()->amr_field(sf->s,n,a,b); };
                     v[n] = (sy==0) ? 0.5*(f(i0,j0) + f(i0+1,j0)) : w[0]*f(i0,j0) + w[1]*f(i0+1,j0) + w[2]*f(i0,j0+1) + w[3]*f(i0+1,j0+1);
                     }
                 },
                 [&](const reefamr_block &B, int key, const double *v)
                 {
                     for(int n=0; n<sediment_exner::amr_nfield; ++n)
                     gsed(B.g)->exner()->amr_field(gsed(B.g)->s,n,B.ic,B.jc) = v[n];
                 });

        // 2. qb around the blocks and in the ring of the native grid (multi-block grid), around the patches
        if(nblk>0)
        fill_run0(sediment_exner::amr_nfield,7958,
                  [&](const reefamr_fill &f, double *v)
                  {
                      for(int n=0; n<sediment_exner::amr_nfield; ++n)
                      v[n] = (f.kind==0) ? gsed(f.g)->exner()->amr_field(gsed(f.g)->s,n,f.si,f.sj) : 0.0;
                  },
                  [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
                  {
                      for(int n=0; n<sediment_exner::amr_nfield; ++n)
                      SP(q)->sed->exner()->amr_field(SP(q)->sed->s,n,f.di,f.dj) = w[n];
                  },
                  [&](const reefamr_fill &f, const double *w)
                  {
                      for(int n=0; n<sediment_exner::amr_nfield; ++n)
                      s0->exner()->amr_field(s0->s,n,f.di,f.dj) = w[n];
                  });

        for(int l=1; l<=maxlev; ++l)
        fill_run(l,sediment_exner::amr_nfield,7960+l,
                 [&](const reefamr_fill &f, double *v)
                 {
                     sediment_f *sf = gsed(f.g);
                     for(int n=0; n<sediment_exner::amr_nfield; ++n)
                     if(f.kind==0)
                     v[n] = sf->exner()->amr_field(sf->s,n,f.si,f.sj);
                     else if(f.kind==1)
                     v[n] = sed_plin((n==0 ? sf->s->qb : sf->exner()->bend_slice(n)),grid(f.g).q,f.si,f.sj,f.ox,f.oy);
                     else
                     v[n] = 0.0;
                 },
                 [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
                 {
                     for(int n=0; n<sediment_exner::amr_nfield; ++n)
                     SP(q)->sed->exner()->amr_field(SP(q)->sed->s,n,f.di,f.dj) = w[n];
                 });
    };

    sync();

    // the lags run on through the interfaces: the relaxations again from the filled cells, until the cells around the
    // grids no longer move them (the change of the first sweep below the tolerance of the relaxation, 1e-6 of the
    // largest equilibrium), at most S 49 passes
    if(lag)
    {
        int pass;
        double res = 0.0;
        for(pass=1; pass<=MAX(p->S49,1); ++pass)
        {
            double d[2] = {0.0,0.0}, f[2] = {0.0,0.0};
            s0->exner()->lag_refresh(p,pgc,s0->s,d,f);
            for(int id=0; id<ng; ++id)
            {
            sflow_amr_patch *c = SP(id);
            sedscope ss(pgc,c->pp);
            c->sed->exner()->lag_refresh(c->pp,pgc,c->sed->s,d,f);
            }
            sync();

            res = 0.0;
            for(int k=0; k<2; ++k)
            {
            d[k] = pgc->globalmax(d[k]);
            f[k] = pgc->globalmax(f[k]);
            if(f[k]>0.0)
            res = MAX(res,d[k]/f[k]);
            }
            if(res<=sediment_exner::lag_tol)
            break;
        }
        if(p->mpirank==0)
        cout<<"lags across the grids: passes "<<MIN(pass,MAX(p->S49,1))<<"  change "<<setprecision(3)<<res<<endl;
    }

    // 3. face fluxes, the box faces of the patches, the matched coarse faces
    vector<std::unique_ptr<slice4>> Fx(ng+1), Fy(ng+1), mask(ng+1);

    for(int g=-1; g<ng; ++g)
    {
        lexer *q = grid(g).q;
        Fx[g+1].reset(new slice4(q));
        Fy[g+1].reset(new slice4(q));
        mask[g+1].reset(new slice4(q));

        for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
        for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
        (*mask[g+1])(ii,jj) = 1.0;

        if(g>=0)
        for(auto &f : SP(g)->fill)
        (*mask[g+1])(f.di,f.dj) = 0.0;

        if(g<0 && nblk>0)
        for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
        for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
        if(nring[lij(q,ii,jj)]==1)
        (*mask[0])(ii,jj) = 0.0;

        if(g<0)
        s0->exner()->amr_flux(p,pgc,s0->s,*Fx[0],*Fy[0]);
        else
        {
        sflow_amr_patch *c = SP(g);
        sedscope ss(pgc,c->pp);
        c->sed->exner()->amr_flux(c->pp,pgc,c->sed->s,*Fx[g+1],*Fy[g+1]);
        }
    }

    for(int id=0; id<ng; ++id)
    {
        sflow_amr_patch &c = *SP(id);
        const int il = EXT-1, ih = EXT+c.nx-1;
        const int jl = EY-1, jh = EY+c.ny-1;
        slice4 &fx = *Fx[id+1], &fy = *Fy[id+1];

        c.sedrec[0].resize(c.ny);
        c.sedrec[1].resize(c.ny);
        c.sedrec[2].resize(c.nx);
        c.sedrec[3].resize(c.nx);

        for(int r=0; r<c.ny; ++r)
        {
        c.sedrec[0][r] = fx(il,EY+r);
        c.sedrec[1][r] = fx(ih,EY+r);
        }
        for(int r=0; r<c.nx; ++r)
        {
        c.sedrec[2][r] = fy(EXT+r,jl);
        c.sedrec[3][r] = fy(EXT+r,jh);
        }
    }

    sedr.resize(rmatch.size());
    for(size_t g=0; g<rmatch.size(); ++g)
    sedr[g].assign(rmatch[g].size(),0.0);

    for(int l=1; l<=maxlev; ++l)
    face_run_at(l,1,7980+l,
                [&](reefamr_patch *q, int side, int r, double *sb)
                {
                    sflow_amr_patch *c = SP(q);
                    sb[0] = 0.5*(c->sedrec[side][r] + c->sedrec[side][r+fr1(side)]);
                },
                [&](int tgt, int idx, const double *v)
                {
                    sedr[tgt][idx] = v[0];
                });

    for(int g=-1; g<ng; ++g)
    {
        const int t = g+1;

        if(t<(int)match.size())
        for(auto &m : match[t])
        {
            sflow_amr_patch &c = *SP(m.child);
            double val = 0.5*(c.sedrec[m.side][m.r] + c.sedrec[m.side][m.r+fr1(m.side)]);
            if(m.sgn<0)
            val = -val;
            if(m.dir==0)
            (*Fx[t])(m.fi,m.fj) = val;
            else
            (*Fy[t])(m.fi,m.fj) = val;
        }

        if(t<(int)rmatch.size())
        for(size_t n=0; n<rmatch[t].size(); ++n)
        {
            const reefamr_match &m = rmatch[t][n];
            const double val = (m.sgn<0) ? -sedr[t][n] : sedr[t][n];
            if(m.dir==0)
            (*Fx[t])(m.fi,m.fj) = val;
            else
            (*Fy[t])(m.fi,m.fj) = val;
        }
    }

    // 4. the covered cells, vz of every grid
    for(int l=maxlev; l>=1; --l)
    block_up(l,1,8000+l,
             [&](reefamr_patch *q, int id, int k, double *v)
             {
                 v[0] = 0.0;
             },
             [&](const reefamr_block &B, int key, const double *v)
             {
                 (*mask[B.g+1])(B.ic,B.jc) = 0.0;
             });

    s0->exner()->amr_vz(p,pgc,s0->s,*Fx[0],*Fy[0],mask[0].get());

    for(int id=0; id<ng; ++id)
    {
        sflow_amr_patch *c = SP(id);
        sedscope ss(pgc,c->pp);
        c->sed->exner()->amr_vz(c->pp,pgc,c->sed->s,*Fx[id+1],*Fy[id+1],mask[id+1].get());
    }

    // the filter (S 100) over all grids, the interfaces of the grids of a level open (reefamr_sed.h)
    reefamr_sedgrids G;
    G.scope = [&](int g, const std::function<void()> &fn)
    {
        if(g<0)
        {
        fn();
        return;
        }
        sflow_amr_patch *c = SP(g);
        sedscope ss(pgc,c->pp);
        fn();
    };
    for(int g=-1; g<ng; ++g)
    {
    G.L.push_back(grid(g).q);
    G.S.push_back(gsed(g));
    G.mask.push_back(mask[g+1].get());
    }

    if(p->S100>0)
    reefamr_sed_filter(*this,G,8300);

    double vzmax = s0->exner()->amr_vzmax(p,s0->s,mask[0].get());
    double dxmin = 1.0e20;
    double brate = s0->exner()->bend_rate;

    for(int id=0; id<ng; ++id)
    {
        sflow_amr_patch *c = SP(id);
        sedscope ss(pgc,c->pp);
        vzmax = MAX(vzmax,c->sed->exner()->amr_vzmax(c->pp,c->sed->s,mask[id+1].get()));
        dxmin = MIN(dxmin,c->sed->exner()->amr_dx(c->pp));
        brate = MAX(brate,c->sed->exner()->bend_rate);
    }

    SLICEBASELOOP
    if((*mask[0])(i,j)<0.5)
    s0->s->vz(i,j) = 0.0;

    // 5. one dtsed for all grids
    vzmax = pgc->globalmax(vzmax);
    dxmin = pgc->globalmin(dxmin);

    s0->exner()->amr_maxvz = vzmax;
    s0->exner()->amr_mindx = (dxmin<1.0e19) ? dxmin : 0.0;
    s0->exner()->amr_bendrate = brate;
    s0->exner()->timestep(p,pgc,s0->s);
    s0->exner()->amr_maxvz = 0.0;
    s0->exner()->amr_mindx = 0.0;
    s0->exner()->amr_bendrate = 0.0;

    for(int id=0; id<ng; ++id)
    SP(id)->pp->dtsed = p->dtsed;

    // 6. bed: Exner, sand slide (reefamr_sed.h), restriction, the flow
    s0->exner()->amr_bed(p,pgc,s0->s);

    for(int id=0; id<ng; ++id)
    {
        sflow_amr_patch *c = SP(id);
        sedscope ss(pgc,c->pp);
        c->sed->exner()->amr_bed(c->pp,pgc,c->sed->s);
    }

    // the sand slide over all grids (reefamr_sed.h)
    if(p->S90>0 && p->sediter%p->S94==0)
    reefamr_sed_slide(*this,pgc,G,8400);

    for(int l=maxlev; l>=1; --l)
    block_up(l,1,8020+l,
             [&](reefamr_patch *q, int id, int k, double *v)
             {
                 slice4 &bz = SP(q)->sed->s->bedzh;
                 const int nby = q->ny>>sy;
                 const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                 double w[4];
                 cw4(q->pp,i0,j0,w);
                 v[0] = (sy==0) ? 0.5*(bz(i0,j0) + bz(i0+1,j0)) : w[0]*bz(i0,j0) + w[1]*bz(i0+1,j0) + w[2]*bz(i0,j0+1) + w[3]*bz(i0+1,j0+1);
             },
             [&](const reefamr_block &B, int key, const double *v)
             {
                 gsed(B.g)->s->bedzh(B.ic,B.jc) = v[0];
             });

    pgc->gcsl_start4(p,s0->s->bedzh,1);
    s0->amr_finish_sflow(p,b,pgc,pflow);

    for(int id=0; id<ng; ++id)
    {
        sflow_amr_patch *c = SP(id);
        sedscope ss(pgc,c->pp);
        c->sed->amr_finish_sflow(c->pp,c->b,pgc,pflow_of(c));
    }

    // the cells around the blocks and the ring, the cells around the patches, coarse to fine
    if(nblk>0)
    sed_fill_bed0(pgc);

    for(int l=1; l<=maxlev; ++l)
    sed_fill_bed(pgc,l);

    bedmemo.clear();

    p->sedtime += p->dtsed;
}

// the state of every grid of this rank (own cells) for the tests: box (level, block, I0, J0, nx, ny), eta, WL,
// bed, U, V, in the format of the regression dump (R3DREG01), to <REEF3D_AMR_DUMP>/amr_r<rank>.bin
void sflow_amr::dump_grids(lexer *p)
{
    static const char *dir = std::getenv("REEF3D_AMR_DUMP");
    if(dir==nullptr)
    return;

    vector<string> names;
    vector<vector<double>> data;

    auto grid2 = [&](int id, fdm2D *d, int lv, int bk, int I0, int J0, int i0, int nx, int ny)
    {
        auto add = [&](const char *v)
        {
            std::ostringstream s;
            s<<"g"<<id<<"_"<<v;
            names.push_back(s.str());
            data.push_back(vector<double>());
            return &data.back();
        };
        vector<double> *a = add("box");
        a->push_back(lv); a->push_back(bk); a->push_back(I0); a->push_back(J0); a->push_back(nx); a->push_back(ny); a->push_back(1);

        slice *f2[5] = {&d->eta,&d->WL,&d->bed,&d->U,&d->V};
        const char *nm[5] = {"eta","WL","bed","U","V"};
        for(int m=0; m<5; ++m)
        {
            vector<double> *f = add(nm[m]);
            for(int ii=i0; ii<i0+nx; ++ii)
            for(int jj=i0; jj<i0+ny; ++jj)
            f->push_back((*f2[m])(ii,jj));
        }
    };

    grid2(0,b0,0,0,O0i,O0j,0,p0->knox,p0->knoy);
    for(int n=0; n<(int)P.size(); ++n)
    grid2(n+1,SP(n)->b,P[n]->lev,P[n]->blk,P[n]->I0,P[n]->J0,EXT,P[n]->nx,P[n]->ny);

    std::ostringstream fn;
    fn<<dir<<"/amr_r"<<p->mpirank<<".bin";
    std::ofstream out(fn.str().c_str(),std::ios::binary);
    const char magic[8] = {'R','3','D','R','E','G','0','1'};
    out.write(magic,8);
    int32_t ival;
    double dval;
    ival=p->mpirank;  out.write((char*)&ival,sizeof(ival));
    ival=p->M10;      out.write((char*)&ival,sizeof(ival));
    ival=p->count;    out.write((char*)&ival,sizeof(ival));
    dval=p->simtime;  out.write((char*)&dval,sizeof(dval));
    ival=int32_t(names.size()); out.write((char*)&ival,sizeof(ival));
    for(size_t q=0; q<names.size(); ++q)
    {
        char name[16];
        std::memset(name,0,16);
        std::strncpy(name,names[q].c_str(),15);
        out.write(name,16);
        int64_t num = int64_t(data[q].size());
        out.write((char*)&num,sizeof(num));
        out.write((char*)data[q].data(),sizeof(double)*data[q].size());
    }
}
