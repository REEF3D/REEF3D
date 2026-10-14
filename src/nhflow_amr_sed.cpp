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
#include"nhflow_amr_fill.h"
#include"sediment_f.h"
#include"sediment_fdm.h"
#include"sediment_exner.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice4.h"
#include"reefamr_sed.h"
#include<memory>

using namespace nhflow_amr_detail;

// Sediment transport (S 10) on the grids of the mesh refinement.
//
// Every grid (level 0 and every patch) has its own sediment_f (the patches: sediment_f with amrgrid 1, no
// log and no probes). A sediment step of level 0 (sediment_f::start_nhflow decides when) runs all grids:
//   1. the physics of every grid: bed slope, roughness, bed shear, bed load and its direction, qb
//   2. qb in the cells around the patches: a patch of the same level, else linear from the coarser level
//   3. the Exner face fluxes of every grid (sediment_exner::amr_flux); the patches record their box faces,
//      a coarse face next to a finer patch takes the mean of its two fine faces (flux matching, as the
//      momentum and continuity fluxes of the flow): the bed volume is conserved over the leaf cells
//   4. vz of every grid; the filter (S 100) over all grids in steps with the exchanges between the grids
//      (reefamr_sed.h): the interfaces of the grids of a level (blocks, the ring of the native grid, patches
//      that meet) are open as the inside of a grid; the covered cells and the cells around a patch with a
//      coarser source take no part
//   5. one dtsed for all grids (the largest |vz| and the smallest cell of all grids)
//   6. the bed of every grid, the sand slide over all grids (reefamr_sed.h: open interfaces of a level, the
//      material slid into the cells around a grid goes to the cell they stand for; covered cells closed), the
//      restriction of the bed to the covered cells (mean per area of the children), the bed to the flow, the
//      bed of the cells around the patches (as at the regrid: linear from the coarser level)
// The flow keeps the water level WL of a column through the bed update; its eta is that of the old bed until
// the first stage (as on the uniform grid). The cells around the patches therefore take the water level from
// the source (fill_stage with S 10), not eta + depth.
// A fresh patch takes the bed of the old patches of its level where they were (regrid_state), else the
// linear (mean preserving) bed of its parent; its sediment starts there (sed_patch_ini).

// the sediment code tests the global domain edges with origin_i/j and gknox/y (filter, sand slide); a patch
// has no global position (origin -1e8): during its sediment calls the edges are those of its cells, the
// cells around it are open or closed by the mask (filter) and by DFBED (sand slide) the driver sets
namespace
{
struct sedscope
{
    lexer *q;
    int oi, oj, gx, gy;
    sedscope(lexer *pp) : q(pp), oi(pp->origin_i), oj(pp->origin_j), gx(pp->gknox), gy(pp->gknoy)
    {
        q->origin_i = q->origin_j = 0;
        q->gknox = q->gknoy = 1<<29;
    }
    ~sedscope()
    {
        q->origin_i = oi; q->origin_j = oj;
        q->gknox = gx; q->gknoy = gy;
    }
};
}

void nhflow_amr::sed_attach(sediment *ps)
{
    if(p0->S10==0)
    return;

    sed0 = dynamic_cast<sediment_f*>(ps);

    if(sed0!=nullptr && (maxlev>0 || nblk>0))
    sed0->grids = this;
}

void nhflow_amr::sed_patch_ini(ghostcell *pgc, nhflow_amr_patch &c)
{
    if(c.sed==nullptr)
    return;

    lexer *pp = c.pp;
    fdm_nhf *d = c.d;
    pscope ps(pgc,d,d0);
    sedscope ss(pp);

    c.sed->ini_nhflow(pp,d,pgc);

    // the bed of the sediment in all cells of the patch, the cells around it included
    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
    c.sed->s->bedzh(ii,jj) = d->bed(ii,jj);
    c.sed->s->bedzh0(ii,jj) = d->bed(ii,jj);
    }
}

// bed of the cells around the patches (coarse to fine) and the depth of the patches

void nhflow_amr::sed_bed_ghosts(ghostcell *pgc)
{
    // multi-block grid: the cells around the blocks and the ring of the native grid
    if(nblk>0)
    {
        fill_bed0(pgc);

        for(int id : lev[0])
        {
            nhflow_amr_patch *c = NP(id);
            lexer *pp = c->pp;
            fdm_nhf *d = c->d;

            for(auto &f : c->fill)
            c->sed->s->bedzh(f.di,f.dj) = d->bed(f.di,f.dj);

            for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
            for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
            {
                d->depth(ii,jj) = pp->wd - d->bed(ii,jj);
                pp->bed[lij(pp,ii,jj)] = d->bed(ii,jj);
            }
        }

        for(int ii=p0->imin; ii<p0->imin+p0->imax; ++ii)
        for(int jj=p0->jmin; jj<p0->jmin+p0->jmax; ++jj)
        if(nring[lij(p0,ii,jj)]==1)
        sed0->s->bedzh(ii,jj) = d0->bed(ii,jj);
    }

    for(int l=1; l<=maxlev; ++l)
    {
        fill_bed(pgc,l);

        for(int id : lev[l])
        {
            nhflow_amr_patch *c = NP(id);
            lexer *pp = c->pp;
            fdm_nhf *d = c->d;

            for(auto &f : c->fill)
            c->sed->s->bedzh(f.di,f.dj) = d->bed(f.di,f.dj);

            for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
            for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
            {
                d->depth(ii,jj) = pp->wd - d->bed(ii,jj);
                pp->bed[lij(pp,ii,jj)] = d->bed(ii,jj);
            }
        }
    }
}

void nhflow_amr::sediment_step_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, sediment_f *s0)
{
    const int ng = (int)P.size();

    // 0. the cells around the patches with the current state (they hold the input of the last stage): two
    //    patches of a level that meet compute the same flux on their common face
    if(nblk>0)
    fill_level0(pgc,-1);

    for(int l=1; l<=maxlev; ++l)
    fill_stage(pgc,l,-1);

    // 1. physics of every grid
    s0->amr_physics_nhflow(p,d,pgc);

    for(int id=0; id<ng; ++id)
    {
        nhflow_amr_patch *c = NP(id);
        pscope ps(pgc,c->d,d0);
        sedscope ss(c->pp);
        c->sed->amr_physics_nhflow(c->pp,c->d,pgc);
    }

    // the covered cells: qb of their children (mean per area): the cells around a patch take qb linear from
    // the coarser level, with the slopes of its neighbours, covered ones among them
    for(int l=maxlev; l>=1; --l)
    block_up(l,1,7880+l,
             [&](reefamr_patch *q, int id, int k, double *v)
             {
                 slice4 &f = NP(q)->sed->s->qb;
                 const int nby = q->ny/2;
                 const int i0 = EXT+2*(k/nby), j0 = EXT+2*(k%nby);
                 double w[4];
                 cw4(q->pp,i0,j0,w);
                 v[0] = w[0]*f(i0,j0) + w[1]*f(i0+1,j0) + w[2]*f(i0,j0+1) + w[3]*f(i0+1,j0+1);
             },
             [&](const reefamr_block &B, int key, const double *v)
             {
                 gsed(B.g)->s->qb(B.ic,B.jc) = v[0];
             });

    // 2. qb around the blocks and in the ring of the native grid (multi-block grid), around the patches
    if(nblk>0)
    fill_run0(1,7898,
              [&](const reefamr_fill &f, double *v)
              {
                  v[0] = (f.kind==0) ? gsed(f.g)->s->qb(f.si,f.sj) : 0.0;
              },
              [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
              {
                  NP(q)->sed->s->qb(f.di,f.dj) = w[0];
              },
              [&](const reefamr_fill &f, const double *w)
              {
                  s0->s->qb(f.di,f.dj) = w[0];
              });

    for(int l=1; l<=maxlev; ++l)
    fill_run(l,1,7900+l,
             [&](const reefamr_fill &f, double *v)
             {
                 if(f.kind==0)
                 v[0] = gsed(f.g)->s->qb(f.si,f.sj);
                 else if(f.kind==1)
                 v[0] = plin(gsed(f.g)->s->qb,glex(f.g),f.si,f.sj,f.ox,f.oy);
                 else
                 v[0] = 0.0;
             },
             [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
             {
                 NP(q)->sed->s->qb(f.di,f.dj) = w[0];
             });

    // 3. face fluxes, the box faces of the patches, the matched coarse faces
    vector<std::unique_ptr<slice4>> Fx(ng+1), Fy(ng+1), mask(ng+1);

    for(int g=-1; g<ng; ++g)
    {
        lexer *q = glex(g);
        Fx[g+1].reset(new slice4(q));
        Fy[g+1].reset(new slice4(q));
        mask[g+1].reset(new slice4(q));

        for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
        for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
        (*mask[g+1])(ii,jj) = 1.0;

        // the cells around a patch take no part in the filter and the sand slide (they are refilled)
        if(g>=0)
        for(auto &f : NP(g)->fill)
        (*mask[g+1])(f.di,f.dj) = 0.0;

        // multi-block grid: the ring of the native grid (copies of the blocks)
        if(g<0 && nblk>0)
        for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
        for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
        if(nring[lij(q,ii,jj)]==1)
        (*mask[0])(ii,jj) = 0.0;

        if(g<0)
        s0->exner()->amr_flux(p,pgc,s0->s,*Fx[0],*Fy[0]);
        else
        {
        nhflow_amr_patch *c = NP(g);
        pscope ps(pgc,c->d,d0);
        sedscope ss(c->pp);
        c->sed->exner()->amr_flux(c->pp,pgc,c->sed->s,*Fx[g+1],*Fy[g+1]);
        }
    }

    for(int id=0; id<ng; ++id)
    {
        nhflow_amr_patch &c = *NP(id);
        const int il = EXT-1, ih = EXT+c.nx-1;
        const int jl = EXT-1, jh = EXT+c.ny-1;
        slice4 &fx = *Fx[id+1], &fy = *Fy[id+1];

        c.sedrec[0].resize(c.ny);
        c.sedrec[1].resize(c.ny);
        c.sedrec[2].resize(c.nx);
        c.sedrec[3].resize(c.nx);

        for(int r=0; r<c.ny; ++r)
        {
        c.sedrec[0][r] = fx(il,EXT+r);
        c.sedrec[1][r] = fx(ih,EXT+r);
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
    face_run_at(l,1,7920+l,
                [&](reefamr_patch *q, int side, int r, double *sb)
                {
                    nhflow_amr_patch *c = NP(q);
                    sb[0] = 0.5*(c->sedrec[side][r] + c->sedrec[side][r+1]);
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
            nhflow_amr_patch &c = *NP(m.child);
            double val = 0.5*(c.sedrec[m.side][m.r] + c.sedrec[m.side][m.r+1]);
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
    block_up(l,1,7940+l,
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
        nhflow_amr_patch *c = NP(id);
        pscope ps(pgc,c->d,d0);
        sedscope ss(c->pp);
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
        nhflow_amr_patch *c = NP(g);
        pscope ps(pgc,c->d,d0);
        sedscope ss(c->pp);
        fn();
    };
    for(int g=-1; g<ng; ++g)
    {
    G.L.push_back(glex(g));
    G.S.push_back(gsed(g));
    G.mask.push_back(mask[g+1].get());
    }

    if(p->S100>0)
    reefamr_sed_filter(*this,G,8100);

    double vzmax = s0->exner()->amr_vzmax(p,s0->s,mask[0].get());
    double dxmin = 1.0e20;

    for(int id=0; id<ng; ++id)
    {
        nhflow_amr_patch *c = NP(id);
        pscope ps(pgc,c->d,d0);
        sedscope ss(c->pp);
        vzmax = MAX(vzmax,c->sed->exner()->amr_vzmax(c->pp,c->sed->s,mask[id+1].get()));
        dxmin = MIN(dxmin,c->sed->exner()->amr_dx(c->pp));
    }

    // level 0: no bed change in the covered cells (their bed is the restriction)
    SLICEBASELOOP
    if((*mask[0])(i,j)<0.5)
    s0->s->vz(i,j) = 0.0;

    // 5. one dtsed for all grids
    vzmax = pgc->globalmax(vzmax);
    dxmin = pgc->globalmin(dxmin);

    s0->exner()->amr_maxvz = vzmax;
    s0->exner()->amr_mindx = (dxmin<1.0e19) ? dxmin : 0.0;
    s0->exner()->timestep(p,pgc,s0->s);
    s0->exner()->amr_maxvz = 0.0;
    s0->exner()->amr_mindx = 0.0;

    for(int id=0; id<ng; ++id)
    NP(id)->pp->dtsed = p->dtsed;

    // 6. bed: Exner, sand slide (reefamr_sed.h), restriction, the flow
    s0->exner()->amr_bed(p,pgc,s0->s);

    for(int id=0; id<ng; ++id)
    {
        nhflow_amr_patch *c = NP(id);
        pscope ps(pgc,c->d,d0);
        sedscope ss(c->pp);
        c->sed->exner()->amr_bed(c->pp,pgc,c->sed->s);
    }

    // the sand slide over all grids (reefamr_sed.h)
    if(p->S90>0 && p->sediter%p->S94==0)
    reefamr_sed_slide(*this,pgc,G,8200);

    for(int l=maxlev; l>=1; --l)
    block_up(l,1,7960+l,
             [&](reefamr_patch *q, int id, int k, double *v)
             {
                 slice4 &b = NP(q)->sed->s->bedzh;
                 const int nby = q->ny/2;
                 const int i0 = EXT+2*(k/nby), j0 = EXT+2*(k%nby);
                 double w[4];
                 cw4(q->pp,i0,j0,w);
                 v[0] = w[0]*b(i0,j0) + w[1]*b(i0+1,j0) + w[2]*b(i0,j0+1) + w[3]*b(i0+1,j0+1);
             },
             [&](const reefamr_block &B, int key, const double *v)
             {
                 gsed(B.g)->s->bedzh(B.ic,B.jc) = v[0];
             });

    pgc->gcsl_start4(p,s0->s->bedzh,1);
    s0->amr_finish_nhflow(p,d,pgc,pflow);

    for(int id=0; id<ng; ++id)
    {
        nhflow_amr_patch *c = NP(id);
        pscope ps(pgc,c->d,d0);
        sedscope ss(c->pp);
        c->sed->amr_finish_nhflow(c->pp,c->d,pgc,c->pflow);
    }

    sed_bed_ghosts(pgc);

    p->sedtime += p->dtsed;
}
