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

#ifndef NHFLOW_AMR_FILL_H_
#define NHFLOW_AMR_FILL_H_

// column fill, restriction and interpolation templates of nhflow_amr (F layout: the pressure
// nodes 0..knoz of a column; with A 281 the nodes of a level are nested in those of the next), as
// fnpf_amr_fill.h

#include"nhflow_amr.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"slice.h"

namespace nhflow_amr_detail
{
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

// restriction of point values to the coarse cell centre: cubic in x and y on the 4x4 fine cells
// around the 2x2 children, the plain average next to solids
inline bool rcubic(const lexer *pp, int i0, int j0)
{
    for(int a=-1; a<=2; ++a)
    for(int b=-1; b<=2; ++b)
    if(pp->flagslice4[(i0+a-pp->imin)*pp->jmax + (j0+b-pp->jmin)]<0)
    return false;
    return true;
}

template<class F>
inline double rc4(F f)
{
    const double w[4] = {-1.0/16.0, 9.0/16.0, 9.0/16.0, -1.0/16.0};
    double r = 0.0;
    for(int a=0; a<4; ++a)
    {
        double ra = 0.0;
        for(int b=0; b<4; ++b)
        ra += w[b]*f(a-1,b-1);
        r += w[a]*ra;
    }
    return r;
}
}

// columns around the level-l patches; sel(g) gives the F-layout array of grid g (from the coarser
// level with A 281: the coarse nodes and the midpoints between them)
template<class SEL>
inline void nhflow_amr::fill_col(int l, int tag, SEL sel)
{
    const int knf = klev(l);
    const int nv = knf+1;

    fill_run(l,nv,tag,
             [&](const reefamr_fill &f, double *v)
             {
                 if(f.kind==0)
                 {
                     lexer *q = glex(f.g);
                     const double *src = sel(f.g);
                     for(int kk=0; kk<=knf; ++kk)
                     v[kk] = src[fidx(q,f.si,f.sj,kk)];
                 }
                 else if(f.kind==1)
                 pcol(f.g,f.si,f.sj,f.ox,f.oy,sel(f.g),knf,v,2);
                 else
                 for(int kk=0; kk<=knf; ++kk)
                 v[kk] = 0.0;
             },
             [&](reefamr_patch *c, int id, const reefamr_fill &f, const double *w)
             {
                 lexer *pp = c->pp;
                 double *dst = sel(id);
                 for(int kk=0; kk<=knf; ++kk)
                 dst[fidx(pp,f.di,f.dj,kk)] = w[kk];
             });
}

// covered columns: restricted from the 2x2 children (cubic where possible), finest first; with
// A 281 coarse node K is fine node 2K
template<class SEL>
inline void nhflow_amr::restrict_col(SEL sel)
{
    using namespace nhflow_amr_detail;

    for(int l=maxlev; l>=1; --l)
    for(int id : lev[l])
    {
        reefamr_patch *c = P[id];
        lexer *pp = c->pp;
        const double *src = sel(id);
        const int nby = c->ny/2;

        for(int bi=0; bi<c->nx/2; ++bi)
        for(int bj=0; bj<nby; ++bj)
        {
            const int k = bi*nby+bj;
            const int g = c->rgrid[k];
            if(g<-1)
            continue;

            lexer *q = glex(g);
            double *dst = sel(g);
            const int i0 = EXT+2*bi, j0 = EXT+2*bj;

            // cubic only on interior blocks: on the edge blocks the 4x4 stencil would read the
            // EXT cells of the patch, which are not unknowns of the composite solve (as fnpf_amr)
            bool hi = (bi>0 && bi<c->nx/2-1 && bj>0 && bj<nby-1) && rcubic(pp,i0,j0);
            const int sI = pp->jmax*pp->kmaxF, sJ = pp->kmaxF;
            const int fz = pp->knoz/q->knoz;

            // A 283: only the children that are unknowns of the pressure (wet and deep); cubic only
            // if the whole 4x4 stencil is; none of them: the coarse column is shallow or dry, P = 0
            int na = 4;
            double wa[4] = {0.25,0.25,0.25,0.25};
            if(shore)
            {
                for(int a=-1; a<=2 && hi; ++a)
                for(int b=-1; b<=2 && hi; ++b)
                if(!wet_at(pp,i0+a,j0+b,2))
                hi = false;

                na = 0;
                for(int a=0; a<2; ++a)
                for(int b=0; b<2; ++b)
                {
                    wa[2*a+b] = wet_at(pp,i0+a,j0+b,2) ? 1.0 : 0.0;
                    na += (int)wa[2*a+b];
                }
                for(int m=0; m<4; ++m)
                wa[m] = (na>0) ? wa[m]/double(na) : 0.0;
            }

            for(int K=0; K<=q->knoz; ++K)
            {
                const int n0 = fidx(pp,i0,j0,fz*K);
                if(hi)
                dst[fidx(q,c->ric[k],c->rjc[k],K)] = rc4([&](int a, int b) { return src[n0+a*sI+b*sJ]; });
                else if(na==4)
                dst[fidx(q,c->ric[k],c->rjc[k],K)] = 0.25*(src[n0] + src[n0+sI] + src[n0+sJ] + src[n0+sI+sJ]);
                else
                dst[fidx(q,c->ric[k],c->rjc[k],K)] = wa[0]*src[n0] + wa[2]*src[n0+sI] + wa[1]*src[n0+sJ] + wa[3]*src[n0+sI+sJ];
            }
        }
    }
}

// interior columns of a patch from its parent grid
template<class SEL>
inline void nhflow_amr::prolong_interior_col(nhflow_amr_patch &c, SEL sel)
{
    int id=-1;
    for(int n=0; n<(int)P.size(); ++n)
    if(P[n]==&c)
    id=n;

    lexer *pp = c.pp;
    const int knf = pp->knoz;
    vector<double> v(knf+1);
    double *dst = sel(id);

    const int nby = c.ny/2;
    for(int bi=0; bi<c.nx/2; ++bi)
    for(int bj=0; bj<nby; ++bj)
    {
        const int k = bi*nby+bj;
        const int g = c.rgrid[k];
        if(g<-1)
        continue;

        for(int a=0; a<2; ++a)
        for(int d=0; d<2; ++d)
        {
            pcol(g,c.ric[k],c.rjc[k],a==0?-1:1,d==0?-1:1,sel(g),knf,&v[0],2);
            for(int kk=0; kk<=knf; ++kk)
            dst[fidx(pp,EXT+2*bi+a,EXT+2*bj+d,kk)] = v[kk];
        }
    }
}

#endif
