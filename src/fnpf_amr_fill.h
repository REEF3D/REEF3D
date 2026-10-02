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

#ifndef FNPF_AMR_FILL_H_
#define FNPF_AMR_FILL_H_

// fill, restriction and interpolation templates of fnpf_amr (fnpf_amr.cpp, fnpf_amr_lap.cpp)

#include"fnpf_amr.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"slice.h"

// --------------------------------------------------------------------- fill, restriction
// cells around the level-l patches; sel(g,m) gives the slice m of grid g
template<class SEL>
inline void fnpf_amr::fill_sl(int l, int ns, int tag, SEL sel)
{
    fill_run(l,ns,tag,
             [&](const reefamr_fill &f, double *v)
             {
                 for(int m=0; m<ns; ++m)
                 {
                     if(f.kind==0)
                     v[m] = sel(f.g,m)(f.si,f.sj);
                     else if(f.kind==1)
                     v[m] = pq(sel(f.g,m),glex(f.g),f.si,f.sj,f.ox,f.oy);
                     else
                     v[m] = 0.0;
                 }
             },
             [&](reefamr_patch *c, int id, const reefamr_fill &f, const double *w)
             {
                 for(int m=0; m<ns; ++m)
                 sel(id,m)(f.di,f.dj) = w[m];
             });
}

// columns around the level-l patches; sel(g) gives the Fi-layout array of grid g
template<class SEL>
inline void fnpf_amr::fill_col(int l, int tag, SEL sel)
{
    const int knf = p0->knoz*((vref==2) ? (1<<l) : 1);
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
                 pcol(f.g,f.si,f.sj,f.ox,f.oy,sel(f.g),knf,v);
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

// restriction of point values to the coarse cell centre between the 2x2 children (i0..i0+1,
// j0..j0+1): cubic in x and y on the 4x4 fine cells around them (4th order), the plain
// average next to solids and in the outermost blocks of the patch (bi, bj): there the 4x4
// stencil reached into the cells around the patch, which are filled after the restriction
// from the coarse values it gives.  The restricted columns then depended on the previous
// fill, and the composite Laplace operator was no fixed linear map of the leaf unknowns:
// the recursive BiCGStab residual converged while the true residual at the patch edges
// stayed orders of magnitude larger, and the psi solves stalled at N 46
inline bool fnpf_amr::rcubic(lexer *pp, int i0, int j0)
{
    if(rorder<4)
    return false;
    for(int a=-1; a<=2; ++a)
    for(int b=-1; b<=2; ++b)
    if(pp->flagslice4[(i0+a-pp->imin)*pp->jmax + (j0+b-pp->jmin)]<0)
    return false;
    return true;
}

template<class F>
inline double fnpf_amr::rc4(F f)
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

// covered cells: restricted from the 2x2 children, finest first
template<class SEL>
inline void fnpf_amr::restrict_sl(int ns, SEL sel)
{
    for(int l=maxlev; l>=1; --l)
    for(int id : lev[l])
    {
        reefamr_patch *c = P[id];
        const int nby = c->ny/2;

        for(int bi=0; bi<c->nx/2; ++bi)
        for(int bj=0; bj<nby; ++bj)
        {
            const int k = bi*nby+bj;
            const int g = c->rgrid[k];
            if(g<-1)
            continue;

            const int i0 = EXT+2*bi, j0 = EXT+2*bj;
            const bool hi = (bi>0 && bi<c->nx/2-1 && bj>0 && bj<nby-1) && rcubic(c->pp,i0,j0);
            for(int m=0; m<ns; ++m)
            {
                slice &f = sel(id,m);
                if(hi)
                sel(g,m)(c->ric[k],c->rjc[k]) = rc4([&](int a, int b) { return f(i0+a,j0+b); });
                else
                sel(g,m)(c->ric[k],c->rjc[k]) = 0.25*(f(i0,j0)+f(i0+1,j0)+f(i0,j0+1)+f(i0+1,j0+1));
            }
        }
    }
}

template<class SEL>
inline void fnpf_amr::restrict_col(SEL sel)
{
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
            const int knc = q->knoz;
            const int fz = pp->knoz/knc;
            double *dst = sel(g);
            const int i0 = EXT+2*bi, j0 = EXT+2*bj;

            const bool hi = (bi>0 && bi<c->nx/2-1 && bj>0 && bj<nby-1) && rcubic(pp,i0,j0);
            const int sI = pp->jmax*pp->kmaxF, sJ = pp->kmaxF;
            for(int K=0; K<=knc; ++K)
            {
                const int kf = fz*K;
                const int n0 = fidx(pp,i0,j0,kf);
                if(hi)
                dst[fidx(q,c->ric[k],c->rjc[k],K)] = rc4([&](int a, int b) { return src[n0+a*sI+b*sJ]; });
                else
                dst[fidx(q,c->ric[k],c->rjc[k],K)] = 0.25*(src[n0] + src[n0+sI] + src[n0+sJ] + src[n0+sI+sJ]);
            }
        }
    }
}

// interior of a patch from its parent grid
template<class SEL>
inline void fnpf_amr::prolong_interior_sl(fnpf_amr_patch &c, int ns, SEL sel)
{
    int id=-1;
    for(int n=0; n<(int)P.size(); ++n)
    if(P[n]==&c)
    id=n;

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
        for(int m=0; m<ns; ++m)
        sel(id,m)(EXT+2*bi+a,EXT+2*bj+d) = pq(sel(g,m),glex(g),c.ric[k],c.rjc[k],a==0?-1:1,d==0?-1:1);
    }
}

template<class SEL>
inline void fnpf_amr::prolong_interior_col(fnpf_amr_patch &c, SEL sel)
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
            pcol(g,c.ric[k],c.rjc[k],a==0?-1:1,d==0?-1:1,sel(g),knf,&v[0]);
            for(int kk=0; kk<=knf; ++kk)
            dst[fidx(pp,EXT+2*bi+a,EXT+2*bj+d,kk)] = v[kk];
        }
    }
}


#endif
