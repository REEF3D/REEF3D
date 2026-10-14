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

#ifndef REEFAMR_SED_H_
#define REEFAMR_SED_H_

#include"sediment_xfilter.h"
#include"sediment_f.h"
#include"sediment_fdm.h"
#include"sandslide.h"
#include"lexer.h"
#include"ghostcell.h"
#include<memory>
#include<functional>
#include<iostream>
#include<cmath>
#include<utility>

// The bed filter (S 100) and the sand slide (S 90) over all grids of the mesh refinement and of a multi-block
// grid (nhflow_amr_sed.cpp, sflow_amr_sed.cpp): the grids run each step in turn and exchange their ghost cells
// between the steps (reefamr::xslice, reefamr::xslice_add), so that the interfaces of the grids of a level
// (blocks, the ring of the native grid, patches that meet) are open as the cells inside a grid. A cell
// around a patch with a source on a coarser level and a covered cell are closed (no filter flux, no slid
// material): the bed volume over the leaf cells stays conserved.

struct reefamr_sedgrids
{
    vector<lexer*> L;               // grid g at g+1 (0: the native grid)
    vector<sediment_f*> S;
    vector<slice*> mask;            // 1: the cells of the grid (not around it, not the ring, not covered)
    std::function<void(int, const std::function<void()>&)> scope;   // runs a function in the scope of grid g
};

template<class R>
void reefamr_sed_filter(R &amr, reefamr_sedgrids &G, int tag)
{
    lexer *p = G.L[0];
    const int n = (int)G.L.size();
    vector<std::unique_ptr<sediment_xfilter>> X(n);

    for(int t=0; t<n; ++t)
    X[t].reset(new sediment_xfilter(G.L[t]));

    auto each = [&](const std::function<void(int)> &fn)
    {
        for(int t=0; t<n; ++t)
        G.scope(t-1,[&]() { fn(t); });
    };
    auto vz = [&](int t) -> slice& { return G.S[t]->s->vz; };
    auto ab = [&](int g, int k) -> slice* { return (k==0) ? (slice*)&X[g+1]->a : (slice*)&X[g+1]->b; };

    // the intermediate fields Sx(.) and Sy(.) belong to the index frame of their grid: a cell around a grid
    // whose source grid is turned against it (its i lines are the j lines here, curvilinear blocks) takes them
    // swapped. The i direction of every cell (from its nodes) goes with the exchange once.
    vector<std::unique_ptr<slice4>> EX(n), EY(n), SW(n);

    auto e1 = [](lexer *q, int ii, int jj, double &ex, double &ey)
    {
        ex = 1.0;
        ey = 0.0;
        if(q->geo_curv==1)
        {
        auto xn = [&](int a, int b) { return q->XN2D[q->nij(a,b)]; };
        auto yn = [&](int a, int b) { return q->YN2D[q->nij(a,b)]; };
        ex = 0.5*(xn(ii+1,jj)+xn(ii+1,jj+1)) - 0.5*(xn(ii,jj)+xn(ii,jj+1));
        ey = 0.5*(yn(ii+1,jj)+yn(ii+1,jj+1)) - 0.5*(yn(ii,jj)+yn(ii,jj+1));
        }
    };

    for(int t=0; t<n; ++t)
    {
        lexer *q = G.L[t];
        EX[t].reset(new slice4(q));
        EY[t].reset(new slice4(q));
        SW[t].reset(new slice4(q));
        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        e1(q,ii,jj,(*EX[t])(ii,jj),(*EY[t])(ii,jj));
    }

    amr.xslice(2,tag,[&](int g, int k) -> slice* { return (k==0) ? EX[g+1].get() : EY[g+1].get(); });

    for(int t=0; t<n; ++t)
    {
        lexer *q = G.L[t];
        double ox, oy;
        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        {
            e1(q,ii,jj,ox,oy);
            const double sx = (*EX[t])(ii,jj), sy = (*EY[t])(ii,jj);
            (*SW[t])(ii,jj) = (fabs(ox*sy-oy*sx)>fabs(ox*sx+oy*sy)) ? 1.0 : 0.0;
        }
    }

    auto xab = [&]()
    {
        amr.xslice(2,tag,ab);

        for(int t=0; t<n; ++t)
        {
            lexer *q = G.L[t];
            slice4 &fa = X[t]->a, &fb = X[t]->b, &sw = *SW[t];
            for(int ii=0; ii<q->knox; ++ii)
            for(int jj=0; jj<q->knoy; ++jj)
            if(sw(ii,jj)>0.5)
            std::swap(fa(ii,jj),fb(ii,jj));
        }
    };

    each([&](int t) { X[t]->mask(G.L[t],G.mask[t]); });
    amr.xslice(1,tag,[&](int g, int k) -> slice* { return &X[g+1]->m; });
    each([&](int t) { X[t]->weights(G.L[t]); });

    for(int qn=0; qn<p->S100; ++qn)
    {
        each([&](int t) { X[t]->copy(G.L[t],vz(t),0); });
        amr.xslice(1,tag,[&](int g, int k) -> slice* { return &X[g+1]->h; });
        each([&](int t) { X[t]->sweep1(G.L[t],1); });
        xab();
        each([&](int t) { X[t]->sweep2(G.L[t],vz(t),1); });

        for(int qqn=0; qqn<p->S101; ++qqn)
        {
            each([&](int t) { X[t]->copy(G.L[t],vz(t),1); });
            amr.xslice(1,tag,[&](int g, int k) -> slice* { return &X[g+1]->dh; });
            each([&](int t) { X[t]->sweep1(G.L[t],2); });
            xab();
            each([&](int t) { X[t]->sweep2(G.L[t],vz(t),0); });
        }
    }
}

template<class R>
void reefamr_sed_slide(R &amr, ghostcell *pgc, reefamr_sedgrids &G, int tag)
{
    lexer *p = G.L[0];
    const int n = (int)G.L.size();

    auto each = [&](const std::function<void(int)> &fn)
    {
        for(int t=0; t<n; ++t)
        G.scope(t-1,[&]() { fn(t); });
    };
    auto cell = [](lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); };

    // the bed (after the Exner step) and the angle of repose around the grids
    amr.xslice(1,tag,[&](int g, int k) -> slice* { return &G.S[g+1]->s->bedzh; });
    amr.xslice(1,tag,[&](int g, int k) -> slice* { return &G.S[g+1]->s->phi; });

    // the open cells: the cells of the grids with their DFBED; around the grids (and in the halo and the ring of
    // the native grid) that of the cell they stand for; closed: covered cells, a source on a coarser level
    vector<std::unique_ptr<slice4>> D(n);
    vector<vector<int>> keep(n);

    for(int t=0; t<n; ++t)
    {
        lexer *q = G.L[t];
        D[t].reset(new slice4(q));
        for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
        for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
        (*D[t])(ii,jj) = ((*G.mask[t])(ii,jj)>0.5) ? double(q->DFBED[cell(q,ii,jj)]) : -1.0;
    }

    amr.xslice(1,tag,[&](int g, int k) -> slice* { return D[g+1].get(); });

    for(int t=0; t<n; ++t)
    {
        lexer *q = G.L[t];
        for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
        for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
        {
            const int c = cell(q,ii,jj);
            const int v = int((*D[t])(ii,jj));
            if(q->DFBED[c]!=v)
            {
            keep[t].push_back(c);
            keep[t].push_back(q->DFBED[c]);
            q->DFBED[c] = v;
            }
        }
    }

    // the iterations of the slide: every grid, the material slid into the cells around the grids to the cells
    // they stand for, the bed, the bed around the grids
    each([&](int t) { G.S[t]->pslide->amr_begin(G.L[t],G.S[t]->s); });

    for(int qn=0; qn<p->S91; ++qn)
    {
        int count = 0;

        each([&](int t)
        {
            sandslide *sl = G.S[t]->pslide;
            sl->amr_donor = G.mask[t];
            count += sl->amr_sweep(G.L[t],pgc,G.S[t]->s);
            sl->amr_donor = nullptr;
        });

        amr.xslice_add(tag,[&](int g) -> slice* { return G.S[g+1]->pslide->amr_fh(); });

        each([&](int t) { G.S[t]->pslide->amr_apply(G.L[t],G.S[t]->s); });

        amr.xslice(1,tag,[&](int g, int k) -> slice* { return &G.S[g+1]->s->bedzh; });

        count = pgc->globalimax(count);
        p->slidecells = count;

        if(count==0)
        break;

        if(p->mpirank==0)
        std::cout<<"sandslide corrections (all grids): "<<count<<std::endl;
    }

    for(int t=0; t<n; ++t)
    for(size_t k=0; k<keep[t].size(); k+=2)
    G.L[t]->DFBED[keep[t][k]] = keep[t][k+1];
}

#endif
