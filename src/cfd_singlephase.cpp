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

#include"cfd_singlephase.h"
#include"cell_index.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"pic_gfm.h"
#include<vector>
#include<cmath>
#include<algorithm>

cfd_singlephase::cfd_singlephase(lexer *p) : km1(p), km2(p), km3(p), pid(p,p->F41==2), pa(p,p->F41==2), lab(p,p->F41==2), bq(p)
{
    fmarg = p->margin;

    // F 41 2: the pressure of the air pockets at the surface (pic_gfm_grad, presscorr)
    if(p->F41==2)
    {
        const int size = p->imax*p->jmax*p->kmax;
        for(int n=0; n<size; ++n)
        {
            pid.V[n] = -2.0;
            pa.V[n] = 0.0;
        }
        pic_gfm_air_field(p,&pa);
    }
}

cfd_singlephase::~cfd_singlephase()
{
    pic_gfm_air_field(nullptr,&pa);
}

// after a projection: the faces next to a cell in the water keep the velocity of the projection (the
// ghost fluid method set them with p = 0 at the surface); the faces between two cells in the air take,
// layer by layer (3 layers), the mean of the neighbour faces of the same component that are known
// (in the water or extended in an earlier layer), and 0 beyond the layers. The same extension as
// pic_f::air_velocity for F 81 1 (the mean, not linear: a linear one grew at the tips of thin layers).
// After each layer the values and the marks of the ghost cells come from the neighbours, so the
// extension does not depend on the decomposition.
void cfd_singlephase::air_velocity(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
    // the patches of the mesh refinement: cfd_amr::single_extend over all patches of a level
    if(p->amr_patch!=0)
    return;

    extend(p,a,pgc,u,0,3);
    extend(p,a,pgc,v,1,3);
    extend(p,a,pgc,w,2,3);

    if(p->F67>0)
    band(p,a,pgc,u,v,w);
}

// F 67: the band of the first F 67 layers of fluid cells in the air next to the water (the cells of
// the ghost fluid method, pic_gfm_wet41_level). In them the potential q with div(u - grad q) = 0, the
// faces between two cells of the band or a cell of the band and one in the air beyond it open, q = 0
// in the air beyond the band; the faces to the water (set by the projection), to a solid and to the
// outside closed, unchanged. Red-black SOR to 1e-7 of the largest velocity in the band (at most 400
// sweeps); level 0 exchanges q between the ranks after every colour and takes the residual over all
// ranks, so that the result does not depend on the decomposition. On a patch of the mesh refinement
// the cells around the patch are closed (each patch alone): the faces of its edges are not changed.
// Then u -= grad q on the open faces.
void cfd_singlephase::band(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
    const int marge = p->margin;
    const int size = p->imax*p->jmax*p->kmax;
    const bool y = p->j_dir==1;
    const bool patch = p->amr_patch!=0;
    const int nl = p->F67;
    const double wl = pic_gfm_wet41_level(p);
    const int sj = p->kmax, si = p->jmax*p->kmax;

    // 0 water, 1 .. nl the layers of the band, nl+1 the air beyond, -1 closed (no fluid cell, outside).
    // Level 0: in the margin as well (the layers of the cells of the neighbours), from the level set
    // there. A patch: the first cells around it as water or as air beyond (q = 0 there, the faces of the
    // patch edges not changed)
    const int lo = patch ? 0 : -(fmarg-1);
    const int ilo = lo, ihi = p->knox-lo, jlo = y ? lo : 0, jhi = y ? p->knoy-lo : 1, klo = lo, khi = p->knoz-lo;
    std::vector<signed char> lay(size,-1);
    for(int i=ilo; i<ihi; ++i)
    for(int j=jlo; j<jhi; ++j)
    for(int k=klo; k<khi; ++k)
    {
        const int q = cell_index(p,i,j,k);
        if(p->flag4[q]>0)
        lay[q] = a->phi(i,j,k)>=wl ? 0 : nl+1;
    }

    const int nb[6][3] = {{-1,0,0},{1,0,0},{0,0,-1},{0,0,1},{0,-1,0},{0,1,0}};
    const int nn = y ? 6 : 4;
    if(patch)
    for(int i=-1; i<=p->knox; ++i)
    for(int j=(y ? -1 : 0); j<=(y ? p->knoy : 0); ++j)
    for(int k=-1; k<=p->knoz; ++k)
    {
        if(i>=0 && i<p->knox && j>=0 && j<p->knoy && k>=0 && k<p->knoz)
        continue;
        const int q = cell_index(p,i,j,k);
        if(p->flag4[q]>0)
        lay[q] = a->phi(i,j,k)>=wl ? 0 : nl+1;
    }

    for(int L=1; L<=nl; ++L)
    {
        std::vector<int> set;
        for(int i=ilo; i<ihi; ++i)
        for(int j=jlo; j<jhi; ++j)
        for(int k=klo; k<khi; ++k)
        {
            const int q = cell_index(p,i,j,k);
            if(lay[q]!=nl+1)
            continue;
            for(int n=0; n<nn; ++n)
            {
                const int ii = i+nb[n][0], jj = j+nb[n][1], kk = k+nb[n][2];
                const int e = patch ? 1 : 0;   // (a patch: the first cells around it as well)
                if(ii<ilo-e || ii>=ihi+e || jj<jlo-(y ? e : 0) || jj>=jhi+(y ? e : 0) || kk<klo-e || kk>=khi+e)
                continue;
                if(lay[cell_index(p,ii,jj,kk)]==L-1)
                {
                    set.push_back(q);
                    break;
                }
            }
        }
        for(int q : set)
        lay[q] = L;
    }

    // a cell of the band takes part only when the band connects it to the air beyond (q = 0 there; on
    // level 0 also a cell of the band of a neighbouring rank): a pocket of air inside the band, closed by
    // water all around, would have no solution (-2: left out, its faces unchanged)
    {
        std::vector<char> con(size,0);
        std::vector<int> front;
        auto band_cell = [&](int q) { return lay[q]>=1 && lay[q]<=nl; };
        for(int i=0; i<p->knox; ++i)
        for(int j=0; j<p->knoy; ++j)
        for(int k=0; k<p->knoz; ++k)
        {
            const int q = cell_index(p,i,j,k);
            if(!band_cell(q))
            continue;
            for(int n=0; n<nn; ++n)
            {
                const int ii = i+nb[n][0], jj = j+nb[n][1], kk = k+nb[n][2];
                const int qn = cell_index(p,ii,jj,kk);
                const bool inner = ii>=0 && ii<p->knox && jj>=0 && jj<p->knoy && kk>=0 && kk<p->knoz;
                if(lay[qn]==nl+1 || (!patch && !inner && band_cell(qn)))
                {
                    con[q] = 1;
                    front.push_back(q);
                    break;
                }
            }
        }
        while(!front.empty())
        {
            std::vector<int> next;
            for(int q : front)
            {
                int i,j,k;
                i = q/si + p->imin;  j = (q%si)/sj + p->jmin;  k = q%sj + p->kmin;
                for(int n=0; n<nn; ++n)
                {
                    const int ii = i+nb[n][0], jj = j+nb[n][1], kk = k+nb[n][2];
                    if(ii<0 || ii>=p->knox || jj<0 || jj>=p->knoy || kk<0 || kk>=p->knoz)
                    continue;
                    const int qn = cell_index(p,ii,jj,kk);
                    if(band_cell(qn) && !con[qn])
                    {
                        con[qn] = 1;
                        next.push_back(qn);
                    }
                }
            }
            front.swap(next);
        }
        for(int i=0; i<p->knox; ++i)
        for(int j=0; j<p->knoy; ++j)
        for(int k=0; k<p->knoz; ++k)
        {
            const int q = cell_index(p,i,j,k);
            if(band_cell(q) && !con[q])
            lay[q] = -2;
        }
    }

    // the cells of the band of this grid (rank, patch), by colour
    std::vector<int> cells[2];
    double umax = 0.0;
    for(int i=0; i<p->knox; ++i)
    for(int j=0; j<p->knoy; ++j)
    for(int k=0; k<p->knoz; ++k)
    {
        const int q = cell_index(p,i,j,k);
        if(lay[q]<1 || lay[q]>nl)
        continue;
        // the colour by the global index, the same on every decomposition
        cells[(i+p->origin_i+j+p->origin_j+k+p->origin_k)&1].push_back(q);
        umax = std::max(umax,std::max(fabs(u.V[q]),fabs(w.V[q])));
        if(y)
        umax = std::max(umax,fabs(v.V[q]));
    }

    int any = (cells[0].size()+cells[1].size())>0 ? 1 : 0;
    if(!patch)
    {
        any = pgc->globalimax(any);
        umax = pgc->globalmax(umax);
    }
    if(any==0)
    return;

    for(int n=0; n<size; ++n)
    bq.V[n] = 0.0;

    auto ijk = [&](int q, int &i, int &j, int &k)
    {
        i = q/si + p->imin;
        j = (q%si)/sj + p->jmin;
        k = q%sj + p->kmin;
    };
    auto qv = [&](int qn) { return (lay[qn]>=1 && lay[qn]<=nl) ? bq.V[qn] : 0.0; };

    const double tol = 1.0e-7*umax + 1.0e-12;
    const double om = 1.6;
    bool conv = false;
    for(int it=0; it<400; ++it)
    {
        double res = 0.0;
        for(int col=0; col<2; ++col)
        {
            for(int q : cells[col])
            {
                int i,j,k;
                ijk(q,i,j,k);
                double diag = 0.0, off = 0.0;
                auto face = [&](int qn, double dxp, double dxn)
                {
                    if(lay[qn]<1)
                    return;
                    const double c = 1.0/(dxp*dxn);
                    diag += c;
                    off += c*qv(qn);
                };
                face(q+si,p->DXP[i+marge],p->DXN[i+marge]);
                face(q-si,p->DXP[i-1+marge],p->DXN[i+marge]);
                face(q+1,p->DZP[k+marge],p->DZN[k+marge]);
                face(q-1,p->DZP[k-1+marge],p->DZN[k+marge]);
                if(y)
                {
                    face(q+sj,p->DYP[j+marge],p->DYN[j+marge]);
                    face(q-sj,p->DYP[j-1+marge],p->DYN[j+marge]);
                }
                if(diag<=0.0)
                continue;

                double d = (u.V[q]-u.V[q-si])/p->DXN[i+marge] + (w.V[q]-w.V[q-1])/p->DZN[k+marge];
                if(y)
                d += (v.V[q]-v.V[q-sj])/p->DYN[j+marge];

                // off - diag q = div
                const double r = off - diag*bq.V[q] - d;
                res = std::max(res,fabs(r)*p->DXN[i+marge]);
                bq.V[q] += om*r/diag;
            }

            if(!patch)
            pgc->gcparax(p,bq,4);
        }

        if(!patch)
        res = pgc->globalmax(res);
        if(res<tol)
        {
            conv = true;
            break;
        }
    }

    // not converged (a part of the band still without a way to the air beyond, across ranks): the
    // extension as it is
    if(!conv)
    return;

    // u -= grad q on the faces between two cells in the air with at least one in the band. The band only
    // takes out the divergence of the extension (of the size of the velocities times the spacing over
    // their length scale). A correction larger than the largest velocity of the band would come from air
    // squeezed through a narrow part of the band (a closing gap, a sheet in a violent flow) and make it
    // fast: then the extension of this grid is kept as it is for this step
    const int ib = patch ? 0 : -1;
    for(int pass=0; pass<2; ++pass)
    {
        double dmax = 0.0;
        auto corr = [&](field &f, int q, int qn, double dxp)
        {
            const int l0 = lay[q], l1 = lay[qn];
            if(l0<1 || l1<1 || (l0>nl && l1>nl))
            return;
            const double d = (qv(qn) - qv(q))/dxp;
            if(pass==0)
            dmax = std::max(dmax,fabs(d));
            else
            f.V[q] -= d;
        };
        for(int i=ib; i<p->knox; ++i)
        for(int j=(y ? ib : 0); j<p->knoy; ++j)
        for(int k=ib; k<p->knoz; ++k)
        {
            const int q = cell_index(p,i,j,k);
            if(j>=0 && k>=0 && i+1<p->knox+(patch ? 0 : 1))
            corr(u,q,q+si,p->DXP[i+marge]);
            if(i>=0 && j>=0 && k+1<p->knoz+(patch ? 0 : 1))
            corr(w,q,q+1,p->DZP[k+marge]);
            if(y && i>=0 && k>=0 && j+1<p->knoy+(patch ? 0 : 1))
            corr(v,q,q+sj,p->DYP[j+marge]);
        }
        if(pass==0)
        {
            if(!patch)
            dmax = pgc->globalmax(dmax);
            if(dmax>umax)
            return;
        }
    }

    if(patch)
    return;

    pgc->start1(p,u,10);
    if(y)
    pgc->start2(p,v,11);
    pgc->start3(p,w,12);
}

void cfd_singlephase::extend(lexer *p, fdm *a, ghostcell *pgc, field &dst, int c, int layers)
{
    const int fm = fmarg;
    const int jlo = p->j_dir==1 ? -fm : 0;
    const int jhi = p->j_dir==1 ? p->knoy+fm-1 : 0;
    const int size = p->imax*p->jmax*p->kmax;

    if(c==1 && p->j_dir==0)
    {
        for(int n=0; n<size; ++n)
        dst.V[n] = 0.0;
        return;
    }

    const int di = c==0, dj = c==1, dk = c==2;

    // the known faces, with the margin (phi is known there): next to a cell in the water
    std::vector<char> known(size,0);

    for(int i=-fm; i<p->knox+fm-di; ++i)
    for(int j=jlo; j<=jhi-(p->j_dir==1 ? dj : 0); ++j)
    for(int k=-fm; k<p->knoz+fm-dk; ++k)
    {
        const double f0 = a->phi(i,j,k), f1 = a->phi(i+di,j+dj,k+dk);
        const double wl = pic_gfm_wet41_level(p);

        // next to a cell of the water (a fluid cell: the level set inside a solid is extended from the
        // water and says nothing about the face; F 66: the cells of the ghost fluid method, pic_gfm_water)
        if((f0>=wl && p->flag4[cell_index(p,i,j,k)]>0) || (f1>=wl && p->flag4[cell_index(p,i+di,j+dj,k+dk)]>0))
        known[cell_index(p,i,j,k)] = 1;
    }

    field &km = c==0 ? static_cast<field&>(km1) : (c==1 ? static_cast<field&>(km2) : static_cast<field&>(km3));

    for(int L=0; L<layers; ++L)
    {
        std::vector<char> next = known;

        for(int i=0; i<p->knox; ++i)
        for(int j=0; j<p->knoy; ++j)
        for(int k=0; k<p->knoz; ++k)
        {
            const int f = cell_index(p,i,j,k);

            if(known[f])
            continue;

            double sum=0.0;
            int cnt=0;

            auto add = [&](int si, int sj, int sk)
            {
                const int g1 = cell_index(p,i+si,j+sj,k+sk);

                if(!known[g1])
                return;

                sum += dst.V[g1];
                ++cnt;
            };

            add(-1,0,0); add(1,0,0);
            add(0,0,-1); add(0,0,1);

            if(p->j_dir==1)
            {
            add(0,-1,0); add(0,1,0);
            }

            if(cnt>0)
            {
            dst.V[f] = sum/double(cnt);
            next[f] = 1;
            }
        }

        known.swap(next);

        for(int n=0; n<size; ++n)
        km.V[n] = known[n] ? 1.0 : 0.0;

        pgc->gcparax(p,km,c+1);

        for(int n=0; n<size; ++n)
        known[n] = km.V[n]>0.5 ? 1 : 0;

        if(c==0)
        pgc->start1(p,dst,10);
        if(c==1)
        pgc->start2(p,dst,11);
        if(c==2)
        pgc->start3(p,dst,12);
    }

    // the air beyond the layers at rest
    for(int i=0; i<p->knox; ++i)
    for(int j=0; j<p->knoy; ++j)
    for(int k=0; k<p->knoz; ++k)
    if(!known[cell_index(p,i,j,k)])
    dst.V[cell_index(p,i,j,k)] = 0.0;

    if(c==0)
    pgc->start1(p,dst,10);
    if(c==1)
    pgc->start2(p,dst,11);
    if(c==2)
    pgc->start3(p,dst,12);
}

// the patches (cfd_amr::single_extend): the known faces of component c in its marks; the faces
// around the patch known
void cfd_singlephase::ext_begin(lexer *p, fdm *a, int c)
{
    field &km = marks(c);
    const int size = p->imax*p->jmax*p->kmax;
    for(int n=0; n<size; ++n)
    km.V[n] = 1.0;

    if(c==1 && p->j_dir==0)
    return;

    // the faces of the patch and those around it next to a cell in the water (a fluid cell); the
    // other faces around the patch between two fluid cells are known (the coarser grid, the
    // siblings), those at the walls and outside are not (their ghost values are boundary values)
    const int fm = fmarg;
    const int di = c==0, dj = c==1, dk = c==2;
    const int jlo = p->j_dir==1 ? -fm : 0;
    const int jhi = p->j_dir==1 ? p->knoy+fm-1 : 0;
    for(int i=-fm; i<p->knox+fm-di; ++i)
    for(int j=jlo; j<=jhi-(p->j_dir==1 ? dj : 0); ++j)
    for(int k=-fm; k<p->knoz+fm-dk; ++k)
    {
        const bool inner = i>=0 && i<p->knox && j>=0 && j<p->knoy && k>=0 && k<p->knoz;
        const int q0 = cell_index(p,i,j,k), q1 = cell_index(p,i+di,j+dj,k+dk);
        const double f0 = a->phi(i,j,k), f1 = a->phi(i+di,j+dj,k+dk);
        const double wl = pic_gfm_wet41_level(p);
        const bool water = (f0>=wl && p->flag4[q0]>0) || (f1>=wl && p->flag4[q1]>0);
        if(inner)
        km.V[q0] = water ? 1.0 : 0.0;
        else
        km.V[q0] = (water || (p->flag4[q0]>0 && p->flag4[q1]>0)) ? 1.0 : 0.0;
    }
}

// one layer: the unknown faces of the patch with known neighbours take their mean
void cfd_singlephase::ext_layer(lexer *p, field &dst, int c)
{
    if(c==1 && p->j_dir==0)
    return;

    field &km = marks(c);
    std::vector<std::pair<int,double>> set;

    for(int i=0; i<p->knox; ++i)
    for(int j=0; j<p->knoy; ++j)
    for(int k=0; k<p->knoz; ++k)
    {
        const int f = cell_index(p,i,j,k);
        if(km.V[f]>0.5)
        continue;

        double sum=0.0;
        int cnt=0;
        auto add = [&](int si, int sj, int sk)
        {
            const int g1 = cell_index(p,i+si,j+sj,k+sk);
            if(km.V[g1]<=0.5)
            return;
            sum += dst.V[g1];
            ++cnt;
        };
        add(-1,0,0); add(1,0,0);
        add(0,0,-1); add(0,0,1);
        if(p->j_dir==1)
        {
        add(0,-1,0); add(0,1,0);
        }

        if(cnt>0)
        set.push_back(std::make_pair(f,sum/double(cnt)));
    }

    for(const auto &q : set)
    {
        dst.V[q.first] = q.second;
        km.V[q.first] = 1.0;
    }
}

// the faces of the patch still unknown: at rest
void cfd_singlephase::ext_end(lexer *p, field &dst, int c)
{
    field &km = marks(c);
    for(int i=0; i<p->knox; ++i)
    for(int j=0; j<p->knoy; ++j)
    for(int k=0; k<p->knoz; ++k)
    if(c==1 && p->j_dir==0)
    dst(i,j,k) = 0.0;
    else if(km(i,j,k)<=0.5)
    dst(i,j,k) = 0.0;
}
