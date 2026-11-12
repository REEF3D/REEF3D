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
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"pic_gfm.h"
#include<vector>

namespace
{
inline int fid(lexer *p, int i, int j, int k)
{
    return (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin;
}
}

cfd_singlephase::cfd_singlephase(lexer *p) : km1(p), km2(p), km3(p), pid(p,p->F41==2), pa(p,p->F41==2), lab(p,p->F41==2)
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

        // next to a cell of the water (a fluid cell: the level set inside a solid is extended from the
        // water and says nothing about the face)
        if((f0>=0.0 && p->flag4[fid(p,i,j,k)]>0) || (f1>=0.0 && p->flag4[fid(p,i+di,j+dj,k+dk)]>0))
        known[fid(p,i,j,k)] = 1;
    }

    field &km = c==0 ? static_cast<field&>(km1) : (c==1 ? static_cast<field&>(km2) : static_cast<field&>(km3));

    for(int L=0; L<layers; ++L)
    {
        std::vector<char> next = known;

        for(int i=0; i<p->knox; ++i)
        for(int j=0; j<p->knoy; ++j)
        for(int k=0; k<p->knoz; ++k)
        {
            const int f = fid(p,i,j,k);

            if(known[f])
            continue;

            double sum=0.0;
            int cnt=0;

            auto add = [&](int si, int sj, int sk)
            {
                const int g1 = fid(p,i+si,j+sj,k+sk);

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
    if(!known[fid(p,i,j,k)])
    dst.V[fid(p,i,j,k)] = 0.0;

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
        const int q0 = fid(p,i,j,k), q1 = fid(p,i+di,j+dj,k+dk);
        const double f0 = a->phi(i,j,k), f1 = a->phi(i+di,j+dj,k+dk);
        const bool water = (f0>=0.0 && p->flag4[q0]>0) || (f1>=0.0 && p->flag4[q1]>0);
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
        const int f = fid(p,i,j,k);
        if(km.V[f]>0.5)
        continue;

        double sum=0.0;
        int cnt=0;
        auto add = [&](int si, int sj, int sk)
        {
            const int g1 = fid(p,i+si,j+sj,k+sk);
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
