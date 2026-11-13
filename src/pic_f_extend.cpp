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

#include"pic_f.h"
#include"cell_index.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"pic_gfm.h"
#include<vector>

// The velocity of the water for the particles (F 81 1, 2; Kelly et al. 2015, section 3.3):
// a velocity face is in the water when phi at the face (the mean of its two cells) is >= 0;
// the faces in the air get, layer by layer (3 layers), the mean over the directions in which the
// neighbour face of the same component is in the water or got a value in an earlier layer: the
// linear extrapolation from that face and the next one, if that has a value too, else the
// neighbour's value (a constant extension moved the surface of a rotating disk by a cell per turn). The particles
// close to the surface then move with the water, not with the light air of the two-phase
// velocity, whose spurious velocities at the interface would otherwise move the surface.
// The layers run over the faces of the subdomain; after each the values (gcval 10, 11, 12) and the
// marks of the faces of the ghost cells come from the neighbours, so the layers do not depend on the
// decomposition (marked from phi in the ghost cells alone, the third layer at a subdomain boundary
// lacked the cells beyond the margin).

void pic_f::water_velocity(lexer *p, fdm *a, ghostcell *pgc)
{
    // single phase: the faces next to the water from the projection, the others as in air_velocity
    extend(p,a,pgc,a->u,*uw,0,3,single);
    extend(p,a,pgc,a->v,*vw,1,3,single);
    extend(p,a,pgc,a->w,*ww,2,3,single);
}

// F 81 1 (single phase): after a projection the faces between two cells in the air take the
// velocity of the water extended into them (3 layers) and 0 beyond; the faces next to a cell in the
// water keep theirs, the projection (p = 0 at the surface, pic_gfm.cpp) set them. The extension is
// the mean of the neighbour faces, not linear: these velocities go back into the water through the
// particles (pic_f::step_begin), and a linear one grows at the tip of a thin layer step by step.
// 2D dam break: the front ran away at 20 m/s with it (with the wall faces as second faces); without
// those, and limited (minmod) or in the first layer only, still 4-9 m/s spikes in the splashes,
// against 3 m/s with the mean. Sloshing tank at 2 cells per wave height: 2.5 % decay per period
// with the mean, 2.0-2.2 % with the linear variants.
void pic_f::air_velocity(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
    if(!single)
    return;

    extend(p,a,pgc,u,u,0,3,true);
    extend(p,a,pgc,v,v,1,3,true);
    extend(p,a,pgc,w,w,2,3,true);
}

// dst = src on the known faces, extended layer by layer into the others (see above). zero (single
// phase): the known faces are those next to a cell in the water (phi >= 0 there), the extension is
// the mean of the neighbours and the rest beyond the layers 0; else the known faces are those with
// phi >= 0 at the face (the mean of its two cells), extended linearly.
void pic_f::extend(lexer *p, fdm *a, ghostcell *pgc, field &src, field &dst, int c, int layers, bool zero)
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

    if(&dst!=&src)
    for(int n=0; n<size; ++n)
    dst.V[n] = src.V[n];

    // the known faces, with the margin (phi is known there)
    std::vector<char> known(size,0);

    for(int i=-fm; i<p->knox+fm-di; ++i)
    for(int j=jlo; j<=jhi-(p->j_dir==1 ? dj : 0); ++j)
    for(int k=-fm; k<p->knoz+fm-dk; ++k)
    {
        const double f0 = a->phi(i,j,k), f1 = a->phi(i+di,j+dj,k+dk);

        // single phase: the cells in the water as for the pressure (pic_gfm_water)
        const double wet = PIC_GFM_WET*p->DXM;

        if(zero ? (f0>=wet || f1>=wet) : (0.5*(f0+f1) >= 0.0))
        known[cell_index(p,i,j,k)] = 1;
    }

    field &km = c==0 ? static_cast<field&>(*km1) : (c==1 ? static_cast<field&>(*km2) : static_cast<field&>(*km3));

    for(int L=0; L<layers; ++L)
    {
        std::vector<char> next = known;

        // the faces of this subdomain; those of the ghost cells, value and mark, from the neighbours
        for(int i=0; i<p->knox; ++i)
        for(int j=0; j<p->knoy; ++j)
        for(int k=0; k<p->knoz; ++k)
        {
            const int f = cell_index(p,i,j,k);

            if(known[f])
            continue;

            double sum=0.0;
            int cnt=0;

            // from the neighbour in a direction: linear from the two faces there, else its value
            auto add = [&](int si, int sj, int sk)
            {
                const int g1 = cell_index(p,i+si,j+sj,k+sk);

                if(!known[g1])
                return;

                const int i2 = i+2*si, j2 = j+2*sj, k2 = k+2*sk;
                const bool inside = i2>=-fm && i2<p->knox+fm && k2>=-fm && k2<p->knoz+fm
                                 && (p->j_dir==0 || (j2>=jlo && j2<=jhi));

                if(!zero && inside && known[cell_index(p,i2,j2,k2)])
                sum += 2.0*dst.V[g1] - dst.V[cell_index(p,i2,j2,k2)];
                else
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

    // single phase: the air beyond the layers at rest
    if(zero)
    {
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
}
