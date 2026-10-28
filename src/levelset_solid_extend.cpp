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

#include"levelset_solid_extend.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"increment.h"
#include<algorithm>
#include<cmath>
#include<vector>

//  The level set inside the solids of the direct forcing.
//
//  With the direct forcing (G 505 0) the cells of a solid are fluid cells: their velocity is forced
//  to zero, and the level set there is not moved.  It keeps the free surface the start (F 60, or
//  the hot start I 30 from FNPF) put through the solid.  A pile through the free surface then holds
//  a water column of its own.  As soon as the free surface outside moves past it (waves, or the
//  surface of the hot start relaxing), the level set has a step at the wall.  The pressure of that
//  column against the water outside drives a jet along the wall through the forcing band, and the
//  air just above the free surface at the wall takes velocities of the order of 1 m/s (a pier with
//  a hot start from FNPF: 0.7 m/s down along the piles where the waves had 0.07 m/s, validation 55).
//
//  Here the level set of the cells inside a solid is extrapolated from the fluid, the way the free
//  surface at rest lies: constant across (x, y), linear along the gravity (z).  The cells of the
//  solid are put in layers by their distance in cells from the fluid (layer 1: a face next to a
//  fluid cell, layer 2: next to layer 1, ...); every cell takes the mean over the directions in
//  which its neighbour lies in a lower layer (towards the fluid): in x and y the value of that
//  neighbour, in z 2 f(+1) - f(+2), the change limited to the height of the cell (|grad f| = 1 of a
//  distance).  The layers come from the sign of the solid level set alone: inside some solids (S 61
//  a wedge) it is not a distance.
//  A free surface at rest (f = h - z) is continued exactly through every wall: at a pile the free
//  surface meets the wall at a right angle and the pile holds the level of the water around it, at a
//  gently sloping bed (a beach) the still water level stays a plane through the bed.  Constant along
//  the normal of the solid would bend the free surface along the shore line of a beach and drive a
//  flow there (0.75 m/s in still water, the case of cfd_2d_beach_still); linear across as well amplifies the slope of
//  the free surface at the wall of a pile a few cells wide into a column of its own (validation 55).
//  The layers are done in order, so that the neighbours towards the fluid are done first: one sweep,
//  a second one after the ghost cells for the solids across a partition edge.  Only solid_band
//  layers: deeper the level set does not reach the fluid.  The solid then holds no water of its own.
//
//  X 48 2 (the default); X 48 0 keeps the level set inside the solids as it is, X 48 1 is the copy
//  into the ghost cells of the flagged solids (G 505 1, ghostcell::solid_forcing_lsm).

namespace
{
const int solid_band = 8;   // layers of cells

// the neighbour after the next one (l2) towards the fluid as well, from the next (l1)
inline bool towards(unsigned char l1, unsigned char l2)
{
    return l2<l1 || l2==0;
}

// linear from the two neighbours towards the fluid (f1 the next at the distance h1, f2 the one after,
// h2 from f1), the change limited to the distance h1
inline double extrapolate(double f1, double f2, double h1, double h2)
{
    return f1 + std::max(-h1, std::min(h1, (f1-f2)*h1/h2));
}
}

void levelset_solid_extend(lexer *p, fdm *a, ghostcell *pgc, field &f, int gcval)
{
    // with a free surface (F 30), solids of the direct forcing (G 505 0)
    if(p->X48!=2 || p->F30==0 || p->solidread==0 || p->G505!=0)
    return;

    int i,j,k;
    const int marge = increment::marge;
    const int jd = p->j_dir==1 ? 1 : 0;
    const unsigned char outside = 255;

    // the layers, on the cells of the partition and its ghost cells across a partition edge or a
    // periodic edge (the fluid there within the ghost cells); fluid 0.  The ghost cells beyond the
    // edges of the domain are left out (outside): the level set there is the boundary condition.
    const bool xm = p->nb1>=0 || p->periodic1==1, xp = p->nb4>=0 || p->periodic1==1;
    const bool ym = p->nb3>=0 || p->periodic2==1, yp = p->nb2>=0 || p->periodic2==1;
    const bool zm = p->nb5>=0 || p->periodic3==1, zp = p->nb6>=0 || p->periodic3==1;
    const int n = p->imax*p->jmax*p->kmax;
    std::vector<unsigned char> layer(n,outside);
    std::vector<int> solid;
    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    for(k=p->kmin; k<p->kmin+p->kmax; ++k)
    {
        if((i<0 && !xm) || (i>=p->knox && !xp) || (k<0 && !zm) || (k>=p->knoz && !zp))
        continue;
        if(jd && ((j<0 && !ym) || (j>=p->knoy && !yp)))
        continue;

        if(a->solid(i,j,k)>=0.0)
        layer[IJK]=0;
        else
        solid.push_back(IJK);
    }

    const int si = p->jmax*p->kmax, sj = p->kmax, sk = 1;
    std::vector<std::vector<int> > band(solid_band+1);
    for(int l=1; l<=solid_band; ++l)
    for(int q : solid)
    {
        if(layer[q]!=outside)
        continue;
        const int ii = q/si, jj = (q/sj)%p->jmax, kk = q%p->kmax;
        bool next = false;
        if(ii>0 && layer[q-si]==l-1) next = true;
        if(ii<p->imax-1 && layer[q+si]==l-1) next = true;
        if(jd && jj>0 && layer[q-sj]==l-1) next = true;
        if(jd && jj<p->jmax-1 && layer[q+sj]==l-1) next = true;
        if(kk>0 && layer[q-sk]==l-1) next = true;
        if(kk<p->kmax-1 && layer[q+sk]==l-1) next = true;
        if(next)
        {
        layer[q]=l;
        band[l].push_back(q);
        }
    }

    // (every rank takes part in the ghost cell exchange)
    for(int pass=0; pass<2; ++pass)
    {
        for(int l=1; l<=solid_band; ++l)
        for(int q : band[l])
        {
            i = q/si + p->imin;
            j = (q/sj)%p->jmax + p->jmin;
            k = q%p->kmax + p->kmin;

            // the cells of the partition (the ghost cells come from the exchange)
            if(i<0 || i>=p->knox || j<0 || j>=p->knoy || k<0 || k>=p->knoz)
            continue;

            double sum = 0.0;
            int count = 0;

            if(layer[q+si]<l)
            {
                sum += f(i+1,j,k);
                ++count;
            }
            if(layer[q-si]<l)
            {
                sum += f(i-1,j,k);
                ++count;
            }
            if(jd && layer[q+sj]<l)
            {
                sum += f(i,j+1,k);
                ++count;
            }
            if(jd && layer[q-sj]<l)
            {
                sum += f(i,j-1,k);
                ++count;
            }
            if(layer[q+sk]<l)
            {
                sum += towards(layer[q+sk],layer[q+2*sk]) ? extrapolate(f(i,j,k+1),f(i,j,k+2),p->DZP[KP],p->DZP[KP1]) : f(i,j,k+1);
                ++count;
            }
            if(layer[q-sk]<l)
            {
                sum += towards(layer[q-sk],layer[q-2*sk]) ? extrapolate(f(i,j,k-1),f(i,j,k-2),p->DZP[KM1],p->DZP[KM2]) : f(i,j,k-1);
                ++count;
            }

            if(count>0)
            f(i,j,k) = sum/double(count);
        }

        pgc->start4(p,f,gcval);
    }
}
