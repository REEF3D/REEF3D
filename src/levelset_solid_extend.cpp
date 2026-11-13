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
//  a second one after the ghost cells for the solids across a partition edge.  The layers go through
//  the whole solid (up to 250 layers): its deep cells are not left to the reinitialisation.
//
//  The fluid below a solid does not take part: a cell next to the fluid only below it is not in
//  layer 1, and the fluid below a cell is not one of its directions.  The fluid under a solid roof (a
//  culvert in a lock wall) is enclosed water: its level set is the distance the reinitialisation
//  finds, to the zero crossing the extension put into the wall above it, not the free surface.
//  Extrapolated upwards into the wall it lowered that crossing, the reinitialisation then the level
//  set of the culvert, and so on: in a 3D culvert under a 0.36 m wall the level set at the culvert
//  roof fell from 0.08 to 0.0 in 4 s; in the full lock the wrong values reached the top of the walls
//  and the air above them (phi -0.04 instead of -0.28 at the wall top after 25 s), and from there
//  the air over the whole lock.  So the solid takes the level set of the free surface around it,
//  sideways and from above.  A solid that touches the fluid only from below (a roof over the whole
//  width at the top of the domain) takes it from below, in layers after the others.  (Until 2026-11
//  the fluid below took part and only 8 layers of level 0 were extended; X 48 4 was this as an
//  option, now the same as 3.)
//
//  X 48 2 (the default); X 48 0 keeps the level set inside the solids as it is, X 48 1 is the copy
//  into the ghost cells of the flagged solids (G 505 1, ghostcell::solid_forcing_lsm).
//
//  X 48 3: as 2, and inside the 6DOF bodies of the direct forcing as well (X 10, the level set a->fb
//  of all bodies at the pose of the stage).  A body between two water levels (a lock gate, its top
//  above both) holds the free surface of the start inside, a step at the face of the high side: the
//  same jet along the face, and the water level of the high side falls (the water goes into the air
//  above, the level set loses it), 0.12 m in 7 s at 0.32 m head (validation 26_cfd_ship_lock).  The
//  layers of the body follow the level set of their own side; with fewer than about four cells across
//  the body the two sides meet in the middle.

namespace
{
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
    // with a free surface (F 30): solids of the direct forcing (G 505 0; X 48 2, 3), the 6DOF bodies
    // (X 10; X 48 3) with the level set a->fb of all bodies
    const bool solids = (p->X48==2 || p->X48==3 || p->X48==4) && p->solidread>0 && p->G505==0;
    const bool bodies = (p->X48==3 || p->X48==4) && p->X10>0;
    
    if(p->F30==0 || (!solids && !bodies))
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

        if(!((solids && a->solid(i,j,k)<0.0) || (bodies && a->fb(i,j,k)<0.0)))
        layer[IJK]=0;
        else
        solid.push_back(IJK);
    }

    const int si = p->jmax*p->kmax, sj = p->kmax, sk = 1;

    // the depth of the layers: through the whole solid, as many layers as the marks hold (254, the
    // mark 255 is outside).  On a patch of the mesh refinement (cfd_amr) the solid cells at the edge of
    // the patch (filled from the coarser level) and those inside it then have the level set of the
    // same extension (8 cells of the patch against 8 of level 0 left a step in the density in a thick
    // wall at a patch edge, and the pressure correction blew the velocities up there)
    const int nband = 250;

    // layer 1: the solid cells next to the fluid; layer l: the solid cells next to layer l-1 (from
    // the cells of layer l-1, the same cells as a sweep over all solid cells per layer)
    std::vector<char> insolid(n,0);
    for(int q : solid)
    insolid[q]=1;

    std::vector<std::vector<int> > band(nband+1);
    for(int q : solid)
    {
        const int ii = q/si, jj = (q/sj)%p->jmax, kk = q%p->kmax;
        bool next = false;
        if(ii>0 && layer[q-si]==0) next = true;
        if(ii<p->imax-1 && layer[q+si]==0) next = true;
        if(jd && jj>0 && layer[q-sj]==0) next = true;
        if(jd && jj<p->jmax-1 && layer[q+sj]==0) next = true;
        if(kk<p->kmax-1 && layer[q+sk]==0) next = true;
        if(next)
        band[1].push_back(q);
    }
    for(int q : band[1])
    layer[q]=1;

    for(int l=2; l<=nband; ++l)
    {
        for(int r : band[l-1])
        {
            const int ii = r/si, jj = (r/sj)%p->jmax, kk = r%p->kmax;
            auto add = [&](int q)
            {
                if(insolid[q] && layer[q]==outside)
                {
                layer[q]=l;
                band[l].push_back(q);
                }
            };
            if(ii<p->imax-1) add(r+si);
            if(ii>0) add(r-si);
            if(jd && jj<p->jmax-1) add(r+sj);
            if(jd && jj>0) add(r-sj);
            if(kk<p->kmax-1) add(r+sk);
            if(kk>0) add(r-sk);
        }
        if(band[l].empty())
        break;
    }

    // a solid that touches the fluid only from below (a roof over the whole width at the top of the
    // domain): its cells next to that fluid start layers after the others, and they take the fluid
    // below them (belowok), so that no solid cell is left to the reinitialisation
    std::vector<char> belowok;
    {
        int lend = 1;
        while(lend<nband && !band[lend+1].empty())
        ++lend;

        std::vector<int> seed;
        for(int q : solid)
        if(layer[q]==outside)
        {
            const int kk = q%p->kmax;
            if(kk>0 && layer[q-sk]==0)
            seed.push_back(q);
        }

        if(!seed.empty() && lend+1<=nband)
        {
            belowok.assign(n,0);
            band[lend+1] = seed;
            for(int q : seed)
            {
            layer[q] = lend+1;
            belowok[q] = 1;
            }

            for(int l=lend+2; l<=nband; ++l)
            {
                for(int r : band[l-1])
                {
                    const int ii = r/si, jj = (r/sj)%p->jmax, kk = r%p->kmax;
                    auto add = [&](int q)
                    {
                        if(insolid[q] && layer[q]==outside)
                        {
                        layer[q]=l;
                        band[l].push_back(q);
                        }
                    };
                    if(ii<p->imax-1) add(r+si);
                    if(ii>0) add(r-si);
                    if(jd && jj<p->jmax-1) add(r+sj);
                    if(jd && jj>0) add(r-sj);
                    if(kk<p->kmax-1) add(r+sk);
                    if(kk>0) add(r-sk);
                }
                if(band[l].empty())
                break;
            }
        }
    }

    // (every rank takes part in the ghost cell exchange)
    for(int pass=0; pass<2; ++pass)
    {
        for(int l=1; l<=nband; ++l)
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
            if(layer[q-sk]<l && !(layer[q-sk]==0 && (belowok.empty() || !belowok[q])))
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
