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


#ifndef SFLOW_DRYWALL_H_
#define SFLOW_DRYWALL_H_

#include"lexer.h"
#include"slice.h"
#include"increment.h"

// Dry land next to deep water (A 221 2): a quay, a breakwater, a cliff.
//
// A wet cell with a dry neighbour is taken out of the non-hydrostatic pressure (sflow_pressure_nh::active,
// sflow_momentum_SI::nh_active), and that hydrostatic cell is a p = 0 neighbour of the non-hydrostatic cells behind
// it. That is the shoreline of a beach: thin water, the pressure hydrostatic. In front of dry land with deep water
// the land took half of the wave height out of the cells in front of it (regular waves normal to a wall of dry cells,
// kh 1.3, A 220 2: 0.68 of the standing wave height in the cell at the land, reflection 0.92 instead of 0.99 of the
// wall at the end of the grid), and with the semi-implicit step and A 220 3 the run did not stay stable.
//
// With A 221 2 a cell deeper than its horizontal cell size keeps its non-hydrostatic pressure when all the dry cells
// around it are land its water cannot reach (ground above the water level of the cell), and a dry cell next to it is
// a wall: dq/dn = 0 in the pressure equation and in the pressure gradients, the normal velocity mirrored in the
// divergence, the depth mirrored in the bed slope terms. Thinner cells and the front of a bore on a dry beach (the dry
// cell ahead lies below the water level) are treated as before. One Cartesian grid (not curvilinear, multi-block,
// mesh refinement).

namespace sflow_drywall
{
inline bool on(lexer *p)
{
    return p->A221==2 && p->geo_curv==0 && p->geo_mb==0 && p->G1==0;
}

// a dry cell of the grid (not a solid, not outside the domain)
inline bool dry(lexer *p, int a, int c)
{
    const int q = (a-p->imin)*p->jmax + (c-p->jmin);
    return p->flagslice4[q]>0 && p->wet[q]==0;
}

// the cell (a,c) next to a non-hydrostatic cell is a wall
inline bool wall(lexer *p, int a, int c)
{
    return on(p) && dry(p,a,c);
}

// the wet cell (i,j) keeps its non-hydrostatic pressure next to dry cells: deeper than its cell size, and every cell
// around it (the diagonals too) that is not wet is dry land above its water level
inline bool cliff(lexer *p, slice &eta, slice &depth, int i, int j)
{
    if(!on(p))
    return false;

    const int marge = increment::marge;

    if(eta(i,j) + depth(i,j) <= MAX(p->DXN[IP], p->j_dir==1 ? p->DYN[JP] : 0.0))
    return false;

    // a cell next to it that is not wet: dry land below the water level of the cell; side: also a cell outside the
    // grid (as before), the corner cells outside the grid do not count (the ghost cells of a row of dry land)
    auto reach = [&](int a, int c, bool side)
    {
        const int q = (i+a-p->imin)*p->jmax + (j+c-p->jmin);

        if(p->wet[q]!=0)
        return false;

        if(p->flagslice4[q]<0)
        return side;

        return eta(i,j) + depth(i+a,j+c) >= 0.0;
    };

    if(reach(1,0,true) || reach(-1,0,true))
    return false;

    if(p->j_dir==1)
    if(reach(0,1,true) || reach(0,-1,true) || reach(1,1,false) || reach(1,-1,false) || reach(-1,1,false) || reach(-1,-1,false))
    return false;

    return true;
}
}

#endif
