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

#ifndef NHFLOW_DRYWALL_H_
#define NHFLOW_DRYWALL_H_

#include"lexer.h"
#include"fdm_nhf.h"
#include"increment.h"

// Dry land next to deep water (A 541 f, f > 0): a quay, a breakwater, a cliff.
//
// A column next to dry cells is taken out of the non-hydrostatic pressure (p->deep = 0, nhflow_fsf_f::wetdry) and its
// dry neighbours are p = 0 boundaries of the Poisson equation. That is the shoreline of a beach: thin water, the
// pressure hydrostatic. In front of dry land with deep water it made the land absorb the waves (regular waves, normal
// to a wall of dry cells: reflection 0.78 at kh 1.3 and 0.60 at kh 3.4 instead of 0.99 and 0.90 of the wall at the end
// of the grid).
//
// With A 541 f a column deeper than f times its horizontal cell size keeps its non-hydrostatic pressure (deep) when
// all the dry columns around it are land its water cannot reach (ground above the water level of the column, no
// solid column), and a dry column next to it is a wall: dp/dn = 0 in the Poisson equation and in the pressure
// gradients, the normal velocity mirrored in the divergence. Thinner columns, the front of a bore on a dry beach
// (the dry cell ahead lies below the water level) and the columns around solids (piers) are treated as before.

namespace nhflow_drywall
{
// the column keeps its non-hydrostatic pressure next to dry cells: deep enough, and every dry column in the band
// of nhflow_fsf_f::wetdry (2 cells in x, 3 in y, the diagonals) is land above its water level
inline bool cliff(lexer *p, fdm_nhf *d, int i, int j, double WL)
{
    const int marge = increment::marge;

    if(p->A541<=0.0 || WL <= p->A541*MAX(p->DXN[IP], p->j_dir==1 ? p->DYN[JP] : 0.0))
    return false;

    const double level = d->bed(i,j) + WL;

    auto reach = [&](int a, int b)
    {
        const int q = (i+a-p->imin)*p->jmax + (j+b-p->jmin);
        return p->wet[q]==0 && (p->flagslice4[q]<0 || d->bed(i+a,j+b) <= level);
    };

    if(reach(1,0) || reach(2,0) || reach(-1,0) || reach(-2,0))
    return false;

    if(p->j_dir==1)
    if(reach(0,1) || reach(0,2) || reach(0,3) || reach(0,-1) || reach(0,-2) || reach(0,-3)
    || reach(1,1) || reach(1,-1) || reach(-1,1) || reach(-1,-1))
    return false;

    return true;
}

// q: a column (slice index) next to a non-hydrostatic one; a dry column of the grid is a wall
inline bool wall(lexer *p, int q)
{
    return p->A541>0.0 && p->flagslice4[q]>0 && p->wet[q]==0;
}

// horizontal gradient of the pressure at the cell centre, second order, dp/dn = 0 at a wall
inline double gradx(lexer *p, double *P, int i, int j, int k)
{
    const int marge = increment::marge;
    auto m = [&](int q) { const int n = (i+q-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin; return 0.5*(P[n]+P[n+1]); };
    auto w = [&](int q) { return wall(p,(i+q-p->imin)*p->jmax + (j-p->jmin)); };

    return ((w(1) ? m(0) : m(1)) - (w(-1) ? m(0) : m(-1)))/(p->DXP[IP]+p->DXP[IM1]);
}

inline double grady(lexer *p, double *P, int i, int j, int k)
{
    const int marge = increment::marge;
    auto m = [&](int q) { const int n = (i-p->imin)*p->jmax*p->kmaxF + (j+q-p->jmin)*p->kmaxF + k-p->kmin; return 0.5*(P[n]+P[n+1]); };
    auto w = [&](int q) { return wall(p,(i-p->imin)*p->jmax + (j+q-p->jmin)); };

    return ((w(1) ? m(0) : m(1)) - (w(-1) ? m(0) : m(-1)))/(p->DYP[JP]+p->DYP[JM1]);
}

// the velocities of the neighbours in the divergence at the pressure node: mirrored at a wall (Uc, Vc: the cell itself)
inline void mirror(lexer *p, int i, int j, double Uc, double Vc, double &U1, double &U2, double &V1, double &V2)
{
    if(p->A541<=0.0)
    return;

    const int q = (i-p->imin)*p->jmax + (j-p->jmin);

    if(wall(p,q-p->jmax)) U1 = -Uc;
    if(wall(p,q+p->jmax)) U2 = -Uc;

    if(p->j_dir==1)
    {
    if(wall(p,q-1)) V1 = -Vc;
    if(wall(p,q+1)) V2 = -Vc;
    }
}
}

#endif
