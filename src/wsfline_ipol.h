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

#ifndef WSFLINE_IPOL_H_
#define WSFLINE_IPOL_H_

#include"lexer.h"
#include"slice.h"
#include"increment.h"

//  The free surface lines (P 52 along x at a y, P 56 along y at an x): the slice value at the
//  coordinate across the line linear between the cell centres around it, as the gauges
//  (ccslipol4); the value of the cell that holds it next to a solid or a domain edge.  The cell
//  that holds a coordinate on a node lies dx/2 off it.

// along x (column i) at y; jn: the cell that holds y
inline double wsfline_at_y(lexer *p, slice &f, int i, int jn, double y)
{
    if(p->j_dir==0)
    return f(i,jn);
    const int a = (y>=p->YP[jn+increment::marge]) ? jn : jn-1;
    auto fl = [&](int jj) { return p->flagslice4[(i-p->imin)*p->jmax + (jj-p->jmin)]>0; };
    if(!fl(a) || !fl(a+1))
    return f(i,jn);
    const double w = (y - p->YP[a+increment::marge])/(p->YP[a+1+increment::marge] - p->YP[a+increment::marge]);
    return (1.0-w)*f(i,a) + w*f(i,a+1);
}

// along y (row j) at x; in: the cell that holds x
inline double wsfline_at_x(lexer *p, slice &f, int in, int j, double x)
{
    const int a = (x>=p->XP[in+increment::marge]) ? in : in-1;
    auto fl = [&](int ii) { return p->flagslice4[(ii-p->imin)*p->jmax + (j-p->jmin)]>0; };
    if(!fl(a) || !fl(a+1))
    return f(in,j);
    const double w = (x - p->XP[a+increment::marge])/(p->XP[a+1+increment::marge] - p->XP[a+increment::marge]);
    return (1.0-w)*f(a,j) + w*f(a+1,j);
}

#endif
