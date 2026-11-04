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
#include"lexer.h"
#include"field.h"
#include<algorithm>

// Trilinear interpolation for the particles on the non-uniform staggered grid, valid in the
// ghost cells of the subdomain too: the stages of the Runge-Kutta step and the particles
// close to a subdomain boundary sample the neighbour's values there, as a rank holding the
// whole domain would. (lexer::posc_i and the ccipol functions are meant for positions inside
// the subdomain: beyond it they clamp to the outermost cell.)
// Positions of the values: u(i,j,k) at (XN[i+1], YP[j], ZP[k]), v at (XP[i], YN[j+1], ZP[k]),
// w at (XP[i], YP[j], ZN[k+1]), cell values (type 4) at (XP[i], YP[j], ZP[k]); indices with
// the margin, i from -marge to knox+marge-1.

namespace
{
// the interval of x among the coordinates c(i) = A[i+base], i = lo ... hi: c(i) <= x < c(i+1),
// clamped to the ends; t the weight of c(i+1), within [0,1] (no extrapolation)
void bracket(const double *A, int base, int lo, int hi, double x, int &i, double &t)
{
    const double *first = A + lo + base;
    const double *last  = A + hi + base + 1;

    int q = int(std::upper_bound(first, last, x) - first) - 1 + lo;
    q = std::max(lo, std::min(q, hi-1));

    const double c0 = A[q+base];
    const double c1 = A[q+1+base];

    i = q;
    t = std::max(0.0, std::min(1.0, (x-c0)/(c1-c0)));
}
}

// the cell (with the margin) that holds x along an axis with the nodes N, n cells
int pic_f::cell(const double *N, int n, double x) const
{
    const double *first = N;
    const double *last  = N + n + 2*marge + 1;

    int q = int(std::upper_bound(first, last, x) - first) - 1 - marge;

    return std::max(-marge, std::min(q, n+marge-1));
}

double pic_f::ipol(lexer *p, field &f, int type, double xp, double yp, double zp) const
{
    int i0,j0=0,k0;
    double tx,ty=0.0,tz;

    // x: faces for u, centres else
    if(type==1)
    bracket(p->XN, 1+marge, -marge, p->knox+marge-2, xp, i0, tx);
    else
    bracket(p->XP, marge, -marge, p->knox+marge-1, xp, i0, tx);

    if(p->j_dir==1)
    {
    if(type==2)
    bracket(p->YN, 1+marge, -marge, p->knoy+marge-2, yp, j0, ty);
    else
    bracket(p->YP, marge, -marge, p->knoy+marge-1, yp, j0, ty);
    }

    if(type==3)
    bracket(p->ZN, 1+marge, -marge, p->knoz+marge-2, zp, k0, tz);
    else
    bracket(p->ZP, marge, -marge, p->knoz+marge-1, zp, k0, tz);

    const double wx[2] = {1.0-tx, tx};
    const double wy[2] = {1.0-ty, ty};
    const double wz[2] = {1.0-tz, tz};
    const int ny = p->j_dir==1 ? 2 : 1;

    double val=0.0;
    for(int a=0; a<2; ++a)
    for(int b=0; b<ny; ++b)
    for(int c=0; c<2; ++c)
    val += wx[a]*wy[b]*wz[c]*f(i0+a, j0+b, k0+c);

    return val;
}
