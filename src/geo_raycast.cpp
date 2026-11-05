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

#include"geo_raycast.h"
#include"lexer.h"

geo_raycast::geo_raycast(lexer *p) : epsi(1.6)
{
}

geo_raycast::~geo_raycast()
{
}

geo_cart geo_raycast::interior(lexer *p)
{
    geo_cart g;

    g.is = 0;
    g.ie = p->knox;
    g.js = 0;
    g.je = p->knoy;
    g.ks = 0;
    g.ke = p->knoz;

    g.interior = true;
    g.clip = true;
    g.raymode = 1;
    g.raystyle = 0;

    g.psi = 1.0e-8*p->DXM;
    g.ext = 10.0*p->DXM;

    return g;
}

geo_cart geo_raycast::extended(lexer *p, double dxm)
{
    geo_cart g;

    // ghost cells of the fields: imin = -margin
    g.is = p->imin;
    g.ie = p->knox - p->imin;
    g.js = p->jmin;
    g.je = p->knoy - p->jmin;
    g.ks = p->kmin;
    g.ke = p->knoz - p->kmin;

    g.interior = false;
    g.clip = false;
    g.raymode = 1;
    g.raystyle = 1;

    g.psi = 1.0e-8*dxm;
    g.ext = 10.0*dxm;

    return g;
}

int geo_raycast::cell_ext(const double *N, int lo, int hi, double x)
{
    // N[c+marge] is node c, cell c spans [N[c],N[c+1])
    const int m = increment::marge;

    if(x < N[lo+m])
    return lo-1;

    if(x >= N[hi+m])
    return hi;

    int a=lo, b=hi;

    while(b-a>1)
    {
        const int c = (a+b)/2;

        if(x < N[c+m])
        b=c;
        else
        a=c;
    }

    return a;
}

int geo_raycast::cell(lexer *p, int dir, double x, const geo_cart &g)
{
    if(g.interior)
    {
        if(dir==0)
        return p->posc_i(x);

        if(dir==1)
        return p->posc_j(x);

        return p->posc_k(x);
    }

    if(dir==0)
    return cell_ext(p->XN,g.is,g.ie,x);

    if(dir==1)
    return cell_ext(p->YN,g.js,g.je,x);

    return cell_ext(p->ZN,g.ks,g.ke,x);
}

bool geo_raycast::checkin(lexer *p, int dir, double Ax, double Ay, double Az,
                          double Bx, double By, double Bz, double Cx, double Cy, double Cz)
{
    // dir 0,2: y ignored in 2D, dir 1: all three directions
    const bool yfree = (dir!=1 && p->j_dir==0);

    bool in = false;

    if(Ax>=p->global_xmin && Ax<=p->global_xmax
    && ((Ay>=p->global_ymin && Ay<=p->global_ymax) || yfree)
    && Az>=p->global_zmin && Az<=p->global_zmax)
    in = true;

    if(Bx>=p->global_xmin && Bx<=p->global_xmax
    && ((By>=p->global_ymin && By<=p->global_ymax) || yfree)
    && Bz>=p->global_zmin && Bz<=p->global_zmax)
    in = true;

    if(Cx>=p->global_xmin && Cx<=p->global_xmax
    && ((Cy>=p->global_ymin && Cy<=p->global_ymax) || yfree)
    && Cz>=p->global_zmin && Cz<=p->global_zmax)
    in = true;

    return in;
}
