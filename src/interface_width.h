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

#ifndef INTERFACE_WIDTH_H_
#define INTERFACE_WIDTH_H_

#include<cmath>
#include<algorithm>
#include"lexer.h"
#include"field.h"
#include"increment.h"
#include"heaviside.h"

// half width of the smoothed interface of the grid in cell (i,j,k) of the level set f:
// F 45 times the local cell size in the direction of the interface normal,
// sqrt(nx^2 DX^2 + ny^2 DY^2 + nz^2 DZ^2) with n = grad f / |grad f|.
// On a uniform grid this is F 45 * dx for every orientation; for a free surface on a grid
// stretched in x it is F 45 * dz. 2D (j_dir = 0) without the y-direction.
inline double interface_width_grid(lexer *p, field &f, int i, int j, int k)
{
    const int marge = increment::marge;

    double nx = (f(i+1,j,k)-f(i-1,j,k))/(p->DXP[IP]+p->DXP[IM1]);
    double ny = p->j_dir==1 ? (f(i,j+1,k)-f(i,j-1,k))/(p->DYP[JP]+p->DYP[JM1]) : 0.0;
    double nz = (f(i,j,k+1)-f(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]);

    double g2 = nx*nx + ny*ny + nz*nz;

    if(g2<1.0e-20)
    {
        if(p->j_dir==0)
        return p->F45*0.5*(p->DXN[IP] + p->DZN[KP]);

        return p->F45*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
    }

    double h2 = (nx*nx*p->DXN[IP]*p->DXN[IP] + ny*ny*p->DYN[JP]*p->DYN[JP] + nz*nz*p->DZN[KP]*p->DZN[KP])/g2;

    return p->F45*sqrt(h2);
}

// half width of the smoothed interface (Heaviside / delta) of the density and the viscosity in
// cell (i,j,k): that of the grid (above), or with G 28 (CFD mesh refinement) one width psi_rho on
// every level, so that a coarse and a fine cell at the same distance from the free surface have
// the same density
inline double interface_width(lexer *p, field &f, int i, int j, int k)
{
    if(p->psi_rho>0.0)
    return p->psi_rho;

    return interface_width_grid(p,f,i,j,k);
}

// the half width of the interface band of the grid p: psi (F 45 times the spacing of the grid), or
// the width psi_rho of G 28 (the band in which the density is smeared)
inline double interface_band(lexer *p)
{
    return p->psi_rho>0.0 ? p->psi_rho : p->psi;
}

// width on the face between cell (i,j,k) and its neighbour (i+aa,j+bb,k+cc): mean of the two cells
// (symmetric, so the face density is the same seen from both cells)
inline double interface_width_face(lexer *p, field &f, int i, int j, int k, int aa, int bb, int cc)
{
    return 0.5*(interface_width(p,f,i,j,k) + interface_width(p,f,i+aa,j+bb,k+cc));
}

// the smoothed step heaviside(x, interface_width_face(..)) of the value x at the face between cell
// (i,j,k) and (i+aa,j+bb,k+cc) (the mean of the level set there).  The width there is at most F 45
// times the largest cell size of the two cells (the mean of the widths of the cells, each a weighted
// mean of the cell sizes; with G 28 the width psi_rho), so away from the interface the step is 0 or
// 1 without computing the width: the same value, for a fraction of the cost (the densities of the
// faces, every assembly of the pressure matrix and every pressure gradient)
inline double heaviside_face(lexer *p, field &f, double x, int i, int j, int k, int aa, int bb, int cc)
{
    const int marge = increment::marge;
    double h = std::max(std::max(p->DXN[IP],p->DZN[KP]),std::max(p->DXN[IP+aa],p->DZN[KP+cc]));
    if(p->j_dir==1)
    h = std::max(h,std::max(p->DYN[JP],p->DYN[JP+bb]));
    if(p->psi_rho>0.0)
    h = p->psi_rho/p->F45;
    if(std::fabs(x) > p->F45*h*(1.0+1.0e-9))
    return x>0.0 ? 1.0 : 0.0;
    return heaviside(x,interface_width_face(p,f,i,j,k,aa,bb,cc));
}

#endif
