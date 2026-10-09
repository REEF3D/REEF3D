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
#include"lexer.h"
#include"field.h"
#include"increment.h"

// half width of the smoothed interface (Heaviside / delta) of the level set f in cell (i,j,k):
// F 45 times the local cell size in the direction of the interface normal,
// sqrt(nx^2 DX^2 + ny^2 DY^2 + nz^2 DZ^2) with n = grad f / |grad f|.
// On a uniform grid this is F 45 * dx for every orientation; for a free surface on a grid
// stretched in x it is F 45 * dz. 2D (j_dir = 0) without the y-direction.
inline double interface_width(lexer *p, field &f, int i, int j, int k)
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

// width on the face between cell (i,j,k) and its neighbour (i+aa,j+bb,k+cc): mean of the two cells
// (symmetric, so the face density is the same seen from both cells)
inline double interface_width_face(lexer *p, field &f, int i, int j, int k, int aa, int bb, int cc)
{
    return 0.5*(interface_width(p,f,i,j,k) + interface_width(p,f,i+aa,j+bb,k+cc));
}

#endif
