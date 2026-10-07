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

#include"vrans_veg.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

// Vegetation turbulence sources (B 295 1), Lopez & Garcia (2001), J. Hydraul. Eng. 127(5):
//   the work of the vegetation drag F_i = 1/2 Cd a |u| u_i (a = N D, frontal area per volume) produces k:
//   S_k = C_fk 1/2 Cd a |u|^3, C_fk = 1;   S_eps = C_feps eps/k S_k, C_feps = 1.33;
//   S_omega = (C_feps - C_fk) omega/k S_k   (from omega = eps/(cmu k))
// Vegetation cells are those with a > 0 (porosity is only reduced with B 308 1).
double vrans_veg::veg_sk(lexer *p, fdm *a)
{
    const double uvel = 0.5*(a->u(i,j,k)+a->u(i-1,j,k));
    const double vvel = 0.5*(a->v(i,j,k)+a->v(i,j-1,k));
    const double wvel = 0.5*(a->w(i,j,k)+a->w(i,j,k-1));
    const double uabs = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);
    
    return 0.5*Cd(i,j,k)*N(i,j,k)*D(i,j,k)*uabs*uabs*uabs;
}

void vrans_veg::kw_source(lexer *p, fdm *a, field &kin)
{
    ke_source(p,a,kin);
}

void vrans_veg::ke_source(lexer *p, fdm *a, field &kin)
{
    int count=0;
    
	if(p->B295==1)
    LOOP
    {
    if(N(i,j,k)*D(i,j,k)>0.0)
    a->rhsvec.V[count] += veg_sk(p,a);
    
    ++count;
    }
}
