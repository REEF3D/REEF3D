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

#include"vrans_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

// Porous-media turbulence closure (Nakayama & Kuwahara 1999), B 295 1:
//   k_inf = 3.7 (1-n) n^1.5 |u|^2,   eps_inf = 39 (1-n)^2.5 n^2 |u|^3 / d50
// k-omega, k-equation: + n beta* k_inf omega_inf = + n eps_inf   (omega_inf = eps_inf/(cmu k_inf))
void vrans_f::kw_source(lexer *p, fdm *a, field &kin)
{
    int count;
    double uvel,vvel,wvel,uu;
    double por;
    double einf;
    
    count=0;
	if(p->B295==1)
    LOOP
    {
    if(a->porosity(i,j,k)<1.0)
    {
        uvel = 0.5*(a->u(i,j,k)+a->u(i-1,j,k));
        vvel = 0.5*(a->v(i,j,k)+a->v(i,j-1,k));
        wvel = 0.5*(a->w(i,j,k)+a->w(i,j,k-1));
        
        uu = uvel*uvel + vvel*vvel + wvel*wvel;
        por = a->porosity(i,j,k);
        
        einf = 39.0*pow(1.0-por,2.5)*pow(por,2.0)*pow(uu,1.5)*(1.0/a->porpart(i,j,k));
        
        a->rhsvec.V[count] += por*MAX(einf,0.0);
    }
    ++count;
    }
}

// k-epsilon, k-equation: + n eps_inf
void vrans_f::ke_source(lexer *p, fdm *a, field &kin)
{
    int count;
    double uvel,vvel,wvel,uu;
    double por;
    double kinf,einf;
    double ke_c_2e=1.92;
    
    count=0;
	if(p->B295==1)
    LOOP
    {
    if(a->porosity(i,j,k)<1.0)
    {
        uvel = 0.5*(a->u(i,j,k)+a->u(i-1,j,k));
        vvel = 0.5*(a->v(i,j,k)+a->v(i,j-1,k));
        wvel = 0.5*(a->w(i,j,k)+a->w(i,j,k-1));
        
        uu = uvel*uvel + vvel*vvel + wvel*wvel;
        por = a->porosity(i,j,k);
        

        kinf = 3.7*(1.0-por)*pow(por,1.5)*uu;
        einf = 39.0*pow(1.0-por,2.5)*pow(por,2.0)*pow(uu,1.5)*(1.0/a->porpart(i,j,k));
        
        a->rhsvec.V[count] += por*MAX(einf,0.0);
    
    }
    ++count;
    }

}
