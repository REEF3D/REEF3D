/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
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

void vrans_f::omega_source(lexer *p, fdm *a, field &kin, field &eps)
{
	int count;
    double uvel,vvel,wvel,uu;
    double por;
    double kinf,winf;
    double ke_c_2e=1.92;
    
    count=0;
	if(p->B295==1)
    LOOP
    if(a->porosity(i,j,k)<1.0)
    {
        uvel = 0.5*(a->u(i,j,k)+a->u(i-1,j,k));
        vvel = 0.5*(a->v(i,j,k)+a->v(i,j-1,k));
        wvel = 0.5*(a->w(i,j,k)+a->w(i,j,k-1));
        
        uu = uvel*uvel + vvel*vvel + wvel*wvel;
        por = a->porosity(i,j,k);
        
        kinf = 3.7*(1.0-por)*pow(por,1.5)*uu;
        winf = 39.0*pow(1.0-por,2.5)*pow(por,2.0)*pow(uu,1.5)*(1.0/porpart(i,j,k))*(p->cmu*(kinf>1.0e-20?kinf:1.0e20));
        
        a->rhsvec.V[count] += por*(winf*winf);
        ++count;  
    }

}
