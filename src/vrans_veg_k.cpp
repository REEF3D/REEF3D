/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"vrans_veg.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void vrans_veg::kw_source(lexer *p, fdm *a, field &kin)
{
    int count;
    double uvel,vvel,wvel,uu;
    double kw;
    double Ckp=1.0;
    
    count=0;
	if(p->B295==1)
    LOOP
    if(a->porosity(i,j,k)<1.0)
    {
        uvel = 0.5*(a->u(i,j,k)+a->u(i-1,j,k));
        vvel = 0.5*(a->v(i,j,k)+a->v(i,j-1,k));
        wvel = 0.5*(a->w(i,j,k)+a->w(i,j,k-1));
        
        uu = uvel*uvel + vvel*vvel + wvel*wvel;
        
        kw = Ckp*Cd(i,j,k)*D(i,j,k)*N(i,j,k)*sqrt(uu*kin(i,j,k));

        a->rhsvec.V[count] += kw;
    
        ++count;  
    }
}

void vrans_veg::ke_source(lexer *p, fdm *a, field &kin)
{
    int count;
    double uvel,vvel,wvel,uu;
    double kw;
    double Ckp=1.0;
    
    count=0;
	if(p->B295==1)
    LOOP
    if(a->porosity(i,j,k)<1.0)
    {
        uvel = 0.5*(a->u(i,j,k)+a->u(i-1,j,k));
        vvel = 0.5*(a->v(i,j,k)+a->v(i,j-1,k));
        wvel = 0.5*(a->w(i,j,k)+a->w(i,j,k-1));
        
        uu = uvel*uvel + vvel*vvel + wvel*wvel;
        
        kw = Ckp*Cd(i,j,k)*D(i,j,k)*N(i,j,k)*sqrt(uu*kin(i,j,k));

        a->rhsvec.V[count] += kw;
    
        ++count;  
    }
}
