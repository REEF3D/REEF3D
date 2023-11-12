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

#include"vrans_veg.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void vrans_veg::u_source(lexer *p, fdm *a)
{	
	// VRANS Vegetation porosity
    count=0;
    if(p->B310>0 || p->B321>0 || p->B322>0)
    ULOOP
	{
        Cd_val = 0.5*(Cd(i,j,k) + Cd(i+1,j,k));
        N_val = 0.5*(N(i,j,k) + N(i+1,j,k));
		D_val = 0.5*(D(i,j,k) + D(i+1,j,k));
        
        Fi = p->B309*0.25*PI*pow(D_val,2.0)*N_val*((a->u(i,j,k) - un(i,j,k))/p->dt);
        
        Fd = 0.5*Cd_val*N_val*D_val*a->u(i,j,k)*fabs(a->u(i,j,k));
        
    	
    a->rhsvec.V[count] += (-Fi -Fd);
	++count;
	}
}

void vrans_veg::v_source(lexer *p, fdm *a)
{
	// VRANS porosity
    count=0;
    if(p->B310>0 || p->B321>0 || p->B322>0)
    VLOOP
	{
        Cd_val = 0.5*(Cd(i,j,k) + Cd(i,j+1,k));
        N_val = 0.5*(N(i,j,k) + N(i,j+1,k));
		D_val = 0.5*(D(i,j,k) + D(i,j+1,k));
        
        Fi = p->B309*0.25*PI*pow(D_val,2.0)*N_val*((a->v(i,j,k) - vn(i,j,k))/p->dt);
        
        Fd = 0.5*Cd_val*N_val*D_val*a->v(i,j,k)*fabs(a->v(i,j,k));

	
    a->rhsvec.V[count] += (-Fi -Fd);
	++count;
	}
}

void vrans_veg::w_source(lexer *p, fdm *a)
{
	// VRANS porosity
    count=0;
    if(p->B310>0 || p->B321>0 || p->B322>0)
    WLOOP
	{
        Cd_val = 0.5*(Cd(i,j,k) + Cd(i,j,k+1));
        N_val = 0.5*(N(i,j,k) + N(i,j,k+1));
		D_val = 0.5*(D(i,j,k) + D(i,j,k+1));
        
        Fi = p->B309*0.25*PI*pow(D_val,2.0)*N_val*((a->w(i,j,k) - wn(i,j,k))/p->dt);
        
        Fd = 0.5*Cd_val*N_val*D_val*a->w(i,j,k)*fabs(a->w(i,j,k));

	
    a->rhsvec.V[count] -= (-Fi -Fd);
	++count;
	}
}

