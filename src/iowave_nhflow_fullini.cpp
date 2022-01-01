/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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
--------------------------------------------------------------------*/

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void iowave::full_initialize_nhflow(lexer *p, fdm*a, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"full NWT initialize"<<endl;
    
    // eta
	SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);

		a->eta(i,j) = wave_eta(p,pgc,xg,0.0);
    }
	
	ULOOP
    {
		xg = xgen1(p);
        yg = ygen1(p);
		dg = distgen(p);
		db = distbeach(p); 
        
        z=p->ZSP[IJK]-p->phimean;

		a->u(i,j,k) = wave_u(p,pgc,xg,yg,z);
	}	
	
	VLOOP
    {
        xg = xgen2(p);
        yg = ygen2(p);
		dg = distgen(p);
		db = distbeach(p); 
        
		z=p->ZSP[IJK]-p->phimean;
		
		a->v(i,j,k) = wave_v(p,pgc,xg,yg,z);
	}
	
	WLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p); 
        
		z=p->ZSN[(i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + (k+1)-p->kmin]-p->phimean;
		
		a->w(i,j,k) = wave_w(p,pgc,xg,yg,z);
	}
	
    if(p->I10==1 || p->I12==1)
	LOOP
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p); 
        
		z=p->ZSP[IJK]-p->phimean;
		
		a->press(i,j,k) = (wave_h(p,pgc,xg,0.0,0.0) - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
	}
	
	
    
}
