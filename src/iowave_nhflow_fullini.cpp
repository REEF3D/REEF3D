/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"iowave.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::full_initialize_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
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

		d->eta(i,j) = wave_eta(p,pgc,xg,0.0);
    }
	
	LOOP
    {
		xg = xgen1(p);
        yg = ygen1(p);
		dg = distgen(p);
		db = distbeach(p); 
        
        z=p->ZSP[IJK]-p->phimean;

		d->U[IJK] = wave_u(p,pgc,xg,yg,z);
	}	
	
	LOOP
    {
        xg = xgen2(p);
        yg = ygen2(p);
		dg = distgen(p);
		db = distbeach(p); 
        
		z=p->ZSP[IJK]-p->phimean;
		
		d->V[IJK] = wave_v(p,pgc,xg,yg,z);
	}
	
	LOOP
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p); 
        
		z=p->ZSP[IJK]-p->phimean;
		
		d->W[IJK] = wave_w(p,pgc,xg,yg,z);
	}
	
    if(p->I10==1 || p->I12==1)
	FLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p); 
		
		d->P[FIJK] = 0.0;//(wave_h(p,pgc,xg,0.0,0.0) - p->ZN[KP])*a->ro(i,j,k)*fabs(p->W22);
	}
	
	
    
}
