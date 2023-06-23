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
    cout<<"full NWT initialize "<<p->dt<<endl;
    
    // eta
	SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
        db = distbeach(p);

		d->eta(i,j) = wave_eta(p,pgc,xg,0.0);
    }
    
    SLICELOOP4
    d->WL(i,j) = MAX(0.0, d->eta(i,j) + p->wd - d->bed(i,j));

    FLOOP
    p->ZSN[FIJK] = p->ZN[KP]*d->WL(i,j) + d->bed(i,j);
    
    
    LOOP
    p->ZSP[IJK]  = p->ZP[KP]*d->WL(i,j) + d->bed(i,j);
	
	LOOP
    {
		xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p); 
        
        z=p->ZSP[IJK]-p->phimean;

		d->U[IJK] = wave_u(p,pgc,xg,yg,z);
        d->UH[IJK] = (d->eta(i,j)+d->depth(i,j))*wave_u(p,pgc,xg,yg,z);
	}	
	
	LOOP
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p); 
        
		z=p->ZSP[IJK]-p->phimean;
		
		d->V[IJK] = wave_v(p,pgc,xg,yg,z);
        d->VH[IJK] = (d->eta(i,j)+d->depth(i,j))*wave_v(p,pgc,xg,yg,z);
	}
	
	LOOP
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p); 
        
		z=p->ZSP[IJK]-p->phimean;
		
		d->W[IJK] = wave_w(p,pgc,xg,yg,z);
        d->WH[IJK] = (d->eta(i,j)+d->depth(i,j))*wave_w(p,pgc,xg,yg,z);
	}
	
    //if(p->I10==1 || p->I12==1)
	FLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p); 
		
		d->P[FIJK] = 0.0;//(wave_h(p,pgc,xg,0.0,0.0) - p->ZN[KP])*a->ro(i,j,k)*fabs(p->W22);
	}
	
	
    
}
