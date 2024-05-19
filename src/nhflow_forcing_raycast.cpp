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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#define WLVL (fabs(d->WL(i,j))>0.00005?d->WL(i,j):1.0e20)

void nhflow_forcing::ray_cast(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    LOOP
    p->ZSP[IJK]  = p->ZP[KP]*d->WL(i,j) + d->bed(i,j);
    
    // sigz
    SLICELOOP4
    {
    if(p->wet[IJ]==0)
    p->sigz[IJ] = 0.0;
    
    if(p->wet[IJ]==1)
    p->sigz[IJ] = 1.0/WLVL;
    }
    
    zmin = 1.0e8;
    zmax = -1.0e8;
    
    LOOP
    {
    zmin = MIN(zmin, p->ZSP[IJK]);
    zmax = MAX(zmax, p->ZSP[IJK]);
    }
    
    //cout<<"ZMIN/ZMAX: "<<zmin<<"  "<<zmax<<endl;
    
    LOOP
	{
    IO[IJK]=1;
	d->SOLID[IJK]=1.0e8;
	}
    
    	
    for(int rayiter=0; rayiter<2; ++rayiter)
    {
        for(int qn=0;qn<entity_sum;++qn)
        {
            if(rayiter==0)
            ray_cast_io(p,d,pgc,tstart[qn],tend[qn]);

            if(rayiter==1)
            {
            pgc->gcparaxintV(p,IO,1);
            
            ray_cast_direct(p,d,pgc,tstart[qn],tend[qn]);
            }
        }
    }
    
    LOOP
    {
        if(IO[IJK]==-1)
        d->SOLID[IJK]=-fabs(d->SOLID[IJK]);
        
        
        if(IO[IJK]==1)
        d->SOLID[IJK]=fabs(d->SOLID[IJK]);
    }
	
	LOOP
	{
		if(d->SOLID[IJK]>100.0*p->DXM)
		d->SOLID[IJK]=100.0*p->DXM;
		
		if(d->SOLID[IJK]<-100.0*p->DXM)
		d->SOLID[IJK]=-100.0*p->DXM;
	}
    
    /*
    LOOP
    {
        if(IO[IJK]==-1)
        d->SOLID[IJK]=-1.0;
        
        
        if(IO[IJK]==1)
        d->SOLID[IJK]=1.0;
    }*/
    
	pgc->start5V(p,d->SOLID,1); 
}