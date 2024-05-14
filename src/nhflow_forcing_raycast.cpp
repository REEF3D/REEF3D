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

void nhflow_forcing::ray_cast(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
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

            if(rayiter==1 && p->X188==2)
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
		if(d->SOLID[IJK]>10.0*p->DXM)
		d->SOLID[IJK]=10.0*p->DXM;
		
		if(d->SOLID[IJK]<-10.0*p->DXM)
		d->SOLID[IJK]=-10.0*p->DXM;
	}
    
	pgc->start5V(p,d->SOLID,1); 
}