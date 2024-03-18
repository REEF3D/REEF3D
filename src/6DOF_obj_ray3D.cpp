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

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fieldint.h"

void sixdof_obj::ray_cast(lexer *p, fdm *a, ghostcell *pgc)
{
	ALOOP
	{
    fbio(i,j,k)=1;
	a->fb(i,j,k)=1.0e9;
	}
	
    for(rayiter=0; rayiter<2; ++rayiter)
    {

        for(int qn=0;qn<entity_sum;++qn)
        {
            if(rayiter==0)
            {
            ray_cast_io_x(p,a,pgc,tstart[qn],tend[qn]);
            ray_cast_io_ycorr(p,a,pgc,tstart[qn],tend[qn]);
            ray_cast_io_zcorr(p,a,pgc,tstart[qn],tend[qn]);
            }
        
            if(rayiter==1 && p->X188==1)
            {
            pgc->gcparaxint(p,fbio,1);
            
            ray_cast_x(p,a,pgc,tstart[qn],tend[qn]);
            ray_cast_y(p,a,pgc,tstart[qn],tend[qn]);
            ray_cast_z(p,a,pgc,tstart[qn],tend[qn]);
            }
            
            if(rayiter==1 && p->X188==2)
            {
            pgc->gcparaxint(p,fbio,1);
            
            ray_cast_direct(p,a,pgc,tstart[qn],tend[qn]);
            }
        }
    }
    
    ALOOP
    {
        if(fbio(i,j,k)==-1)
        a->fb(i,j,k)=-fabs(a->fb(i,j,k));
        
        
        if(fbio(i,j,k)==1)
        a->fb(i,j,k)=fabs(a->fb(i,j,k));
    }
    
	
	ALOOP
	{
		if(a->fb(i,j,k)>10.0*p->DXM)
		a->fb(i,j,k)=10.0*p->DXM;
		
		if(a->fb(i,j,k)<-10.0*p->DXM)
		a->fb(i,j,k)=-10.0*p->DXM;
	}
    
	pgc->start4a(p,a->fb,50);

}





