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

#include"nhflow_geometry.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#define WLVL (fabs(d->WL(i,j))>0.00005?d->WL(i,j):1.0e20)

void nhflow_geometry::ray_cast(lexer *p, fdm_nhf *d, ghostcell *pgc, double *LS)
{
    zmin = 1.0e8;
    zmax = -1.0e8;
    

    LOOP
    WETDRY
    {
    zmin = MIN(zmin, p->ZSP[IJK]);
    zmax = MAX(zmax, p->ZSP[IJK]);
    }
    
    LOOP
	{
    IO[IJK]=1;
	LS[IJK]=1.0e8;
	}
    
    	
    // geometry core kernels (geo_raycast)
    for(int rayiter=0; rayiter<2; ++rayiter)
    {
        for(int qn=0;qn<entity_sum;++qn)
        {
            if(rayiter==0)
            {
            const int mode = (qn<int(ent_raymode.size())) ? ent_raymode[qn] : 1;
            
            georay.sigma_io(p,tri_x,tri_y,tri_z,tstart[qn],tend[qn],DSM,mode,IO,CL,CR);
            
            if(qn<int(ent_invert.size()) && ent_invert[qn]==1)
            LOOP
            IO[IJK] = -IO[IJK];
            }
            
            if(rayiter==1)
            {
            pgc->startintV(p,IO,1);
            
            georay.sigma_band(p,tri_x,tri_y,tri_z,tstart[qn],tend[qn],NB*DSM,LS);
            }
        }
    }
    
    LOOP
    WETDRY
    {
        if(IO[IJK]==-1)
        LS[IJK]=-fabs(LS[IJK]);
        
        
        if(IO[IJK]==1)
        LS[IJK]=fabs(LS[IJK]);
    }
	
	LOOP
    WETDRY
	{
		if(LS[IJK]>100.0*DSM)
		LS[IJK]=100.0*DSM;
		
		if(LS[IJK]<-100.0*DSM)
		LS[IJK]=-100.0*DSM;
	}
    
    LOOP
    if(p->wet[IJ]==0)
    LS[IJK]=100.0*DSM;
    
    
	pgc->start5V(p,LS,1); 
}