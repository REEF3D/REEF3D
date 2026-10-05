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

#include"6DOF_obj_nhflow.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#define WLVL (fabs(d->WL(i,j))>0.00005?d->WL(i,j):1.0e20)

void sixdof_obj_nhflow::ray_cast(lexer *p, fdm_nhf *d, ghostcell *pgc)
{    
    zmin = 1.0e8;
    zmax = -1.0e8;
    
    
    LOOP
    WETDRY
    {
    zmin = MIN(zmin, p->ZSP[IJK]);
    zmax = MAX(zmax, p->ZSP[IJK]);
    }
    
    ray_cast_nhflow_grid(p,d,pgc,IO,CL,CR,DSM);
}

void sixdof_obj_nhflow::ray_cast_nhflow_grid(lexer *p, fdm_nhf *d, ghostcell *pgc, int *IO, int *CL, int *CR, double DSM)
{
    LOOP
	{
    IO[IJK]=1;
	d->FB[IJK]=1.0e8;
	}
    pgc->start5V(p,d->FB,1); 
    
    // geometry core kernels (geo_raycast)
    for(int rayiter=0; rayiter<2; ++rayiter)
    {
        for(int qn=0;qn<entity_sum;++qn)
        {
            if(rayiter==0)
            georay.sigma_io(p,tri_x,tri_y,tri_z,tstart[qn],tend[qn],DSM,1,IO,CL,CR);
            
            if(rayiter==1)
            {
            pgc->startintV(p,IO,1);
            
            georay.sigma_band(p,tri_x,tri_y,tri_z,tstart[qn],tend[qn],NB*DSM,d->FB);
            }
        }
    }
    
    int nband=0;

    // cell centres on the surface (distance at round-off level, e.g. a box face through a row of
    // cell centres) get FB = 0 exactly: the parity test perturbs the ray position in a fixed
    // direction, so such cells would come out inside on one side of the body and outside on the
    // other, and the sign-dependent forcing (X 15) would make a symmetric body asymmetric
    const double fb_eps = 1.0e-10*DSM;
    
    LOOP
    WETDRY
    {
        if(IO[IJK]==-1)
        d->FB[IJK]=-fabs(d->FB[IJK]);
        
        if(IO[IJK]==1)
        d->FB[IJK]=fabs(d->FB[IJK]);
        
        if(fabs(d->FB[IJK])<fb_eps)
        d->FB[IJK]=0.0;
        
        d->test[IJK] = IO[IJK];
    }
	
    LOOP
    WETDRY
	{
		if(d->FB[IJK] >  NB*DSM) d->FB[IJK] =  NB*DSM;
		if(d->FB[IJK] < -NB*DSM) d->FB[IJK] = -NB*DSM;
	}
    
    LOOP
    if(p->wet[IJ]==0)
    d->FB[IJK]=100.0*p->DXM;
    
        
	pgc->start5V(p,d->FB,1); 
}