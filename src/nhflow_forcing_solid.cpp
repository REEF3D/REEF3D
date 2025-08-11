/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_forcing::solid_forcing(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double alpha, double *U, double *V, double *W, slice &WL)
{

// update Heaviside
    pgc->start5V(p,d->FHB,1);
    
    uf=vf=wf=0.0;
    
     if(p->A521==0)
    LOOP
    {
        H = Hsolidface(p,d,0,0,0);
       
        FX[IJK] += H*(uf - U[IJK])/(alpha*p->dt);
        FY[IJK] += H*(vf - V[IJK])/(alpha*p->dt);
        FZ[IJK] += H*(wf - W[IJK])/(alpha*p->dt);
        
        d->FHB[IJK] = min(d->FHB[IJK] + H, 1.0); 
    }
    
     if(p->A521==1)
    LOOP
    {
        H = Hsolidface(p,d,0,0,0);
        
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
    // Normal vectors calculation 
		nx = -(d->SOLID[Ip1JK] - d->SOLID[Im1JK])/(p->DXN[IP] + p->DXN[IM1]);
		ny = -(d->SOLID[IJp1K] - d->SOLID[IJm1K])/(p->DYN[JP] + p->DYN[JM1]);
		nz = -(d->SOLID[IJKp1] - d->SOLID[IJKm1])/(p->DZN[KP]*WL(i,j) + p->DZN[KM1]*WL(i,j));

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;
        
        
        if(d->SOLID[IJK]<=0.0)
        {
        FX[IJK] += H*(uf - U[IJK])/(alpha*p->dt);
        FY[IJK] += H*(vf - V[IJK])/(alpha*p->dt);
        FZ[IJK] += H*(wf - W[IJK])/(alpha*p->dt);
        }

        if(d->SOLID[IJK]>0.0)
        {
        FX[IJK] += fabs(nx)*H*(uf - U[IJK])/(alpha*p->dt);
        FY[IJK] += fabs(ny)*H*(vf - V[IJK])/(alpha*p->dt);
        FZ[IJK] += fabs(nz)*H*(wf - W[IJK])/(alpha*p->dt);
        }
    
    }
    
    pgc->start5V(p,d->FHB,50);
}