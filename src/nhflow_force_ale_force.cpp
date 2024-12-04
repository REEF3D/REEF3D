/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"nhflow_force_ale.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include <math.h>

void nhflow_force_ale::force_ale_force(lexer* p, fdm_nhf *d, ghostcell *pgc)
{	
    double ztot=0; // check for strip total
    
    double uvel,vvel,wvel;

	Fx=Fy=0;
	
    for(k=0; k<p->knoz; ++k)
	{
        uvel = d->U[IJK];
        vvel = d->V[IJK];
        wvel = d->W[IJK];
        
        dudsig_= dudsig(p,d,pgc); 
        dvdsig_= dvdsig(p,d,pgc); 
        
        // Term 1 from eqn (9) of Pakozdi et al (2021) MS
        ax1= (uvel - un[k])/(p->dt);  
        ay1= (vvel - vn[k])/ (p->dt);
        
        // Term 2
        ax2 = uvel*(dudxi(p,d,pgc) + (dudsig_*p->sigx[FIJK]));
        ay2 = vvel*(dvdxi(p,d,pgc) + (dvdsig_*p->sigy[FIJK]));
        
        // Term 3
        ax3 = (wvel - (p->sig[FIJK]*dndt(p,d,pgc)))* dudsig_*p->sigz[IJ];
        ay3 = (wvel - (p->sig[FIJK]*dndt(p,d,pgc)))* dvdsig_*p->sigz[IJ];
      
        // Sum up acceleration
        ax = ax1 + ax2 + ax3;
        ay = ay1 + ay2 + ay3;
        
        // Force on current strip
        Fx1 = d->WL(i,j)*((cm*ax*p->W1*PI*rc*rc*p->DZN[KP]) + (cd*uvel*fabs(uvel)*0.5*p->W1*2.0*rc*p->DZN[KP]));
        Fy1 = d->WL(i,j)*((cm*ay*p->W1*PI*rc*rc*p->DZN[KP]) + (cd*vvel*fabs(vvel)*0.5*p->W1*2.0*rc*p->DZN[KP]));
        
        // Sum up forces
        Fx += Fx1;
        Fy += Fy1;
        ztot += p->DZN[KP]; // checking total dz=1
        
        // Storing current time step information for next time step gradient calculation
        un[k] = uvel; 
        vn[k] = vvel;
	}
	
	etan=d->eta(i,j);
}

double nhflow_force_ale::dndt(lexer *p, fdm_nhf *d, ghostcell *pgc) 
{
    double dndt = (d->eta(i,j) - etan)/ p->dt;
		 
    return dndt;
}

double nhflow_force_ale::dudsig(lexer *p, fdm_nhf *d, ghostcell *pgc) 	
{
    double dudsig_ = 0;
    
    dudsig_ = (d->U[IJKp1] - d->U[IJKm1])/(p->DZN[KP]+p->DZN[KM1]);

	return dudsig_;        
}

double nhflow_force_ale::dvdsig(lexer *p, fdm_nhf *d, ghostcell *pgc) 	
{
	double dvdsig_ = 0;

    dvdsig_ = (d->V[IJKp1] - d->V[IJKm1])/(p->DZN[KP]+p->DZN[KM1]);

	return dvdsig_;        
}

double nhflow_force_ale::dudxi(lexer *p, fdm_nhf *d, ghostcell *pgc) 
{
    return (d->U[Ip1JK] - d->U[Im1JK])/(p->DXN[IP1] + p->DXN[IM1]); 
}

double nhflow_force_ale::dvdxi(lexer *p, fdm_nhf *d, ghostcell *pgc) 
{
    return (d->V[IJp1K] - d->V[IJm1K])/(p->DYN[JP1] + p->DYN[JM1]); 
}





