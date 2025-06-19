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
Authors: Arun Kamath, Hans Bihs
--------------------------------------------------------------------*/

#include"fnpf_force_ale.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include <math.h>

void fnpf_force_ale::force_ale_force(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{	
    double ztot=0; // check for strip total
    
    double uvel,vvel,wvel;

	Fx=Fy=0;
	
    for(k=0; k<p->knoz; ++k)
	{
        uvel = 0.5*(c->U[FIJK]+c->U[FIJKp1]);
        vvel = 0.5*(c->V[FIJK]+c->V[FIJKp1]);
        wvel = 0.5*(c->W[FIJK]+c->W[FIJKp1]);
        
        dudsig_= dudsig(p,c,pgc); 
        dvdsig_= dvdsig(p,c,pgc); 
        
        // Term 1 from eqn (9) of Pakozdi et al (2021) MS
        ax1= (uvel - un[k])/(p->dt);  
        ay1= (vvel - vn[k])/ (p->dt);
        
        // Term 2
        ax2 = uvel*(dudxi(p,c,pgc) + (dudsig_*p->sigx[FIJK]));
        ay2 = vvel*(dvdxi(p,c,pgc) + (dvdsig_*p->sigy[FIJK]));
        
        // Term 3
        ax3 = (wvel - (p->sig[FIJK]*dndt(p,c,pgc)))* dudsig_*p->sigz[IJ];
        ay3 = (wvel - (p->sig[FIJK]*dndt(p,c,pgc)))* dvdsig_*p->sigz[IJ];
      
        // Sum up acceleration
        ax = ax1 + ax2 + ax3;
        ay = ay1 + ay2 + ay3;
        
        // Force on current strip
        Fx1 = c->WL(i,j)*((cm*ax*p->W1*PI*rc*rc*p->DZN[KP]) + (cd*uvel*fabs(uvel)*0.5*p->W1*2.0*rc*p->DZN[KP]));
        Fy1 = c->WL(i,j)*((cm*ay*p->W1*PI*rc*rc*p->DZN[KP]) + (cd*vvel*fabs(vvel)*0.5*p->W1*2.0*rc*p->DZN[KP]));
        
        // Sum up forces
        Fx += Fx1;
        Fy += Fy1;
        ztot += p->DZN[KP]; // checking total dz=1
        
        // Storing current time step information for next time step gradient calculation
        un[k] = uvel; 
        vn[k] = vvel;
	}
	
	etan=c->eta(i,j);
}

double fnpf_force_ale::dndt(lexer *p, fdm_fnpf *c, ghostcell *pgc) 
{
    double dndt = (c->eta(i,j) - etan)/ p->dt;
		 
    return dndt;
}

double fnpf_force_ale::dudsig(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	
{
    double dudsig_ = 0;
    
    dudsig_ = (c->U[FIJKp1] - c->U[FIJK])/(p->DZN[KP]);

	return dudsig_;        
}

double fnpf_force_ale::dvdsig(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	
{
	double dvdsig_ = 0;

    dvdsig_ = (c->V[FIJKp1] - c->V[FIJK])/(p->DZN[KP]);

	return dvdsig_;        
}

double fnpf_force_ale::dudxi(lexer *p, fdm_fnpf *c, ghostcell *pgc) 
{
    return (0.5*(c->U[FIp1JK]+c->U[FIp1JKp1]) - 0.5*(c->U[FIm1JK]+c->U[FIm1JKp1]))/(p->DXN[IP1] + p->DXN[IM1]); 
}

double fnpf_force_ale::dvdxi(lexer *p, fdm_fnpf *c, ghostcell *pgc) 
{
    return (0.5*(c->V[FIJp1K]+c->V[FIJm1Kp1]) - 0.5*(c->V[FIJm1K]+c->V[FIJm1Kp1]))/(p->DYN[JP1] + p->DYN[JM1]); 
}





