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

#include"fnpf_print_kinematics.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include <math.h>

void fnpf_print_kinematics::kinematics_calc(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{	
    for(k=0; k<p->knoz+1; ++k)
	{
        dudsig_= dudsig(p,c,pgc); 
        dvdsig_= dvdsig(p,c,pgc); 
        
        // Term 1 from eqn (9) of Pakozdi et al (2021) MS
        ax1= (c->U[FIJK] - un[k])/(p->dt);  
        ay1= (c->V[FIJK] - vn[k])/ (p->dt);
        
        // Term 2
        ax2 = c->U[FIJK]*(dudxi(p,c,pgc) + (dudsig_*p->sigx[FIJK]));
        ay2 = c->V[FIJK]*(dvdxi(p,c,pgc) + (dvdsig_*p->sigy[FIJK]));
        
        // Term 3
        ax3 = (c->W[FIJK] - (p->sig[FIJK]*dndt(p,c,pgc)))* dudsig_*p->sigz[IJ];
        ay3 = (c->W[FIJK] - (p->sig[FIJK]*dndt(p,c,pgc)))* dvdsig_*p->sigz[IJ];
      
        // Sum up acceleration
        ax[k] = ax1 + ax2 + ax3;
        ay[k] = ay1 + ay2 + ay3;
        
        
        // Storing current time step information for next time step gradient calculation
        un[k] = c->U[FIJK]; 
        vn[k] = c->V[FIJK];
	}
	
	etan=c->eta(i,j);
}

double fnpf_print_kinematics::dndt(lexer *p, fdm_fnpf *c, ghostcell *pgc) // to calculate dn dt for ax3
{
    double dndt = (c->eta(i,j) - etan)/ p->dt;
		 
    return dndt;
}

double fnpf_print_kinematics::dudsig(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	// getting dudsig for ax2 and 3
{
    double dudsig_ = 0;
	double dudsig2_ = 0;
    
    if(k<p->knoz)
    {
        dudsig_ =  (c->U[FIJKp1] - c->U[FIJKm1])/(p->DZN[KP1] + p->DZN[KM1]);
    }

    if(k==p->knoz)
    { 
        dudsig_ = (c->U[FIJK] - c->U[FIJKm1])/(p->DZN[KM1]);
    }

	return dudsig_;        
}

double fnpf_print_kinematics::dvdsig(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	// getting dvdsig for ax2 and 3
{
	double dvdsig_ = 0;

    if(k<p->knoz)
    dvdsig_ =  (c->V[FIJKp1] - c->V[FIJKm1])/(p->DZN[KP1] + p->DZN[KM1]);

    if(k==p->knoz)
    dvdsig_ = (c->V[FIJK] - c->V[FIJKm1])/(p->DZN[KM1]);

    return dvdsig_;        
}

double fnpf_print_kinematics::dudxi(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	// getting dudxi
{
    return (c->U[FIp1JK] - c->U[FIm1JK])/(p->DXN[IP1] + p->DXN[IM1]); 
}

double fnpf_print_kinematics::dvdxi(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	// getting dvdxi
{
    return (c->V[FIJp1K] - c->V[FIJm1K])/(p->DYN[JP1] + p->DYN[JM1]); 
}





