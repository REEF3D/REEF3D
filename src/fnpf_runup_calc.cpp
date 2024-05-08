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
Author: Arun Kamath
--------------------------------------------------------------------*/

#include"fnpf_runup.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include <math.h>

void fnpf_runup::fnpf_runup_calc(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{	
    // -> here go your runup formulas....
    
    // R1 = ....?
    
    
    
    
    
    
    
    
    
    
    
    
}

double fnpf_runup::dndt(lexer *p, fdm_fnpf *c, ghostcell *pgc) // to calculate dn dt for ax3
{
    double dndt = (c->eta(i,j) - etan)/ p->dt;
 
    return dndt;
}

double fnpf_runup::dudsig(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	// getting dudsig for ax2 and 3
{
    double dudsig_ = 0;
	double dudsig2_ = 0;
    
    if(k<p->knoz)
    dudsig_ =  (c->U[FIJKp1] - c->U[FIJKm1])/(p->DZN[KP1] + p->DZN[KM1]);


    if(k==p->knoz)
    dudsig_ = (c->U[FIJK] - c->U[FIJKm1])/(p->DZN[KM1]);

	return dudsig_;        
}

double fnpf_runup::dvdsig(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	// getting dvdsig for ax2 and 3
{
	double dvdsig_ = 0;

    if(k<p->knoz)
    dvdsig_ =  (c->V[FIJKp1] - c->V[FIJKm1])/(p->DZN[KP1] + p->DZN[KM1]);


    if(k==p->knoz)
    dvdsig_ = (c->V[FIJK] - c->V[FIJKm1])/(p->DZN[KM1]);

	return dvdsig_;        
}


double fnpf_runup::dudxi(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	// getting dudxi
{
    return (c->U[FIp1JK] - c->U[FIm1JK])/(p->DXN[IP1] + p->DXN[IM1]); 
}

double fnpf_runup::dvdxi(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	// getting dvdxi
{
    return (c->V[FIJp1K] - c->V[FIJm1K])/(p->DYN[JP1] + p->DYN[JM1]); 
}






