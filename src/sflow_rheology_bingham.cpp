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

#include"sflow_rheology_f.h"
#include"lexer.h"
#include"fdm2D.h"

double sflow_rheology_f::bingham(lexer *p, fdm2D *b, double vel, double u_abs, double press, double HIJ)
{
        if(p->W101==0)  // HB
        tau0=p->W96;
        
        if(p->W101==1)  // HB-C dry sand
        tau0 = MAX(0.0,tanphi*press + p->W102_c);//*(1.0-exp(-p->W103*vel));
    
    /*
    denom =  (p->W1*fabs(p->W22)*((b->eta(i+1,j)-b->eta(i-1,j))/p->DXM + (b->eta(i,j+1)-b->eta(i,j-1))/p->DXM));
    
    denom = fabs(denom)>1.0e-20?denom:1.0e20;
    
    hc = tau0/denom;
    
    // if h <hc , tau_b = tau0 otherwise 1.5tau0+3Ku/h
    if(hc < b->hp(i,j))
    val = tau0;
    
    if(hc > b->hp(i,j))*/
    val = (1.5*tau0 + 3.0*p->W97*(1.0/HIJ)*u_abs)*(vel/(fabs(u_abs)>1.0e-20?u_abs:1.0e20));
    
    
    return val;
}