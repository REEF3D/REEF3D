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

#include"vrans_nhflow_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"

void vrans_nhflow_f::u_source(lexer *p, fdm_nhf *d, slice &WL)
{
	// VRANS porosity
    if(p->B200>0)
    LOOP
	{
        H = Hporface(p,d,0,0,0);  
        
        porval   = d->POR[IJK];
        partval  = d->PORPART[IJK];
        alphaval = APOR[IJK];
        betaval  = BPOR[IJK];
        viscval  = d->VISC[IJK];
		
        
        Aporval = Apor(porval,partval,alphaval,viscval);
        Bporval = Bpor(porval,partval,betaval);
            

        porousterm = Aporval*d->U[IJK] + Bporval*d->U[IJK]*fabs(d->U[IJK]); 
    	
    d->F[IJK] -= H*WL(i,j)*porousterm;
	}
}

void vrans_nhflow_f::v_source(lexer *p, fdm_nhf *d, slice &WL)
{
	// VRANS porosity
    if(p->B200>0)
    LOOP
	{
        H = Hporface(p,d,0,0,0);  
        
        porval   = d->POR[IJK];
        partval  = d->PORPART[IJK];
        alphaval = APOR[IJK];
        betaval  = BPOR[IJK];
        viscval  = d->VISC[IJK];
		
        
        Aporval = Apor(porval,partval,alphaval,viscval);
        Bporval = Bpor(porval,partval,betaval);
            

        porousterm = Aporval*d->V[IJK] + Bporval*d->V[IJK]*fabs(d->V[IJK]); 
    	
    d->G[IJK] -= H*WL(i,j)*porousterm;
	}
}

void vrans_nhflow_f::w_source(lexer *p, fdm_nhf *d, slice &WL)
{
	// VRANS porosity
    if(p->B200>0)
    LOOP
	{
        H = Hporface(p,d,0,0,0);  
        
        porval   = d->POR[IJK];
        partval  = d->PORPART[IJK];
        alphaval = APOR[IJK];
        betaval  = BPOR[IJK];
        viscval  = d->VISC[IJK];
		
        
        Aporval = Apor(porval,partval,alphaval,viscval);
        Bporval = Bpor(porval,partval,betaval);
            

        porousterm = Aporval*d->W[IJK] + Bporval*d->W[IJK]*fabs(d->W[IJK]); 
    	
    d->H[IJK] -= H*WL(i,j)*porousterm;
	}
}

double vrans_nhflow_f::Apor(double por, double part, double alpha, double visc)
{
	val = alpha*(pow(1.0-por,2.0)/pow(por,3.0))*(viscval/pow(part,2.0));
    
    if(val!=val)
    val=0.0;
	
	return val;
}

double vrans_nhflow_f::Bpor(double por, double part, double beta)
{
	val = beta*(1.0 + 7.5/Cval)*((1.0-por)/pow(por,3.0))*(1.0/part);
    
    if(val!=val)
    val=0.0;
	
	return val;
}
