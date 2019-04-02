/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"sflow_sediment_f.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"


void sflow_sediment_f::bedload(lexer *p, fdm2D *b, ghostcell *pgc)
{
    bedload_vanRijn(p,b,pgc);
}

void sflow_sediment_f::bedload_vanRijn(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double Ti,r;
    double shearvel_eff, shearvel_crit;
    double d50,Rstar,Ds;
    
    d50=p->S20;
    Rstar=(p->S22-p->W1)/p->W1;
    Ds= d50*pow((Rstar*9.81)/(p->W2*p->W2),1.0/3.0);
    
	SLICELOOP4
    {
        //cout<<tau(i,j)<<" "<<taucr(i,j)<<endl;

        shearvel_eff = sqrt(tau(i,j)/p->W1);
        shearvel_crit = sqrt(taucr(i,j)/p->W1);

        Ti=MAX((shearvel_eff*shearvel_eff-shearvel_crit*shearvel_crit)/(shearvel_crit*shearvel_crit),0.0);

        if(shearvel_eff>shearvel_crit)
        b->qb(i,j) =(0.053*pow(d50,1.5)*sqrt(9.81*Rstar)*pow(Ti,2.1))/pow(Ds,0.3)  ;

        if(shearvel_eff<=shearvel_crit)
        b->qb(i,j) = 0.0;
	}
    
    pgc->gcsl_start4(p,b->qb,1);
    pgc->dgcslpol(p,b->qb,p->dgcsl4,p->dgcsl4_count,14);
    b->qb.ggcpol(p);
}