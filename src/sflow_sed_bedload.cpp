/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
    if(p->S11==1)
    bedload_vanRijn(p,b,pgc);
    
    if(p->S11==2)
    bedload_MPM(p,b,pgc);
    
    if(p->S11==3)
    bedload_EF(p,b,pgc);
    
    if(p->S11==4)
    bedload_vanRijn_C(p,b,pgc);
}

void sflow_sediment_f::bedload_vanRijn(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double Ti,r;
    double tau_eff, tau_crit;
    double d50,Rstar,Ds;
    
    d50=p->S20;
    Rstar=(p->S22-p->W1)/p->W1;
    Ds= d50*pow((Rstar*9.81)/(p->W2*p->W2),1.0/3.0);
    
	SLICELOOP4
    {
        b->test(i,j) = 0.0;
        
        tau_eff = tau(i,j);
        tau_crit = taucr(i,j);
        
        if(tau_eff>tau_crit)
        b->test(i,j) = 1.0;

        Ti=MAX((tau_eff-tau_crit)/tau_crit,0.0);

        if(tau_eff>tau_crit)
        b->qb(i,j)  = (0.053*pow(d50,1.5)*sqrt(9.81*Rstar)*pow(Ti,2.1))/pow(Ds,0.3)  ;

        if(tau_eff<=tau_crit)
        b->qb(i,j) = 0.0;
	}
    
    pgc->gcsl_start4(p,b->qb,1);
}

void sflow_sediment_f::bedload_MPM(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double shearvel_eff, shearvel_crit;
     
	SLICELOOP4
    {
        shearvel_eff = sqrt(tau(i,j)/p->W1);
        shearvel_crit = sqrt(taucr(i,j)/p->W1);
        
        b->test(i,j) = 0.0;
        
        if(shearvel_eff>shearvel_crit)
        b->test(i,j) = 1.0;

        if(shearvel_eff>shearvel_crit)
        b->qb(i,j)  = 8.0*pow(MAX(shearvel_eff - shearvel_crit,0.0),1.5)* p->S20*sqrt(((p->S22-p->W1)/p->W1)*fabs(p->W22)*p->S20);

        if(shearvel_eff<=shearvel_crit)
        b->qb(i,j) = 0.0;
	}
    
    pgc->gcsl_start4(p,b->qb,1);
}


void sflow_sediment_f::bedload_EF(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double Ti,r;
    double shearvel_eff, shearvel_crit;
    double d50;
    
    d50=p->S20;

	SLICELOOP4
    {
        shearvel_eff = sqrt(tau(i,j)/p->W1);
        shearvel_crit = sqrt(taucr(i,j)/p->W1);
        
        b->test(i,j) = 0.0;
        
        if(shearvel_eff>shearvel_crit)
        b->test(i,j) = 1.0;

        
        if(shearvel_eff>shearvel_crit)
        b->qb(i,j)  = d50*sqrt((p->S22/p->W1-1.0)*9.81*d50) * 11.6* (shearvel_eff-shearvel_crit)*(sqrt(shearvel_eff) - 0.7*sqrt(shearvel_crit));

        if(shearvel_eff<=shearvel_crit)
        b->qb(i,j) = 0.0;
	}
    
    pgc->gcsl_start4(p,b->qb,1);
}

void sflow_sediment_f::bedload_vanRijn_C(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double Ti,r;
    double tau_eff, tau_crit;
    double d50,Rstar,Ds;
    
    d50=p->S20;
    Rstar=(p->S22-p->W1)/p->W1;
    Ds= d50*pow((Rstar*9.81)/(p->W2*p->W2),1.0/3.0);
    
	SLICELOOP4
    {
        b->test(i,j) = 0.0;
        
        tau_eff = tau(i,j);
        tau_crit = taucr(i,j);
        
        if(tau_eff>tau_crit)
        b->test(i,j) = 1.0;

        Ti=MAX((tau_eff-tau_crit)/tau_crit,0.0);

        if(tau_eff>tau_crit)
        b->qb(i,j)  = (0.053*pow(d50,1.5)*sqrt(9.81*Rstar)*pow(Ti,2.1))/pow(Ds,0.3)  ;

        if(tau_eff<=tau_crit)
        b->qb(i,j) = 0.0;
	}
    
    pgc->gcsl_start4(p,b->qb,1);
}
