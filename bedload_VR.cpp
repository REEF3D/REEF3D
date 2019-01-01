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

#include"bedload_VR.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

bedload_VR::bedload_VR(lexer *p, turbulence *pturb) : bedshear(p,pturb), epsi(1.6*p->DXM)
{
    rhosed=p->S22;
    rhowat=p->W1;
    g=9.81;
    d50=p->S20;
    visc=p->W2;
    kappa=0.4;
    ks=p->S21*d50;
    Rstar=(rhosed-rhowat)/rhowat;
    Ds= d50*pow((Rstar*g)/(visc*visc),1.0/3.0);
}

bedload_VR::~bedload_VR()
{
}

void bedload_VR::start(lexer* p, fdm* a, ghostcell* pgc)
{
    double Ti,r;
	double qb;
	
	SLICELOOP4
    {

        taubed(p,a,pgc,tau_eff,shearvel_eff,shields_eff);
        taucritbed(p,a,pgc,tau_crit,shearvel_crit,shields_crit);

        Ti=MAX((shearvel_eff*shearvel_eff-shearvel_crit*shearvel_crit)/(shearvel_crit*shearvel_crit),0.0);

        if(shearvel_eff>shearvel_crit)
        qb =(0.053*pow(d50,1.5)*sqrt(g*Rstar)*pow(Ti,2.1))/pow(Ds,0.3)  ;

        if(shearvel_eff<=shearvel_crit)
        qb=0.0;
		
		a->bedload(i,j) = qb;
        
	}
    
    SLICELOOP1
    {
        taubedx(p,a,pgc,tau_eff,shearvel_eff,shields_eff);
        taucritbed(p,a,pgc,tau_crit,shearvel_crit,shields_crit);

        Ti=MAX((shearvel_eff*shearvel_eff-shearvel_crit*shearvel_crit)/(shearvel_crit*shearvel_crit),0.0);

        if(shearvel_eff>shearvel_crit)
        qb =(0.053*pow(d50,1.5)*sqrt(g*Rstar)*pow(Ti,2.1))/pow(Ds,0.3)  ;

        if(shearvel_eff<=shearvel_crit)
        qb=0.0;
		
		a->qbx(i,j) = qb;
    }
    
    SLICELOOP2
    {
        taubedy(p,a,pgc,tau_eff,shearvel_eff,shields_eff);
        taucritbed(p,a,pgc,tau_crit,shearvel_crit,shields_crit);

        Ti=MAX((shearvel_eff*shearvel_eff-shearvel_crit*shearvel_crit)/(shearvel_crit*shearvel_crit),0.0);

        if(shearvel_eff>shearvel_crit)
        qb =(0.053*pow(d50,1.5)*sqrt(g*Rstar)*pow(Ti,2.1))/pow(Ds,0.3)  ;

        if(shearvel_eff<=shearvel_crit)
        qb=0.0;
		
		a->qby(i,j) = qb;
        
    }
    
    pgc->gcsl_start4(p,a->bedload,1);
    pgc->dgcslpol(p,a->bedload,p->dgcsl4,p->dgcsl4_count,14);
    a->bedload.ggcpol(p);
    
    pgc->gcsl_start1(p,a->qbx,1);
    pgc->dgcslpol(p,a->qbx,p->dgcsl1,p->dgcsl1_count,11);
    a->qbx.ggcpol(p);
    
    pgc->gcsl_start2(p,a->qby,1);
    pgc->dgcslpol(p,a->qby,p->dgcsl2,p->dgcsl2_count,12);
    a->qby.ggcpol(p);
}
