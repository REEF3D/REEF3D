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

#include"bedload_VR.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

bedload_VR::bedload_VR(lexer *p, turbulence *pturb) : bedshear(p,pturb), bedload_noneq(p), epsi(1.6*p->DXM)
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
    
    // noneq ini
	
	SLICELOOP4
    {

        taubed(p,a,pgc,tau_eff,shearvel_eff,shields_eff);
        taucritbed(p,a,pgc,tau_crit,shearvel_crit,shields_crit);

        //Ti=MAX((shearvel_eff*shearvel_eff-shearvel_crit*shearvel_crit)/(shearvel_crit*shearvel_crit),0.0);
        Ti=MAX((shields_eff-shields_crit)/(shields_crit),0.0);
        
        if(shearvel_eff>shearvel_crit)
        qb =(0.053*pow(d50,1.5)*sqrt(g*Rstar)*pow(Ti,2.1))/pow(Ds,0.3)  ;

        if(shearvel_eff<=shearvel_crit)
        qb=0.0;
		
		a->bedload(i,j) = qb;
        
	}
    
    pgc->gcsl_start4a(p,a->bedload,1);    
    

    /*slice4 tt(p);
    
    SLICELOOP4
    {
    taubed(p,a,pgc,tau_eff,shearvel_eff,shields_eff);
    taucritbed(p,a,pgc,tau_crit,shearvel_crit,shields_crit);
    
    
    tt(i,j) = shields_crit;
    }
    
    ALOOP
    {
    a->test(i,j,k) = tt(i,j);
    }
    pgc->start4a(p,a->test,1);*/
}
