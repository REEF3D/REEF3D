/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"bedload_EF.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

bedload_EF::bedload_EF(lexer *p, turbulence *pturb) : bedshear(p,pturb), epsi(1.6*p->DXM)
{
    rhosed=p->S22;
    rhowat=p->W1;
    g=9.81;
    d50=p->S20;
    visc=p->W2;
    kappa=0.4;
    ks=p->S21*d50;
    repose=p->S25* (3.14159/180.0);
    Rstar=(rhosed-rhowat)/rhowat;
    Ds= d50*pow((Rstar*g)/(visc*visc),1.0/3.0);
}

bedload_EF::~bedload_EF()
{
}

void bedload_EF::start(lexer* p, fdm* a, ghostcell* pgc)
{

	double qb,qbx,qby,Ts,Tb;
	
	SLICELOOP4
    {
        taubed(p,a,pgc,tau_eff,shearvel_eff,shields_eff);
        taucritbed(p,a,pgc,tau_crit,shearvel_crit,shields_crit);

        Ts = shields_crit;
	    Tb = shields_eff;

        if(Tb>Ts)
        qb =  d50*sqrt((rhosed/rhowat-1.0)*g*d50) * 11.6* (Tb-Ts)*(sqrt(Tb) - 0.7*sqrt(Ts));

        if(Tb<=Ts)
        qb=0.0;
	
        a->bedload(i,j) = qb;
	}
    
    pgc->gcsl_start4(p,a->bedload,1);
}
