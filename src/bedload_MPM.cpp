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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"bedload_MPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"

bedload_MPM::bedload_MPM(lexer* p, turbulence *pturb) : bedshear(p,pturb), epsi(1.6*p->DXM)
{
    rhosed=p->S22;
    rhowat=p->W1;
    g=9.81;
    d50=p->S20;
    shields=p->S30;
    eta=0.053;
    visc=p->W2;
    kappa=0.4;
    ks=2.5*d50;
    repose=p->S25*(PI/180.0);

}

bedload_MPM::~bedload_MPM()
{
}

void bedload_MPM::start(lexer* p, fdm* a, ghostcell* pgc)
{
    double qb;

	SLICELOOP4
    {

        taubed(p,a,pgc,tau_eff,shearvel_eff,shields_eff);
        taucritbed(p,a,pgc,tau_crit,shearvel_crit,shields_crit);


        if(shields_eff>shields_crit)
        qb = 8.0*pow(MAX(shields_eff - shields_crit,0.0),1.5)* p->S20*sqrt(((p->S22-p->W1)/p->W1)*fabs(p->W22)*p->S20);

        if(shields_eff<=shields_crit)
        qb=0.0;
		
        a->bedload(i,j) = qb;
	}
    
    pgc->gcsl_start4(p,a->bedload,1);
}
