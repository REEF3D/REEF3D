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

#include"bedconc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"

bedconc::bedconc(lexer *p, turbulence *pturb) : bedshear(p,pturb)
{
    rhosed=p->S22;
    rhowat=p->W1;
    g=9.81;
    d50=p->S20;
    shields=p->S30;
    visc=p->W2;
    kappa=0.4;
    ks=2.5*d50;
    Rstar=(rhosed-rhowat)/rhowat;
}

bedconc::~bedconc()
{
}

double bedconc::cbed(lexer* p,fdm* a,ghostcell *pgc, field& topo)
{
	double adist=2.0*d50;
	
    taubed(p,a,pgc,tau_eff);
    taucritbed(p,a,pgc,tau_crit);

    Ti=MAX((tau_eff-tau_crit)/tau_crit,0.0);

    Ds= d50*pow((Rstar*g)/(visc*visc),1.0/3.0);

    val = (0.015*d50*pow(Ti,1.5))/(pow(Ds,0.3)*adist);

    return val;
}

