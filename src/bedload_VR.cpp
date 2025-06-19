/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"bedload_VR.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

bedload_VR::bedload_VR(lexer *p)
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

void bedload_VR::start(lexer* p, ghostcell* pgc, sediment_fdm *s)
{
    double Ti,r,f,Ts,Tb;
	double qb;
	
	SLICELOOP4
    {
        Ts = s->shields_crit(i,j);
	    Tb = s->shields_eff(i,j);

        Ti=MAX((Tb-Ts)/(Ts),0.0);
        
        if(s->active(i,j)==1 && Tb>=Ts)
        qb = (0.053*pow(d50,1.5)*sqrt(g*Rstar)*pow(Ti,2.1))/pow(Ds,0.3);


        if(s->active(i,j)==0 || Tb<Ts)
        qb=0.0;
		
		s->qbe(i,j) = qb;
	}
    
    pgc->gcsl_start4a(p,s->qbe,1);    
    
}
