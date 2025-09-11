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

#include"bedload_MPM.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

bedload_MPM::bedload_MPM(lexer* p) 
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
}

bedload_MPM::~bedload_MPM()
{
}

void bedload_MPM::start(lexer* p, ghostcell* pgc, sediment_fdm *s)
{
    double qb,Ts,Tb;

	SEDSLICELOOP
    {
        Ts = s->shields_crit(i,j);
	    Tb = s->shields_eff(i,j);

        if(s->active(i,j)==1 && Tb>=Ts)
        qb = 8.0*pow(MAX(Tb - Ts,0.0),1.5)* p->S20*sqrt(((p->S22-p->W1)/p->W1)*fabs(p->W22)*p->S20);

        if(s->active(i,j)==0 || Tb<Ts)
        qb=0.0;
		
        s->qbe(i,j) = qb;
	}
    
    pgc->gcsl_start4(p,s->qbe,1);
}
