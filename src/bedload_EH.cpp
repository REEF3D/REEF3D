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

#include"bedload_EH.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

bedload_EH::bedload_EH(lexer *p)
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
    fh= 1.0e-8;
}

bedload_EH::~bedload_EH()
{
}

void bedload_EH::start(lexer* p, ghostcell* pgc, sediment_fdm *s)
{
	double qb,qbx,qby,Ts,Tb;
	
	SLICELOOP4
    {
        Ts = s->shields_crit(i,j);
	    Tb = s->tau_eff(i,j);

        if(Tb>Ts)
        if(s->active(i,j)==1)
        qb =  (0.1/fh)*pow((s->tau_eff(i,j)*rhowat)/(g*d50*(rhosed-rhowat)),2.5);

        if(Tb<=Ts || s->active(i,j)==0)
        qb=0.0;
	
        s->qbe(i,j) = qb;
	}
    
    pgc->gcsl_start4(p,s->qbe,1);
    
}
