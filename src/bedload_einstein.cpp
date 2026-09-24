/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"bedload_einstein.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

bedload_einstein::bedload_einstein(lexer* p)
{
    rhosed=p->S22;
    g=fabs(p->W22);
    d50=p->S20;
}

bedload_einstein::~bedload_einstein()
{
}

void bedload_einstein::start(lexer* p, ghostcell* pgc, sediment_fdm *s)
{
    // Einstein-Brown:  Phi = 2.15*exp(-0.391/theta)  for theta < 0.182
    //                  Phi = 40*theta^3              for theta >= 0.182
    double qb,Tb,Rstar;

	SEDSLICELOOP
    {
        rhowat = s->ro(i,j);
        Rstar = (rhosed-rhowat)/rhowat;

        Tb = s->shields_eff(i,j);

        qb=0.0;

        if(s->active(i,j)==1 && Tb>1.0e-10)
        {
            if(Tb<0.182)
            qb = 2.15*exp(-0.391/Tb);

            if(Tb>=0.182)
            qb = 40.0*Tb*Tb*Tb;

        qb *= sqrt(Rstar*g*d50*d50*d50);
        }

        s->qbe(i,j) = qb;
	}

    pgc->gcsl_start4(p,s->qbe,1);
}
