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

#include"bedload_EH.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

bedload_EH::bedload_EH(lexer *p)
{
    rhosed=p->S22;
    rhowat=p->W1;
    g=fabs(p->W22);
    d50=p->S20;
}

bedload_EH::~bedload_EH()
{
}

void bedload_EH::start(lexer* p, ghostcell* pgc, sediment_fdm *s)
{
    // Engelund-Hansen (1967), total load:  Phi = 0.1*theta^2.5/f,  f = 2*u*^2/V^2
    // SFLOW: P,Q are depth-averaged -> f from u* and V directly
    // CFD/NHFLOW: P,Q are near-bed -> f from Engelund's resistance law V/u* = 6 + 2.5 ln(h/ks)
	double qb,Ts,Tb,Rstar,fc,h,ks,uvel,vvel,u_abs,ustar;

	SEDSLICELOOP
    {
        rhowat = s->ro(i,j);
        Rstar = (rhosed-rhowat)/rhowat;

        Ts = s->shields_crit(i,j);
        Tb = s->shields_eff(i,j);

        fc = 0.0;

        if(p->A10==2)
        {
        uvel = 0.5*(s->P(i,j)+s->P(i-1,j));
        vvel = 0.5*(s->Q(i,j)+s->Q(i,j-1));
        u_abs = sqrt(uvel*uvel + vvel*vvel);
        ustar = s->shearvel_eff(i,j);

        if(u_abs>1.0e-10)
        fc = 2.0*(ustar*ustar)/(u_abs*u_abs);
        }

        if(p->A10!=2)
        {
        h  = s->waterlevel(i,j);
        ks = s->ks_eff(i,j);

        if(h>ks && ks>1.0e-20)
        fc = 2.0/pow(6.0 + 2.5*log(h/ks),2.0);
        }

        qb=0.0;

        if(s->active(i,j)==1 && Tb>=Ts && fc>1.0e-10)
        qb = (0.1/fc)*pow(Tb,2.5)*sqrt(Rstar*g*d50*d50*d50);

        s->qbe(i,j) = qb;
	}

    pgc->gcsl_start4(p,s->qbe,1);
}
