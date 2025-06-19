/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field4a.h"
#include"sediment_fdm.h"
#include"turbulence.h"

partres::partres(lexer *p, ghostcell *pgc) : P(p,pgc), bedch(p), Tau(p), Ts(p), cellSum(p), irand(100000), drand(100000.0),
                                               dPx(p),dPy(p),dPz(p),dTx(p),dTy(p),dTz(p)
{
    p->Darray(betaQ73,p->Q73);
	p->Darray(tan_betaQ73,p->Q73);
	p->Darray(dist_Q73,p->Q73);


	for(n=0;n<p->Q73;++n)
	betaQ73[n] = (p->Q73_b[n]+90.0)*(PI/180.0);

	for(n=0;n<p->Q73;++n)
	tan_betaQ73[n] = tan(betaQ73[n]);
    
    relax_ini(p);
    
    printcount=0;
    timestep_ini=0;
}

partres::~partres()
{

}








