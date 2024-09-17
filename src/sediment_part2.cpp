/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"sediment_part2.h"
#include"sediment_fdm.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"ioflow.h"
#include"topo.h"
#include"bedshear.h"
#include"patchBC_interface.h"
#include"bedslope.h"

sediment_part2::sediment_part2(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb, patchBC_interface *ppBC)
{

    pBC = ppBC;
    
    sediment_logic(p,a,pgc,pturb);

	p->gcin4a_count=p->gcin_count;
	p->gcout4a_count=p->gcout_count;
    
    pslope = new bedslope(p);
	
    
    volume_token=0;
    
    if(p->F50==1)
	gcval_eta = 51;
    
    if(p->F50==2)
	gcval_eta = 52;
    
    if(p->F50==3)
	gcval_eta = 53;
    
    if(p->F50==4)
	gcval_eta = 54;
}

sediment_part2::~sediment_part2()
{
}

void sediment_part2::start_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reinitopo *preto, solver *psolv)
{
    sedcalc=0;
    
	if((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47))
	{

		sediment_algorithm_cfd(p,a,pgc,pflow,preto,psolv);
		

    
    sedcalc=1;
	}
    
    if(sedcalc==0)
    {
    fill_bedk(p,a,pgc);
    waterlevel(p,a,pgc);
    pbedshear->taubed(p,a,pgc,s);
    }
}