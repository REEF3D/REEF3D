/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
#include"sediment_f.h"
#include"sediment_fdm.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"ioflow.h"
#include"topo.h"
#include"bedshear.h"
#include"patchBC_interface.h"

sediment_f::sediment_f(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb, patchBC_interface *ppBC): bedslope(p)
{

    pBC = ppBC;
    
    sediment_logic(p,a,pgc,pturb);

	p->gcin4a_count=p->gcin_count;
	p->gcout4a_count=p->gcout_count;
	
    
    volume_token=0;
}

sediment_f::~sediment_f()
{
}

void sediment_f::start_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reinitopo *preto, solver *psolv)
{
    // bedshear stress
    sedcalc=0;
    
	if((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47))
	{
		if(p->S42==1 && p->count%p->S44==0)
		sediment_algorithm_cfd(p,a,pgc,pflow,preto,psolv);
		
		if(p->S42==2 && p->simtime>=p->sedsimtime)
		{
		sediment_algorithm_cfd(p,a,pgc,pflow,preto,psolv);
		p->sedsimtime = p->simtime + p->S46;
		}
		
		if(p->S42==3  && p->simtime/p->wT>=p->sedwavetime)
		{            
		sediment_algorithm_cfd(p,a,pgc,pflow,preto,psolv);
		p->sedwavetime = p->simtime/p->wT + p->S48;
		}
    
    sedcalc=1;
	}
    
    if(sedcalc==0)
    {
    fill_bedk(p,a,pgc);
    pbedshear->taubed(p,a,pgc,s);
    }
}

void sediment_f::start_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow, slice &P, slice &Q)
{
    
    if((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47 && p->count>0))
	{
		if(p->S42==1 && p->count%p->S44==0)
		sediment_algorithm_sflow(p,b,pgc,pflow,P,Q);
		
		if(p->S42==2 && p->simtime>=p->sedsimtime)
		{
		sediment_algorithm_sflow(p,b,pgc,pflow,P,Q);
		p->sedsimtime = p->simtime + p->S46;
		}
		
		if(p->S42==3  && p->simtime/p->wT>=p->sedwavetime )
		{
		sediment_algorithm_sflow(p,b,pgc,pflow,P,Q);
		p->sedwavetime = p->simtime/p->wT + p->S48;
		}
	}
}




