/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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
#include"sediment_exner.h"
#include"reinitopo.h"
#include"suspended.h"
#include"bedload_VR.h"
#include"bedload_einstein.h"
#include"bedload_MPM.h"
#include"bedload_MPM.h"
#include"bedload_EF.h"
#include"bedload_void.h"
#include"bedshear.h"
#include"sandslide_f.h"
#include"sandslide_f2.h"
#include"sandslide_f3.h"
#include"sandslide_pde.h"
#include"sandslide_v.h"
#include"topo_relax.h"
#include"vrans_v.h"
#include"vrans_f.h"
#include"reduction_void.h"
#include"reduction_parker.h"
#include"reduction_deyemp.h"
#include"reduction_deyana.h"
#include"reduction_FD.h"

sediment_f::sediment_f(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb): bedslope(p)
{
    s = new sediment_fdm(p);
    
    
    if(p->S11==0)
    pbed = new bedload_void();

    if(p->S11==1)
    pbed = new bedload_VR(p);

    if(p->S11==2)
    pbed = new bedload_MPM(p);
	
	if(p->S11==3)
    pbed = new bedload_EF(p);
    
    if(p->S11==4)
    pbed = new bedload_einstein(p);
    
    
    if(p->S90==0)
    pslide=new sandslide_v(p);   
    
    if(p->S90==1)
    pslide=new sandslide_f(p);
    
    if(p->S90==2)
    pslide=new sandslide_f2(p);
    
    if(p->S90==3)
    pslide=new sandslide_f3(p);
    
    if(p->S90==4)
    pslide=new sandslide_pde(p);
    
    if(p->S10!=2 && p->A10==6)
	pvrans = new vrans_v(p,a,pgc);
	
	if(p->S10==2 && p->A10==6)
	pvrans = new vrans_f(p,a,pgc);
    
    
    if(p->S80==0)
    preduce=new reduction_void(p);

    if(p->S80==1)
    preduce=new reduction_parker(p);

    if(p->S80==2)
    preduce=new reduction_deyemp(p);

    if(p->S80==3)
    preduce=new reduction_deyana(p);
	
	if(p->S80==4)
    preduce=new reduction_FD(p);
    
    ptopo = new sediment_exner(p,pgc);
    
	
	p->gcin4a_count=p->gcin_count;
	p->gcout4a_count=p->gcout_count;
	
    
    prelax = new topo_relax(p);
	
	pbedshear  = new bedshear(p,pturb);
    
    volume_token=0;
}

sediment_f::~sediment_f()
{
}

void sediment_f::start_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow,
                                    reinitopo *preto, suspended *psusp)
{
    // bedshear stress
    sedcalc=0;
    
	if((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47))
	{
		if(p->S42==1 && p->count%p->S44==0)
		sediment_algorithm_cfd(p,a,pgc,pflow,preto,psusp);
		
		if(p->S42==2 && p->simtime>=p->sedsimtime)
		{
		sediment_algorithm_cfd(p,a,pgc,pflow,preto,psusp);
		p->sedsimtime = p->simtime + p->S46;
		}
		
		if(p->S42==3  && p->simtime/p->wT>=p->sedwavetime)
		{
		sediment_algorithm_cfd(p,a,pgc,pflow,preto,psusp);
		p->sedwavetime = p->simtime/p->wT + p->S48;
		}
    
    sedcalc=1;
	}
    
    if(sedcalc==0)
    pbedshear->taubed(p,a,pgc,s);
}

void sediment_f::start_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow, slice &P, slice &Q)
{
    
    if((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47))
	{
		if(p->S42==1 && p->count%p->S44==0)
		sediment_algorithm_sflow(p,b,pgc,pflow,P,Q);
		
		if(p->S42==2 && p->simtime>=p->sedsimtime)
		{
		sediment_algorithm_sflow(p,b,pgc,pflow,P,Q);
		p->sedsimtime = p->simtime + p->S46;
		}
		
		if(p->S42==3  && p->simtime/p->wT>=p->sedwavetime)
		{
		sediment_algorithm_sflow(p,b,pgc,pflow,P,Q);
		p->sedwavetime = p->simtime/p->wT + p->S48;
		}
	}
}




