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
#include"reinitopo.h"
#include"suspended.h"
#include"bedload.h"
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
    
    if(p->S10!=2)
	pvrans = new vrans_v(p,a,pgc);
	
	if(p->S10==2)
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
    
	
	p->gcin4a_count=p->gcin_count;
	p->gcout4a_count=p->gcout_count;
	
    
    prelax = new topo_relax(p);
	
	pbedshear  = new bedshear(p,pturb);
    
    volume_token=0;
}

sediment_f::~sediment_f()
{
}

void sediment_f::start(lexer *p, fdm *a, convection *pconvec, ghostcell *pgc, ioflow *pflow,
                                    topo *ptopo, reinitopo *preto, suspended *psusp, bedload *pbed)
{
	if((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47))
	{
		if(p->S42==1 && p->count%p->S44==0)
		sediment_algorithm(p,a,pconvec,pgc,pflow,ptopo,preto,psusp,pbed);
		
		if(p->S42==2 && p->simtime>=p->sedsimtime)
		{
		sediment_algorithm(p,a,pconvec,pgc,pflow,ptopo,preto,psusp,pbed);
		p->sedsimtime = p->simtime + p->S46;
		}
		
		if(p->S42==3  && p->simtime/p->wT>=p->sedwavetime)
		{
		sediment_algorithm(p,a,pconvec,pgc,pflow,ptopo,preto,psusp,pbed);
		p->sedwavetime = p->simtime/p->wT + p->S48;
		}
	}
}

void sediment_f::sediment_algorithm(lexer *p, fdm *a, convection *pconvec, ghostcell *pgc, ioflow *pflow,
                                    topo *ptopo, reinitopo *preto, suspended *psusp, bedload *pbed)
{
    starttime=pgc->timer();
    
    pgc->start1(p,a->u,14);
	pgc->start2(p,a->v,15);
	pgc->start3(p,a->w,16);
    
    // find bedk
    fill_bedk(p,a,pgc);
    
    // bedslope
    if(p->S83==2)
    slope_cds(p,a,pgc,s);
    
    if(p->S83==5)
    slope_weno(p,a,pgc,s);
    
    // bedslope reduction 
    preduce->start(p,a,pgc,s);
    
    // bedshear stress
	pbedshear->taubed(p,a,pgc,s);
    pbedshear->taucritbed(p,a,pgc,s);

    // bedload
    pbed->start(p,a,pgc,s);
	
    // Exner
    ptopo->start(a,p,pconvec,pgc,preto,s);
    
    // sandslide
    pslide->start(p,a,pgc,s);
    
	prelax->start(p,a,pgc);
	
    
    // filter bedzh
	if(p->S100>0)
	filter(p,a,pgc,a->bedzh,p->S100,p->S101);
	
	topo_zh_update(p,a,pgc);
    preto->start(a,p,a->topo,pconvec,pgc);

    volume_calc(p,a,pgc);
    
    pgc->start1(p,a->u,10);
	pgc->start2(p,a->v,11);
	pgc->start3(p,a->w,12);
    
    if(p->mpirank==0)
    cout<<"Topo: update grid..."<<endl;
    
    update(p,a,pgc,pflow);
    bedlevel(p,a,pgc); 
	
	pgc->start4(p,a->conc,40);
    
    if(p->mpirank==0 && p->count>0)
    cout<<"Sediment Timestep: "<<p->dtsed<<" Sediment Total Timestep: "<<p->dtsed<<"  Total Time: "<<setprecision(7)<<p->sedtime<<endl;

	if(p->mpirank==0)
    cout<<"Sediment CompTime: "<<setprecision(5)<<pgc->timer()-starttime<<endl<<endl;

}







