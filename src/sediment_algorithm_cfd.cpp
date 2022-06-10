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
#include"ioflow.h"
#include"topo.h"
#include"reinitopo.h"
#include"suspended.h"
#include"bedload.h"
#include"bedshear.h"
#include"sandslide.h"
#include"topo_relax.h"
#include"bedshear_reduction.h"

void sediment_f::sediment_algorithm_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow,
                                    topo *ptopo, reinitopo *preto, suspended *psusp)
{
    starttime=pgc->timer();
    
    // Prep CFD ---------------
    // vel prep --------
    pgc->start1(p,a->u,14);
	pgc->start2(p,a->v,15);
	pgc->start3(p,a->w,16);
    
    // find bedk -------
    fill_bedk(p,a,pgc);
    
    // fill P/Q --------
    //-------------------------
    
    
    
    // Prep SFLOW ++++++++
    // fill P/Q ++++++++++
    // +++++++++++++++++++
    
    
    
    // bedslope cds ******
    if(p->S83==2)
    slope_cds(p,pgc,s);
    
    // bedslope weno -------
    if(p->S83==5)
    slope_weno(p,pgc,s,a->topo);
    
    // bedslope reduction ******
    preduce->start(p,pgc,s);
    
    // bedshear stress -------
	pbedshear->taubed(p,a,pgc,s);
    pbedshear->taucritbed(p,a,pgc,s);

    // bedload *******
    pbed->start(p,pgc,s);
	
    // Exner *******
    ptopo->start(p,pgc,s);
    
    // sandslide ********
    pslide->start(p,pgc,s);
    
    // relax bedzh *******
	prelax->start(p,pgc,s);
	
    // filter bedzh *******
	if(p->S100>0)
	filter(p,a,pgc,s->bedzh,p->S100,p->S101);
    
    // -------------------
    // update cfd -------
	topo_zh_update(p,a,pgc,s);
    preto->start(p,a,pgc,a->topo);

    volume_calc(p,a,pgc);
    
    pgc->start1(p,a->u,10);
	pgc->start2(p,a->v,11);
	pgc->start3(p,a->w,12);
    
    if(p->mpirank==0)
    cout<<"Topo: update grid..."<<endl;
    
    update_cfd(p,a,pgc,pflow);
    bedlevel(p,a,pgc); 
	
	pgc->start4(p,a->conc,40);
    // ---------------------
    
    if(p->mpirank==0 && p->count>0)
    cout<<"Sediment Timestep: "<<p->dtsed<<" Sediment Total Timestep: "<<p->dtsed<<"  Total Time: "<<setprecision(7)<<p->sedtime<<endl;

	if(p->mpirank==0)
    cout<<"Sediment CompTime: "<<setprecision(5)<<pgc->timer()-starttime<<endl<<endl;
}







