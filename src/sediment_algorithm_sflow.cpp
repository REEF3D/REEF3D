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
#include"sandslide.h"
#include"topo_relax.h"
#include"bedshear_reduction.h"

void sediment_f::sediment_algorithm_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow, slice &P, slice &Q)
{
    starttime=pgc->timer();
    
    ++p->sediter;
    
    // prep SFLOW
    prep_sflow(p,b,pgc,P,Q);
    
    // bedslope cds ******
    slope_cds(p,pgc,s);
    
    // bedslope reduction ******
    preduce->start(p,pgc,s);
    
    // bedshear stress -------
	pbedshear->taubed(p,b,pgc,s);
    pbedshear->taucritbed(p,b,pgc,s);

    // bedload *******
    pbed->start(p,pgc,s);
	
    // relax *******
	prelax->start(p,pgc,s);
    
    // Exner *******
    ptopo->start(p,pgc,s);
    
    // sandslide ********
    pslide->start(p,pgc,s);
    
    // relax *******
	prelax->start(p,pgc,s);
	
    // filter bedzh *******
	if(p->S100>0)
	filter(p,pgc,s->bedzh,p->S100,p->S101);
    
    // update sflow  --------
    update_sflow(p,b,pgc,pflow);
    
    // sediment log
    sedimentlog(p);
    
    if(p->mpirank==0 && p->count>0)
    cout<<"Sediment Iter: "<<p->sediter<<" Sediment Timestep: "<<p->dtsed<<"  Total Time: "<<setprecision(7)<<p->sedtime<<endl;

	if(p->mpirank==0)
    cout<<"Sediment CompTime: "<<setprecision(5)<<pgc->timer()-starttime<<endl<<endl;

}