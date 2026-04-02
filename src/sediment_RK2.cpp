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

#include"sediment_RK2.h"
#include"sediment_fdm.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"topo.h"
#include"reinitopo.h"
#include"nhflow_suspended.h"
#include"nhflow_diffusion.h"
#include"nhflow_scalar_convection.h"
#include"bedload.h"
#include"bedconc_VR.h"
#include"bedshear.h"
#include"sandslide.h"
#include"topo_relax.h"
#include"bedslope.h"
#include"bedshear_reduction.h"
#include"bedload_direction.h"

sediment_RK2::sediment_RK2(lexer *p, ghostcell *pgc, turbulence *pturb, patchBC_interface *ppBC) : sediment_f(p,pgc,pturb,ppBC), bedzh_n(p)
{

    pBC = ppBC;
    

    
    psed = this;
}

sediment_RK2::~sediment_RK2()
{
}

void sediment_RK2::RK2_step1_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow)
{
    starttime=pgc->timer();
    
    ++p->sediter;
    
    // RK Step 1
    SLICELOOP4
    bedzh_n(i,j) = s->bedzh(i,j);
    
    // prep NHFLOW -------
    prep_nhflow(p,d,pgc);
    
    // bedslope cds ******
    pslope->slope_cds(p,pgc,s);
    
    // bedslope reduction ******
    preduce->start(p,pgc,s);
    
    // bedshear stress -------
	pbedshear->taubed(p,d,pgc,s);
    pbedshear->taucritbed(p,d,pgc,s);
    
    // bedload *******
    pbed->start(p,pgc,s);
    
    // bedload_direction *******
    pbeddir->start(p,pgc,s);
    
    // suspended load -------
    pcbed->start(p,pgc,s);
	
    // relax *******
	prelax->start(p,pgc,s);
    
    // Exner *******
    ptopo->start_RK(p,pgc,s);
    p->sedtime+=p->dtsed;
    
    // RK Step 1
    SLICELOOP4
    s->bedzh(i,j) += p->dtsed*s->vz(i,j);
    
    // sandslide ********
    if(p->sediter%p->S94==0)
    pslide->start(p,pgc,s);
    
    // relax *******
	prelax->start(p,pgc,s);
    
    // update sflow  --------
    update_nhflow(p,d,pgc,pflow);
    
    // sediment print
    print_probes(p,pgc,s,pflow);
    
    // sediment log
    sedimentlog(p);
    
    
    if(p->mpirank==0 && p->count>0)
    cout<<"Sediment Iter RK1: "<<p->sediter<<" Sediment Timestep: "<<p->dtsed<<"  Total Time: "<<setprecision(7)<<p->sedtime<<endl;

	if(p->mpirank==0)
    cout<<"Sediment CompTime: "<<setprecision(5)<<pgc->timer()-starttime<<endl;
}

void sediment_RK2::RK2_step2_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow)
{
    starttime=pgc->timer();
    
    ++p->sediter;
    
    // prep NHFLOW -------
    prep_nhflow(p,d,pgc);
    
    // bedslope cds ******
    pslope->slope_cds(p,pgc,s);
    
    // bedslope reduction ******
    preduce->start(p,pgc,s);
    
    // bedshear stress -------
	pbedshear->taubed(p,d,pgc,s);
    pbedshear->taucritbed(p,d,pgc,s);
    
    // bedload *******
    pbed->start(p,pgc,s);
    
    // bedload_direction *******
    pbeddir->start(p,pgc,s);
    
    // suspended load -------
    pcbed->start(p,pgc,s);
	
    // relax *******
	prelax->start(p,pgc,s);
    
    // Exner *******
    ptopo->start_RK(p,pgc,s);
    p->sedtime+=p->dtsed;
    
    // RK Step 2
    SLICELOOP4
    s->bedzh(i,j) = 0.5*bedzh_n(i,j) +  0.5*s->bedzh(i,j) + 0.5*p->dtsed*s->vz(i,j);
    
    // sandslide ********
    if(p->sediter%p->S94==0)
    pslide->start(p,pgc,s);
    
    // relax *******
	prelax->start(p,pgc,s);
    
    // update sflow  --------
    update_nhflow(p,d,pgc,pflow);
    
    // sediment print
    print_probes(p,pgc,s,pflow);
    
    // sediment log
    sedimentlog(p);
    
    
    if(p->mpirank==0 && p->count>0)
    cout<<"Sediment Iter RK2: "<<p->sediter<<" Sediment Timestep: "<<p->dtsed<<"  Total Time: "<<setprecision(7)<<p->sedtime<<endl;

	if(p->mpirank==0)
    cout<<"Sediment CompTime: "<<setprecision(5)<<pgc->timer()-starttime<<endl;
    
    
}