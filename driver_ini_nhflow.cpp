/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"driver.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"freesurface_header.h"
#include"turbulence_header.h"
#include"momentum_header.h"
#include"pressure_header.h"
#include"fdm_header.h"
#include"sediment_header.h"
#include"convection_header.h"
#include"solver_header.h"
#include"field_header.h"
#include"heat_header.h"
#include"concentration_header.h"
#include"benchmark_header.h"
#include"6DOF_header.h"
#include"waves_header.h"
#include"lexer.h"
#include"cart1.h"
#include"cart2.h"
#include"cart3.h"
#include"cart4.h"
#include"cart4a.h"
#include<sys/stat.h>
#include<sys/types.h>

#define WLVL (fabs(a->WL(i,j))>1.0e-20?a->WL(i,j):1.0-20)

void driver::driver_ini_nhflow()
{

    p->count=0;

    p->cellnumtot=pgc->globalisum(p->cellnum);
    p->pointnumtot=pgc->globalisum(p->pointnum);

	log_ini();
    
    if(p->mpirank==0)
    cout<<"number of cells: "<<p->cellnumtot<<endl;

    if(p->mpirank==0)
    cout<<"starting driver_ini_NHFLOW"<<endl;
    
    // SIGMA grid
     // bed ini
    SLICELOOP4
	a->bed(i,j) = p->bed[IJ];
    
    pgc->gcsl_start4(p,a->bed,50);
    
    // eta ini
	SLICELOOP4
	a->eta(i,j) = 0.0;

    pgc->gcsl_start4(p,a->eta,50);
    
    SLICELOOP4
    a->WL(i,j) = MAX(0.0,a->eta(i,j) + p->wd - a->bed(i,j));
    
    // sigma ini
    p->sigma_ini(p,a,pgc,a->eta);
    p->sigma_update(p,a,pgc,a->eta);

    
    //ioflow ini
    pflow->ini_nhflow(p,a,pgc);
    
    pflow->eta_relax(p,pgc,a->eta);
    pgc->gcsl_start4(p,a->eta,50);
    
    if(p->P150==0)
	pdata = new data_void(p,a,pgc);
	
	if(p->P150>0)
	pdata = new data_f(p,a,pgc);
	
	pdata->start(p,a,pgc);
	
    pheat->heat_ini(p,a,pgc,pheat);
	pconc->ini(p,a,pgc,pconc);

    ptstep->ini(a,p,pgc);
	pflow->gcio_update(p,a,pgc);
	pflow->pressure_io(p,a,pgc);
    
    
    
    // inflow ini
	pflow->discharge(p,a,pgc);
	pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	potflow->start(p,a,psolv,pgc);
    pflow->wavegen_precalc(p,pgc);
    
	if(p->I12>=1)
	pini->hydrostatic(p,a,pgc);

	if(p->I11==1)
	ptstep->start(a,p,pgc,pturb);
    
    if(p->I13==1)
    pturb->ini(p,a,pgc);
	
	pflow->pressure_io(p,a,pgc);

	pgc->start1(p,a->u,10);
	pgc->start2(p,a->v,11);
	pgc->start3(p,a->w,12);

    pgc->start4(p,a->press,40);
    pgc->dgcpol(p,a->topo,p->dgc4,p->dgc4_count,14);
	a->topo.ggcpol(p);
	
	if(p->I40==1)
	pini->stateini(p,a,pgc,pturb);
    
	pgc->start4(p,a->press,40);
	
    p->sigma_update(p,a,pgc,a->eta);
    pprint->start(a,p,pgc,pturb,pheat,pflow,psolv,pdata,pconc,psed);


// ini variables
    for(int qn=0; qn<2; ++qn)
    {
    pturb->ktimesave(p,a,pgc);
    pturb->etimesave(p,a,pgc);
    }

    p->gctime=0.0;
    p->xtime=0.0;
	p->reinitime=0.0;
	p->wavetime=0.0;
	p->field4time=0.0;
}




