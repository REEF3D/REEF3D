/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"driver.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm_nhf.h"
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

#define WLVL (fabs(a->WL(i,j))>1.0e-20?a->WL(i,j):1.0e20)

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
	d->bed(i,j) = p->bed[IJ];
    
    pgc->gcsl_start4(p,d->bed,50);
    
    // eta ini
	SLICELOOP4
    {
	d->eta(i,j) = 0.0;
    p->wet[IJ] = 1;
    }
    
    pgc->gcsl_start4(p,d->eta,50);
    
    SLICELOOP4
    d->WL(i,j) = MAX(0.0,d->eta(i,j) + p->wd - d->bed(i,j));
    
    // sigma ini
    p->sigma_ini(p,d,pgc,d->eta);
    p->sigma_update(p,d,pgc,d->eta,d->eta,1.0);

    //ioflow ini
    pflow->ini_nhflow(p,d,pgc); // replace a with d

    pflow->eta_relax(p,pgc,d->eta);
    pgc->gcsl_start4(p,d->eta,50);

    if(p->P150==0)
	pdata = new data_void(p,a,pgc);
	
	if(p->P150>0)
	pdata = new data_f(p,a,pgc);
	
	pdata->start(p,a,pgc);

    pnhfstep->ini(p,d,pgc);
 
	pflow->gcio_update(p,a,pgc); 
	pflow->pressure_io(p,a,pgc);
     
    // inflow ini
	pflow->discharge_nhflow(p,d,pgc);

    pflow->wavegen_precalc(p,pgc);

	if(p->I11==1)
	ptstep->start(a,p,pgc,pturb);
    
    if(p->I13==1)
    pturb->ini(p,a,pgc);
	
	pflow->pressure_io(p,a,pgc);

    
    SLICELOOP4
    d->eta_n(i,j) = d->eta(i,j);

	pgc->start4V(p,d->U,d->bc,10);
    pgc->start4V(p,d->V,d->bc,11);
    pgc->start4V(p,d->W,d->bc,12);
    pgc->start4V(p,d->P,d->bc,540);
    
    pnhf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta,d->eta_n,1.0);
    p->sigma_update(p,d,pgc,d->eta,d->eta,1.0);

    SLICELOOP4
    d->WL(i,j) = MAX(0.0, d->eta(i,j) + p->wd - d->bed(i,j));
    
    //pprint->start(a,p,pgc,pturb,pheat,pflow,psolv,pdata,pconc,pmp,psed);

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




