/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
#include"nhflow_header.h"
#include"lexer.h"
#include<sys/stat.h>
#include<sys/types.h>

#define WLVL (fabs(a->WL(i,j))>1.0e-20?a->WL(i,j):1.0e20)

void driver::driver_ini_nhflow()
{
    
    pnhf->ini(p,d,pgc,pflow);

	log_ini();
    
    if(p->mpirank==0)
    cout<<"starting driver_ini_NHFLOW"<<endl;
    

    // sigma ini
    pnhfmom->inidisc(p,d,pgc,pnhfsf);
    
    
    //ioflow ini
    pflow->ini_nhflow(p,d,pgc);
    pnhfsf->wetdry(p,d,pgc,d->U,d->V,d->W,d->WL); 
    
    // sigma ini
    pnhfmom->inidisc(p,d,pgc,pnhfsf);
    
    for(int qn=0;qn<20;++qn)
    {
    pflow->eta_relax(p,pgc,d->eta);
    pflow->WL_relax(p,pgc,d->WL,d->depth);
    }
    pgc->gcsl_start4(p,d->eta,50);

    pnhfstep->ini(p,d,pgc);
 
	pflow->gcio_update_nhflow(p,d,pgc); 

    // inflow ini
	pflow->discharge_nhflow(p,d,pgc);
    pflow->wavegen_precalc_nhflow(p,d,pgc);
    
    SLICELOOP4
    d->WL(i,j) = d->eta(i,j) + d->depth(i,j);
    
    SLICELOOP4
    d->eta_n(i,j) = d->eta(i,j);
    
    LOOP
    {
    d->RO[IJK] = p->W1;
    d->VISC[IJK] = p->W2;
    }
    
    SLICELOOP4
    d->ks(i,j) = p->B50;
    
    pgc->gcsl_start4(p,d->ks,50);
    
    pgc->start4V(p,d->RO,1);
    pgc->start4V(p,d->VISC,1);

	pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->W,12);
    pgc->start5V(p,d->P,540);
    
    // forcing ini
    pnhfdf->forcing_ini(p,d,pgc);

    pnhfsf->ini(p,d,pgc,pflow,d->U,d->V,d->W);
    pnhfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    
    pflow->gcio_update_nhflow(p,d,pgc); 
    
    // potential ini
    pnhfpot->start(p,d,ppoissonsolv,pgc);
    
    pflow->discharge_nhflow(p,d,pgc);
    pflow->inflow_nhflow(p,d,pgc,d->U,d->V,d->W,d->UH,d->VH,d->WH);
    
    // turbulence ini
    pnhfturb->ini(p, d, pgc);
    
    //sediment ini
    psed->ini_nhflow(p,d,pgc);
    
    //6DOF ini
    p6dof->initialize(p, d, pgc, pnet);
    
    pnhfprint->start(p,d,pgc,pflow,pnhfturb,psed);

// ini variables
    for(int qn=0; qn<2; ++qn)
    {
    pnhfturb->ktimesave(p,d,pgc);
    pnhfturb->etimesave(p,d,pgc);
    }

    p->gctime=0.0;
    p->xtime=0.0;
	p->reinitime=0.0;
	p->wavecalctime=0.0;
	p->field4time=0.0;
}




