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

#include"driver.h"
#include"ghostcell.h"
#include"lexer.h"
#include"nhflow_amr.h"
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
#include"river_ini.h"
#include"6DOF_header.h"
#include"vrans_header.h"
#include"nhflow_header.h"
#include"nhflow_thinbody.h"
#include"seastate_nhflow.h"
#include<sys/stat.h>
#include<sys/types.h>

#define WLVL (fabs(a->WL(i,j))>1.0e-20?a->WL(i,j):1.0e20)

void driver::driver_ini_nhflow()
{
    
    pnhf->ini(p,d,pgc,pflow,pnhfsf);

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
    
    // sigma ini
    pnhfmom->inidisc(p,d,pgc,pnhfsf);
    
    LOOP
    {
    d->RO[IJK] = p->W1;
    d->VISC[IJK] = p->W2;
    }
    
    SLICELOOP4
    d->ks(i,j) = p->B50;
    
    // a roughness without bed friction: B 50 alone acts only in the wall functions of a turbulence model
    if(p->mpirank==0 && p->A519==0 && p->B50!=0.0001)
    cout<<"!!! B 50 "<<p->B50<<" is set but A 519 0: no bed friction on the momentum (A 519 1 for the bed)"<<endl;
    
    if(p->S10>0 && p->S28==0)
    SLICELOOP4
    d->ks(i,j) = p->S20;
    
    if(p->S10>0 && p->S28==1)
    SLICELOOP4
    d->ks(i,j) = p->S21*p->S20;
    
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
    
    // steady 1D river flow (F 65 2, river_ini.h): depth-averaged velocities of the wet cells on all layers
    if(p->F65==2)
    {
    river_ini river;
        if(river.start(p,pgc,d->bed,0,1)==1)
        {
        double uc,vc;
        
            LOOP
            {
            river.velocity(p,i,j,d->bed(i,j),uc,vc);
            
            if(p->wet[IJ]==0)
            uc=vc=0.0;
            
            d->U[IJK] = uc;
            d->V[IJK] = vc;
            d->W[IJK] = 0.0;
            d->UH[IJK] = d->WL(i,j)*uc;
            d->VH[IJK] = d->WL(i,j)*vc;
            d->WH[IJK] = 0.0;
            }
        
        pgc->start4V(p,d->U,10);
        pgc->start4V(p,d->V,11);
        pgc->start4V(p,d->W,12);
        pgc->start4V(p,d->UH,14);
        pgc->start4V(p,d->VH,15);
        pgc->start4V(p,d->WH,16);
        }
    }
    
    pflow->discharge_nhflow(p,d,pgc);
    pflow->inflow_nhflow(p,d,pgc,d->U,d->V,d->W,d->UH,d->VH,d->WH,d->WL);
    
    // REEF3D::SEASTATE coupling (A 750 1 in-process, A 750 2 remote waves on their own grid)
    if(p->A750==1 || p->A750==2)
    {
    d->wave = new seastate_nhflow(p,d,pgc);
    d->wave->ini(p,d,pgc);
    }
    
    // turbulence ini
    pnhfturb->ini(p, d, pgc);
    
    // sediment ini
    psed->ini_nhflow(p,d,pgc);
    
    // 6DOF ini
    p6dof->initialize(p, d, pgc);
    
    // sharp thin bodies (membranes X 330 'mobility sharp'): wall fluxes through the flux hook of the convection
    if(d->thinbody!=nullptr)
    pnhfconvec->set_hook(d->thinbody,0);
    
    // VRANS ini
    pnhfvrans->initialize(p,d,pgc);
    
    // mesh refinement ini: patches from the initial state (A 560: the turbulence of level 0, whose fields the blocks
    // and the patches exchange with)
    if(pnhfamr!=nullptr && p->A560!=0)
    pnhfamr->turb_attach(pnhfturb);
    if(pnhfamr!=nullptr)
    pnhfamr->ini(p,d,pgc);

    // the wave gauges and probes of the printer on the finest grid
    if(pnhfamr!=nullptr)
    pprint->amr_columns(p,pnhfamr);
    
    // S 10: the sediment step runs on all grids (nhflow_amr_sed.cpp)
    if(pnhfamr!=nullptr)
    pnhfamr->sed_attach(psed);
    
    // the first output (not on a restart, I 45: the outputs continue)
    if(p->I45!=1)
    {
    pprint->start(p,d,pgc,pflow,pnhfturb,psed);
    
    if(pnhfamr!=nullptr)
    pnhfamr->print(p,d,pgc);
    }

    // Lagrangian particles ini
    pnhfpart->ini(p,d,pgc);

// ini variables
    for(int qn=0; qn<2; ++qn)
    {
    pnhfturb->ktimesave(p,d,pgc);
    pnhfturb->etimesave(p,d,pgc);
    }

    // restart from the latest checkpoint (I 45), on any decomposition
    if(p->I45==1)
    ckpt_read();

    p->gctime=0.0;
    p->xtime=0.0;
	p->reinitime=0.0;
	p->wavecalctime=0.0;
	p->field4time=0.0;
}




