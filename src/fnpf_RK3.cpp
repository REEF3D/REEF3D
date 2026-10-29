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

#include"fnpf_RK3.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"convection.h"
#include"ioflow.h"
#include"solver.h"
#include"reini.h"
#include"fnpf_laplace_cds2.h"
#include"fnpf_fsfbc.h"
#include"fnpf_fsfbc_wd.h"
#include"fnpf_body.h"
#include"fnpf_amr.h"
#include"fnpf_prof.h"

fnpf_RK3::fnpf_RK3(lexer *p, fdm_fnpf *c, ghostcell *pgc) : fnpf_ini(p,c,pgc),fnpf_sigma(p,c,pgc),
                                                      erk1(p),erk2(p),frk1(p),frk2(p),ek(p),fk(p),en(p)
{
    gcval=250;
    if(p->j_dir==0)
    gcval=150;
   
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    // 3D
    gcval_eta = 55;
    gcval_fifsf = 60;
    
    // 2D
    if(p->j_dir==0)
    {
    gcval_eta = 155;
    gcval_fifsf = 160;
    }
    
    plap0 = new fnpf_laplace_cds2(p);
    plap = plap0;
        
    if(p->A343==0)
    pf = new fnpf_fsfbc(p,c,pgc);
    
    if(p->A343>=1)
    pf = new fnpf_fsfbc_wd(p,c,pgc);
    
    // resolved bodies (X 10): no-op unless present
    pbody = fnpf_body::create(p,c,pgc);
    plap = pbody->laplace(plap);
}

fnpf_RK3::~fnpf_RK3()
{
}

void fnpf_RK3::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, convection *pconvec, ioflow *pflow, reini *preini)
{	   
    SLICELOOP4
    en(i,j) = c->eta(i,j);
    
// Step 1
    // fsf eta
    {fnpf_prof::scope ps_(fnpf_prof::fsfbc); pf->kfsfbc(p,c,pgc);}
    SLICELOOP4
    ek(i,j) = c->K(i,j);
    
    // fsf Fi (dfsfbc only uses the eta passed to it)
    {fnpf_prof::scope ps_(fnpf_prof::fsfbc); pf->dfsfbc(p,c,pgc,c->eta);}
    SLICELOOP4
    fk(i,j) = c->K(i,j);
    
    // mesh refinement: tendencies of the patches (the body loads need them)
    if(pamr!=nullptr)
    pamr->stage_tendency(p,c,pgc,0);
    
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->stage(p,c,pgc,psolv,pf,ek,fk,0);}
    
    SLICELOOP4
	erk1(i,j) = c->eta(i,j) + p->dt*ek(i,j);
    
    SLICELOOP4
	frk1(i,j) = c->Fifsf(i,j) + p->dt*fk(i,j);
    
    {fnpf_prof::scope ps_(fnpf_prof::damping); pf->damping(p,c,pgc,erk1,gcval_eta,1.0);}
    {fnpf_prof::scope ps_(fnpf_prof::damping); pf->damping(p,c,pgc,frk1,gcval_fifsf,1.0);}
   
    // wavegen and coastline
    {fnpf_prof::scope ps_(fnpf_prof::relax); pflow->eta_relax(p,pgc,erk1);}
    {fnpf_prof::scope ps_(fnpf_prof::wetdry); pf->wetdry(p,c,pgc,erk1,frk1);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,erk1,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_eta(p,c,pgc,erk1);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_fi(p,c,pgc,frk1);}
    {fnpf_prof::scope ps_(fnpf_prof::relax); pflow->fifsf_relax(p,pgc,frk1);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,frk1,gcval_fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->surface(p,c,pgc,erk1,frk1,gcval_eta,gcval_fifsf);}
    
    // mesh refinement: stage values of the patches, restriction, cells around the patches
    if(pamr!=nullptr)
    pamr->stage_surface(p,c,pgc,erk1,frk1,0);
    
    // fsfdisc and sigma update
    {fnpf_prof::scope ps_(fnpf_prof::breaking); pf->breaking_zones(p,pflow);}
    {fnpf_prof::scope ps_(fnpf_prof::breaking); pf->breaking(p,c,pgc,erk1,en,frk1,1.0);}
    // A 352: the filter changed eta and Fifsf in place after their halo exchange
    if((p->A350==1 || p->A350==4) && p->A352>0)
    {
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,erk1,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,frk1,gcval_fifsf);}
    }
    {fnpf_prof::scope ps_(fnpf_prof::inflow); pflow->inflow_fnpf(p,c,pgc,c->Fi,c->Uin,frk1,erk1);}
    {fnpf_prof::scope ps_(fnpf_prof::fsfdisc); pf->fsfdisc(p,c,pgc,erk1,frk1);}
    {fnpf_prof::scope ps_(fnpf_prof::sigma); sigma_update(p,c,pgc,pf,erk1);}
  
    // Set Boundary Conditions Fi
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); fsfbc_sig(p,c,pgc,frk1,c->Fi);}
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); bedbc_sig(p,c,pgc,c->Fi,pf);}
    
    // solve Fi
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    {fnpf_prof::scope ps_(fnpf_prof::laplace); plap->start(p,c,pgc,psolv,pf,c->Fi,frk1);}
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    {fnpf_prof::scope ps_(fnpf_prof::fsfwvel); pf->fsfwvel(p,c,pgc,erk1,frk1);}

// Step 2
    // fsf eta
    {fnpf_prof::scope ps_(fnpf_prof::fsfbc); pf->kfsfbc(p,c,pgc);}
    SLICELOOP4
    ek(i,j) = c->K(i,j);
    
    // fsf Fi (dfsfbc only uses the eta passed to it)
    {fnpf_prof::scope ps_(fnpf_prof::fsfbc); pf->dfsfbc(p,c,pgc,erk1);}
    SLICELOOP4
    fk(i,j) = c->K(i,j);
    
    // mesh refinement: tendencies of the patches (the body loads need them)
    if(pamr!=nullptr)
    pamr->stage_tendency(p,c,pgc,1);
    
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->stage(p,c,pgc,psolv,pf,ek,fk,1);}
    
    SLICELOOP4
	erk2(i,j) = 0.75*c->eta(i,j) + 0.25*erk1(i,j) + 0.25*p->dt*ek(i,j);
    
    SLICELOOP4
	frk2(i,j) = 0.75*c->Fifsf(i,j) + 0.25*frk1(i,j) + 0.25*p->dt*fk(i,j);
    
    {fnpf_prof::scope ps_(fnpf_prof::damping); pf->damping(p,c,pgc,erk2,gcval_eta,0.25);}
    {fnpf_prof::scope ps_(fnpf_prof::damping); pf->damping(p,c,pgc,frk2,gcval_fifsf,0.25);}
    
    // wavegen and coastline
    {fnpf_prof::scope ps_(fnpf_prof::relax); pflow->eta_relax(p,pgc,erk2);}
    {fnpf_prof::scope ps_(fnpf_prof::wetdry); pf->wetdry(p,c,pgc,erk2,frk2);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,erk2,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_eta(p,c,pgc,erk2);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_fi(p,c,pgc,frk2);}
    {fnpf_prof::scope ps_(fnpf_prof::relax); pflow->fifsf_relax(p,pgc,frk2);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,frk2,gcval_fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->surface(p,c,pgc,erk2,frk2,gcval_eta,gcval_fifsf);}
    
    if(pamr!=nullptr)
    pamr->stage_surface(p,c,pgc,erk2,frk2,1);
    
    // fsfdisc and sigma update
    {fnpf_prof::scope ps_(fnpf_prof::breaking); pf->breaking(p,c,pgc,erk2,en,frk2,0.5);}
    // A 352: the filter changed eta and Fifsf in place after their halo exchange
    if((p->A350==1 || p->A350==4) && p->A352>0)
    {
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,erk2,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,frk2,gcval_fifsf);}
    }
    {fnpf_prof::scope ps_(fnpf_prof::inflow); pflow->inflow_fnpf(p,c,pgc,c->Fi,c->Uin,frk2,erk2);}
    {fnpf_prof::scope ps_(fnpf_prof::fsfdisc); pf->fsfdisc(p,c,pgc,erk2,frk2);}
    {fnpf_prof::scope ps_(fnpf_prof::sigma); sigma_update(p,c,pgc,pf,erk2);}
    
    // Set Boundary Conditions Fi
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); fsfbc_sig(p,c,pgc,frk2,c->Fi);}
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); bedbc_sig(p,c,pgc,c->Fi,pf);}
    
    // solve Fi
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    {fnpf_prof::scope ps_(fnpf_prof::laplace); plap->start(p,c,pgc,psolv,pf,c->Fi,frk2);}
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    {fnpf_prof::scope ps_(fnpf_prof::fsfwvel); pf->fsfwvel(p,c,pgc,erk2,frk2);}

// Step 3 
    // fsf eta
    {fnpf_prof::scope ps_(fnpf_prof::fsfbc); pf->kfsfbc(p,c,pgc);}
    SLICELOOP4
    ek(i,j) = c->K(i,j);
    
    // fsf Fi (dfsfbc only uses the eta passed to it)
    {fnpf_prof::scope ps_(fnpf_prof::fsfbc); pf->dfsfbc(p,c,pgc,erk2);}
    SLICELOOP4
    fk(i,j) = c->K(i,j);
    
    // mesh refinement: tendencies of the patches (the body loads need them)
    if(pamr!=nullptr)
    pamr->stage_tendency(p,c,pgc,2);
    
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->stage(p,c,pgc,psolv,pf,ek,fk,2);}
    
    SLICELOOP4
	c->eta(i,j) = (1.0/3.0)*c->eta(i,j) + (2.0/3.0)*erk2(i,j) + (2.0/3.0)*p->dt*ek(i,j);
    
    SLICELOOP4
	c->Fifsf(i,j) = (1.0/3.0)*c->Fifsf(i,j) + (2.0/3.0)*frk2(i,j) + (2.0/3.0)*p->dt*fk(i,j);
    
    {fnpf_prof::scope ps_(fnpf_prof::damping); pf->damping(p,c,pgc,c->eta,gcval_eta,2.0/3.0);}
    {fnpf_prof::scope ps_(fnpf_prof::damping); pf->damping(p,c,pgc,c->Fifsf,gcval_fifsf,2.0/3.0);}
    
    // wavegen and coastline
    {fnpf_prof::scope ps_(fnpf_prof::relax); pflow->eta_relax(p,pgc,c->eta);}
    {fnpf_prof::scope ps_(fnpf_prof::wetdry); pf->wetdry(p,c,pgc,c->eta,c->Fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->eta,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_eta(p,c,pgc,c->eta);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_fi(p,c,pgc,c->Fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::relax); pflow->fifsf_relax(p,pgc,c->Fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->surface(p,c,pgc,c->eta,c->Fifsf,gcval_eta,gcval_fifsf);}
    
    if(pamr!=nullptr)
    pamr->stage_surface(p,c,pgc,c->eta,c->Fifsf,2);
    
    // fsfdisc and sigma update
    {fnpf_prof::scope ps_(fnpf_prof::breaking); pf->breaking(p,c,pgc,c->eta,en,c->Fifsf,1.0);}
    // A 352: the filter changed eta and Fifsf in place after their halo exchange
    if((p->A350==1 || p->A350==4) && p->A352>0)
    {
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->eta,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);}
    }
    {fnpf_prof::scope ps_(fnpf_prof::inflow); pflow->inflow_fnpf(p,c,pgc,c->Fi,c->Uin,c->Fifsf,c->eta);}
    {fnpf_prof::scope ps_(fnpf_prof::fsfdisc); pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::sigma); sigma_update(p,c,pgc,pf,c->eta);}
    
    // Set Boundary Conditions Fi
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);}
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); bedbc_sig(p,c,pgc,c->Fi,pf);}
    
    // solve Fi
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    {fnpf_prof::scope ps_(fnpf_prof::laplace); plap->start(p,c,pgc,psolv,pf,c->Fi,c->Fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    {fnpf_prof::scope ps_(fnpf_prof::fsfwvel); pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);}

    //---------------------------------

    {fnpf_prof::scope ps_(fnpf_prof::bc3d); bedbc_sig(p,c,pgc,c->Fi,pf);}
    {fnpf_prof::scope ps_(fnpf_prof::velcalc); velcalc_sig(p,c,pgc,c->Fi);}
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->velocity(p,c,pgc);}
    
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_vel(p,c,pgc,c->U);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_vel(p,c,pgc,c->V);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_vel(p,c,pgc,c->W);}
    
    if(pamr!=nullptr)
    pamr->step_end(p,c,pgc);
}

// mesh refinement (G 1): the Laplace equation is solved on all grids together
void fnpf_RK3::attach_amr(fnpf_amr *a)
{
    pamr = a;
    a->attach_level0(pf,this,this,this);
    a->attach_body(pbody);
    pbody->amr_attach(a);
    plap = a->laplace(plap,plap0);
}

void fnpf_RK3::inidisc_step1(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, solver *psolv)
{
    SLICELOOP4
	c->eta(i,j) = 0.0;
    
    pgc->gcsl_start4(p,c->eta,50);
    
    etaloc_sig(p,c,pgc);
    sigma_ini(p,c,pgc,pf,c->eta);
}

void fnpf_RK3::inidisc_step2(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, solver *psolv)
{	
    for(int qn=0;qn<20;++qn)
    {
    pf->coastline_eta(p,c,pgc,c->eta);
    pf->coastline_fi(p,c,pgc,c->Fifsf);
    pf->coastline_fi_ini(p,c,pgc,c->Fifsf);
    }
    
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
    
    pf->fsfdisc_ini(p,c,pgc,c->eta,c->Fifsf);
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);   
    pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);
    sigma_update(p,c,pgc,pf,c->eta);
    pbody->initialize(p,c,pgc);
    
    fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    pgc->start7V(p,c->Fi,c->bc,gcval);
    
    // velocity 
    if(p->I30>=1)
    {
    pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);
    velcalc_sig(p,c,pgc,c->Fi);
    pbody->velocity(p,c,pgc);
    }
}

void fnpf_RK3::ini_wetdry(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{	
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);   // coastline ini
}

