/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
#include"fnpf_RK4.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"convection.h"
#include"ioflow.h"
#include"solver.h"
#include"reini.h"
#include"fnpf_laplace_cds2.h"
#include"fnpf_laplace_cds2_v2.h"
#include"fnpf_laplace_cds4.h"
#include"fnpf_laplace_cds4_bc2.h"
#include"onephase.h"
#include"fnpf_fsfbc.h"
#include"fnpf_fsfbc_wd.h"

fnpf_RK4::fnpf_RK4(lexer *p, fdm_fnpf *c, ghostcell *pgc) : fnpf_ini(p,c,pgc),fnpf_sigma(p,c,pgc),
                                                      erk1(p),erk2(p),erk3(p),erk(p),en(p),
                                                      frk1(p),frk2(p),frk3(p),frk(p)
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
    
    if(p->A320==1)
    plap = new fnpf_laplace_cds2(p);
    
    if(p->A320==2)
    plap = new fnpf_laplace_cds4(p);
    
    if(p->A320==3)
    plap = new fnpf_laplace_cds4_bc2(p);
    
    if(p->A320==5)
    plap = new fnpf_laplace_cds2_v2(p,pgc);
    
    
    if(p->A343==0)
    pf = new fnpf_fsfbc(p,c,pgc);
    
    if(p->A343>=1)
    pf = new fnpf_fsfbc_wd(p,c,pgc);
}

fnpf_RK4::~fnpf_RK4()
{
}

void fnpf_RK4::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, convection *pconvec, ioflow *pflow, reini *preini, onephase* poneph)
{	
    
// Step 1
    // fsf eta
    pf->kfsfbc(p,c,pgc);
    pf->damping(p,c,pgc,c->eta,gcval_eta,0.5);

    SLICELOOP4
    {
	erk1(i,j) = c->K(i,j);
    erk(i,j) = c->eta(i,j) + 0.5*p->dt*c->K(i,j);
    }
    
    // fsf Fi
    pf->dfsfbc(p,c,pgc,c->eta);
    pf->damping(p,c,pgc,c->Fifsf,gcval_fifsf,0.5);

	SLICELOOP4
    {
	frk1(i,j) = c->K(i,j);
    frk(i,j)  = c->Fifsf(i,j) + 0.5*p->dt*c->K(i,j);
    }
    
    pflow->eta_relax(p,pgc,erk);
    pgc->gcsl_start4(p,erk,gcval_eta);
    pf->coastline_eta(p,c,pgc,erk);
    pf->coastline_fi(p,c,pgc,frk);
    pflow->fifsf_relax(p,pgc,frk);
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->breaking(p, c, pgc, en, c->eta, frk,0.5);
    pf->wetdry(p,c,pgc,erk,frk);
    pflow->inflow_fnpf(p,c,pgc,c->Fi,c->Uin,c->Fifsf,c->eta);
    pf->fsfdisc(p,c,pgc,erk,frk);
    sigma_update(p,c,pgc,pf,erk);
    
    // Set Boundary Conditions
    pflow->fivec_relax(p,pgc,c->Fi);
    fsfbc_sig(p,c,pgc,frk,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pgc->start7V(p,c->Fi,c->bc,gcval);
    plap->start(p,c,pgc,psolv,pf,c->Fi,frk1);
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    pf->fsfwvel(p,c,pgc,erk,frk);
    
    SLICELOOP4
    en(i,j) = erk(i,j);

// Step 2
    // fsf eta
    pf->kfsfbc(p,c,pgc);
    pf->damping(p,c,pgc,erk,gcval_eta,0.5);
    
    SLICELOOP4
    {
	erk2(i,j) = c->K(i,j);
    erk(i,j)  = c->eta(i,j) + 0.5*p->dt*c->K(i,j);
    }
    
    // fsf Fi
    pf->dfsfbc(p,c,pgc,en);
    pf->damping(p,c,pgc,frk,gcval_fifsf,0.5);
    
    SLICELOOP4
    {
	frk2(i,j) = c->K(i,j);
    frk(i,j)  = c->Fifsf(i,j) + 0.5*p->dt*c->K(i,j);
    }
    
    pflow->eta_relax(p,pgc,erk);
    pgc->gcsl_start4(p,erk,gcval_eta);
    pf->coastline_eta(p,c,pgc,erk);
    pf->coastline_fi(p,c,pgc,frk);
    pflow->fifsf_relax(p,pgc,frk);
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->breaking(p, c, pgc, en, c->eta, frk,0.5);
    pf->wetdry(p,c,pgc,erk,frk);
    pflow->inflow_fnpf(p,c,pgc,c->Fi,c->Uin,frk,erk);
    pf->fsfdisc(p,c,pgc,erk,frk);
    sigma_update(p,c,pgc,pf,erk);
    
    // Set Boundary Conditions
    pflow->fivec_relax(p,pgc,c->Fi);
    fsfbc_sig(p,c,pgc,frk,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pgc->start7V(p,c->Fi,c->bc,gcval);
    plap->start(p,c,pgc,psolv,pf,c->Fi,frk2);
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    pf->fsfwvel(p,c,pgc,erk,frk);
    
    SLICELOOP4
    en(i,j) = erk(i,j);
    
// Step 3
    // fsf eta
    pf->kfsfbc(p,c,pgc);
    pf->damping(p,c,pgc,erk,gcval_eta,1.0);
    
    SLICELOOP4
    {
	erk3(i,j) = c->K(i,j);
    erk(i,j)  = c->eta(i,j) + p->dt*c->K(i,j);
    }
    
    // fsf Fi
    pf->dfsfbc(p,c,pgc,en);
    pf->damping(p,c,pgc,frk,gcval_fifsf,1.0);
    
    SLICELOOP4
    {
	frk3(i,j) = c->K(i,j);
    frk(i,j)  = c->Fifsf(i,j) + p->dt*c->K(i,j);
    }
    
    pflow->eta_relax(p,pgc,erk);
    pgc->gcsl_start4(p,erk,gcval_eta);
    pf->coastline_eta(p,c,pgc,erk);
    pf->coastline_fi(p,c,pgc,frk);
    pflow->fifsf_relax(p,pgc,frk);
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->breaking(p,c,pgc,en,c->eta,frk,1.0);
    pf->wetdry(p,c,pgc,erk,frk);
    pflow->inflow_fnpf(p,c,pgc,c->Fi,c->Uin,frk,erk);
    pf->fsfdisc(p,c,pgc,erk,frk);
    sigma_update(p,c,pgc,pf,erk);
    
    // Set Boundary Conditions
    pflow->fivec_relax(p,pgc,c->Fi);
    fsfbc_sig(p,c,pgc,frk,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pgc->start7V(p,c->Fi,c->bc,gcval);
    plap->start(p,c,pgc,psolv,pf,c->Fi,frk3);
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    pf->fsfwvel(p,c,pgc,erk,frk);
    
    SLICELOOP4
    en(i,j) = erk(i,j);

// Step 4 
    // fsf eta
    pf->kfsfbc(p,c,pgc);
    pf->damping(p,c,pgc,erk,gcval_eta,1.0);
    
    SLICELOOP4
    c->eta(i,j) = c->eta(i,j) + p->dt*(1.0/6.0)*(erk1(i,j) + 2.0*erk2(i,j) + 2.0*erk3(i,j) + c->K(i,j));
    
    // fsf Fi
    pf->dfsfbc(p,c,pgc,en);
    pf->damping(p,c,pgc,frk,gcval_fifsf,1.0);
    
    SLICELOOP4
	c->Fifsf(i,j) = c->Fifsf(i,j) + p->dt*(1.0/6.0)*(frk1(i,j) + 2.0*frk2(i,j) + 2.0*frk3(i,j) + c->K(i,j));
    
    pflow->eta_relax(p,pgc,c->eta);
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pf->coastline_eta(p,c,pgc,c->eta);
    pf->coastline_fi(p,c,pgc,c->Fifsf);
    pflow->fifsf_relax(p,pgc,c->Fifsf);
    pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->breaking(p, c, pgc, en, c->eta, c->Fifsf,1.0);
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);
    pflow->inflow_fnpf(p,c,pgc,c->Fi,c->Uin,frk,erk);
    pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);
    sigma_update(p,c,pgc,pf,c->eta);
    
    // Set Boundary Conditions
    pflow->fivec_relax(p,pgc,c->Fi);
    fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pgc->start7V(p,c->Fi,c->bc,gcval);
    plap->start(p,c,pgc,psolv,pf,c->Fi,c->Fifsf);
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);


    bedbc_sig(p,c,pgc,c->Fi,pf);
    velcalc_sig(p,c,pgc,c->Fi);
}

void fnpf_RK4::inidisc(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, solver *psolv)
{	
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    etaloc_sig(p,c,pgc);
    fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);
    sigma_ini(p,c,pgc,pf,c->eta);
    pf->fsfdisc_ini(p,c,pgc,c->eta,c->Fifsf);
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);   // coastline ini
    pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);
    sigma_update(p,c,pgc,pf,c->eta);
    
    pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);

    
    for(int qn=0; qn<10; ++qn)
    {
    pf->coastline_eta(p,c,pgc,c->eta);
    pf->coastline_fi(p,c,pgc,c->Fifsf);
    }
    
    
    velcalc_sig(p,c,pgc,c->Fi);
    
    pgc->start7V(p,c->U,c->bc,210);
    pgc->start7V(p,c->V,c->bc,210);
    pgc->start7V(p,c->W,c->bc,210);
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
    
    
    if(p->I40==1)
    {
    fnpf_restart(p,c,pgc);
    
    
    sigma_update(p,c,pgc,pf,c->eta);
    
  
    for(int qn=0;qn<0;++qn)
    {
    fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pgc->start7V(p,c->Fi,c->bc,gcval);
    plap->start(p,c,pgc,psolv,pf,c->Fi,c->Fifsf);
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);
    }
    }
}

void fnpf_RK4::ini_wetdry(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{	
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);   // coastline ini

    pf->coastline_eta(p,c,pgc,c->eta);
    pf->coastline_fi(p,c,pgc,c->Fifsf);
}

