/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
#include"fnpf_sg_RK4.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"convection.h"
#include"ioflow.h"
#include"solver.h"
#include"reini.h"
#include"fnpf_sg_laplace_cds2.h"
#include"fnpf_sg_laplace_cds4.h"
#include"fnpf_sg_laplace_cds24.h"
#include"onephase.h"
#include"fnpf_sg_fsfbc.h"
#include"fnpf_sg_fsfbc_wd.h"

fnpf_sg_RK4::fnpf_sg_RK4(lexer *p, fdm_fnpf *c, ghostcell *pgc) : fnpf_sg_ini(p,c,pgc),fnpf_sigma(p,c,pgc),
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
    plap = new fnpf_sg_laplace_cds2(p);
    
    if(p->A320==2)
    plap = new fnpf_sg_laplace_cds4(p);
    
    if(p->A320==3)
    plap = new fnpf_sg_laplace_cds24(p);
    
    
    if(p->A343==0)
    pf = new fnpf_sg_fsfbc(p,c,pgc);
    
    if(p->A343>=1)
    pf = new fnpf_sg_fsfbc_wd(p,c,pgc);
}

fnpf_sg_RK4::~fnpf_sg_RK4()
{
}

void fnpf_sg_RK4::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, convection *pconvec, ioflow *pflow, reini *preini, onephase* poneph)
{	
    
// Step 1
    pflow->inflow_fnpf(p,pgc,c->Fi,c->Uin,c->Fifsf,c->eta);
    
    // fsf eta
    pf->kfsfbc(p,c,pgc);

    SLICELOOP4
    {
	erk1(i,j) = p->dt*c->K(i,j);
    erk(i,j) = c->eta(i,j) + 0.5*erk1(i,j);
    }
    
    // fsf Fi
    pf->dfsfbc(p,c,pgc,c->eta);

	SLICELOOP4
    {
	frk1(i,j) = p->dt*c->K(i,j);
    frk(i,j)  = c->Fifsf(i,j) + 0.5*frk1(i,j);
    }
    
    pflow->eta_relax(p,pgc,erk);
    pgc->gcsl_start4(p,erk,gcval_eta);
    pf->coastline(p,c,pgc,erk);
    pflow->fifsf_relax(p,pgc,frk);
    pgc->gcsl_start4(p,frk1,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->breaking(p, c, pgc, erk, c->eta, frk,0.5);
    pf->wetdry(p,c,pgc,erk,frk);
    pf->fsfdisc(p,c,pgc,erk,frk);
    sigma_update(p,c,pgc,pf,erk);
    
    
    // Set Boundary Conditions
    pflow->fivec_relax(p,pgc,c->Fi);
    fsfbc_sig(p,c,pgc,frk,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pgc->start7V(p,c->Fi,c->bc,250);
    plap->start(p,c,pgc,psolv,pf,c->Fi);
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,250);
    pf->fsfwvel(p,c,pgc,erk,frk);
    
    SLICELOOP4
    en(i,j) = erk(i,j);

// Step 2
    pflow->inflow_fnpf(p,pgc,c->Fi,c->Uin,frk,erk);
    
    // fsf eta
    pf->kfsfbc(p,c,pgc);
    
    SLICELOOP4
    {
	erk2(i,j) = p->dt*c->K(i,j);
    erk(i,j)  = c->eta(i,j) + 0.5*erk2(i,j);
    }
    
    // fsf Fi
    pf->dfsfbc(p,c,pgc,en);
    
    SLICELOOP4
    {
	frk2(i,j) = p->dt*c->K(i,j);
    frk(i,j)  = c->Fifsf(i,j) + 0.5*frk2(i,j);
    }
    
    pflow->eta_relax(p,pgc,erk);
    pgc->gcsl_start4(p,erk,gcval_eta);
    pf->coastline(p,c,pgc,erk);
    pflow->fifsf_relax(p,pgc,frk);
    pgc->gcsl_start4(p,frk1,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->breaking(p, c, pgc, erk, c->eta, frk,0.5);
    pf->wetdry(p,c,pgc,erk,frk);
    pf->fsfdisc(p,c,pgc,erk,frk);
    sigma_update(p,c,pgc,pf,erk);
    
    // Set Boundary Conditions
    pflow->fivec_relax(p,pgc,c->Fi);
    fsfbc_sig(p,c,pgc,frk,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pgc->start7V(p,c->Fi,c->bc,250);
    plap->start(p,c,pgc,psolv,pf,c->Fi);
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,250);
    pf->fsfwvel(p,c,pgc,erk,frk);
    
    SLICELOOP4
    en(i,j) = erk(i,j);
    
// Step 3
    pflow->inflow_fnpf(p,pgc,c->Fi,c->Uin,frk,erk);
    
    // fsf eta
    pf->kfsfbc(p,c,pgc);
    
    SLICELOOP4
    {
	erk3(i,j) = p->dt*c->K(i,j);
    erk(i,j)  = c->eta(i,j) + erk3(i,j);
    }
    
    // fsf Fi
    pf->dfsfbc(p,c,pgc,en);
    
    SLICELOOP4
    {
	frk3(i,j) = p->dt*c->K(i,j);
    frk(i,j)  = c->Fifsf(i,j) + frk3(i,j);
    }
    
    pflow->eta_relax(p,pgc,erk);
    pgc->gcsl_start4(p,erk,gcval_eta);
    pf->coastline(p,c,pgc,erk);
    pflow->fifsf_relax(p,pgc,frk);
    pgc->gcsl_start4(p,frk1,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->breaking(p, c, pgc, erk, c->eta, frk,0.5);
    pf->wetdry(p,c,pgc,erk,frk);
    pf->fsfdisc(p,c,pgc,erk,frk);
    sigma_update(p,c,pgc,pf,erk);
    
    // Set Boundary Conditions
    pflow->fivec_relax(p,pgc,c->Fi);
    fsfbc_sig(p,c,pgc,frk,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pgc->start7V(p,c->Fi,c->bc,250);
    plap->start(p,c,pgc,psolv,pf,c->Fi);
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,250);
    pf->fsfwvel(p,c,pgc,erk,frk);

// Step 4 
    pflow->inflow_fnpf(p,pgc,c->Fi,c->Uin,frk,erk);
    
    // fsf eta
    pf->kfsfbc(p,c,pgc);
    
    SLICELOOP4
    c->eta(i,j) = c->eta(i,j) + (1.0/6.0)*(erk1(i,j) + 2.0*erk2(i,j) + 2.0*erk3(i,j) + p->dt*c->K(i,j));
    
    // fsf Fi
    pf->dfsfbc(p,c,pgc,erk);
    
    SLICELOOP4
	c->Fifsf(i,j) = c->Fifsf(i,j) + (1.0/6.0)*(frk1(i,j) + 2.0*frk2(i,j) + 2.0*frk3(i,j) + p->dt*c->K(i,j));
    
    pflow->eta_relax(p,pgc,c->eta);
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pf->coastline(p,c,pgc,c->eta);
    pflow->fifsf_relax(p,pgc,c->Fifsf);
    pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->breaking(p, c, pgc, erk, c->eta, c->Fifsf,0.5);
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);
    pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);
    sigma_update(p,c,pgc,pf,c->eta);
    
    // Set Boundary Conditions
    pflow->fivec_relax(p,pgc,c->Fi);
    fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pgc->start7V(p,c->Fi,c->bc,250);
    plap->start(p,c,pgc,psolv,pf,c->Fi);
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,250);
    pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);

    LOOP
    c->test(i,j,k) = c->Fifsf(i,j);

    velcalc_sig(p,c,pgc,c->Fi);
}

void fnpf_sg_RK4::inidisc(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{	
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    etaloc_sig(p,c,pgc);
    fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);
    sigma_ini(p,c,pgc,pf,c->eta);
    pf->fsfdisc_ini(p,c,pgc,c->eta,c->Fifsf);
    pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);
    sigma_update(p,c,pgc,pf,c->eta);
    
    pgc->start7V(p,c->Fi,c->bc,250);
    pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);
    

    LOOP
    c->test(i,j,k) = c->Fz(i,j);
    
    pgc->start4(p,c->test,50);
    
    bedbc_sig(p,c,pgc,c->Fi,pf);
    velcalc_sig(p,c,pgc,c->Fi);

}

