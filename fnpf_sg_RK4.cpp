/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
#include"field.h"
#include"convection.h"
#include"ioflow.h"
#include"solver.h"
#include"reini.h"
#include"fnpf_sg_laplace_cds2.h"
#include"fnpf_sg_laplace_cds4.h"
#include"onephase.h"
#include"fnpf_sg_fsfbc.h"

fnpf_sg_RK4::fnpf_sg_RK4(lexer *p, fdm_fnpf *c, ghostcell *pgc) : fnpf_sg_ini(p,c,pgc),fnpf_sigma(p,c,pgc),
                                                      erk1(p),erk2(p),erk3(p),erk(p),en(p),
                                                      frk1(p),frk2(p),frk3(p),frk(p)
{
    gcval=250;
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    if(p->F50==1)
	gcval_eta = 51;
    
    if(p->F50==2)
	gcval_eta = 52;
    
    if(p->F50==3)
	gcval_eta = 53;
    
    if(p->F50==4)
	gcval_eta = 54;
    
    gcval_eta = 50;
    gcval_fifsf = 50;
    
    if(p->A320==1)
    plap = new fnpf_sg_laplace_cds2(p);
    
    if(p->A320==2)
    plap = new fnpf_sg_laplace_cds4;
    
    pf = new fnpf_sg_fsfbc(p,c,pgc);
}

fnpf_sg_RK4::~fnpf_sg_RK4()
{
}

void fnpf_sg_RK4::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, convection *pconvec, ioflow *pflow, reini *preini, onephase* poneph)
{	
// Step 1
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
    pflow->fifsf_relax(p,pgc,frk);
    pgc->gcsl_start4(p,frk1,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->fsfdisc(p,c,pgc,erk,frk);
    sigma_update(p,c,pgc,pf,erk);
    
    // Set Boundary Conditions
    fsfbc_sig(p,c,pgc,frk,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,250);
    plap->start(p,c,pgc,psolv,pf,c->Fi);
    pgc->start7V(p,c->Fi,250);
    pf->fsfwvel(p,c,pgc,erk,frk);
    
    SLICELOOP4
    en(i,j) = erk(i,j);

// Step 2
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
    pflow->fifsf_relax(p,pgc,frk);
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->fsfdisc(p,c,pgc,erk,frk);
    sigma_update(p,c,pgc,pf,erk);
    
    // Set Boundary Conditions
    fsfbc_sig(p,c,pgc,frk,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,250);
    plap->start(p,c,pgc,psolv,pf,c->Fi);
    pgc->start7V(p,c->Fi,250);
    pf->fsfwvel(p,c,pgc,erk,frk);
    
    SLICELOOP4
    en(i,j) = erk(i,j);
    
// Step 3
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
    pflow->fifsf_relax(p,pgc,frk);
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->fsfdisc(p,c,pgc,erk,frk);
    sigma_update(p,c,pgc,pf,erk);
    
    // Set Boundary Conditions
    fsfbc_sig(p,c,pgc,frk,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,250);
    plap->start(p,c,pgc,psolv,pf,c->Fi);
    pgc->start7V(p,c->Fi,250);
    pf->fsfwvel(p,c,pgc,erk,frk);

// Step 4 
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
    pflow->fifsf_relax(p,pgc,c->Fifsf);
    pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);
    sigma_update(p,c,pgc,pf,c->eta);
    
    // Set Boundary Conditions
    fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,250);
    plap->start(p,c,pgc,psolv,pf,c->Fi);
    pgc->start7V(p,c->Fi,250);
    pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);

    LOOP
    c->test(i,j,k) = c->Fifsf(i,j);

    pflow->eta_relax(p,pgc,c->eta);
    pflow->fifsf_relax(p,pgc,c->Fifsf);
    
    velcalc_sig(p,c,pgc,c->Fi);
}

void fnpf_sg_RK4::inidisc(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{	
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    etaloc_sig(p,c,pgc);
    fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);
    sigma_ini(p,c,pgc,pf,c->eta);
    pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);
    sigma_update(p,c,pgc,pf,c->eta);
    
    pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);
    

    LOOP
    c->test(i,j,k) = c->Fz(i,j);
    
    pgc->start4(p,c->test,50);

}

