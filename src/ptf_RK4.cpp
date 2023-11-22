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

#include"ptf_RK4.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"
#include"field4.h"
#include"convection.h"
#include"ioflow.h"
#include"solver.h"
#include"reini.h"
#include"ptf_laplace_cds2.h"
#include"ptf_laplace_cds4.h"
#include"onephase.h"
#include"ptf_fsf_update.h"
#include"ptf_bed_update.h"

ptf_RK4::ptf_RK4(lexer *p, fdm_ptf *e, ghostcell *pgc) : ptf_fsfbc(p,e,pgc),   
                                                      erk1(p),erk2(p),erk3(p),erk(p),
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
    plap = new ptf_laplace_cds2(p,e,pgc);
    
    if(p->A320==2)
    plap = new ptf_laplace_cds4;
    
    pfsfupdate = new ptf_fsf_update(p,e,pgc);
    
    pbedupdate = new ptf_bed_update(p,e,pgc);

}

ptf_RK4::~ptf_RK4()
{
}

void ptf_RK4::start(lexer *p, fdm_ptf *e, ghostcell *pgc, solver_ptf *psolv, convection *pconvec, ioflow *pflow, reini *preini, onephase* poneph_ptf)
{	
    //pflow->inflow(p,e,pgc,e->u,e->v,e->w);

// Step 1

    // fsf eta
    kfsfbc(p,e,pgc);

    SLICELOOP4
    {
	erk1(i,j) = e->K(i,j);
    erk(i,j)  = e->eta(i,j) + 0.5*p->dt*erk1(i,j);
    }
    
    pgc->gcsl_start4(p,erk,gcval_eta);

    // fsf Fi
    dfsfbc(p,e,pgc,e->eta);

    SLICELOOP4
    {
	frk1(i,j) = e->K(i,j);
    frk(i,j)  = e->Fifsf(i,j) + 0.5*p->dt*frk1(i,j);
    }
    
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    
    // Set Boundary Conditions    
    pflow->eta_relax(p,pgc,erk);
    pflow->fifsf_relax(p,pgc,frk);
    pfsfupdate->fsfupdate(p,e,pgc,pflow,poneph_ptf,erk);
    pfsfupdate->etaloc(p,e,pgc);
    pfsfupdate->fsfbc(p,e,pgc,frk,e->Fi,erk);
    pbedupdate->waterdepth(p,e,pgc);
    pbedupdate->bedbc(p,e,pgc,e->Fi);
    
    // fsfdisc
    fsfdisc(p,e,pgc,erk,frk,e->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,e->Fi,e->phi);
    pgc->start4(p,e->Fi,gcval);
    plap->start(p,e,pgc,psolv,e->Fi,frk,erk);
    pfsfupdate->fsfbc(p,e,pgc,frk,e->Fi,erk);
    pgc->start4(p,e->Fi,gcval);
    fsfwvel(p,e,pgc,erk,frk);
     
// Step 2
    
    // fsf eta
    kfsfbc(p,e,pgc);
    
    SLICELOOP4
    {
	erk2(i,j) = e->K(i,j);
    erk(i,j)  = e->eta(i,j) + 0.5*p->dt*erk2(i,j);
    }
    
    pgc->gcsl_start4(p,erk,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,e,pgc,erk);
    
    SLICELOOP4
    {
	frk2(i,j) = e->K(i,j);
    frk(i,j)  = e->Fifsf(i,j) + 0.5*p->dt*frk2(i,j);
    }
    
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    
    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,erk);
    pflow->fifsf_relax(p,pgc,frk);
    pfsfupdate->fsfupdate(p,e,pgc,pflow,poneph_ptf,erk);
    pfsfupdate->etaloc(p,e,pgc);
    pfsfupdate->fsfbc(p,e,pgc,frk,e->Fi,erk);
    pbedupdate->waterdepth(p,e,pgc);
    pbedupdate->bedbc(p,e,pgc,e->Fi);
    
    // fsfdisc
    fsfdisc(p,e,pgc,erk,frk,e->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,e->Fi,e->phi);
    pgc->start4(p,e->Fi,gcval);
    plap->start(p,e,pgc,psolv,e->Fi,frk,erk);
    pfsfupdate->fsfbc(p,e,pgc,frk,e->Fi,erk);
    pgc->start4(p,e->Fi,gcval);
    fsfwvel(p,e,pgc,erk,frk);
   
// Step 3
    
    // fsf eta
    kfsfbc(p,e,pgc);
    
    SLICELOOP4
    {
	erk3(i,j) = e->K(i,j);
    erk(i,j)  = e->eta(i,j) + p->dt*erk3(i,j);
    }
    
    pgc->gcsl_start4(p,erk,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,e,pgc,erk);
    
    SLICELOOP4
    {
	frk3(i,j) = e->K(i,j);
    frk(i,j)  = e->Fifsf(i,j) + p->dt*frk3(i,j);
    }
    
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    
    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,erk);
    pflow->fifsf_relax(p,pgc,frk);
    pfsfupdate->fsfupdate(p,e,pgc,pflow,poneph_ptf,erk);
    pfsfupdate->etaloc(p,e,pgc);
    pfsfupdate->fsfbc(p,e,pgc,frk,e->Fi,erk);
    pbedupdate->waterdepth(p,e,pgc);
    pbedupdate->bedbc(p,e,pgc,e->Fi);
    
    // fsfdisc
    fsfdisc(p,e,pgc,erk,frk,e->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,e->Fi,e->phi);
    pgc->start4(p,e->Fi,gcval);
    plap->start(p,e,pgc,psolv,e->Fi,frk,erk);
    pfsfupdate->fsfbc(p,e,pgc,frk,e->Fi,erk);
    pgc->start4(p,e->Fi,gcval);
    fsfwvel(p,e,pgc,erk,frk);

// Step 4
    
    // fsf eta
    kfsfbc(p,e,pgc);
    
    SLICELOOP4
    e->eta(i,j) = e->eta(i,j) + (1.0/6.0)*p->dt*(erk1(i,j) + 2.0*erk2(i,j) + 2.0*erk3(i,j) + e->K(i,j));
    
    pgc->gcsl_start4(p,e->eta,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,e,pgc,erk);
    
    SLICELOOP4
	e->Fifsf(i,j) = e->Fifsf(i,j) + (1.0/6.0)*p->dt*(frk1(i,j) + 2.0*frk2(i,j) + 2.0*frk3(i,j) + e->K(i,j));
    
    pgc->gcsl_start4(p,e->Fifsf,gcval_fifsf);

    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,e->eta);
    pflow->fifsf_relax(p,pgc,e->Fifsf);
    pfsfupdate->fsfupdate(p,e,pgc,pflow,poneph_ptf,e->eta);
    pfsfupdate->etaloc(p,e,pgc);
    pfsfupdate->fsfbc(p,e,pgc,e->Fifsf,e->Fi,e->eta);
    pbedupdate->waterdepth(p,e,pgc);
    pbedupdate->bedbc(p,e,pgc,e->Fi);
    
    // fsfdisc
    fsfdisc(p,e,pgc,e->eta,e->Fifsf,e->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,e->Fi,e->phi);
    pgc->start4(p,e->Fi,gcval);
    plap->start(p,e,pgc,psolv,e->Fi,e->Fifsf,e->eta);
    pfsfupdate->fsfbc(p,e,pgc,e->Fifsf,e->Fi,e->eta);
    pgc->start4(p,e->Fi,gcval);
    fsfwvel(p,e,pgc,e->eta,e->Fifsf);
    
    FLUIDLOOP
    e->test(i,j,k) = e->Fifsf(i,j);

    pfsfupdate->velcalc(p,e,pgc,e->Fi);
}

void ptf_RK4::ini(lexer *p, fdm_ptf *e, ghostcell *pgc, ioflow *pflow, reini *preini, onephase *poneph_ptf)
{	
    pfsfupdate->fsfupdate(p,e,pgc,pflow,poneph_ptf,e->eta);
    pfsfupdate->etaloc(p,e,pgc);
    
    pbedupdate->waterdepth(p,e,pgc);
    
    // potential ini
    //pflow->fi_relax(p,pgc,e->Fi,e->phi);
    pflow->fifsf_relax(p,pgc,e->Fifsf);
    pgc->start4(p,e->Fi,250);
    pgc->gcsl_start4(p,e->Fifsf,50);
    
    pfsfupdate->fsfbc(p,e,pgc,e->Fifsf,e->Fi,e->eta);
    
    pfsfupdate->fsfepol(p,e,pgc,e->eta,e->Fi);
    
    
    pgc->gcsl_start4(p,e->eta,50);
    
    FLUIDLOOP
    e->test(i,j,k) = e->Fifsf(i,j);
    
    pfsfupdate->velcalc(p,e,pgc,e->Fi);
}

void ptf_RK4::inidisc(lexer *p, fdm_ptf *e, ghostcell *pgc)
{	
}


