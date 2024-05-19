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

#include"ptf_RK4.h"
#include"lexer.h"
#include"fdm.h"
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

ptf_RK4::ptf_RK4(lexer *p, fdm *a, ghostcell *pgc) : ptf_fsfbc(p,a,pgc),   
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
    plap = new ptf_laplace_cds2(p,a,pgc);
    
    if(p->A320==2)
    plap = new ptf_laplace_cds4;
    
    pfsfupdate = new ptf_fsf_update(p,a,pgc);
    
    pbedupdate = new ptf_bed_update(p,a,pgc);

}

ptf_RK4::~ptf_RK4()
{
}

void ptf_RK4::start(lexer *p, fdm *a, ghostcell *pgc, solver *psolv, convection *pconvec, ioflow *pflow, reini *preini, onephase* poneph)
{	
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);

// Step 1

    // fsf eta
    kfsfbc(p,a,pgc);

    SLICELOOP4
    {
	erk1(i,j) = a->K(i,j);
    erk(i,j)  = a->eta(i,j) + 0.5*p->dt*erk1(i,j);
    }
    
    pgc->gcsl_start4(p,erk,gcval_eta);

    // fsf Fi
    dfsfbc(p,a,pgc,a->eta);

    SLICELOOP4
    {
	frk1(i,j) = a->K(i,j);
    frk(i,j)  = a->Fifsf(i,j) + 0.5*p->dt*frk1(i,j);
    }
    
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    
    // Set Boundary Conditions    
    pflow->eta_relax(p,pgc,erk);
    pflow->fifsf_relax(p,pgc,frk);
    pfsfupdate->fsfupdate(p,a,pgc,pflow,poneph,erk);
    pfsfupdate->etaloc(p,a,pgc);
    pfsfupdate->fsfbc(p,a,pgc,frk,a->Fi);
    pbedupdate->waterdepth(p,a,pgc);
    pbedupdate->bedbc(p,a,pgc,a->Fi);
    
    // fsfdisc
    fsfdisc(p,a,pgc,erk,frk,a->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi,frk);
    pfsfupdate->fsfbc(p,a,pgc,frk,a->Fi);
    pgc->start4(p,a->Fi,gcval);
    fsfwvel(p,a,pgc,erk,frk);
     
// Step 2
    
    // fsf eta
    kfsfbc(p,a,pgc);
    
    SLICELOOP4
    {
	erk2(i,j) = a->K(i,j);
    erk(i,j)  = a->eta(i,j) + 0.5*p->dt*erk2(i,j);
    }
    
    pgc->gcsl_start4(p,erk,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,a,pgc,erk);
    
    SLICELOOP4
    {
	frk2(i,j) = a->K(i,j);
    frk(i,j)  = a->Fifsf(i,j) + 0.5*p->dt*frk2(i,j);
    }
    
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    
    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,erk);
    pflow->fifsf_relax(p,pgc,frk);
    pfsfupdate->fsfupdate(p,a,pgc,pflow,poneph,erk);
    pfsfupdate->etaloc(p,a,pgc);
    pfsfupdate->fsfbc(p,a,pgc,frk,a->Fi);
    pbedupdate->waterdepth(p,a,pgc);
    pbedupdate->bedbc(p,a,pgc,a->Fi);
    
    // fsfdisc
    fsfdisc(p,a,pgc,erk,frk,a->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi,frk);
    pfsfupdate->fsfbc(p,a,pgc,frk,a->Fi);
    pgc->start4(p,a->Fi,gcval);
    fsfwvel(p,a,pgc,erk,frk);
   
// Step 3
    
    // fsf eta
    kfsfbc(p,a,pgc);
    
    SLICELOOP4
    {
	erk3(i,j) = a->K(i,j);
    erk(i,j)  = a->eta(i,j) + p->dt*erk3(i,j);
    }
    
    pgc->gcsl_start4(p,erk,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,a,pgc,erk);
    
    SLICELOOP4
    {
	frk3(i,j) = a->K(i,j);
    frk(i,j)  = a->Fifsf(i,j) + p->dt*frk3(i,j);
    }
    
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    
    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,erk);
    pflow->fifsf_relax(p,pgc,frk);
    pfsfupdate->fsfupdate(p,a,pgc,pflow,poneph,erk);
    pfsfupdate->etaloc(p,a,pgc);
    pfsfupdate->fsfbc(p,a,pgc,frk,a->Fi);
    pbedupdate->waterdepth(p,a,pgc);
    pbedupdate->bedbc(p,a,pgc,a->Fi);
    
    // fsfdisc
    fsfdisc(p,a,pgc,erk,frk,a->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi,frk);
    pfsfupdate->fsfbc(p,a,pgc,frk,a->Fi);
    pgc->start4(p,a->Fi,gcval);
    fsfwvel(p,a,pgc,erk,frk);

// Step 4
    
    // fsf eta
    kfsfbc(p,a,pgc);
    
    SLICELOOP4
    a->eta(i,j) = a->eta(i,j) + (1.0/6.0)*p->dt*(erk1(i,j) + 2.0*erk2(i,j) + 2.0*erk3(i,j) + a->K(i,j));
    
    pgc->gcsl_start4(p,a->eta,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,a,pgc,erk);
    
    SLICELOOP4
	a->Fifsf(i,j) = a->Fifsf(i,j) + (1.0/6.0)*p->dt*(frk1(i,j) + 2.0*frk2(i,j) + 2.0*frk3(i,j) + a->K(i,j));
    
    pgc->gcsl_start4(p,a->Fifsf,gcval_fifsf);

    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,a->eta);
    pflow->fifsf_relax(p,pgc,a->Fifsf);
    pfsfupdate->fsfupdate(p,a,pgc,pflow,poneph,a->eta);
    pfsfupdate->etaloc(p,a,pgc);
    pfsfupdate->fsfbc(p,a,pgc,a->Fifsf,a->Fi);
    pbedupdate->waterdepth(p,a,pgc);
    pbedupdate->bedbc(p,a,pgc,a->Fi);
    
    // fsfdisc
    fsfdisc(p,a,pgc,a->eta,a->Fifsf,a->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi,a->Fifsf);
    pfsfupdate->fsfbc(p,a,pgc,a->Fifsf,a->Fi);
    pgc->start4(p,a->Fi,gcval);
    fsfwvel(p,a,pgc,a->eta,a->Fifsf);
    
    FLUIDLOOP
    a->test(i,j,k) = a->Fifsf(i,j);

    pfsfupdate->velcalc(p,a,pgc,a->Fi);
}

void ptf_RK4::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reini *preini, onephase *poneph)
{	
    pfsfupdate->fsfupdate(p,a,pgc,pflow,poneph,a->eta);
    pfsfupdate->etaloc(p,a,pgc);
    
    pbedupdate->waterdepth(p,a,pgc);
    
    // potential ini
    //pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pflow->fifsf_relax(p,pgc,a->Fifsf);
    pgc->start4(p,a->Fi,250);
    pgc->gcsl_start4(p,a->Fifsf,50);
    
    pfsfupdate->fsfbc(p,a,pgc,a->Fifsf,a->Fi);
    
    pfsfupdate->fsfepol(p,a,pgc,a->eta,a->Fi);
    
    
    pgc->gcsl_start4(p,a->eta,50);
    
    FLUIDLOOP
    a->test(i,j,k) = a->Fifsf(i,j);
    
    pfsfupdate->velcalc(p,a,pgc,a->Fi);
}

void ptf_RK4::inidisc(lexer *p, fdm *a, ghostcell *pgc)
{	
}


