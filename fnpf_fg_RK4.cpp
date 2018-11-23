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

#include"fnpf_fg_RK4.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field4.h"
#include"convection.h"
#include"ioflow.h"
#include"solver.h"
#include"reini.h"
#include"fnpf_fg_laplace_cds2.h"
#include"fnpf_fg_laplace_cds4.h"
#include"onephase.h"

fnpf_fg_RK4::fnpf_fg_RK4(lexer *p, fdm *a, ghostcell *pgc) : fnpf_fg_fsfbc(p,a,pgc), fnpf_fg_fsf_update(p,a,pgc),
                                                      fnpf_fg_bed_update(p,a,pgc),  
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
    plap = new fnpf_fg_laplace_cds2;
    
    if(p->A320==2)
    plap = new fnpf_fg_laplace_cds4;

}

fnpf_fg_RK4::~fnpf_fg_RK4()
{
}

void fnpf_fg_RK4::start(lexer *p, fdm *a, ghostcell *pgc, solver *psolv, convection *pconvec, ioflow *pflow, reini *preini, onephase* poneph)
{	

// Step 1
    fsfdisc(p,a,pgc,a->eta,a->Fifsf,a->Fi);
    
    // fsf eta
    kfsfbc(p,a,pgc);

    SLICELOOP4
    {
	erk1(i,j) = p->dt*a->K(i,j);
    erk(i,j)  = a->eta(i,j) + 0.5*erk1(i,j);
    }
    
    pgc->gcsl_start4(p,erk,gcval_eta);

    // fsf Fi
    dfsfbc(p,a,pgc,a->eta);

    SLICELOOP4
    {
	frk1(i,j) = p->dt*a->K(i,j);
    frk(i,j)  = a->Fifsf(i,j) + 0.5*frk1(i,j);
    }
    
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    

    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,erk);
    pflow->fifsf_relax(p,pgc,frk);
    fsfupdate(p,a,pgc,pflow,poneph,erk);
    etaloc(p,a,pgc);
    fsfbc(p,a,pgc,frk,a->Fi);
    waterdepth(p,a,pgc);
    bedbc(p,a,pgc,a->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi);

     
// Step 2

    fsfdisc(p,a,pgc,erk,frk,a->Fi);
    
    // fsf eta
    kfsfbc(p,a,pgc);
    
    SLICELOOP4
    {
	erk2(i,j) = p->dt*a->K(i,j);
    erk(i,j)  = a->eta(i,j) + 0.5*erk2(i,j);
    }
    
    pgc->gcsl_start4(p,erk,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,a,pgc,erk);
    
    SLICELOOP4
    {
	frk2(i,j) = p->dt*a->K(i,j);
    frk(i,j)  = a->Fifsf(i,j) + 0.5*frk2(i,j);
    }
    
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    
    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,erk);
    pflow->fifsf_relax(p,pgc,frk);
    fsfupdate(p,a,pgc,pflow,poneph,erk);
    etaloc(p,a,pgc);
    fsfbc(p,a,pgc,frk,a->Fi);
    waterdepth(p,a,pgc);
    bedbc(p,a,pgc,a->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi);
   
// Step 3
    
    fsfdisc(p,a,pgc,erk,frk,a->Fi);
    
    // fsf eta
    kfsfbc(p,a,pgc);
    
    SLICELOOP4
    {
	erk3(i,j) = p->dt*a->K(i,j);
    erk(i,j)  = a->eta(i,j) + erk3(i,j);
    }
    
    pgc->gcsl_start4(p,erk,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,a,pgc,erk);
    
    SLICELOOP4
    {
	frk3(i,j) = p->dt*a->K(i,j);
    frk(i,j)  = a->Fifsf(i,j) + frk3(i,j);
    }
    
    pgc->gcsl_start4(p,frk,gcval_fifsf);
    
    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,erk);
    pflow->fifsf_relax(p,pgc,frk);
    fsfupdate(p,a,pgc,pflow,poneph,erk);
    etaloc(p,a,pgc);
    fsfbc(p,a,pgc,frk,a->Fi);
    waterdepth(p,a,pgc);
    bedbc(p,a,pgc,a->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi);

// Step 4
    
    fsfdisc(p,a,pgc,erk,frk,a->Fi);
    
    // fsf eta
    kfsfbc(p,a,pgc);
    
    SLICELOOP4
    a->eta(i,j) = a->eta(i,j) + (1.0/6.0)*(erk1(i,j) + 2.0*erk2(i,j) + 2.0*erk3(i,j) + p->dt*a->K(i,j));
    
    pgc->gcsl_start4(p,a->eta,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,a,pgc,erk);
    
    SLICELOOP4
	a->Fifsf(i,j) = a->Fifsf(i,j) + (1.0/6.0)*(frk1(i,j) + 2.0*frk2(i,j) + 2.0*frk3(i,j) + p->dt*a->K(i,j));
    
    pgc->gcsl_start4(p,a->Fifsf,gcval_fifsf);

    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,a->eta);
    pflow->fifsf_relax(p,pgc,a->Fifsf);
    fsfupdate(p,a,pgc,pflow,poneph,a->eta);
    etaloc(p,a,pgc);
    fsfbc(p,a,pgc,a->Fifsf,a->Fi);
    waterdepth(p,a,pgc);
    bedbc(p,a,pgc,a->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi);
    
    FLUIDLOOP
    a->test(i,j,k) = a->Fifsf(i,j);

    velcalc(p,a,pgc,a->Fi);
    
}

void fnpf_fg_RK4::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reini *preini, onephase *poneph)
{	
    fsfupdate(p,a,pgc,pflow,poneph,a->eta);
    etaloc(p,a,pgc);
    
    // potential ini
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pflow->fifsf_relax(p,pgc,a->Fifsf);
    pgc->start4(p,a->Fi,250);
    
    
    fsfbc(p,a,pgc,a->Fifsf,a->Fi);
    
    fsfepol(p,a,pgc,a->eta,a->Fi);
    
    pgc->gcsl_start4(p,a->Fifsf,50);
    pgc->gcsl_start4(p,a->eta,50);
}

void fnpf_fg_RK4::inidisc(lexer *p, fdm *a, ghostcell *pgc)
{	
}


