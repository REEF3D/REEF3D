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

#include"ptf_RK3.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"
#include"field4.h"
#include"convection.h"
#include"ioflow.h"
#include"solver_ptf.h"
#include"reini.h"
#include"ptf_laplace_cds2.h"
#include"ptf_laplace_cds4.h"
#include"onephase_ptf.h"
#include"onephase_ptf_f.h"
#include"onephase_ptf_v.h"
#include"ptf_fsf_update.h"
#include"ptf_bed_update.h"
#include"ptf_breaking.h"

ptf_RK3::ptf_RK3(lexer *p, fdm_ptf *e, ghostcell *pgc) : ptf_fsfbc(p,e,pgc),erk1(p),erk2(p),frk1(p),frk2(p),Firk_0(p),Firk_1(p),Firk_2(p),Firk_3(p)
{
    gcval=250;
    
    if(p->F50==1)
	gcval_eta = 51;
    
    if(p->F50==2)
	gcval_eta = 52;
    
    if(p->F50==3)
	gcval_eta = 53;
    
    if(p->F50==4)
	gcval_eta = 54;
    
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
    plap = new ptf_laplace_cds2(p,e,pgc);
    
    if(p->A320==2)
    plap = new ptf_laplace_cds4;
    
    pfsfupdate = new ptf_fsf_update(p,e,pgc);
    
    pbedupdate = new ptf_bed_update(p,e,pgc);
    
}

ptf_RK3::~ptf_RK3()
{
}

void ptf_RK3::start(lexer *p, fdm_ptf *e, ghostcell *pgc, solver_ptf *psolv, convection *pconvec, ioflow *pflow, reini *preini, onephase_ptf* poneph_ptf)
{	
   // pflow->inflow(p,e,pgc,e->u,e->v,e->w);
    
    
// Step
    Firk_0=e->Fi;
    // fsf eta
    kfsfbc(p,e,pgc);

    SLICELOOP4
	erk1(i,j) = e->eta(i,j) + p->dt*e->K(i,j);
    
    pgc->gcsl_start4(p,erk1,gcval_eta);

    // fsf Fi
    dfsfbc(p,e,pgc,e->eta);

    SLICELOOP4
	frk1(i,j) = e->Fifsf(i,j) + p->dt*e->K(i,j);
    
    pgc->gcsl_start4(p,frk1,gcval_fifsf);

    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,erk1);
    pflow->fifsf_relax(p,pgc,frk1);
    if(p->A350>0)
    //breaking(p, a, pgc, erk1, e->eta, frk1,1.0);
    pfsfupdate->fsfupdate(p,e,pgc,pflow,poneph_ptf,erk1);
    pfsfupdate->etaloc(p,e,pgc);
    pfsfupdate->fsfbc(p,e,pgc,frk1,e->Fi,erk1);
    pbedupdate->waterdepth(p,e,pgc);
    pbedupdate->bedbc(p,e,pgc,e->Fi);
    
    // fsfdisc
    fsfdisc(p,e,pgc,erk1,frk1,e->Fi);
    
    // solve Fi
    // pflow->fivec_relax(p,pgc,e->Fi);
    pflow->fi_relax(p,pgc,e->Fi,e->phi);
    pgc->start4(p,e->Fi,gcval);
    plap->start(p,e,pgc,psolv,e->Fi,frk1,erk1);
    pfsfupdate->fsfbc(p,e,pgc,frk1,e->Fi,erk1);
    pgc->start4(p,e->Fi,gcval);
    fsfwvel(p,e,pgc,erk1,frk1);
    
// Step 2
    
    Firk_1=e->Fi;
    // fsf eta
    kfsfbc(p,e,pgc);
    
    SLICELOOP4
	erk2(i,j) = 0.75*e->eta(i,j) + 0.25*erk1(i,j) + 0.25*p->dt*e->K(i,j);
    
    pgc->gcsl_start4(p,erk2,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,e,pgc,erk1);
    
    SLICELOOP4
	frk2(i,j) = 0.75*e->Fifsf(i,j) + 0.25*frk1(i,j) + 0.25*p->dt*e->K(i,j);
    
    pgc->gcsl_start4(p,frk2,gcval_fifsf);
    
    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,erk2);
    pflow->fifsf_relax(p,pgc,frk2);
    if(p->A350>0)
    //breaking(p, a, pgc, erk2, erk1, frk2, 0.25);
    pfsfupdate->fsfupdate(p,e,pgc,pflow,poneph_ptf,erk2);
    pfsfupdate->etaloc(p,e,pgc);
    pfsfupdate->fsfbc(p,e,pgc,frk2,e->Fi,erk2);
    pbedupdate->waterdepth(p,e,pgc);
    pbedupdate->bedbc(p,e,pgc,e->Fi);
    
    // fsfdisc
    fsfdisc(p,e,pgc,erk2,frk2,e->Fi);
    
    // solve Fi
    // pflow->fivec_relax(p,pgc,e->Fi);
    pflow->fi_relax(p,pgc,e->Fi,e->phi);
    pgc->start4(p,e->Fi,gcval);
    plap->start(p,e,pgc,psolv,e->Fi,frk2,erk2);
   pfsfupdate->fsfbc(p,e,pgc,frk2,e->Fi,erk2);
   pgc->start4(p,e->Fi,gcval);
    fsfwvel(p,e,pgc,erk2,frk2);
    
// Step 3  
    
    Firk_2=e->Fi;
    // fsf eta
    kfsfbc(p,e,pgc);
    
    SLICELOOP4
	e->eta(i,j) = (1.0/3.0)*e->eta(i,j) + (2.0/3.0)*erk2(i,j) + (2.0/3.0)*p->dt*e->K(i,j);
    
    pgc->gcsl_start4(p,e->eta,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,e,pgc,erk2);
    
    SLICELOOP4
	e->Fifsf(i,j) = (1.0/3.0)*e->Fifsf(i,j) + (2.0/3.0)*frk2(i,j) + (2.0/3.0)*p->dt*e->K(i,j);
    
    pgc->gcsl_start4(p,e->Fifsf,gcval_fifsf);
    
    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,e->eta);
    pflow->fifsf_relax(p,pgc,e->Fifsf);
    if(p->A350>0)
    //breaking(p, a, pgc, e->eta, erk2, e->Fifsf, 2.0/3.0);
    pfsfupdate->fsfupdate(p,e,pgc,pflow,poneph_ptf,e->eta);
    pfsfupdate->etaloc(p,e,pgc);
    pfsfupdate->fsfbc(p,e,pgc,e->Fifsf,e->Fi,e->eta);
    pbedupdate->waterdepth(p,e,pgc);
    pbedupdate->bedbc(p,e,pgc,e->Fi);
    
    // fsfdisc
    fsfdisc(p,e,pgc,e->eta,e->Fifsf,e->Fi);
    
    // solve Fi
   // pflow->fivec_relax(p,pgc,e->Fi);
    pflow->fi_relax(p,pgc,e->Fi,e->phi);
    pgc->start4(p,e->Fi,gcval);
    plap->start(p,e,pgc,psolv,e->Fi,e->Fifsf,e->eta);
    pfsfupdate->fsfbc(p,e,pgc,e->Fifsf,e->Fi,e->eta);
    pgc->start4(p,e->Fi,gcval);
    fsfwvel(p,e,pgc,e->eta,e->Fifsf);
    Firk_3=e->Fi;
    
    FLUIDLOOP
    {
    if(p->flag4[IJK]<0)
    e->test(i,j,k) = 0.0;
    
    if(p->flag4[IJK]>0)
    e->test(i,j,k) = 1.0;
    
    }
    
    pfsfupdate->velcalc(p,e,pgc,e->Fi);
    
    if(p->A401>0)
        pfsfupdate->presscalc(p,e,pgc,Firk_0,Firk_1,Firk_2,Firk_3,e->eta);
    
    
}

void ptf_RK3::ini(lexer *p, fdm_ptf *e, ghostcell *pgc, ioflow *pflow, reini *preini, onephase_ptf *poneph_ptf)
{	
    pfsfupdate->fsfupdate(p,e,pgc,pflow,poneph_ptf,e->eta);
    pfsfupdate->etaloc(p,e,pgc);
    
    pbedupdate->waterdepth(p,e,pgc);
    
    // potential ini
    //pflow->fivec_relax(p,pgc,e->Fi);
    pflow->fi_relax(p,pgc,e->Fi,e->phi);
    pflow->fifsf_relax(p,pgc,e->Fifsf);
    pgc->start4(p,e->Fi,250);
    pgc->gcsl_start4(p,e->Fifsf,50);
    
    pfsfupdate->fsfbc(p,e,pgc,e->Fifsf,e->Fi,e->eta);
    
    pfsfupdate->fsfepol(p,e,pgc,e->eta,e->Fi);
    
    
    pgc->gcsl_start4(p,e->eta,50);
    
    FLUIDLOOP
    e->test(i,j,k) = e->Fifsf(i,j);
    
    pfsfupdate->velcalc(p,e,pgc,e->Fi);
    if(p->A401>0)
        pfsfupdate->presscalc(p,e,pgc,Firk_0,Firk_1,Firk_2,Firk_3,e->eta);
}

void ptf_RK3::inidisc(lexer *p, fdm_ptf *e, ghostcell *pgc)
{	
}
