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
#include"iowave.h"
#include"solver.h"
#include"reini.h"
#include"ptf_laplace_cds2.h"
#include"ptf_laplace_cds4.h"
#include"onephase.h"
#include"ptf_fsf_update.h"
#include"ptf_bed_update.h"

ptf_RK3::ptf_RK3(lexer *p, fdm_ptf *a, ghostcell *pgc) : ptf_fsfbc(p,a,pgc),erk1(p),erk2(p),frk1(p),frk2(p)
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
    plap = new ptf_laplace_cds2(p,a,pgc);
    
    if(p->A320==2)
    plap = new ptf_laplace_cds4;
    
    pfsfupdate = new ptf_fsf_update(p,a,pgc);
    
    pbedupdate = new ptf_bed_update(p,a,pgc);
}

ptf_RK3::~ptf_RK3()
{
}

void ptf_RK3::start(lexer *p, fdm_ptf *a, ghostcell *pgc, solver *psolv, convection *pconvec, iowave *pwave, reini *preini, onephase* poneph)
{	
    //pwave->inflow(p,a,pgc,a->u,a->v,a->w);
   // pwave->inflow_ptf(p,a,pgc,a->Fi_,a->Uin_,a->Fifsf,a->eta);
// Step 1
    
    // fsf eta
    kfsfbc(p,a,pgc);

    SLICELOOP4
	erk1(i,j) = a->eta(i,j) + p->dt*a->K(i,j);
    
    pgc->gcsl_start4(p,erk1,gcval_eta);

    // fsf Fi
    dfsfbc(p,a,pgc,a->eta);

    SLICELOOP4
	frk1(i,j) = a->Fifsf(i,j) + p->dt*a->K(i,j);
    
    pgc->gcsl_start4(p,frk1,gcval_fifsf);

    // Set Boundary Conditions
   // pwave->fifsf_relax(p,pgc,frk1);
   cout<<"up to fiexchange works"<<endl;
   pgc->start4(p,a->Fi,50);
   FLUIDLOOP
        a->Fi_[IJK]=a->Fi(i,j,k);
    pgc->start4V(p,a->Fi_,50);
    cout<<"fiexchange  works"<<endl;
    pwave->inflow_ptf(p,a,pgc,a->Fi_,a->Uin_,frk1,erk1);
    pgc->start4V(p,a->Fi_,50);
    cout<<"inflow works"<<endl;
    FLUIDLOOP
        a->Fi(i,j,k)=a->Fi_[IJK];
    pgc->start4(p,a->Fi,50);
    cout<<"fiexchange 2 works"<<endl;
    cout<<"1 works"<<endl;
    pfsfupdate->fsfupdate(p,a,pgc,pwave,poneph,erk1);
    pfsfupdate->etaloc(p,a,pgc);
    pfsfupdate->fsfbc(p,a,pgc,frk1,a->Fi,erk1);
    pbedupdate->waterdepth(p,a,pgc);
    pbedupdate->bedbc(p,a,pgc,a->Fi);
    
    // fsfdisc
    fsfdisc(p,a,pgc,erk1,frk1,a->Fi);
    
    // solve Fi
   // pwave->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi,frk1,erk1);
    pfsfupdate->fsfbc(p,a,pgc,frk1,a->Fi,erk1);
    pgc->start4(p,a->Fi,gcval);
    fsfwvel(p,a,pgc,erk1,frk1);

// Step 2
    
    // fsf eta
    kfsfbc(p,a,pgc);
    
    SLICELOOP4
	erk2(i,j) = 0.75*a->eta(i,j) + 0.25*erk1(i,j) + 0.25*p->dt*a->K(i,j);
    
    pgc->gcsl_start4(p,erk2,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,a,pgc,erk1);
    
    SLICELOOP4
	frk2(i,j) = 0.75*a->Fifsf(i,j) + 0.25*frk1(i,j) + 0.25*p->dt*a->K(i,j);
    
    pgc->gcsl_start4(p,frk2,gcval_fifsf);
    
    // Set Boundary Conditions
  //  pwave->eta_relax(p,pgc,erk2);
  //  pwave->fifsf_relax(p,pgc,frk2);
    pfsfupdate->fsfupdate(p,a,pgc,pwave,poneph,erk2);
    pfsfupdate->etaloc(p,a,pgc);
    cout<<"2 test"<<endl;
    FLUIDLOOP
        a->Fi_[IJK]=a->Fi(i,j,k);
    pwave->inflow_ptf(p,a,pgc,a->Fi_,a->Uin_,frk2,erk2);
    FLUIDLOOP
        a->Fi(i,j,k)=a->Fi_[IJK];
    cout<<"2 works"<<endl;
    pfsfupdate->fsfbc(p,a,pgc,frk2,a->Fi,erk2);
    pbedupdate->waterdepth(p,a,pgc);
    pbedupdate->bedbc(p,a,pgc,a->Fi);
    
    // fsfdisc
    fsfdisc(p,a,pgc,erk2,frk2,a->Fi);
    
    // solve Fi
  //  pwave->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi,frk2,erk2);
    pfsfupdate->fsfbc(p,a,pgc,frk2,a->Fi,erk2);
    pgc->start4(p,a->Fi,gcval);
    fsfwvel(p,a,pgc,erk2,frk2);

// Step 3  

    // fsf eta
    kfsfbc(p,a,pgc);
    
    SLICELOOP4
	a->eta(i,j) = (1.0/3.0)*a->eta(i,j) + (2.0/3.0)*erk2(i,j) + (2.0/3.0)*p->dt*a->K(i,j);
    
    pgc->gcsl_start4(p,a->eta,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,a,pgc,erk2);
    
    SLICELOOP4
	a->Fifsf(i,j) = (1.0/3.0)*a->Fifsf(i,j) + (2.0/3.0)*frk2(i,j) + (2.0/3.0)*p->dt*a->K(i,j);
    
    pgc->gcsl_start4(p,a->Fifsf,gcval_fifsf);
    
    // Set Boundary Conditions
  //  pwave->eta_relax(p,pgc,a->eta);
  //  pwave->fifsf_relax(p,pgc,a->Fifsf);
    pfsfupdate->fsfupdate(p,a,pgc,pwave,poneph,a->eta);
    pfsfupdate->etaloc(p,a,pgc);
    FLUIDLOOP
        a->Fi_[IJK]=a->Fi(i,j,k);
    pwave->inflow_ptf(p,a,pgc,a->Fi_,a->Uin_,a->Fifsf,a->eta);
    FLUIDLOOP
        a->Fi(i,j,k)=a->Fi_[IJK];
    pfsfupdate->fsfbc(p,a,pgc,a->Fifsf,a->Fi,a->eta);
    pbedupdate->waterdepth(p,a,pgc);
    pbedupdate->bedbc(p,a,pgc,a->Fi);
    
    // fsfdisc
    fsfdisc(p,a,pgc,a->eta,a->Fifsf,a->Fi);
    
    // solve Fi
  //  pwave->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi,a->Fifsf,a->eta);
    pfsfupdate->fsfbc(p,a,pgc,a->Fifsf,a->Fi,a->eta);
    pgc->start4(p,a->Fi,gcval);
    fsfwvel(p,a,pgc,a->eta,a->Fifsf);
    
    FLUIDLOOP
    {
    if(p->flag4[IJK]<0)
    a->test(i,j,k) = 0.0;
    
    if(p->flag4[IJK]>0)
    a->test(i,j,k) = 1.0;
    
    }
    
    pfsfupdate->velcalc(p,a,pgc,a->Fi);
}

void ptf_RK3::ini(lexer *p, fdm_ptf *a, ghostcell *pgc, iowave *pwave, reini *preini, onephase *poneph)
{	
    pfsfupdate->fsfupdate(p,a,pgc,pwave,poneph,a->eta);
    pfsfupdate->etaloc(p,a,pgc);
    
    pbedupdate->waterdepth(p,a,pgc);
    
    // potential ini
    //pwave->fi_relax(p,pgc,a->Fi,a->phi);
  //  pwave->fifsf_relax(p,pgc,a->Fifsf);
    pgc->start4(p,a->Fi,250);
    pgc->gcsl_start4(p,a->Fifsf,50);
    
    pfsfupdate->fsfbc(p,a,pgc,a->Fifsf,a->Fi,a->eta);
    
    pfsfupdate->fsfepol(p,a,pgc,a->eta,a->Fi);
    
    
    pgc->gcsl_start4(p,a->eta,50);
    
    FLUIDLOOP
    a->test(i,j,k) = a->Fifsf(i,j);
    
    pfsfupdate->velcalc(p,a,pgc,a->Fi);
    cout<<"ini works"<<endl;
}

void ptf_RK3::inidisc(lexer *p, fdm_ptf *a, ghostcell *pgc)
{	
}
