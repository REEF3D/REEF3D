/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"sflow_momentum_RK2_SC.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sflow_convection.h"
#include"sflow_pressure.h"
#include"sflow_diffusion.h"
#include"sflow_fsf.h"
#include"sflow_forcing.h"
#include"ioflow.h"
#include"solver2D.h"
#include"sflow_rough_manning.h"
#include"sflow_rough_void.h"
#include"sflow_rheology_f.h"
#include"sflow_rheology_v.h"
#include"6DOF.h"

sflow_momentum_RK2_SC::sflow_momentum_RK2_SC(lexer *p, fdm2D *b, sflow_convection *pconvection, sflow_diffusion *ppdiff, sflow_pressure* ppressure,
                                                    solver2D *psolver, solver2D *ppoissonsolver, ioflow *pioflow, sflow_fsf *pfreesurf,sflow_forcing *ppsfdf,
                                                    sixdof *pp6dof)
                                                    :UHRK1(p),VHRK1(p),WHRK1(p),WLRK1(p)
{
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

	pconvec=pconvection;
	pdiff=ppdiff;
	ppress=ppressure;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
	pfsf=pfreesurf;
    p6dof=pp6dof;
    psfdf=ppsfdf;
    
    if(p->A218==0)
    prough = new sflow_rough_void(p);
    
    if(p->A218==1)
    prough = new sflow_rough_manning(p);
    
    
    if(p->W90==0)
    prheo = new sflow_rheology_v(p);
    
    if(p->W90==1)
    prheo = new sflow_rheology_f(p);
}

sflow_momentum_RK2_SC::~sflow_momentum_RK2_SC()
{
}

void sflow_momentum_RK2_SC::start(lexer *p, fdm2D* b, ghostcell* pgc)
{	/*        
    pflow->discharge2D(p,b,pgc);
    pflow->inflow2D(p,b,pgc,b->P,b->Q,b->bed,b->eta);
    pflow->inflow2D(p,b,pgc,UHRK1,VHRK1,b->bed,b->eta);
	
//Step 1
//--------------------------------------------------------	
    reconstruct(p,b,pgc,pfsf,pss,precon,b->WL,b->U,b->V,b->W,b->UH,b->VH,b->WH);

    // FSF
    starttime=pgc->timer();
    pconvec->start(p,b,4,b->WL);

    pfsf->rk2_step1(p, d, pgc, pflow, b->UH, b->VH, b->WH, WLRK1, WLRK1, 1.0);
    p->fsftime+=pgc->timer()-starttime;

    // U
	starttime=pgc->timer();
	pflow->isource2D(p,b,pgc); 
	ppress->upgrad(p,b,etark1,b->eta);
	irhs(p,b,pgc,b->P,1.0);
    //prough->u_source(p,b,b->P);
    //prheo->u_source(p,b,b->P,b->Q);
    //p6dof->isource2D(p,b,pgc);
	pconvec->start(p,b,1,WLRK1);
	//pdiff->diff_u(p,b,pgc,psolv,b->UH,b->VH,1.0);
                
    SLICELOOP4
	UHRK1[IJK] = b->UH(i,j)
				+ p->dt*CPORNH*b->F(i,j);
	
    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pflow->jsource2D(p,b,pgc);
	ppress->vpgrad(p,b,etark1,b->eta);
	jrhs(p,b,pgc,b->Q,1.0);
    //prough->v_source(p,b,b->Q);
    //prheo->v_source(p,b,b->P,b->Q);
    //p6dof->jsource2D(p,b,pgc);
	pconvec->start(p,b,2,WLRK1);
	//pdiff->diff_v(p,b,pgc,psolv,b->P,b->Q,1.0);

	SLICELOOP4
	UHRK1(i,j) = b->V(i,j)
			  + p->dt*b->G(i,j);
	
    p->vtime=pgc->timer()-starttime;
	
    // W
    ppress->wpgrad(p,b,etark1,b->eta);
    pconvec->start(p,b,3,d->WLRK1);
    //pdiff->diff_w(p,b,pgc,psolv,b->P,b->Q,b->ws,1.0);
    
    SLICELOOP4
	WHRK1(i,j) = b->WH(i,j)
			  + p->dt*b->L(i,j);
              
    pgc->gcsl_start4(p,WHRK1,12);
    
    psfdf->forcing(p,b,pgc,p6dof,0,1.0,UHRK1,VHRK1,WHRK1,etark1,b->hp,0);
    
    velcalc(p,b,pgc,UHRK1,VHRK1,WHRK1,WLRK1);

	// press
    ppress->start(p,b,pgc,ppoissonsolv,pflow, UHRK1, VHRK1, b->P, b->Q, WHRK1, etark1, 1.0);
    velcalc(p,b,pgc,UHRK1,VHRK1,WHRK1,WLRK1);

	pflow->pm_relax(p,pgc,b->press);

	pflow->um_relax(p,pgc,UHRK1,b->bed,b->eta);
	pflow->vm_relax(p,pgc,VHRK1,b->bed,b->eta);
    pflow->wm_relax(p,pgc,WHRK1,b->bed,b->eta);

	pgc->gcsl_start1(p,UHRK1,gcval_u);
	pgc->gcsl_start2(p,VHRK1,gcval_v);
    pgc->gcsl_start4(p,WHRK1,gcval_w);
        
//Step 2
//--------------------------------------------------------
    
    reconstruct(p,b,pgc,pfsf,pss,precon,WLRK1,d->U,d->V,d->W,UHRK1,VHRK1,WHRK1);
    
    // FSF
    starttime=pgc->timer();
    
    pconvec->start(p,b,4,WLRK1);
    pfsf->rk2_step2(p, d, pgc, pflow, UHRK1,VHRK1,WHRK1, WLRK1, WLRK1, 0.5);
    
    p->fsftime+=pgc->timer()-starttime;
    
    
	// U
	starttime=pgc->timer();

	pflow->isource2D(p,b,pgc);
	ppress->upgrad(p,b,b->eta,etark1);
	irhs(p,b,pgc,UHRK1,0.5);
    //prough->u_source(p,b,UHRK1);
    //prheo->u_source(p,b,UHRK1,VHRK1);
    //p6dof->isource2D(p,b,pgc);
	pconvec->start(p,b,1,d->WL);
	//pdiff->diff_u(p,b,pgc,psolv,UHRK1,VHRK1,0.5);

	SLICELOOP4
	b->UH(i,j) = 0.5*b->UH(i,j) + 0.5*UHRK1(i,j)
				+ 0.5*p->dt*b->F(i,j);
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pflow->jsource2D(p,b,pgc);
	ppress->vpgrad(p,b,b->eta,etark1);
	jrhs(p,b,pgc,VHRK1,0.5);
    //prough->v_source(p,b,VHRK1);
    //prheo->v_source(p,b,UHRK1,VHRK1);
    //p6dof->jsource2D(p,b,pgc);
	pconvec->start(p,b,2,d->WL);
	//pdiff->diff_v(p,b,pgc,psolv,UHRK1,VHRK1,0.5);

	SLICELOOP4
	b->VH(i,j) = 0.5*b->VH(i,j) + 0.5*VHRK1(i,j)
			  + 0.5*p->dt*b->G(i,j);
	
    p->vtime+=pgc->timer()-starttime;
	
    // W     
    ppress->wpgrad(p,b,b->eta,etark1);
    pconvec->start(p,b,3,d->WL);
    //pdiff->diff_w(p,b,pgc,psolv,UHRK1,VHRK1,WHRK1,0.5);
    
    SLICELOOP4
	b->WH(i,j) = 0.5*b->WH(i,j) + 0.5*WHRK1(i,j)
			  + 0.5*p->dt*b->L(i,j);
              
    pgc->gcsl_start4(p,b->ws,12);
    
    psfdf->forcing(p,b,pgc,p6dof,1,0.5,b->P,b->Q,b->ws,b->eta,b->hp,1);
              
	//--------------------------------------------------------
	// pressure
	ppress->start(p,b,pgc,ppoissonsolv,pflow, b->P, b->Q, UHRK1, VHRK1, b->ws, b->eta, 0.5);
	
    pflow->pm_relax(p,pgc,b->press);
    
	pflow->um_relax(p,pgc,b->P,b->bed,b->eta);
	pflow->vm_relax(p,pgc,b->Q,b->bed,b->eta);
    pflow->wm_relax(p,pgc,b->ws,b->bed,b->eta);

	pgc->gcsl_start1(p,b->P,gcval_u);
	pgc->gcsl_start2(p,b->Q,gcval_v);
    pgc->gcsl_start4(p,b->ws,gcval_w);

    pfsf->breaking_persist(p,b,pgc,b->eta,b->eta_n,1.0);
    
    SLICELOOP4
    b->eta_n(i,j) = b->eta(i,j);

    pgc->gcsl_start4(p,b->eta_n,gcval_eta);*/
}

void sflow_momentum_RK2_SC::irhs(lexer *p, fdm2D *b, ghostcell *pgc, slice &f, double alpha)
{
	n=0;
	if(p->D20<4)
	SLICELOOP1
	{
	b->maxF=MAX(fabs(b->F(i,j)),b->maxF);
	b->rhsvec.V[n]=0.0;
	++n;
	}
}

void sflow_momentum_RK2_SC::jrhs(lexer *p, fdm2D *b, ghostcell *pgc, slice &f, double alpha)
{
	n=0;
	if(p->D20<4)
	SLICELOOP2
	{
	b->maxG=MAX(fabs(b->G(i,j)),b->maxG);
	b->rhsvec.V[n]=0.0;
	++n;
	}
}

