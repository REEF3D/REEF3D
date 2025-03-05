/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MEFCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Fabian Knoblauch
--------------------------------------------------------------------*/

#include"momentum_FC3_PLIC.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"convection.h"
#include"diffusion.h"
#include"pressure.h"
#include"poisson.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"fluid_update_vof.h"
#include"VOF_PLIC.h"
#include"nhflow.h"
#include"heat.h"
#include"concentration.h"

momentum_FC3_PLIC::momentum_FC3_PLIC(lexer *p, fdm *a, ghostcell *pgc, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, ioflow *pioflow,
                                                    heat *&pheat, concentration *&pconc,
                                                    fsi *ppfsi, reini *ppreini)
                                                    :momentum_forcing(p),bcmom(p),udiff(p),vdiff(p),wdiff(p),urk1(p),urk2(p),vrk1(p),
                                                    vrk2(p),wrk1(p),wrk2(p),ls(p),frk1(p),frk2(p),fx(p),fy(p),fz(p),lu(p),lv(p),lw(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_vof=1;
    gcval_phi=1;

	pconvec=pconvection;
	pdiff=pdiffusion;
	ppress=ppressure;
	ppois=ppoisson;
	pturb=pturbulence;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
    pfsi=ppfsi;
    pupdate=new fluid_update_vof(p,a,pgc);
    pplic= new VOF_PLIC(p,a,pgc,pheat);
    preini=ppreini;
    
}

momentum_FC3_PLIC::~momentum_FC3_PLIC()
{
}

void momentum_FC3_PLIC::start(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof, vector<net*>& pnet)
{	
    
    LOOP
    {
	a->L(i,j,k)=0.0;
    ls(i,j,k)=a->vof(i,j,k);
    }
    
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,lu,lv,lw);
	pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
	pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);
    
		
//Step 1
//--------------------------------------------------------

    // FSF
    
    pplic->RKcalcL(a,p,pgc,a->u,a->v,a->w);
    
	LOOP
    {
        frk1(i,j,k) = ls(i,j,k) + a->L(i,j,k);
        
        if(frk1(i,j,k)<0.0)
            frk1(i,j,k)=0.0;
        if(frk1(i,j,k)>1.0)
            frk1(i,j,k)=1.0;
    }
    pgc->start4(p,frk1,gcval_vof);
   // pflow->vof_relax(p,a,pgc,frk1);
	pgc->start4(p,frk1,gcval_vof);
    
    LOOP
        a->vof(i,j,k) = frk1(i,j,k);
    pgc->start4(p,a->vof,gcval_vof);
    pplic->RK_redistance(a,p,pgc);
    pgc->start4(p,a->phi,gcval_phi);
    
    preini->start(a,p,a->phi,pgc,pflow);
    pgc->start4(p,a->phi,gcval_phi);
    preini->start(a,p,a->phi,pgc,pflow);
    pgc->start4(p,a->phi,gcval_phi);
    preini->start(a,p,a->phi,pgc,pflow);
    pgc->start4(p,a->phi,gcval_phi);
    
    pupdate->start(p,a,pgc);

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans); 
	bcmom_start(a,p,pgc,pturb,lu,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,lu,lu,lv,lw,1.0);
	pconvec->start(p,a,lu,1,lu,lv,lw);
	pdiff->diff_u(p,a,pgc,psolv,udiff,lu,lu,lv,lw,1.0);

	ULOOP
	urk1(i,j,k) = udiff(i,j,k)
				+ p->dt*CPOR1*a->F(i,j,k);

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,lv,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,lv,lu,lv,lw,1.0);
	pconvec->start(p,a,lv,2,lu,lv,lw);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,lv,lu,lv,lw,1.0);

	VLOOP
	vrk1(i,j,k) = vdiff(i,j,k)
				+ p->dt*CPOR2*a->G(i,j,k);

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,lw,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,lw,lu,lv,lw,1.0);
	pconvec->start(p,a,lw,3,lu,lv,lw);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,lw,lu,lv,lw,1.0);

	WLOOP
	wrk1(i,j,k) = wdiff(i,j,k)
				+ p->dt*CPOR3*a->H(i,j,k);
	
    p->wtime=pgc->timer()-starttime;
    
    pgc->start1(p,urk1,gcval_u);
	pgc->start2(p,vrk1,gcval_v);
    pgc->start3(p,wrk1,gcval_w);
    
    momentum_forcing_start(a, p, pgc, p6dof, pvrans, pnet, pfsi,
                           urk1, vrk1, wrk1, fx, fy, fz, 0, 1.0, false);
                           
    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk1, vrk1, wrk1, 1.0);
	
	pflow->u_relax(p,a,pgc,urk1);
	pflow->v_relax(p,a,pgc,vrk1);
	pflow->w_relax(p,a,pgc,wrk1);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,urk1,gcval_u);
	pgc->start2(p,vrk1,gcval_v);
	pgc->start3(p,wrk1,gcval_w);
    
    pupdate->start(p,a,pgc);
	
//Step 2
//--------------------------------------------------------
	
    // FSF
    LOOP
        a->L(i,j,k)=0.0;

	pplic->RKcalcL(a,p,pgc,urk1,vrk1,wrk1);

	LOOP
    {
        frk2(i,j,k) = 0.75*ls(i,j,k)
                + 0.25*frk1(i,j,k)
                + 0.25*a->L(i,j,k);
                
        if(frk1(i,j,k)<0.0)
            frk1(i,j,k)=0.0;
        if(frk1(i,j,k)>1.0)
            frk1(i,j,k)=1.0;
	}			
	pgc->start4(p,frk2,gcval_vof);
    //pflow->vof_relax(p,a,pgc,frk2);
	pgc->start4(p,frk2,gcval_vof);
    
    LOOP
        a->vof(i,j,k) =  frk2(i,j,k);
        
    pgc->start4(p,a->vof,gcval_vof);
    pplic->RK_redistance(a,p,pgc);
    pgc->start4(p,a->phi,gcval_phi);
    
    preini->start(a,p,a->phi,pgc,pflow);
    pgc->start4(p,a->phi,gcval_phi);
    preini->start(a,p,a->phi,pgc,pflow);
    pgc->start4(p,a->phi,gcval_phi);
    preini->start(a,p,a->phi,pgc,pflow);
    pgc->start4(p,a->phi,gcval_phi);
    
    pupdate->start(p,a,pgc);
    
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,lu,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,urk1,urk1,vrk1,wrk1,0.25);
	pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
	pdiff->diff_u(p,a,pgc,psolv,udiff,urk1,urk1,vrk1,wrk1,1.0);

	ULOOP
	urk2(i,j,k) = 0.75*lu(i,j,k) + 0.25*udiff(i,j,k)
				+ 0.25*p->dt*CPOR1*a->F(i,j,k);
                
    p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,lv,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,vrk1,urk1,vrk1,wrk1,0.25);
	pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,vrk1,urk1,vrk1,wrk1,1.0);

	VLOOP
	vrk2(i,j,k) = 0.75*lv(i,j,k) + 0.25*vdiff(i,j,k)
				+ 0.25*p->dt*CPOR2*a->G(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,lw,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,wrk1,urk1,vrk1,wrk1,0.25);
	pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,wrk1,urk1,vrk1,wrk1,1.0);

	WLOOP
	wrk2(i,j,k) = 0.75*lw(i,j,k) + 0.25*wdiff(i,j,k)
				+ 0.25*p->dt*CPOR3*a->H(i,j,k);

    p->wtime+=pgc->timer()-starttime;
    
    pgc->start1(p,urk2,gcval_u);
	pgc->start2(p,vrk2,gcval_v);
    pgc->start3(p,wrk2,gcval_w);
    
    momentum_forcing_start(a, p, pgc, p6dof, pvrans, pnet, pfsi,
                           urk2, vrk2, wrk2, fx, fy, fz, 1, 0.25, false);
                           
    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk2, vrk2, wrk2, 0.25);
	
	pflow->u_relax(p,a,pgc,urk2);
	pflow->v_relax(p,a,pgc,vrk2);
	pflow->w_relax(p,a,pgc,wrk2);
	pflow->p_relax(p,a,pgc,a->press);
	
	pgc->start1(p,urk2,gcval_u);
	pgc->start2(p,vrk2,gcval_v);
	pgc->start3(p,wrk2,gcval_w);

    pupdate->start(p,a,pgc);

//Step 3
//--------------------------------------------------------
    
    // FSF
    LOOP
        a->L(i,j,k)=0.0;

	pplic->RKcalcL(a,p,pgc,urk2,vrk2,wrk2);

	LOOP
    {
        ls(i,j,k) =   (1.0/3.0)*ls(i,j,k)
                + (2.0/3.0)*frk2(i,j,k)
                + (2.0/3.0)*a->L(i,j,k);
        
        if(ls(i,j,k)<0.0)
            ls(i,j,k)=0.0;
        if(ls(i,j,k)>1.0)
            ls(i,j,k)=1.0;
    }

    pgc->start4(p,ls,gcval_vof);
   // pflow->vof_relax(p,a,pgc,ls);
	pgc->start4(p,ls,gcval_vof);
    
    LOOP
        a->vof(i,j,k) =  ls(i,j,k);
    
    pgc->start4(p,a->vof,gcval_vof);
    pplic->RK_redistance(a,p,pgc);
    pgc->start4(p,a->phi,gcval_phi);
    
    preini->start(a,p,a->phi,pgc,pflow);
    pgc->start4(p,a->phi,gcval_phi);
    preini->start(a,p,a->phi,pgc,pflow);
    pgc->start4(p,a->phi,gcval_phi);
    preini->start(a,p,a->phi,pgc,pflow);
    pgc->start4(p,a->phi,gcval_phi);
    
    pupdate->start(p,a,pgc);
    
    
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,lu,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,urk2,urk2,vrk2,wrk2,2.0/3.0);
	pconvec->start(p,a,urk2,1,urk2,vrk2,wrk2);
	pdiff->diff_u(p,a,pgc,psolv,udiff,urk2,urk2,vrk2,wrk2,1.0);

	ULOOP
	a->u(i,j,k) = (1.0/3.0)*lu(i,j,k) + (2.0/3.0)*udiff(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR1*a->F(i,j,k);
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,lv,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,vrk2,urk2,vrk2,wrk2,2.0/3.0);
	pconvec->start(p,a,vrk2,2,urk2,vrk2,wrk2);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,vrk2,urk2,vrk2,wrk2,1.0);

	VLOOP
	a->v(i,j,k) = (1.0/3.0)*lv(i,j,k) + (2.0/3.0)*vdiff(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR2*a->G(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,lw,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,wrk2,urk2,vrk2,wrk2,2.0/3.0);
	pconvec->start(p,a,wrk2,3,urk2,vrk2,wrk2);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,wrk2,urk2,vrk2,wrk2,1.0);

	WLOOP
	a->w(i,j,k) = (1.0/3.0)*lw(i,j,k) + (2.0/3.0)*wdiff(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR3*a->H(i,j,k);
	
    p->wtime+=pgc->timer()-starttime;
    
    pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
    
    momentum_forcing_start(a, p, pgc, p6dof, pvrans, pnet, pfsi,
                           a->u, a->v, a->w, fx, fy, fz, 2, 2.0/3.0, true);

	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, a->u, a->v,a->w,2.0/3.0);
	
	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
    
    pupdate->start(p,a,pgc);
   
}

void momentum_FC3_PLIC::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	ULOOP
	{
    a->maxF=MAX(fabs(a->rhsvec.V[n] + a->gi),a->maxF);
	a->F(i,j,k) += (a->rhsvec.V[n] + a->gi + p->W29_x)*PORVAL1;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void momentum_FC3_PLIC::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	VLOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n] + a->gj),a->maxG);
	a->G(i,j,k) += (a->rhsvec.V[n] + a->gj + p->W29_y)*PORVAL2;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void momentum_FC3_PLIC::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	WLOOP
	{
    a->maxH=MAX(fabs(a->rhsvec.V[n] + a->gk),a->maxH);
	a->H(i,j,k) += (a->rhsvec.V[n] + a->gk + p->W29_z)*PORVAL3;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void momentum_FC3_PLIC::utimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FC3_PLIC::vtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FC3_PLIC::wtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}
