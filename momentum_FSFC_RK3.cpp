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
for more detaia->phi.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"momentum_FSFC_RK3.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"discrete.h"
#include"diffusion.h"
#include"pressure.h"
#include"poisson.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"reini.h"
#include"picard.h"
#include"fluid_update_fsf.h"
#include"fluid_update_fsf_heat.h"
#include"fluid_update_fsf_comp.h"
#include"fluid_update_void.h"
#include"fluid_update_fsf_concentration.h"
#include"fluid_update_fsf_entrain.h"
#include"fluid_update_rheology.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"heat.h"
#include"concentration.h"

momentum_FSFC_RK3::momentum_FSFC_RK3(lexer *p, fdm *a, ghostcell *pgc, discrete *pdiscrete, diffusion *pdiffusion, pressure* ppressure, 
                                                    poisson* ppoisson, turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, 
                                                    ioflow *pioflow, discrete *ppfsfdisc, reini *ppreini, heat *&pheat, concentration *&pconc)
                                                    :bcmom(p),urk1(p),urk2(p),vrk1(p),vrk2(p),wrk1(p),wrk2(p),ark1(p),ark2(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
	
	gcval_urk=20;
	gcval_vrk=21;
	gcval_wrk=22;

	pdisc=pdiscrete;
	pdiff=pdiffusion;
	ppress=ppressure;
	ppois=ppoisson;
	pturb=pturbulence;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
    pfsfdisc=ppfsfdisc;
    preini=ppreini;
    
    
    // fsf
    if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;

	if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90==0)
	pupdate = new fluid_update_fsf(p,a,pgc);
	
	if(p->F30>0 && p->H10==0 && p->W30==1 && p->F300==0 && p->W90==0)
	pupdate = new fluid_update_fsf_comp(p,a,pgc);
	
	if(p->F30>0 && p->H10>0 && p->F300==0 && p->W90==0)
	pupdate = new fluid_update_fsf_heat(p,a,pgc,pheat);
	
	if(p->F30>0 && p->C10>0 && p->F300==0 && p->W90==0)
	pupdate = new fluid_update_fsf_concentration(p,a,pgc,pconc);
	
	if(p->F30>0 && p->F101>0 && p->F300==0 && p->W90==0)
	pupdate = new fluid_update_fsf_entrain(p,a,pgc,pconc);
	
	if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90>0)
	pupdate = new fluid_update_rheology(p,a,pgc);
	
	if(p->F300>0)
	pupdate = new fluid_update_void();

	if(p->F46==2)
	ppicard = new picard_f(p);

	if(p->F46==3)
	ppicard = new picard_lsm(p);

	if(p->F46!=2 && p->F46!=3)
	ppicard = new picard_void(p);
    
    //cout<<p->mpirank<<" FSFC 004"<<endl;
}

momentum_FSFC_RK3::~momentum_FSFC_RK3()
{
}

void momentum_FSFC_RK3::start(lexer *p, fdm* a, ghostcell* pgc, momentum *pmom)
{	
	
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
	pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);
    
    pflow->fsfinflow(p,a,pgc);
    pflow->fsfrkin(p,a,pgc,ark1);
    pflow->fsfrkin(p,a,pgc,ark2);
    pflow->fsfrkout(p,a,pgc,ark1);
    pflow->fsfrkout(p,a,pgc,ark2);
    ppicard->volcalc(p,a,pgc,a->phi);
	
	pflow->phi_relax(p,pgc,a->phi);

		
//Step 1
//--------------------------------------------------------

    // FSF
    /*starttime=pgc->timer();

    FLUIDLOOP
	a->L(i,j,k)=0.0;

	pfsfdisc->start(p,a,a->phi,4,a->u,a->v,a->w);
	

	FLUIDLOOP
	ark1(i,j,k) = a->phi(i,j,k)
				+ p->dt*a->L(i,j,k);
	
	pflow->phi_relax(p,pgc,ark1,1.0);
	
	pgc->start4(p,ark1,gcval_phi);
    
	if(p->F48==1)
	preini->start(a,p,ark1, pgc, pflow);
    
    pupdate->start(p,a,pgc);*/

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc); 
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a,pgc,a->u,a->u,a->v,a->w,1.0);
	pdisc->start(p,a,a->u,1,a->u,a->v,a->w);
	pdiff->diff_u(p,a,pgc,psolv,a->u,a->v,a->w,1.0);

	ULOOP
	urk1(i,j,k) = a->u(i,j,k)
				+ p->dt*CPOR1*a->F(i,j,k);

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a,pgc,a->v,a->u,a->v,a->w,1.0);
	pdisc->start(p,a,a->v,2,a->u,a->v,a->w);
	pdiff->diff_v(p,a,pgc,psolv,a->u,a->v,a->w,1.0);

	VLOOP
	vrk1(i,j,k) = a->v(i,j,k)
				+ p->dt*CPOR2*a->G(i,j,k);

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a,pgc,a->w,a->u,a->v,a->w,1.0);
	pdisc->start(p,a,a->w,3,a->u,a->v,a->w);
	pdiff->diff_w(p,a,pgc,psolv,a->u,a->v,a->w,1.0);

	WLOOP
	wrk1(i,j,k) = a->w(i,j,k)
				+ p->dt*CPOR3*a->H(i,j,k);
	
    p->wtime=pgc->timer()-starttime;

	pgc->start1(p,urk1,gcval_urk);
	pgc->start2(p,vrk1,gcval_vrk);
	pgc->start3(p,wrk1,gcval_wrk);
	
	urk1.ggcpol(p);
	vrk1.ggcpol(p);
	wrk1.ggcpol(p);
	
    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pmom,pflow, urk1, vrk1, wrk1, 1.0);
	
	pflow->u_relax(p,a,pgc,urk1);
	pflow->v_relax(p,a,pgc,vrk1);
	pflow->w_relax(p,a,pgc,wrk1);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,urk1,gcval_urk);
	pgc->start2(p,vrk1,gcval_vrk);
	pgc->start3(p,wrk1,gcval_wrk);
    
	
    
    // FSF
    starttime=pgc->timer();

    FLUIDLOOP
	a->L(i,j,k)=0.0;

	pfsfdisc->start(p,a,a->phi,4,urk1,vrk1,wrk1);
	

	FLUIDLOOP
	ark1(i,j,k) = a->phi(i,j,k)
				+ p->dt*a->L(i,j,k);
	
	pflow->phi_relax(p,pgc,ark1);
	
	pgc->start4(p,ark1,gcval_phi);
    
	if(p->F48==1)
	preini->start(a,p,ark1, pgc, pflow);
    
    pupdate->start(p,a,pgc);
    

//Step 2
//--------------------------------------------------------
	
    // FSF
    /*FLUIDLOOP
	a->L(i,j,k)=0.0;

	pdisc->start(p,a,ark1,4,urk1,vrk1,wrk1);

	FLUIDLOOP
	ark2(i,j,k) = 0.75*a->phi(i,j,k)
				   + 0.25*ark1(i,j,k)
				   + 0.25*p->dt*a->L(i,j,k);
				
	pflow->phi_relax(p,pgc,ark2,0.25);
	
	pgc->start4(p,ark2,gcval_phi);
	if(p->F48==1)
	preini->start(a,p,ark2, pgc, pflow);
    
    pupdate->start(p,a,pgc);*/
	
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a,pgc,urk1,urk1,vrk1,wrk1,0.25);
	pdisc->start(p,a,urk1,1,urk1,vrk1,wrk1);
	pdiff->diff_u(p,a,pgc,psolv,urk1,vrk1,wrk1,0.25);

	ULOOP
	urk2(i,j,k) = 0.75*a->u(i,j,k) + 0.25*urk1(i,j,k)
				+ 0.25*p->dt*CPOR1*a->F(i,j,k);
                
    p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a,pgc,vrk1,urk1,vrk1,wrk1,0.25);
	pdisc->start(p,a,vrk1,2,urk1,vrk1,wrk1);
	pdiff->diff_v(p,a,pgc,psolv,urk1,vrk1,wrk1,0.25);

	VLOOP
	vrk2(i,j,k) = 0.75*a->v(i,j,k) + 0.25*vrk1(i,j,k)
				+ 0.25*p->dt*CPOR2*a->G(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a,pgc,wrk1,urk1,vrk1,wrk1,0.25);
	pdisc->start(p,a,wrk1,3,urk1,vrk1,wrk1);
	pdiff->diff_w(p,a,pgc,psolv,urk1,vrk1,wrk1,0.25);

	WLOOP
	wrk2(i,j,k) = 0.75*a->w(i,j,k) + 0.25*wrk1(i,j,k)
				+ 0.25*p->dt*CPOR3*a->H(i,j,k);

    p->wtime+=pgc->timer()-starttime;

	pgc->start1(p,urk2,gcval_urk);
	pgc->start2(p,vrk2,gcval_vrk);
	pgc->start3(p,wrk2,gcval_wrk);
	
	urk2.ggcpol(p);
	vrk2.ggcpol(p);
	wrk2.ggcpol(p);

    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pmom,pflow, urk2, vrk2, wrk2, 0.25);
	
	pflow->u_relax(p,a,pgc,urk2);
	pflow->v_relax(p,a,pgc,vrk2);
	pflow->w_relax(p,a,pgc,wrk2);
	pflow->p_relax(p,a,pgc,a->press);
	
	pgc->start1(p,urk2,gcval_urk);
	pgc->start2(p,vrk2,gcval_vrk);
	pgc->start3(p,wrk2,gcval_wrk);
    
    
    // FSF
    FLUIDLOOP
	a->L(i,j,k)=0.0;

	pdisc->start(p,a,ark1,4,urk2,vrk2,wrk2);

	FLUIDLOOP
	ark2(i,j,k) = 0.75*a->phi(i,j,k)
				   + 0.25*ark1(i,j,k)
				   + 0.25*p->dt*a->L(i,j,k);
				
	pflow->phi_relax(p,pgc,ark2);
	
	pgc->start4(p,ark2,gcval_phi);
	if(p->F48==1)
	preini->start(a,p,ark2, pgc, pflow);
    
    pupdate->start(p,a,pgc);

//Step 3
//--------------------------------------------------------
    
    
    // FSF
    /*FLUIDLOOP
	a->L(i,j,k)=0.0;

	pdisc->start(p,a,ark2,4,urk2,vrk2,wrk2);

	FLUIDLOOP
	a->phi(i,j,k) =     (1.0/3.0)*a->phi(i,j,k)
				  + (2.0/3.0)*ark2(i,j,k)
				  + (2.0/3.0)*p->dt*a->L(i,j,k);

    pflow->phi_relax(p,pgc,a->phi,2.0/3.0);
	pgc->start4(p,a->phi,gcval_phi);

    //ppart->start(p,a,pgc,pflow);
    pgc->start4(p,a->phi,gcval_phi);
	
	p->lsmtime=pgc->timer()-starttime;

    if(p->count%p->F41==0)
	preini->start(a,p,a->phi pgc, pflow);
	

    ppicard->correct_ls(p,a,pgc,a->phi);
	//ppart->picardmove(p,a,pgc);

	pflow->periodic(a->phi,p);

	pupdate->start(p,a,pgc);*/


	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a,pgc,urk2,urk2,vrk2,wrk2,2.0/3.0);
	pdisc->start(p,a,urk2,1,urk2,vrk2,wrk2);
	pdiff->diff_u(p,a,pgc,psolv,urk2,vrk2,wrk2,2.0/3.0);

	ULOOP
	a->u(i,j,k) = (1.0/3.0)*a->u(i,j,k) + (2.0/3.0)*urk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR1*a->F(i,j,k);
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a,pgc,vrk2,urk2,vrk2,wrk2,2.0/3.0);
	pdisc->start(p,a,vrk2,2,urk2,vrk2,wrk2);
	pdiff->diff_v(p,a,pgc,psolv,urk2,vrk2,wrk2,2.0/3.0);

	VLOOP
	a->v(i,j,k) = (1.0/3.0)*a->v(i,j,k) + (2.0/3.0)*vrk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR2*a->G(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a,pgc,wrk2,urk2,vrk2,wrk2,2.0/3.0);
	pdisc->start(p,a,wrk2,3,urk2,vrk2,wrk2);
	pdiff->diff_w(p,a,pgc,psolv,urk2,vrk2,wrk2,2.0/3.0);

	WLOOP
	a->w(i,j,k) = (1.0/3.0)*a->w(i,j,k) + (2.0/3.0)*wrk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR3*a->H(i,j,k);
	
    p->wtime+=pgc->timer()-starttime;

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
	
	a->u.ggcpol(p);
	a->v.ggcpol(p);
	a->w.ggcpol(p);

	//--------------------------------------------------------
	// pressure
	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pmom,pflow, a->u, a->v,a->w,2.0/3.0);
	
	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);

	pflow->periodic(a->u,p);
	pflow->periodic(a->v,p);
	pflow->periodic(a->w,p);
    
    
    // FSF
    FLUIDLOOP
	a->L(i,j,k)=0.0;

	pdisc->start(p,a,ark2,4,a->u, a->v,a->w);

	FLUIDLOOP
	a->phi(i,j,k) =     (1.0/3.0)*a->phi(i,j,k)
				  + (2.0/3.0)*ark2(i,j,k)
				  + (2.0/3.0)*p->dt*a->L(i,j,k);

    pflow->phi_relax(p,pgc,a->phi);
	pgc->start4(p,a->phi,gcval_phi);

    //ppart->start(p,a,pgc,pflow);
    pgc->start4(p,a->phi,gcval_phi);
	
	p->lsmtime=pgc->timer()-starttime;

    if(p->count%p->F41==0)
	preini->start(a,p,a->phi, pgc, pflow);
	

    ppicard->correct_ls(p,a,pgc,a->phi);
	//ppart->picardmove(p,a,pgc);

	pflow->periodic(a->phi,p);

	pupdate->start(p,a,pgc);

}

void momentum_FSFC_RK3::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{

    pgc->forcing1(p,a,f,uvel,vvel,wvel,alpha);
    
	n=0;
	if(p->D20<3)
	ULOOP
	{
    a->maxF=MAX(fabs(a->rhsvec.V[n]),a->maxF);
	a->F(i,j,k) += (a->rhsvec.V[n] + a->gi)*PORVAL1;
	a->rhsvec.V[n]=0.0;
	++n;
	}
	
	n=0;
	if(p->D20==3)
	ULOOP
	{
	a->rhsvec.V[n]+=a->gi;
	++n;
	}
}

void momentum_FSFC_RK3::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
    pgc->forcing2(p,a,f,uvel,vvel,wvel,alpha);
    
	n=0;
	if(p->D20<3)
	VLOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n]),a->maxG);
	a->G(i,j,k) += (a->rhsvec.V[n] + a->gj)*PORVAL2;
	a->rhsvec.V[n]=0.0;
	++n;
	}
	
	n=0;
	if(p->D20==3)
	VLOOP
	{
	a->rhsvec.V[n]+=a->gj;
	++n;
	}
}

void momentum_FSFC_RK3::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
    pgc->forcing3(p,a,f,uvel,vvel,wvel,alpha);
    
	n=0;
	if(p->D20<3)
	WLOOP
	{
    a->maxH=MAX(fabs(a->rhsvec.V[n]),a->maxH);
	a->H(i,j,k) += (a->rhsvec.V[n] + a->gk)*PORVAL3;
	a->rhsvec.V[n]=0.0;
	++n;
	}
	
	n=0;
	if(p->D20==3)
	WLOOP
	{
	a->rhsvec.V[n]+=a->gk;
	++n;
	}
}

void momentum_FSFC_RK3::utimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FSFC_RK3::vtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FSFC_RK3::wtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FSFC_RK3::fillaij1(lexer *p, fdm *a, ghostcell* pgc, solver *psolv)
{
}

void momentum_FSFC_RK3::fillaij2(lexer *p, fdm *a, ghostcell* pgc, solver *psolv)
{
}

void momentum_FSFC_RK3::fillaij3(lexer *p, fdm *a, ghostcell* pgc, solver *psolv)
{
}
