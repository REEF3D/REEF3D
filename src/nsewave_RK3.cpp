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
#include"nsewave_RK3.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_fsf.h"
#include"fluid_update_fsf_heat.h"
#include"fluid_update_fsf_comp.h"
#include"fluid_update_fsf_concentration.h"
#include"fluid_update_rheology.h"
#include"fluid_update_void.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"heat.h"
#include"concentration.h"
#include"turbulence.h"
#include"pressure.h"
#include"convection.h"
#include"diffusion.h"

nsewave_RK3::nsewave_RK3(lexer *p, fdm *a, ghostcell *pgc, heat *&pheat, concentration *&pconc) : eta(p),
                etark1(p),etark2(p),L(p),P(p),Q(p),epsi(1.6*p->DXM),bcmom(p),
                urk1(p),urk2(p),vrk1(p),vrk2(p),wrk1(p),wrk2(p)
{
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;

    
    if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;


    pupdate = new fluid_update_fsf(p,a,pgc);

	ppicard = new picard_void(p);
    
    p->phimean=p->F60;
    
    SLICELOOP4
    eta(i,j) = 0.0;

    
    pgc->gcsl_start4(p,eta,gcval_phi);
    
    LOOP
    a->phi(i,j,k) = eta(i,j) + p->phimean - p->pos_z();
    
    pgc->start4(p,a->phi,gcval_phi);
}

nsewave_RK3::~nsewave_RK3()
{
}


void nsewave_RK3::start(lexer* p, fdm* a, ghostcell* pgc, momentum *pmom, diffusion *pdiff, turbulence *pturb,
                      convection* pconvec, pressure *ppress, poisson *ppois, solver *ppoissonsolv, solver *psolv, 
                      ioflow* pflow, vrans* pvrans, sixdof *p6dof, vector<net*>& pnet)
{
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
	pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);
		
//Step 1
//--------------------------------------------------------

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans); 
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,a->u,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
	pdiff->diff_u(p,a,pgc,psolv,a->u,a->v,a->w,1.0);

	ULOOP
	urk1(i,j,k) = a->u(i,j,k)
				+ p->dt*CPOR1*a->F(i,j,k);

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,a->v,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
	pdiff->diff_v(p,a,pgc,psolv,a->u,a->v,a->w,1.0);

	VLOOP
	vrk1(i,j,k) = a->v(i,j,k)
				+ p->dt*CPOR2*a->G(i,j,k);

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,a->w,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
	pdiff->diff_w(p,a,pgc,psolv,a->u,a->v,a->w,1.0);

	WLOOP
	wrk1(i,j,k) = a->w(i,j,k)
				+ p->dt*CPOR3*a->H(i,j,k);
	
    p->wtime=pgc->timer()-starttime;

	pgc->start1(p,urk1,gcval_u);
	pgc->start2(p,vrk1,gcval_v);
	pgc->start3(p,wrk1,gcval_w);

	
    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk1, vrk1, wrk1, 1.0);
	
	pflow->u_relax(p,a,pgc,urk1);
	pflow->v_relax(p,a,pgc,vrk1);
	pflow->w_relax(p,a,pgc,wrk1);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,urk1,gcval_u);
	pgc->start2(p,vrk1,gcval_v);
	pgc->start3(p,wrk1,gcval_w);
    
    // Calculate Eta
    eta_disc(p,a,pgc,urk1,vrk1);
    SLICELOOP4
    etark1(i,j) = eta(i,j) + p->dt*L(i,j);
    
    pflow->eta_relax(p,pgc,etark1);
    pgc->gcsl_start4(p,etark1,gcval_phi);
    
    FLUIDLOOP
    a->phi(i,j,k) = etark1(i,j) + p->phimean - p->pos_z();
    
    pgc->start4(p,a->phi,gcval_phi);
    pupdate->start(p,a,pgc);
    
//Step 2
//--------------------------------------------------------
	
	
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,urk1,urk1,vrk1,wrk1,0.25);
	pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
	pdiff->diff_u(p,a,pgc,psolv,urk1,vrk1,wrk1,0.25);

	ULOOP
	urk2(i,j,k) = 0.75*a->u(i,j,k) + 0.25*urk1(i,j,k)
				+ 0.25*p->dt*CPOR1*a->F(i,j,k);
                
    p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,vrk1,urk1,vrk1,wrk1,0.25);
	pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
	pdiff->diff_v(p,a,pgc,psolv,urk1,vrk1,wrk1,80.25);

	VLOOP
	vrk2(i,j,k) = 0.75*a->v(i,j,k) + 0.25*vrk1(i,j,k)
				+ 0.25*p->dt*CPOR2*a->G(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,wrk1,urk1,vrk1,wrk1,0.25);
	pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
	pdiff->diff_w(p,a,pgc,psolv,urk1,vrk1,wrk1,0.25);

	WLOOP
	wrk2(i,j,k) = 0.75*a->w(i,j,k) + 0.25*wrk1(i,j,k)
				+ 0.25*p->dt*CPOR3*a->H(i,j,k);

    p->wtime+=pgc->timer()-starttime;

	pgc->start1(p,urk2,gcval_u);
	pgc->start2(p,vrk2,gcval_v);
	pgc->start3(p,wrk2,gcval_w);
	
    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk2, vrk2, wrk2,0.25);
	
	pflow->u_relax(p,a,pgc,urk2);
	pflow->v_relax(p,a,pgc,vrk2);
	pflow->w_relax(p,a,pgc,wrk2);
	pflow->p_relax(p,a,pgc,a->press);
	
	pgc->start1(p,urk2,gcval_u);
	pgc->start2(p,vrk2,gcval_v);
	pgc->start3(p,wrk2,gcval_w);
    
    // Calculate Eta
    eta_disc(p,a,pgc,urk2,vrk2);
    SLICELOOP4
    etark2(i,j) = 0.75*eta(i,j) + 0.25*etark1(i,j) + 0.25*p->dt*L(i,j);
    
    pflow->eta_relax(p,pgc,etark2);
    pgc->gcsl_start4(p,etark2,gcval_phi);
    
    FLUIDLOOP
    a->phi(i,j,k) = etark2(i,j) + p->phimean - p->pos_z();
    
    pgc->start4(p,a->phi,gcval_phi);
    pupdate->start(p,a,pgc);
    
//Step 3
//--------------------------------------------------------

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,urk2,urk2,vrk2,wrk2,2.0/3.0);
	pconvec->start(p,a,urk2,1,urk2,vrk2,wrk2);
	pdiff->diff_u(p,a,pgc,psolv,urk2,vrk2,wrk2,2.0/3.0);

	ULOOP
	a->u(i,j,k) = (1.0/3.0)*a->u(i,j,k) + (2.0/3.0)*urk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR1*a->F(i,j,k);
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,vrk2,urk2,vrk2,wrk2,2.0/3.0);
	pconvec->start(p,a,vrk2,2,urk2,vrk2,wrk2);
	pdiff->diff_v(p,a,pgc,psolv,urk2,vrk2,wrk2,2.0/3.0);

	VLOOP
	a->v(i,j,k) = (1.0/3.0)*a->v(i,j,k) + (2.0/3.0)*vrk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR2*a->G(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,wrk2,urk2,vrk2,wrk2,2.0/3.0);
	pconvec->start(p,a,wrk2,3,urk2,vrk2,wrk2);
	pdiff->diff_w(p,a,pgc,psolv,urk2,vrk2,wrk2,2.0/3.0);

	WLOOP
	a->w(i,j,k) = (1.0/3.0)*a->w(i,j,k) + (2.0/3.0)*wrk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR3*a->H(i,j,k);
	
    p->wtime+=pgc->timer()-starttime;

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);

	//--------------------------------------------------------
	// pressure
	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, a->u, a->v,a->w,2.0/3.0);
	
	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
    
    // Calculate Eta
    eta_disc(p,a,pgc,a->u,a->v);
    SLICELOOP4
    eta(i,j) = (1.0/3.0)*eta(i,j) + (2.0/3.0)*etark2(i,j) + (2.0/3.0)*p->dt*L(i,j);

    
    pflow->eta_relax(p,pgc,eta);
    pgc->gcsl_start4(p,eta,gcval_phi);
    
    FLUIDLOOP
    a->phi(i,j,k) = eta(i,j) + p->phimean - p->pos_z();
    
    pgc->start4(p,a->phi,gcval_phi);
    pupdate->start(p,a,pgc);
    
}

void nsewave_RK3::eta_disc(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v)
{
    SLICELOOP1
    P(i,j)=0.0;
    
    SLICELOOP2
    Q(i,j)=0.0;

    ULOOP
    {
    phival = 0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));
    
        if(phival>epsi)
		H=1.0;

		if(phival<-epsi)
		H=0.0;

		if(fabs(phival)<=epsi)
		H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
        
        P(i,j) += u(i,j,k)*p->DXM*H;
    }
    
    VLOOP
    {
    phival = 0.5*(a->phi(i,j,k)+a->phi(i,j+1,k));
    
        if(phival>epsi)
		H=1.0;

		if(phival<-epsi)
		H=0.0;

		if(fabs(phival)<=epsi)
		H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
        
        Q(i,j) += v(i,j,k)*p->DXM*H;
    }
    
    pgc->gcsl_start1(p,P,10);
    pgc->gcsl_start2(p,Q,11);
    
    
    SLICELOOP4
    L(i,j) = -(P(i,j)-P(i-1,j) + Q(i,j)-Q(i,j-1))/p->DXM;
    
}

void nsewave_RK3::update(lexer *p, fdm *a, ghostcell *pgc, slice &f)
{
    pupdate->start(p,a,pgc);
}

void nsewave_RK3::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{
    
    p->phimean=p->F60;
    
    SLICELOOP4
    eta(i,j) = 0.0;
    
    pflow->eta_relax(p,pgc,eta);
    
    pgc->gcsl_start4(p,eta,gcval_phi);
    
    LOOP
    a->phi(i,j,k) = eta(i,j) + p->phimean - p->pos_z();
    
    pgc->start4(p,a->phi,gcval_phi);
    
    pupdate->start(p,a,pgc);
}


void nsewave_RK3::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{ 
	n=0;
	ULOOP
	{
    a->maxF=MAX(fabs(a->rhsvec.V[n] + a->gi),a->maxF);
	a->F(i,j,k) += (a->rhsvec.V[n] + a->gi)*PORVAL1;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void nsewave_RK3::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	VLOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n] + a->gj),a->maxG);
	a->G(i,j,k) += (a->rhsvec.V[n] + a->gj)*PORVAL2;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void nsewave_RK3::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	if(p->D20<4)
	WLOOP
	{
    a->maxH=MAX(fabs(a->rhsvec.V[n] + a->gk),a->maxH);
	a->H(i,j,k) += (a->rhsvec.V[n] + a->gk)*PORVAL3;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}
