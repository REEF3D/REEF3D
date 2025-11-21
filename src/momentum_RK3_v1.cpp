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

#include"momentum_RK3_v1.h"
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
#include"fluid_update_rheology.h"
#include"fluid_update_void.h"
#include"nhflow.h"

momentum_RK3_v1::momentum_RK3_v1(lexer *p, fdm *a, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, 
                                                    ioflow *pioflow, fsi *ppfsi)
                                                    :momentum_forcing(p),bcmom(p),udiff(p),vdiff(p),wdiff(p),urk(p),vrk(p),wrk(p),fx(p),fy(p),fz(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;

	pconvec=pconvection;
	pdiff=pdiffusion;
	ppress=ppressure;
	ppois=ppoisson;
	pturb=pturbulence;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;    
    pfsi=ppfsi;

    if(p->W90==0  || p->F300>0)
	pupdate = new fluid_update_void();
    
    if(p->W90>0 && p->F300==0)
    pupdate = new fluid_update_rheology(p);
    
    alpha[0] = 0.0;
    alpha[1] = 0.75;
    alpha[2] = 1.0/3.0;
    
    beta[0] = 1.0;
    beta[1] = 0.25;
    beta[2] = 2.0/3.0;
}

momentum_RK3_v1::~momentum_RK3_v1()
{
}

void momentum_RK3_v1::start(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof)
{	
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	pflow->rkinflow(p,a,pgc,urk,vrk,wrk);
	pflow->rkinflow(p,a,pgc,urk,vrk,wrk);
		
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
	pdiff->diff_u(p,a,pgc,psolv,udiff,a->u,a->u,a->v,a->w,1.0);

	ULOOP
	urk(i,j,k) = udiff(i,j,k) + p->dt*CPOR1*a->F(i,j,k);

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,a->v,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,a->v,a->u,a->v,a->w,1.0);

	VLOOP
	vrk(i,j,k) = vdiff(i,j,k) + p->dt*CPOR2*a->G(i,j,k);

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,a->w,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,a->w,a->u,a->v,a->w,1.0);

	WLOOP
	wrk(i,j,k) = wdiff(i,j,k) + p->dt*CPOR3*a->H(i,j,k);
	
    p->wtime=pgc->timer()-starttime;
    
    momentum_forcing_start(a, p, pgc, p6dof, pfsi,
                           urk, vrk, wrk, fx, fy, fz, 0, 1.0, false);
    
    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk, vrk, wrk, 1.0);
	
	pflow->u_relax(p,a,pgc,urk);
	pflow->v_relax(p,a,pgc,vrk);
	pflow->w_relax(p,a,pgc,wrk);
	pflow->p_relax(p,a,pgc,a->press);
    
    ULOOP
    urk(i,j,k) = alpha[0]*a->u(i,j,k) + beta[0]*urk(i,j,k);
    
    VLOOP
    vrk(i,j,k) = alpha[0]*a->v(i,j,k) + beta[0]*vrk(i,j,k);
    
    WLOOP
    wrk(i,j,k) = alpha[0]*a->w(i,j,k) + beta[0]*wrk(i,j,k);

	pgc->start1(p,urk,gcval_u);
	pgc->start2(p,vrk,gcval_v);
	pgc->start3(p,wrk,gcval_w);
    
//Step 2
//--------------------------------------------------------
	
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,urk,urk,vrk,wrk,1.0);
	pconvec->start(p,a,urk,1,urk,vrk,wrk);
	pdiff->diff_u(p,a,pgc,psolv,udiff,urk,urk,vrk,wrk,1.0);

	ULOOP
	urk(i,j,k) = udiff(i,j,k) + p->dt*CPOR1*a->F(i,j,k);
                
    p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,vrk,urk,vrk,wrk,1.0);
	pconvec->start(p,a,vrk,2,urk,vrk,wrk);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,vrk,urk,vrk,wrk,1.0);

	VLOOP
	vrk(i,j,k) = vdiff(i,j,k) + p->dt*CPOR2*a->G(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,wrk,urk,vrk,wrk,1.0);
	pconvec->start(p,a,wrk,3,urk,vrk,wrk);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,wrk,urk,vrk,wrk,1.0);

	WLOOP
	wrk(i,j,k) = wdiff(i,j,k) + p->dt*CPOR3*a->H(i,j,k);

    p->wtime+=pgc->timer()-starttime;

    momentum_forcing_start(a, p, pgc, p6dof, pfsi,
                           urk, vrk, wrk, fx, fy, fz, 1, 1.0, false);

    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk, vrk, wrk, 1.0);
	
	pflow->u_relax(p,a,pgc,urk);
	pflow->v_relax(p,a,pgc,vrk);
	pflow->w_relax(p,a,pgc,wrk);
	pflow->p_relax(p,a,pgc,a->press);
	
    ULOOP
    urk(i,j,k) = alpha[1]*a->u(i,j,k) + beta[1]*urk(i,j,k);
    
    VLOOP
    vrk(i,j,k) = alpha[1]*a->v(i,j,k) + beta[1]*vrk(i,j,k);
    
    WLOOP
    wrk(i,j,k) = alpha[1]*a->w(i,j,k) + beta[1]*wrk(i,j,k);
    
	pgc->start1(p,urk,gcval_u);
	pgc->start2(p,vrk,gcval_v);
	pgc->start3(p,wrk,gcval_w);

//Step 3
//--------------------------------------------------------

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,urk,urk,vrk,wrk,1.0);
	pconvec->start(p,a,urk,1,urk,vrk,wrk);
	pdiff->diff_u(p,a,pgc,psolv,udiff,urk,urk,vrk,wrk,1.0);

	ULOOP
	a->u(i,j,k) = udiff(i,j,k) + p->dt*CPOR1*a->F(i,j,k);
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,vrk,urk,vrk,wrk,1.0);
	pconvec->start(p,a,vrk,2,urk,vrk,wrk);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,vrk,urk,vrk,wrk,1.0);

	VLOOP
	a->v(i,j,k) = vdiff(i,j,k) + p->dt*CPOR2*a->G(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,wrk,urk,vrk,wrk,1.0);
	pconvec->start(p,a,wrk,3,urk,vrk,wrk);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,wrk,urk,vrk,wrk,1.0);

	WLOOP
	wrk(i,j,k) = wdiff(i,j,k) + p->dt*CPOR3*a->H(i,j,k);
	
    p->wtime+=pgc->timer()-starttime;

    momentum_forcing_start(a, p, pgc, p6dof, pfsi,
                           a->u, a->v, a->w, fx, fy, fz, 2, 1.0, true);

	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow,a->u,a->v,a->w,1.0);
	
	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);
    
    ULOOP
    a->u(i,j,k) = alpha[2]*a->u(i,j,k) + beta[2]*urk(i,j,k);
    
    VLOOP
    a->v(i,j,k) = alpha[2]*a->v(i,j,k) + beta[2]*vrk(i,j,k);
    
    WLOOP
    a->w(i,j,k) = alpha[2]*a->w(i,j,k) + beta[2]*wrk(i,j,k);

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
}

void momentum_RK3_v1::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	ULOOP
	{
    a->maxF=MAX(fabs(a->rhsvec.V[n] + a->gi),a->maxF);
	a->F(i,j,k) += (a->rhsvec.V[n] + a->gi + p->W29_x + a->Fext(i,j,k))*PORVAL1;
    
	a->rhsvec.V[n]=0.0;
    a->Fext(i,j,k)=0.0;
	++n;
	}
}

void momentum_RK3_v1::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	VLOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n] + a->gj),a->maxG);
	a->G(i,j,k) += (a->rhsvec.V[n] + a->gj + p->W29_y + a->Gext(i,j,k))*PORVAL2;
    
	a->rhsvec.V[n] = 0.0;
    a->Gext(i,j,k) = 0.0;
	++n;
	}
}

void momentum_RK3_v1::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	WLOOP
	{
    a->maxH=MAX(fabs(a->rhsvec.V[n] + a->gk),a->maxH);
	a->H(i,j,k) += (a->rhsvec.V[n] + a->gk + p->W29_z + a->Hext(i,j,k))*PORVAL3;
    
	a->rhsvec.V[n] = 0.0;
    a->Hext(i,j,k) = 0.0;
	++n;
	}
}

