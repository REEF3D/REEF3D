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

#include"momentum_RK3_gen.h"
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

momentum_RK3_gen::momentum_RK3_gen(lexer *p, fdm *a, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, 
                                                    ioflow *pioflow, fsi *ppfsi)
                                                    :momentum_forcing(p),bcmom(p),udiff(p),vdiff(p),wdiff(p),urk(p),urk1(p),urk2(p),vrk(p),vrk1(p),
                                                    vrk2(p),wrk(p),wrk1(p),wrk2(p),fx(p),fy(p),fz(p)
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
    
    
    alpha[0] = 8.0/15.0;
    alpha[1] = 1.0/15.0;
    alpha[2] = 1.0/6.0;
    
    beta[0] = 8.0/15.0;
    beta[1] = 1.0/15.0;
    beta[2] = 1.0/6.0;
    
    /*
    alpha[0] = 1.0;
    alpha[1] = 0.25;
    alpha[2] = 2.0/3.0;
    
    beta[0] = 1.0;
    beta[1] = 0.75;
    beta[2] = 1.0/3.0;*/
}

momentum_RK3_gen::~momentum_RK3_gen()
{
}

void momentum_RK3_gen::start(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof)
{	
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
	pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);
		
//Step 1
//--------------------------------------------------------
    
    loop=0;
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans); 
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,a->u,a->u,a->v,a->w,alpha[loop]);
	pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
	pdiff->diff_u(p,a,pgc,psolv,urk,a->u,a->u,a->v,a->w,alpha[loop]);

	ULOOP
    {
    urk1(i,j,k) = a->F(i,j,k);
                
	urk(i,j,k) = beta[loop]*urk(i,j,k)
				+ p->dt*alpha[loop]*CPOR1*a->F(i,j,k);
    }

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,a->v,a->u,a->v,a->w,alpha[loop]);
	pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
	pdiff->diff_v(p,a,pgc,psolv,vrk,a->v,a->u,a->v,a->w,alpha[loop]);

	VLOOP
    {
	vrk1(i,j,k) = a->G(i,j,k);
                
    vrk(i,j,k) = beta[loop]*vrk(i,j,k)
				+ p->dt*alpha[loop]*CPOR2*a->G(i,j,k);
    }

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,a->w,a->u,a->v,a->w,alpha[loop]);
	pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
	pdiff->diff_w(p,a,pgc,psolv,wrk,a->w,a->u,a->v,a->w,alpha[loop]);

	WLOOP
    {
	wrk1(i,j,k) = a->H(i,j,k);
                
    wrk(i,j,k) = beta[loop]*wrk(i,j,k)
				+ p->dt*alpha[loop]*CPOR3*a->H(i,j,k);
    }
	
    p->wtime=pgc->timer()-starttime;
    
    momentum_forcing_start(a, p, pgc, p6dof, pfsi,
                           urk1, vrk1, wrk1, fx, fy, fz, 0, alpha[loop], false);
    
    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk1, vrk1, wrk1, alpha[loop]);
	
	pflow->u_relax(p,a,pgc,urk1);
	pflow->v_relax(p,a,pgc,vrk1);
	pflow->w_relax(p,a,pgc,wrk1);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,urk1,gcval_u);
	pgc->start2(p,vrk1,gcval_v);
	pgc->start3(p,wrk1,gcval_w);
    
//Step 2
//--------------------------------------------------------
	
    loop=1;
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,urk1,urk1,vrk1,wrk1,alpha[loop]);
	pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
	pdiff->diff_u(p,a,pgc,psolv,urk,a->u,urk1,vrk1,wrk1,alpha[loop]);

	ULOOP
    {
	urk2(i,j,k) = a->F(i,j,k);
                
    urk(i,j,k) = beta[loop]*a->u(i,j,k) + alpha[loop]*udiff(i,j,k)
				+ alpha[loop]*p->dt*CPOR1*a->F(i,j,k);
    }
                
    p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,vrk1,urk1,vrk1,wrk1,alpha[loop]);
	pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
	pdiff->diff_v(p,a,pgc,psolv,vrk,a->v,urk1,vrk1,wrk1,alpha[loop]);

	VLOOP
	vrk2(i,j,k) = beta[loop]*a->v(i,j,k) + alpha[loop]*vdiff(i,j,k)
				+ alpha[loop]*p->dt*CPOR2*a->G(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,wrk1,urk1,vrk1,wrk1,alpha[loop]);
	pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
	pdiff->diff_w(p,a,pgc,psolv,wrk,a->w,urk1,vrk1,wrk1,alpha[loop]);

	WLOOP
	wrk2(i,j,k) = beta[loop]*a->w(i,j,k) + alpha[loop]*wdiff(i,j,k)
				+ alpha[loop]*p->dt*CPOR3*a->H(i,j,k);

    p->wtime+=pgc->timer()-starttime;

    momentum_forcing_start(a, p, pgc, p6dof, pfsi,
                           urk2, vrk2, wrk2, fx, fy, fz, 1, alpha[loop], false);

    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk2, vrk2, wrk2, alpha[loop]);
	
	pflow->u_relax(p,a,pgc,urk2);
	pflow->v_relax(p,a,pgc,vrk2);
	pflow->w_relax(p,a,pgc,wrk2);
	pflow->p_relax(p,a,pgc,a->press);
	
	pgc->start1(p,urk2,gcval_u);
	pgc->start2(p,vrk2,gcval_v);
	pgc->start3(p,wrk2,gcval_w);

//Step 3
//--------------------------------------------------------
    
    loop=2;
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,urk2,urk2,vrk2,wrk2,alpha[loop]);
	pconvec->start(p,a,urk2,1,urk2,vrk2,wrk2);
	pdiff->diff_u(p,a,pgc,psolv,urk,a->u,urk2,vrk2,wrk2,alpha[loop]);

	ULOOP
	a->u(i,j,k) = beta[loop]*a->u(i,j,k) + alpha[loop]*udiff(i,j,k)
				+ alpha[loop]*p->dt*CPOR1*a->F(i,j,k);
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,vrk2,urk2,vrk2,wrk2,alpha[loop]);
	pconvec->start(p,a,vrk2,2,urk2,vrk2,wrk2);
	pdiff->diff_v(p,a,pgc,psolv,vrk,a->v,urk2,vrk2,wrk2,alpha[loop]);

	VLOOP
	a->v(i,j,k) = beta[loop]*a->v(i,j,k) + alpha[loop]*vdiff(i,j,k)
				+ alpha[loop]*p->dt*CPOR2*a->G(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,wrk2,urk2,vrk2,wrk2,alpha[loop]);
	pconvec->start(p,a,wrk2,3,urk2,vrk2,wrk2);
	pdiff->diff_w(p,a,pgc,psolv,wrk,a->w,urk2,vrk2,wrk2,alpha[loop]);

	WLOOP
	a->w(i,j,k) = beta[loop]*a->w(i,j,k) + alpha[loop]*wdiff(i,j,k)
				+ alpha[loop]*p->dt*CPOR3*a->H(i,j,k);
	
    p->wtime+=pgc->timer()-starttime;

    momentum_forcing_start(a, p, pgc, p6dof, pfsi,
                           a->u, a->v, a->w, fx, fy, fz, 2, alpha[loop], true);

	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow,a->u,a->v,a->w,alpha[loop]);
	
	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
}

void momentum_RK3_gen::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
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

void momentum_RK3_gen::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
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

void momentum_RK3_gen::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
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

