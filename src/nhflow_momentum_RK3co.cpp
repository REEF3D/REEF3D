/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"nhflow_momentum_RK3co.h"
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
#include"nhflow_fsf.h"

nhflow_momentum_RK3co::nhflow_momentum_RK3co(lexer *p, fdm *a, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, ioflow *pioflow,
                                                    nhflow *ppnh, nhflow_fsf *ppnhfsf)
                                                    :bcmom(p),udiff(p),vdiff(p),wdiff(p),urk1(p),urk2(p),vrk1(p),vrk2(p),wrk1(p),wrk2(p),
                                                    etark1(p),etark2(p)
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
    pnh=ppnh;
    pnhfsf=ppnhfsf;
    
    if(p->W90>0)
	pupdate = new fluid_update_rheology(p,a);
    
    if(p->W90==0)
	pupdate = new fluid_update_void();
}

nhflow_momentum_RK3co::~nhflow_momentum_RK3co()
{
}

void nhflow_momentum_RK3co::start(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{	
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->U,a->V,a->W);
	pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
	pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);
		
//Step 1
//--------------------------------------------------------

    pnhfsf->step1(p, a, pgc, pflow, a->U, a->V, a->W, etark1, etark2, 1.0);

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans); 
	bcmom_start(a,p,pgc,pturb,a->U,gcval_u);
	ppress->upgrad(p,a,etark1,a->eta);
	irhs(p,a,pgc,a->U,a->U,a->V,a->W,1.0);
	pconvec->start(p,a,a->U,1,a->U,a->V,a->W);
	pdiff->diff_u(p,a,pgc,psolv,udiff,a->U,a->V,a->W,1.0);

	LOOP
	urk1(i,j,k) = a->U(i,j,k)
				+ p->dt*CPOR1*a->L(i,j,k);

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->V,gcval_v);
	ppress->vpgrad(p,a,etark1,a->eta);
	jrhs(p,a,pgc,a->V,a->U,a->V,a->W,1.0);
	pconvec->start(p,a,a->V,2,a->U,a->V,a->W);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,a->U,a->V,a->W,1.0);

	LOOP
	vrk1(i,j,k) = a->V(i,j,k)
				+ p->dt*CPOR2*a->L(i,j,k);

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->W,gcval_w);
	ppress->wpgrad(p,a,etark1,a->eta);
	krhs(p,a,pgc,a->W,a->U,a->V,a->W,1.0);
	pconvec->start(p,a,a->W,3,a->U,a->V,a->W);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,a->U,a->V,a->W,1.0);

	LOOP
	wrk1(i,j,k) = a->W(i,j,k)
				+ p->dt*CPOR3*a->L(i,j,k);
	
    p->wtime=pgc->timer()-starttime;
    
    pgc->start1(p,urk1,gcval_u);
	pgc->start2(p,vrk1,gcval_v);
    pgc->start3(p,wrk1,gcval_w);
    
    pnh->kinematic_fsf(p,a,a->U,a->V,wrk1,etark1,a->eta,1.0);
    p->omega_update(p,a,pgc,a->U,a->V,wrk1,etark1,etark1,1.0);

    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk1, vrk1, wrk1, 1.0);
	
	pflow->u_relax(p,a,pgc,urk1);
	pflow->v_relax(p,a,pgc,vrk1);
	pflow->w_relax(p,a,pgc,wrk1);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,urk1,gcval_u);
	pgc->start2(p,vrk1,gcval_v);
	pgc->start3(p,wrk1,gcval_w);
    
    pnh->kinematic_fsf(p,a,urk1,vrk1,wrk1,etark1,a->eta,1.0);
    p->omega_update(p,a,pgc,urk1,vrk1,wrk1,etark1,etark1,1.0);
    
    pupdate->start(p,a,pgc);
    
//Step 2
//--------------------------------------------------------
	
    pnhfsf->step2(p, a, pgc, pflow, urk1, vrk1, wrk1, etark1, etark2, 0.25);
    
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->U,gcval_u);
	ppress->upgrad(p,a,etark2,etark1);
	irhs(p,a,pgc,urk1,urk1,vrk1,wrk1,0.25);
	pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
	pdiff->diff_u(p,a,pgc,psolv,udiff,urk1,vrk1,wrk1,1.0);

	LOOP
	urk2(i,j,k) = 0.75*a->U(i,j,k) + 0.25*urk1(i,j,k)
				+ 0.25*p->dt*CPOR1*a->L(i,j,k);
                
    p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->V,gcval_v);
	ppress->vpgrad(p,a,etark2,etark1);
	jrhs(p,a,pgc,vrk1,urk1,vrk1,wrk1,0.25);
	pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,urk1,vrk1,wrk1,1.0);

	LOOP
	vrk2(i,j,k) = 0.75*a->V(i,j,k) + 0.25*vrk1(i,j,k)
				+ 0.25*p->dt*CPOR2*a->L(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->W,gcval_w);
	ppress->wpgrad(p,a,etark2,etark1);
	krhs(p,a,pgc,wrk1,urk1,vrk1,wrk1,0.25);
	pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,urk1,vrk1,wrk1,1.0);

	LOOP
	wrk2(i,j,k) = 0.75*a->W(i,j,k) + 0.25*wrk1(i,j,k)
				+ 0.25*p->dt*CPOR3*a->L(i,j,k);

    p->wtime+=pgc->timer()-starttime;
    
    pgc->start1(p,urk2,gcval_u);
	pgc->start2(p,vrk2,gcval_v);
    pgc->start3(p,wrk2,gcval_w);
    
    pnh->kinematic_fsf(p,a,urk1,vrk1,wrk2,etark2,etark1,0.25);
    p->omega_update(p,a,pgc,urk1,vrk1,wrk2,etark2,etark1,0.25);

    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk2, vrk2, wrk2, 0.25);
	
	pflow->u_relax(p,a,pgc,urk2);
	pflow->v_relax(p,a,pgc,vrk2);
	pflow->w_relax(p,a,pgc,wrk2);
	pflow->p_relax(p,a,pgc,a->press);
	
	pgc->start1(p,urk2,gcval_u);
	pgc->start2(p,vrk2,gcval_v);
	pgc->start3(p,wrk2,gcval_w);
    
    pnh->kinematic_fsf(p,a,urk2,vrk2,wrk2,etark2,etark1,0.25);
    p->omega_update(p,a,pgc,urk2,vrk2,wrk2,etark2,etark1,0.25);
    
    pupdate->start(p,a,pgc);

//Step 3
//--------------------------------------------------------
    
    pnhfsf->step3(p, a, pgc, pflow, urk2, vrk2, wrk2, etark1, etark2, 2.0/3.0);
    
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->U,gcval_u);
	ppress->upgrad(p,a,a->eta,etark2);
	irhs(p,a,pgc,urk2,urk2,vrk2,wrk2,2.0/3.0);
	pconvec->start(p,a,urk2,1,urk2,vrk2,wrk2);
	pdiff->diff_u(p,a,pgc,psolv,udiff,urk2,vrk2,wrk2,1.0);

	LOOP
	a->U(i,j,k) = (1.0/3.0)*a->U(i,j,k) + (2.0/3.0)*urk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR1*a->L(i,j,k);
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->V,gcval_v);
	ppress->vpgrad(p,a,a->eta,etark2);
	jrhs(p,a,pgc,vrk2,urk2,vrk2,wrk2,2.0/3.0);
	pconvec->start(p,a,vrk2,2,urk2,vrk2,wrk2);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,urk2,vrk2,wrk2,1.0);

	LOOP
	a->V(i,j,k) = (1.0/3.0)*a->V(i,j,k) + (2.0/3.0)*vrk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR2*a->L(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->W,gcval_w);
	ppress->wpgrad(p,a,a->eta,etark2);
	krhs(p,a,pgc,wrk2,urk2,vrk2,wrk2,2.0/3.0);
	pconvec->start(p,a,wrk2,3,urk2,vrk2,wrk2);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,urk2,vrk2,wrk2,1.0);

	LOOP
	a->W(i,j,k) = (1.0/3.0)*a->W(i,j,k) + (2.0/3.0)*wrk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR3*a->L(i,j,k);
	
    p->wtime+=pgc->timer()-starttime;
    
    pgc->start1(p,a->U,gcval_u);
	pgc->start2(p,a->V,gcval_v);
    pgc->start3(p,a->W,gcval_w);
    
    pnh->kinematic_fsf(p,a,urk2,vrk2,a->W,a->eta,etark2,2.0/3.0);
    p->omega_update(p,a,pgc,urk2,vrk2,a->W,a->eta,etark2,2.0/3.0);

	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, a->U, a->V, a->W, 2.0/3.0);
	
	pflow->u_relax(p,a,pgc,a->U);
	pflow->v_relax(p,a,pgc,a->V);
	pflow->w_relax(p,a,pgc,a->W);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->U,gcval_u);
	pgc->start2(p,a->V,gcval_v);
	pgc->start3(p,a->W,gcval_w);
    
    pnh->kinematic_fsf(p,a,a->U,a->V,a->W,a->eta,etark2,2.0/3.0);
    p->omega_update(p,a,pgc,a->U,a->V,a->W,a->eta,etark2,2.0/3.0);
    
    pupdate->start(p,a,pgc);
}

void nhflow_momentum_RK3co::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	LOOP
	{
    a->maxF=MAX(fabs(a->rhsvec.V[n]),a->maxF);
	a->L(i,j,k) += (a->rhsvec.V[n])*PORVAL1;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void nhflow_momentum_RK3co::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	LOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n]),a->maxG);
	a->L(i,j,k) += (a->rhsvec.V[n])*PORVAL2;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void nhflow_momentum_RK3co::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	LOOP
	{
    a->maxH=MAX(fabs(a->rhsvec.V[n] + a->gk),a->maxH);
    
    if(p->D38==0)
    a->L(i,j,k) += (a->rhsvec.V[n] + a->gk)*PORVAL3;
    
    if(p->D38>0)
	a->L(i,j,k) += (a->rhsvec.V[n])*PORVAL3;
    
    
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void nhflow_momentum_RK3co::timecheck(lexer *p,fdm *a,ghostcell *pgc,field &u,field &v,field &w)
{
}

void nhflow_momentum_RK3co::utimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void nhflow_momentum_RK3co::vtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void nhflow_momentum_RK3co::wtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

