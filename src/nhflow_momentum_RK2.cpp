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

#include"nhflow_momentum_RK2.h"
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

nhflow_momentum_RK2::nhflow_momentum_RK2(lexer *p, fdm *a, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, ioflow *pioflow,
                                                    nhflow *ppnh, nhflow_fsf *ppnhfsf)
                                                    :bcmom(p),udiff(p),vdiff(p),wdiff(p),urk1(p),vrk1(p),wrk1(p),
                                                    etark1(p)
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

nhflow_momentum_RK2::~nhflow_momentum_RK2()
{
}

void nhflow_momentum_RK2::start(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{	/*
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
	pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
		
//Step 1
//--------------------------------------------------------

    pnhfsf->step1(p, a, pgc, pflow, a->u, a->v, a->w, etark1, etark1, 1.0);

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans); 
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,etark1,a->eta);
	irhs(p,a,pgc,a->u,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
	pdiff->diff_u(p,a,pgc,psolv,udiff,a->u,a->v,a->w,1.0);

	ULOOP
	urk1(i,j,k) = a->u(i,j,k)
				+ p->dt*CPOR1*a->F(i,j,k);

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,etark1,a->eta);
	jrhs(p,a,pgc,a->v,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,a->u,a->v,a->w,1.0);

	VLOOP
	vrk1(i,j,k) = a->v(i,j,k)
				+ p->dt*CPOR2*a->G(i,j,k);

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,etark1,a->eta);
	krhs(p,a,pgc,a->w,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,a->u,a->v,a->w,1.0);

	WLOOP
	wrk1(i,j,k) = a->w(i,j,k)
				+ p->dt*CPOR3*a->H(i,j,k);
	
    p->wtime=pgc->timer()-starttime;
    
    pgc->start1(p,urk1,gcval_u);
	pgc->start2(p,vrk1,gcval_v);
    pgc->start3(p,wrk1,gcval_w);
    
    pnh->kinematic_fsf(p,a,a->u,a->v,wrk1,etark1,a->eta,1.0);
    p->omega_update(p,a,pgc,a->u,a->v,wrk1,etark1,etark1,1.0);

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
    
    pnhfsf->step3(p, a, pgc, pflow, urk1, vrk1, wrk1, etark1, etark1, 0.5);
    
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,etark1);
	irhs(p,a,pgc,urk1,urk1,vrk1,wrk1,0.5);
	pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
	pdiff->diff_u(p,a,pgc,psolv,udiff,urk1,vrk1,wrk1,1.0);

	ULOOP
	a->u(i,j,k) = 0.5*a->u(i,j,k) + 0.5*urk1(i,j,k)
				+ 0.5*p->dt*CPOR1*a->F(i,j,k);
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,etark1);
	jrhs(p,a,pgc,vrk1,urk1,vrk1,wrk1,0.5);
	pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,urk1,vrk1,wrk1,1.0);

	VLOOP
	a->v(i,j,k) = 0.5*a->v(i,j,k) + 0.5*vrk1(i,j,k)
				+ 0.5*p->dt*CPOR2*a->G(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,etark1);
	krhs(p,a,pgc,wrk1,urk1,vrk1,wrk1,0.5);
	pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,urk1,vrk1,wrk1,1.0);

	WLOOP
	a->w(i,j,k) = 0.5*a->w(i,j,k) + 0.5*wrk1(i,j,k)
				+ 0.5*p->dt*CPOR3*a->H(i,j,k);
	
    p->wtime+=pgc->timer()-starttime;
    
    pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
    pgc->start3(p,a->w,gcval_w);
    
    pnh->kinematic_fsf(p,a,urk1,vrk1,a->w,a->eta,etark1,0.5);
    p->omega_update(p,a,pgc,urk1,vrk1,a->w,a->eta,etark1,0.5);

	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, a->u, a->v, a->w, 0.5);
	
	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
    
    pnh->kinematic_fsf(p,a,a->u,a->v,a->w,a->eta,etark1,0.5);
    p->omega_update(p,a,pgc,a->u,a->v,a->w,a->eta,etark1,0.5);
    
    pupdate->start(p,a,pgc);*/
}

void nhflow_momentum_RK2::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	ULOOP
	{
    a->maxF=MAX(fabs(a->rhsvec.V[n]),a->maxF);
	a->F(i,j,k) += (a->rhsvec.V[n])*PORVAL1;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void nhflow_momentum_RK2::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	VLOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n]),a->maxG);
	a->G(i,j,k) += (a->rhsvec.V[n])*PORVAL2;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void nhflow_momentum_RK2::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	WLOOP
	{
    a->maxH=MAX(fabs(a->rhsvec.V[n] + a->gk),a->maxH);
    
    if(p->D38==0)
    a->H(i,j,k) += (a->rhsvec.V[n] + a->gk)*PORVAL3;
    
    if(p->D38>0)
	a->H(i,j,k) += (a->rhsvec.V[n])*PORVAL3;
    
    
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void nhflow_momentum_RK2::timecheck(lexer *p,fdm *a,ghostcell *pgc,field &u,field &v,field &w)
{
}

void nhflow_momentum_RK2::utimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void nhflow_momentum_RK2::vtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void nhflow_momentum_RK2::wtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

