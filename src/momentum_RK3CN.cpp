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
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#include"momentum_RK3CN.h"
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

momentum_RK3CN::momentum_RK3CN(lexer *p, fdm *a, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, 
                                                    ioflow *pioflow, fsi *ppfsi)
                                                    :momentum_forcing(p),bcmom(p),udiff(p),vdiff(p),wdiff(p),urk1(p),urk2(p),vrk1(p),
                                                    vrk2(p),wrk1(p),wrk2(p),fx(p),fy(p),fz(p)
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
    
    if(p->W90==1 && p->F300==0)
	pupdate = new fluid_update_rheology(p,a);
    

}

momentum_RK3CN::~momentum_RK3CN()
{
}

void momentum_RK3CN::start(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof, vector<net*>& pnet)
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
        addirhs(p,a,pgc,a->u,a->u,a->v,a->w,1.0);
	pdiff->diff_u(p,a,pgc,psolv,urk1,a->u,a->u,a->v,a->w,8.0/15.0);

        p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	
        jrhs(p,a,pgc,a->v,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
        addjrhs(p,a,pgc,a->v,a->u,a->v,a->w,1.0);        
        pdiff->diff_v(p,a,pgc,psolv,vrk1,a->v,a->u,a->v,a->w,8.0/15.0);

        p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);

	krhs(p,a,pgc,a->w,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
        addkrhs(p,a,pgc,a->w,a->u,a->v,a->w,1.0);
	pdiff->diff_w(p,a,pgc,psolv,wrk1,a->w,a->u,a->v,a->w,8.0/15.0);
	
        p->wtime=pgc->timer()-starttime;
    
        pgc->start1(p,urk1,gcval_u);
	pgc->start2(p,vrk1,gcval_v);
        pgc->start3(p,wrk1,gcval_w);

    
    momentum_forcing_start(a, p, pgc, p6dof, pvrans, pnet, pfsi,
                           urk1, vrk1, wrk1, fx, fy, fz, 0, 1.0, false);
    
    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk1, vrk1, wrk1, 8.0/15.0);

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
	
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);

	irhs(p,a,pgc,urk1,urk1,vrk1,wrk1,1.0);
        addirhs(p,a,pgc,urk1,urk1,vrk1,wrk1,1.0);
	pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
        addirhs(p,a,pgc,urk1,urk1,vrk1,wrk1,25.0/8.0);
        pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
        addirhs(p,a,pgc,a->u,a->u,a->v,a->w,-17.0/8.0);
	pdiff->diff_u(p,a,pgc,psolv,urk2,urk1,urk1,vrk1,wrk1,2.0/15.0);
                
        p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);

	jrhs(p,a,pgc,vrk1,urk1,vrk1,wrk1,1.0);
        addjrhs(p,a,pgc,vrk1,urk1,vrk1,wrk1,1.0);
	pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
        addjrhs(p,a,pgc,vrk1,urk1,vrk1,wrk1,25.0/8.0);
        pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
        addjrhs(p,a,pgc,a->v,a->u,a->v,a->w,-17.0/8.0);
	pdiff->diff_v(p,a,pgc,psolv,vrk2,vrk1,urk1,vrk1,wrk1,2.0/15.0);
	
        p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
 
	krhs(p,a,pgc,wrk1,urk1,vrk1,wrk1,1.0);
        addkrhs(p,a,pgc,wrk1,urk1,vrk1,wrk1,1.0);
	pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
        addkrhs(p,a,pgc,wrk1,urk1,vrk1,wrk1,25.0/8.0);
        pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
        addkrhs(p,a,pgc,a->w,a->u,a->v,a->w,-17/8.0);
	pdiff->diff_w(p,a,pgc,psolv,wrk2,wrk1,urk1,vrk1,wrk1,2.0/15.0);

        p->wtime+=pgc->timer()-starttime;
    
        pgc->start1(p,urk2,gcval_u);
	pgc->start2(p,vrk2,gcval_v);
            pgc->start3(p,wrk2,gcval_w);

    momentum_forcing_start(a, p, pgc, p6dof, pvrans, pnet, pfsi,
                           urk2, vrk2, wrk2, fx, fy, fz, 1, 0.25, false);

    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk2, vrk2, wrk2, 2.0/15.0);
	
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

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);

	irhs(p,a,pgc,urk2,urk2,vrk2,wrk2,1.0);
        addirhs(p,a,pgc,urk2,urk2,vrk2,wrk2,1.0);
	pconvec->start(p,a,urk2,1,urk2,vrk2,wrk2);
        addirhs(p,a,pgc,urk2,urk2,vrk2,wrk2,9.0/4.0);        
        pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
        addirhs(p,a,pgc,urk1,urk1,vrk1,wrk1,-5.0/4.0);
	pdiff->diff_u(p,a,pgc,psolv,a->u,urk2,urk2,vrk2,wrk2,1.0/3.0);
	
        p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);

	jrhs(p,a,pgc,vrk2,urk2,vrk2,wrk2,1.0);
        addjrhs(p,a,pgc,vrk2,urk2,vrk2,wrk2,1.0);
	pconvec->start(p,a,vrk2,2,urk2,vrk2,wrk2);
        addjrhs(p,a,pgc,vrk2,urk2,vrk2,wrk2,9.0/4.0);
        pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
        addjrhs(p,a,pgc,vrk1,urk1,vrk1,wrk1,-5.0/4.0);
	pdiff->diff_v(p,a,pgc,psolv,a->v,vrk2,urk2,vrk2,wrk2,1.0/3.0);
	
        p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);

	krhs(p,a,pgc,wrk2,urk2,vrk2,wrk2,1.0);
        addkrhs(p,a,pgc,wrk2,urk2,vrk2,wrk2,1.0);
	pconvec->start(p,a,wrk2,3,urk2,vrk2,wrk2);
        addkrhs(p,a,pgc,wrk2,urk2,vrk2,wrk2,9.0/4.0);
        pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
        addkrhs(p,a,pgc,wrk1,urk1,vrk1,wrk1,-5.0/4.0);
	pdiff->diff_w(p,a,pgc,psolv,a->w,wrk2,urk2,vrk2,wrk2,1.0/3.0);
	
        p->wtime+=pgc->timer()-starttime;

        pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);

    momentum_forcing_start(a, p, pgc, p6dof, pvrans, pnet, pfsi,
                           a->u, a->v, a->w, fx, fy, fz, 2, 2.0/3.0, true);

	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, a->u, a->v,a->w,1.0/3.0);
	
	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);

        pupdate->start(p,a,pgc);
}

void momentum_RK3CN::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
        double dens;
        n = 0;
	ULOOP
	{
        a->maxF=MAX(fabs(a->rhsvec.V[n] + a->gi),a->maxF);
        if (p->H10>0 && p->W90==0 && p->H3==2){
            dens = ((a->dro(i+1,j,k)+a->dro(i,j,k))/(a->ro(i+1,j,k)+a->ro(i,j,k)));}
        else {
            dens = 1.0;
        }
	a->F(i,j,k) += (a->rhsvec.V[n] + a->gi*dens + p->W29_x)*PORVAL1;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void momentum_RK3CN::addirhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
       n = 0;
       ULOOP
       {
       a->rhsvec.V[n] += alpha*a->F(i,j,k);
       a->F(i,j,k) = 0.0;
       ++n;
       }
}

void momentum_RK3CN::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	double dens;
	n = 0;
	VLOOP
	{
        a->maxG=MAX(fabs(a->rhsvec.V[n] + a->gj),a->maxG);
        if (p->H10>0 && p->W90==0 && p->H3==2){
            dens = ((a->dro(i,j+1,k)+a->dro(i,j,k))/(a->ro(i,j+1,k)+a->ro(i,j,k)));}
        else {
            dens = 1.0;
        }
	a->G(i,j,k) += (a->rhsvec.V[n] + a->gj*dens + p->W29_y)*PORVAL2;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void momentum_RK3CN::addjrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
       n = 0;
       VLOOP
       {
       a->rhsvec.V[n] += alpha*a->G(i,j,k);
       a->G(i,j,k) = 0.0;
       ++n;
       }
}

void momentum_RK3CN::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	double dens;
	n = 0;
	WLOOP
	{
        a->maxH=MAX(fabs(a->rhsvec.V[n] + a->gk),a->maxH);
        if (p->H10>0 && p->W90==0 && p->H3==2){
            dens = ((a->dro(i,j,k+1)+a->dro(i,j,k))/(a->ro(i,j,k+1)+a->ro(i,j,k)));}
        else {
            dens = 1.0;
        }
	a->H(i,j,k) += (a->rhsvec.V[n] + a->gk*dens + p->W29_z)*PORVAL3;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void momentum_RK3CN::addkrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
       n = 0;
       WLOOP
       {
       a->rhsvec.V[n] += alpha*a->H(i,j,k);
       a->H(i,j,k) = 0.0;
       ++n;
       }
}

void momentum_RK3CN::utimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_RK3CN::vtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_RK3CN::wtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

