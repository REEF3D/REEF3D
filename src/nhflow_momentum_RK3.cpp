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

#include"nhflow_momentum_RK3.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"nhflow_convection.h"
#include"diffusion.h"
#include"nhflow_pressure.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"nhflow.h"
#include"nhflow.h"
#include"nhflow_fsf.h"
#include"vrans.h"

nhflow_momentum_RK3::nhflow_momentum_RK3(lexer *p, fdm_nhf *d, ghostcell *pgc)
                                                    :bcmom(p), etark1(p), etark2(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    p->Darray(URK1,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VRK1,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WRK1,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(URK2,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VRK2,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WRK2,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(UDIFF,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VDIFF,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WDIFF,p->imax*p->jmax*(p->kmax+2));
}

nhflow_momentum_RK3::~nhflow_momentum_RK3()
{
}

void nhflow_momentum_RK3::start(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, nhflow_convection *pconvec, diffusion *pdiff, 
                                     nhflow_pressure *ppress, solver *psolv, nhflow *pnhf, nhflow_fsf *pfsf, vrans *pvrans)
{	
    pflow->discharge_nhflow(p,d,pgc);
    pflow->inflow_nhflow(p,d,pgc,d->U,d->V,d->W);
    pflow->rkinflow_nhflow(p,d,pgc,URK1,VRK1,WRK1);
    pflow->rkinflow_nhflow(p,d,pgc,URK2,VRK2,WRK2);
		
//Step 1
//--------------------------------------------------------

    pfsf->step1(p, d, pgc, pflow, d->U, d->V, d->W, etark1, etark2, 1.0);
    

	// U
	starttime=pgc->timer();

	//pturb->isource(p,a);
	pflow->isource_nhflow(p,d,pgc,pvrans); 
	//bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,d,etark1,d->eta);
	irhs(p,d,pgc);
	pconvec->start(p,d,d->U,1,d->U,d->V,d->W);
	//pdiff->diff_u(p,a,pgc,psolv,udiff,a->u,a->v,a->w,1.0);

	LOOP
	URK1[IJK] = d->U[IJK]
				+ p->dt*CPORNH*d->F[IJK];

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();
/*
	//pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,etark1,a->eta);
	jrhs(p,d,pgc);
	pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
	//pdiff->diff_v(p,a,pgc,psolv,vdiff,a->u,a->v,a->w,1.0);

	LOOP
	VRK1[IJK] = d->V[IJK]
				+ p->dt*CPORNH*d->G[IJK];

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	//pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,etark1,a->eta);
	krhs(p,d,pgc);
	pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
	//pdiff->diff_w(p,a,pgc,psolv,wdiff,a->u,a->v,a->w,1.0);

	LOOP
	WRK1[IJK] = d->W[IJK]
				+ p->dt*CPORNH*d->H[IJK];
	
    p->wtime=pgc->timer()-starttime;
    
    pgc->start1(p,urk1,gcval_u);
	pgc->start2(p,vrk1,gcval_v);
    pgc->start3(p,wrk1,gcval_w);
    
    pnh->kinematic_fsf(p,a,a->u,a->v,wrk1,etark1,a->eta,1.0);
    p->omega_update(p,a,pgc,a->u,a->v,wrk1,etark1,etark1,1.0);

    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk1, vrk1, wrk1, 1.0);
    
    start(lexer *p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow,
                        double *U, double *V, double *W, double alpha)
	
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

	//pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,etark2,etark1);
	irhs(p,d,pgc);
	pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
	//pdiff->diff_u(p,a,pgc,psolv,udiff,urk1,vrk1,wrk1,1.0);

	LOOP
	URK2[IJK] = 0.75*d->U[IJK] + 0.25*URK1[IJK]
				+ 0.25*p->dt*CPORNH*d->F[IJK];
                
    p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	//pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,etark2,etark1);
	jrhs(p,d,pgc);
	pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
	//pdiff->diff_v(p,a,pgc,psolv,vdiff,urk1,vrk1,wrk1,1.0);

	LOOP
	VRK2[IJK] = 0.75*d->V[IJK] + 0.25*VRK1[IJK]
				+ 0.25*p->dt*CPORNH*d->G[IJK];
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	//pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,etark2,etark1);
	krhs(p,d,pgc);
	pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
	//pdiff->diff_w(p,a,pgc,psolv,wdiff,urk1,vrk1,wrk1,1.0);

	LOOP
	WRK2[IJK] = 0.75*d->W[IJK] + 0.25*WRK1[IJK]
				+ 0.25*p->dt*CPORNH*d->H[IJK];

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

	//pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,etark2);
	irhs(p,d,pgc);
	pconvec->start(p,a,urk2,1,urk2,vrk2,wrk2);
	//pdiff->diff_u(p,a,pgc,psolv,udiff,urk2,vrk2,wrk2,1.0);

	LOOP
	d->U[IJK] = (1.0/3.0)*d->U[IJK] + (2.0/3.0)*URK2[IJK]
				+ (2.0/3.0)*p->dt*CPORNH*d->F[IJK];
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	//pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,etark2);
	jrhs(p,d,pgc);
	pconvec->start(p,a,vrk2,2,urk2,vrk2,wrk2);
	//pdiff->diff_v(p,a,pgc,psolv,vdiff,urk2,vrk2,wrk2,1.0);

	LOOP
	d->V[IJK] = (1.0/3.0)*d->V[IJK] + (2.0/3.0)*VRK2[IJK]
				+ (2.0/3.0)*p->dt*CPORNH*d->G[IJK];
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	//pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,etark2);
	krhs(p,d,pgc);
	pconvec->start(p,a,wrk2,3,urk2,vrk2,wrk2);
	//pdiff->diff_w(p,a,pgc,psolv,wdiff,urk2,vrk2,wrk2,1.0);

	LOOP
	d->W[IJK] = (1.0/3.0)*d->W[IJK] + (2.0/3.0)*WRK2[IJK]
				+ (2.0/3.0)*p->dt*CPORNH*d->H[IJK];
	
    p->wtime+=pgc->timer()-starttime;
    
    pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
    pgc->start3(p,a->w,gcval_w);
    
    pnh->kinematic_fsf(p,a,urk2,vrk2,a->w,a->eta,etark2,2.0/3.0);
    p->omega_update(p,a,pgc,urk2,vrk2,a->w,a->eta,etark2,2.0/3.0);

	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, a->u, a->v, a->w, 2.0/3.0);
	
	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
    
    pnh->kinematic_fsf(p,a,a->u,a->v,a->w,a->eta,etark2,2.0/3.0);
    p->omega_update(p,a,pgc,a->u,a->v,a->w,a->eta,etark2,2.0/3.0);
    
    pupdate->start(p,a,pgc);*/
}

void nhflow_momentum_RK3::irhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
	n=0;
	LOOP
	{
    d->maxF=MAX(fabs(d->rhsvec.V[n]),d->maxF);
	d->F[IJK] += (d->rhsvec.V[n])*PORVALNH;
	d->rhsvec.V[n]=0.0;
	++n;
	}
}

void nhflow_momentum_RK3::jrhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    /*
	n=0;
	VLOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n]),a->maxG);
	a->G[IJK] += (a->rhsvec.V[n])*PORVAL2;
	a->rhsvec.V[n]=0.0;
	++n;
	}*/
}

void nhflow_momentum_RK3::krhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    /*
	n=0;
	WLOOP
	{
    a->maxH=MAX(fabs(a->rhsvec.V[n] + a->gk),a->maxH);
    
    if(p->D38==0)
    a->H[IJK] += (a->rhsvec.V[n] + a->gk)*PORVAL3;
    
    if(p->D38>0)
	a->H[IJK] += (a->rhsvec.V[n])*PORVAL3;
    
    
	a->rhsvec.V[n]=0.0;
	++n;
	}
    */
}



