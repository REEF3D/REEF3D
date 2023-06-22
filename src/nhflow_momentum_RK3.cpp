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
#include"nhflow_reconstruct.h"
#include"nhflow_fsf_reconstruct.h"
#include"nhflow_convection.h"
#include"nhflow_signal_speed.h"
#include"diffusion.h"
#include"nhflow_pressure.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"nhflow.h"
#include"nhflow.h"
#include"nhflow_fsf.h"
#include"nhflow_turbulence.h"
#include"vrans.h"

nhflow_momentum_RK3::nhflow_momentum_RK3(lexer *p, fdm_nhf *d, ghostcell *pgc)
                                                    : bcmom(p), nhflow_sigma(p), etark1(p), etark2(p)
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

void nhflow_momentum_RK3::start(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, nhflow_signal_speed *pss, nhflow_fsf_reconstruct *pfsfrecon, 
                                     nhflow_reconstruct *precon, nhflow_convection *pconvec, diffusion *pdiff, 
                                     nhflow_pressure *ppress, solver *psolv, nhflow *pnhf, nhflow_fsf *pfsf, nhflow_turbulence *pnhfturb, vrans *pvrans)
{	
    pflow->discharge_nhflow(p,d,pgc);
    pflow->inflow_nhflow(p,d,pgc,d->U,d->V,d->W);
    pflow->rkinflow_nhflow(p,d,pgc,URK1,VRK1,WRK1);
    pflow->rkinflow_nhflow(p,d,pgc,URK2,VRK2,WRK2);
		
//Step 1
//--------------------------------------------------------
    // reconstruct eta
    pfsfrecon->reconstruct_2D(p, pgc, d, d->eta, d->ETAs, d->ETAn, d->ETAe, d->ETAw);
    
    // reconstruct U 
    precon->reconstruct_3D_x(p, pgc, d, d->U, d->Us, d->Un);
    precon->reconstruct_3D_y(p, pgc, d, d->U, d->Ue, d->Uw);
    precon->reconstruct_3D_z(p, pgc, d, d->U, d->Ub, d->Ut);
    
    // reconstruct  V
    precon->reconstruct_3D_x(p, pgc, d, d->V, d->Vs, d->Vn);
    precon->reconstruct_3D_y(p, pgc, d, d->V, d->Ve, d->Vw);
    precon->reconstruct_3D_z(p, pgc, d, d->V, d->Vb, d->Vt);
    
    // reconstruct  W
    precon->reconstruct_3D_x(p, pgc, d, d->W, d->Ws, d->Wn);
    precon->reconstruct_3D_y(p, pgc, d, d->W, d->We, d->Ww);
    precon->reconstruct_3D_z(p, pgc, d, d->W, d->Wb, d->Wt);
    
    // reconstruct UH
    precon->reconstruct_3D_x(p, pgc, d, d->UH, d->UHs, d->UHn);
    precon->reconstruct_3D_y(p, pgc, d, d->UH, d->UHe, d->UHw);
    
    // reconstruct  VH
    precon->reconstruct_3D_x(p, pgc, d, d->VH, d->VHs, d->VHn);
    precon->reconstruct_3D_y(p, pgc, d, d->VH, d->VHe, d->VHw);
    
    // reconstruct  WH
    precon->reconstruct_3D_x(p, pgc, d, d->WH, d->WHs, d->WHn);
    precon->reconstruct_3D_y(p, pgc, d, d->WH, d->WHe, d->WHw);
    
    pss->signal_speed_update(p, pgc, d, d->Us, d->Un, d->Ve, d->Vw, d->Ds, d->Dn, d->De, d->Dw);
    
    
    // FSF
    pfsf->step1(p, d, pgc, pflow, d->U, d->V, d->W, etark1, etark2, 1.0);
    sigma_update(p,d,pgc,etark1,d->eta,1.0);
    omega_update(p,d,pgc,d->U,d->V,d->W,etark1,d->eta,1.0);
    
    pconvec->precalc(p,d,d->U,1,d->U,d->V,d->W,d->eta);
    
	// U
	starttime=pgc->timer();

	pnhfturb->isource(p,d);
	pflow->isource_nhflow(p,d,pgc,pvrans); 
	//bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,d,etark1,d->eta);
	irhs(p,d,pgc);
	pconvec->start(p,d,d->U,1,d->U,d->V,d->W,etark1);
	//pdiff->diff_u(p,a,pgc,psolv,udiff,a->u,a->v,a->w,1.0);

	LOOP
	URK1[IJK] = d->U[IJK]
				+ p->dt*CPORNH*d->F[IJK];

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pnhfturb->jsource(p,d);
	pflow->jsource_nhflow(p,d,pgc,pvrans); 
	//bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
    ppress->vpgrad(p,d,etark1,d->eta);
	jrhs(p,d,pgc);
    pconvec->start(p,d,d->V,2,d->U,d->V,d->W,etark1);
	//pdiff->diff_v(p,a,pgc,psolv,vdiff,a->u,a->v,a->w,1.0);

	LOOP
	VRK1[IJK] = d->V[IJK]
				+ p->dt*CPORNH*d->G[IJK];

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pnhfturb->ksource(p,d);
	//pflow->ksource_nhflow(p,d,pgc,pvrans); 
	//bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,d,etark1,d->eta);
	krhs(p,d,pgc);
	pconvec->start(p,d,d->W,3,d->U,d->V,d->W,etark1);
	//pdiff->diff_w(p,a,pgc,psolv,wdiff,a->u,a->v,a->w,1.0);

	LOOP
	WRK1[IJK] = d->W[IJK]
				+ p->dt*CPORNH*d->H[IJK];
	
    p->wtime=pgc->timer()-starttime;
    
    pgc->start1V(p, URK1, gcval_u);
    pgc->start2V(p, VRK1, gcval_v);
    pgc->start3V(p, WRK1, gcval_w);
    

    pfsf->kinematic_fsf(p,d,d->U,d->V,WRK1,etark1,d->eta,1.0);
    omega_update(p,d,pgc,d->U,d->V,WRK1,etark1,d->eta,1.0);

    //pflow->pressure_io(p,a,pgc);
	ppress->start(p,d,psolv,pgc,pflow,URK1,VRK1,WRK1,1.0);
	

    pflow->U_relax(p,pgc,URK1);
    pflow->V_relax(p,pgc,VRK1);
    pflow->W_relax(p,pgc,WRK1);

	pflow->P_relax(p,pgc,d->P);

	pgc->start1V(p,URK1,gcval_u);
    pgc->start2V(p,VRK1,gcval_v);
    pgc->start3V(p,WRK1,gcval_w);
    clearrhs(p,d,pgc);
    
    pfsf->kinematic_fsf(p,d,URK1,VRK1,WRK1,etark1,d->eta,1.0);    
    omega_update(p,d,pgc,URK1,VRK1,WRK1,etark1,d->eta,1.0);
    
    //pfsf->step1(p, d, pgc, pflow, URK1, VRK1, WRK1, etark1, etark2, 1.0);
    
    //pupdate->start(p,a,pgc);
    
//Step 2
//--------------------------------------------------------
	
    pfsf->step2(p, d, pgc, pflow, URK1, VRK1, WRK1, etark1, etark2, 0.25);
    sigma_update(p,d,pgc,etark2,etark1,0.25);
    omega_update(p,d,pgc,URK1,VRK1,WRK1,etark2,etark1,0.25);
    
    pconvec->precalc(p,d,d->U,1,URK1,VRK1,WRK1,etark1);
    
	// U
	starttime=pgc->timer();

	pnhfturb->isource(p,d);
	//pflow->isource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,d,etark2,etark1);
	irhs(p,d,pgc);
    pconvec->start(p,d,URK1,1,URK1,VRK1,WRK1,etark2);
	//pdiff->diff_u(p,a,pgc,psolv,udiff,urk1,vrk1,wrk1,1.0);

	LOOP
	URK2[IJK] = 0.75*d->U[IJK] + 0.25*URK1[IJK]
				+ 0.25*p->dt*CPORNH*d->F[IJK];
                
    p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	pnhfturb->jsource(p,d);
	//pflow->jsource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,d,etark2,etark1);
	jrhs(p,d,pgc);
	pconvec->start(p,d,VRK1,2,URK1,VRK1,WRK1,etark2);
	//pdiff->diff_v(p,a,pgc,psolv,vdiff,urk1,vrk1,wrk1,1.0);

	LOOP
	VRK2[IJK] = 0.75*d->V[IJK] + 0.25*VRK1[IJK]
				+ 0.25*p->dt*CPORNH*d->G[IJK];
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pnhfturb->ksource(p,d);
	//pflow->ksource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,d,etark2,etark1);
	krhs(p,d,pgc);
	pconvec->start(p,d,WRK1,2,URK1,VRK1,WRK1,etark2);
	//pdiff->diff_w(p,a,pgc,psolv,wdiff,urk1,vrk1,wrk1,1.0);

	LOOP
	WRK2[IJK] = 0.75*d->W[IJK] + 0.25*WRK1[IJK]
				+ 0.25*p->dt*CPORNH*d->H[IJK];

    p->wtime+=pgc->timer()-starttime;
    
    pgc->start1V(p,URK2,gcval_u);
    pgc->start2V(p,VRK2,gcval_v);
    pgc->start3V(p,WRK2,gcval_w);
    
    pfsf->kinematic_fsf(p,d,URK1,VRK1,WRK2,etark2,etark1,0.25); 
    omega_update(p,d,pgc,URK1,VRK1,WRK2,etark2,etark1,0.25);

    //pflow->pressure_io(p,a,pgc);
	ppress->start(p,d,psolv,pgc,pflow,URK2,VRK2,WRK2,0.25);
    
	
	pflow->U_relax(p,pgc,URK2);
    pflow->V_relax(p,pgc,VRK2);
    pflow->W_relax(p,pgc,WRK2);

	pflow->P_relax(p,pgc,d->P);

	pgc->start1V(p,URK2,gcval_u);
    pgc->start2V(p,VRK2,gcval_v);
    pgc->start3V(p,WRK2,gcval_w);
    clearrhs(p,d,pgc);
    
    pfsf->kinematic_fsf(p,d,URK2,VRK2,WRK2,etark2,etark1,0.25); 
    omega_update(p,d,pgc,URK2,VRK2,WRK2,etark2,etark1,0.25);
    
    //pfsf->step2(p, d, pgc, pflow, URK2, VRK2, WRK2, etark1, etark2, 0.25);
    
    //pupdate->start(p,a,pgc);

//Step 3
//--------------------------------------------------------
    
    pfsf->step3(p, d, pgc, pflow, URK2, VRK2, WRK2, etark1, etark2, 2.0/3.0);
    sigma_update(p,d,pgc,d->eta,etark2,2.0/3.0);
    omega_update(p,d,pgc,URK2,VRK2,WRK2,d->eta,etark2,2.0/3.0);
    
    pconvec->precalc(p,d,d->U,1,URK2,VRK2,WRK2,etark2);
    
	// U
	starttime=pgc->timer();

	pnhfturb->isource(p,d);
	//pflow->isource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,d,d->eta,etark2);
	irhs(p,d,pgc);
    pconvec->start(p,d,URK2,1,URK2,VRK2,WRK2,d->eta);
	//pdiff->diff_u(p,a,pgc,psolv,udiff,urk2,vrk2,wrk2,1.0);

	LOOP
	d->U[IJK] = (1.0/3.0)*d->U[IJK] + (2.0/3.0)*URK2[IJK]
				+ (2.0/3.0)*p->dt*CPORNH*d->F[IJK];
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pnhfturb->jsource(p,d);
	//pflow->jsource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,d,d->eta,etark2);
	jrhs(p,d,pgc);
	pconvec->start(p,d,VRK2,2,URK2,VRK2,WRK2,d->eta);
	//pdiff->diff_v(p,a,pgc,psolv,vdiff,urk2,vrk2,wrk2,1.0);

	LOOP
	d->V[IJK] = (1.0/3.0)*d->V[IJK] + (2.0/3.0)*VRK2[IJK]
				+ (2.0/3.0)*p->dt*CPORNH*d->G[IJK];
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pnhfturb->ksource(p,d);
	//pflow->ksource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,d,d->eta,etark2);
	krhs(p,d,pgc);
	pconvec->start(p,d,WRK2,3,URK2,VRK2,WRK2,d->eta);
	//pdiff->diff_w(p,a,pgc,psolv,wdiff,urk2,vrk2,wrk2,1.0);

	LOOP
	d->W[IJK] = (1.0/3.0)*d->W[IJK] + (2.0/3.0)*WRK2[IJK]
				+ (2.0/3.0)*p->dt*CPORNH*d->H[IJK];
	
    p->wtime+=pgc->timer()-starttime;
    
    pgc->start1V(p,d->U,gcval_u);
    pgc->start2V(p,d->V,gcval_v);
    pgc->start3V(p,d->W,gcval_w);
    
    pfsf->kinematic_fsf(p,d,URK2,VRK2,d->W,d->eta,etark2,2.0/3.0);
    omega_update(p,d,pgc,URK2,VRK2,d->W,d->eta,etark2,2.0/3.0);

	//pflow->pressure_io(p,a,pgc);
    ppress->start(p,d,psolv,pgc,pflow,d->U,d->V,d->W,2.0/3.0);
	
	pflow->U_relax(p,pgc,d->U);
    pflow->V_relax(p,pgc,d->V);
    pflow->W_relax(p,pgc,d->W);

	pflow->P_relax(p,pgc,d->P);

	pgc->start1V(p,d->U,gcval_u);
    pgc->start2V(p,d->V,gcval_v);
    pgc->start3V(p,d->W,gcval_w);
    clearrhs(p,d,pgc);
    
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta,etark2,2.0/3.0);
    omega_update(p,d,pgc,d->U,d->V,d->W,d->eta,etark2,2.0/3.0);
    
    //pfsf->step3(p, d, pgc, pflow, d->U, d->V, d->W, etark1, etark2, 2.0/3.0);
    //pupdate->start(p,a,pgc);
}

void nhflow_momentum_RK3::irhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    /*
	n=0;
	LOOP
	{
    d->maxF=MAX(fabs(d->rhsvec.V[n]),d->maxF);
	d->F[IJK] += (d->rhsvec.V[n])*PORVALNH;
	d->rhsvec.V[n]=0.0;
	++n;
	}*/
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
    a->maxH=MAX(fabs(a->rhsvec.V[n]),a->maxH);
    
    if(p->D38==0)
    a->H[IJK] += (a->rhsvec.V[n])*PORVAL3;
    
    if(p->D38>0)
	a->H[IJK] += (a->rhsvec.V[n])*PORVAL3;
    
    
	a->rhsvec.V[n]=0.0;
	++n;
	}
    */
}

void nhflow_momentum_RK3::clearrhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    
	n=0;
	LOOP
	{
	d->rhsvec.V[n]=0.0;
	++n;
	}
    
}

void nhflow_momentum_RK3::inidisc(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_fsf *pfsf)
{
    sigma_ini(p,d,pgc,d->eta);
    sigma_update(p,d,pgc,d->eta,d->eta,1.0);
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta,d->eta,1.0);
}
     



