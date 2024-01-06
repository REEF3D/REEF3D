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

#include"nhflow_momentum_RK2.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"nhflow_reconstruct.h"
#include"nhflow_convection.h"
#include"nhflow_signal_speed.h"
#include"nhflow_diffusion.h"
#include"nhflow_pressure.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"nhflow_fsf.h"
#include"nhflow_turbulence.h"
#include"vrans.h"
#include"6DOF.h"

#define WLVL (fabs(WL(i,j))>(1.0*p->A544)?WL(i,j):1.0e20)

nhflow_momentum_RK2::nhflow_momentum_RK2(lexer *p, fdm_nhf *d, ghostcell *pgc, sixdof *pp6dof)
                                                    : bcmom(p), nhflow_sigma(p), WLRK1(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_uh=20;
	gcval_vh=21;
	gcval_wh=22;
    
    p->Darray(UHRK1,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VHRK1,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WHRK1,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(UHDIFF,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VHDIFF,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WHDIFF,p->imax*p->jmax*(p->kmax+2));
    
    sigma_ini(p,d,pgc,d->eta);
    
    p6dof = pp6dof;
}

nhflow_momentum_RK2::~nhflow_momentum_RK2()
{
}

void nhflow_momentum_RK2::start(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, nhflow_signal_speed *pss, 
                                     nhflow_reconstruct *precon, nhflow_convection *pconvec, nhflow_diffusion *pnhfdiff, 
                                     nhflow_pressure *ppress, solver *psolv, nhflow *pnhf, nhflow_fsf *pfsf, nhflow_turbulence *pnhfturb, vrans *pvrans)
{	
    pflow->discharge_nhflow(p,d,pgc);
    pflow->inflow_nhflow(p,d,pgc,d->U,d->V,d->W,d->UH,d->VH,d->WH);
    pflow->rkinflow_nhflow(p,d,pgc,d->U,d->V,d->W,UHRK1,VHRK1,WHRK1);
		
//Step 1
//--------------------------------------------------------

    sigma_update(p,d,pgc,d->WL);
    reconstruct(p,d,pgc,pfsf,pss,precon,d->WL,d->U,d->V,d->W,d->UH,d->VH,d->WH);
    
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);   
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);

    // FSF
    starttime=pgc->timer();
    pconvec->start(p,d,4,d->WL);

    pfsf->rk2_step1(p, d, pgc, pflow, d->UH, d->VH, d->WH, WLRK1, WLRK1, 1.0);
    //sigma_update(p,d,pgc,WLRK1);
    omega_update(p,d,pgc,WLRK1,d->U,d->V,d->W);
    p->fsftime+=pgc->timer()-starttime;
    
	// U
	starttime=pgc->timer();

	pnhfturb->isource(p,d);
	pflow->isource_nhflow(p,d,pgc,pvrans); 
	//bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,d,WLRK1);
    p6dof->isource(p,d,pgc);
	irhs(p,d,pgc);
	pconvec->start(p,d,1,WLRK1);
	pnhfdiff->diff_u(p,d,pgc,psolv,UHDIFF,d->UH,d->UH,d->VH,d->WH,1.0);

	LOOP
	UHRK1[IJK] = UHDIFF[IJK]
				+ p->dt*CPORNH*d->F[IJK];

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pnhfturb->jsource(p,d);
	pflow->jsource_nhflow(p,d,pgc,pvrans); 
	//bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
    ppress->vpgrad(p,d,WLRK1);
    p6dof->jsource(p,d,pgc);
	jrhs(p,d,pgc);
    pconvec->start(p,d,2,WLRK1);
	pnhfdiff->diff_v(p,d,pgc,psolv,VHDIFF,d->VH,d->UH,d->VH,d->WH,1.0);

	LOOP
	VHRK1[IJK] = VHDIFF[IJK]
				+ p->dt*CPORNH*d->G[IJK];

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pnhfturb->ksource(p,d);
	//pflow->ksource_nhflow(p,d,pgc,pvrans); 
	//bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,d,WLRK1);
	krhs(p,d,pgc);
	pconvec->start(p,d,3,WLRK1);
	pnhfdiff->diff_w(p,d,pgc,psolv,WHDIFF,d->WH,d->UH,d->VH,d->WH,1.0);

	LOOP
	WHRK1[IJK] = WHDIFF[IJK]
				+ p->dt*CPORNH*d->H[IJK];
	
    p->wtime=pgc->timer()-starttime;

    // fsfcorr
    pfsf->ucorr(p,d,UHRK1,WLRK1,1.0);
    pfsf->vcorr(p,d,VHRK1,WLRK1,1.0);

    velcalc(p,d,pgc,UHRK1,VHRK1,WHRK1,WLRK1);
    
    //pflow->pressure_io(p,a,pgc);
	ppress->start(p,d,psolv,pgc,pflow,WLRK1,UHRK1,VHRK1,WHRK1,1.0);
    velcalc(p,d,pgc,UHRK1,VHRK1,WHRK1,WLRK1);

    pflow->U_relax(p,pgc,d->U,UHRK1);
    pflow->V_relax(p,pgc,d->V,VHRK1);
    pflow->W_relax(p,pgc,d->W,WHRK1);
	pflow->P_relax(p,pgc,d->P);

	pgc->start4V(p,UHRK1,gcval_uh);
    pgc->start4V(p,VHRK1,gcval_vh);
    pgc->start4V(p,WHRK1,gcval_wh);

    clearrhs(p,d,pgc);
    
//Step 2
//--------------------------------------------------------
 
    sigma_update(p,d,pgc,WLRK1);
    reconstruct(p,d,pgc,pfsf,pss,precon,WLRK1,d->U,d->V,d->W,UHRK1,VHRK1,WHRK1);
    
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
    
    // FSF
    starttime=pgc->timer();
    
    pconvec->start(p,d,4,WLRK1);
    pfsf->rk2_step2(p, d, pgc, pflow, UHRK1,VHRK1,WHRK1, WLRK1, WLRK1, 0.5);
    //sigma_update(p,d,pgc,d->WL);
    omega_update(p,d,pgc,d->WL,d->U,d->V,d->W);
    
    p->fsftime+=pgc->timer()-starttime;
    
	// U
	starttime=pgc->timer();

	pnhfturb->isource(p,d);
	//pflow->isource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,d,d->WL);
    p6dof->isource(p,d,pgc);
	irhs(p,d,pgc);
    pconvec->start(p,d,1,d->WL);
	pnhfdiff->diff_u(p,d,pgc,psolv,UHDIFF,UHRK1,UHRK1,VHRK1,WHRK1,0.5);

	LOOP
	d->UH[IJK] = 0.5*d->UH[IJK] + 0.5*UHDIFF[IJK]
				+ 0.5*p->dt*CPORNH*d->F[IJK];
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pnhfturb->jsource(p,d);
	//pflow->jsource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,d,d->WL);
    p6dof->jsource(p,d,pgc);
	jrhs(p,d,pgc);
	pconvec->start(p,d,2,d->WL);
	pnhfdiff->diff_v(p,d,pgc,psolv,VHDIFF,VHRK1,UHRK1,VHRK1,WHRK1,0.5);

	LOOP
	d->VH[IJK] = 0.5*d->VH[IJK] + 0.5*VHDIFF[IJK]
				+ 0.5*p->dt*CPORNH*d->G[IJK];
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pnhfturb->ksource(p,d);
	//pflow->ksource(p,a,pgc,pvrans);
	//bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,d,d->WL);
	krhs(p,d,pgc);
	pconvec->start(p,d,3,d->WL);
	pnhfdiff->diff_w(p,d,pgc,psolv,WHDIFF,WHRK1,UHRK1,VHRK1,WHRK1,0.5);

	LOOP
	d->WH[IJK] = 0.5*d->WH[IJK] + 0.5*WHDIFF[IJK]
				+ 0.5*p->dt*CPORNH*d->H[IJK];
	
    p->wtime+=pgc->timer()-starttime;
    
    // fsfcorr
    pfsf->ucorr(p,d,d->UH,d->WL,0.5);
    pfsf->vcorr(p,d,d->VH,d->WL,0.5);
    
    velcalc(p,d,pgc,d->UH,d->VH,d->WH,d->WL);
    
	//pflow->pressure_io(p,a,pgc);
    ppress->start(p,d,psolv,pgc,pflow,d->WL,d->UH,d->VH,d->WH,0.5);
    
    velcalc(p,d,pgc,d->UH,d->VH,d->WH,d->WL);
	
	pflow->U_relax(p,pgc,d->U,d->UH);
    pflow->V_relax(p,pgc,d->V,d->VH);
    pflow->W_relax(p,pgc,d->W,d->WH);

	pflow->P_relax(p,pgc,d->P);

	pgc->start4V(p,d->UH,gcval_uh);
    pgc->start4V(p,d->VH,gcval_vh);
    pgc->start4V(p,d->WH,gcval_wh);
    
    clearrhs(p,d,pgc);
}

void nhflow_momentum_RK2::reconstruct(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_fsf *pfsf, nhflow_signal_speed *pss, 
                                     nhflow_reconstruct *precon, slice &WL, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    starttime=pgc->timer();
    
    // reconstruct eta
    precon->reconstruct_2D_x(p, pgc, d, d->eta, d->ETAs, d->ETAn);
    precon->reconstruct_2D_y(p, pgc, d, d->eta, d->ETAe, d->ETAw);
    precon->reconstruct_2D_WL(p, pgc, d);

    // reconstruct U 
    precon->reconstruct_3D_x(p, pgc, d, U, d->Us, d->Un);
    precon->reconstruct_3D_y(p, pgc, d, U, d->Ue, d->Uw);
    precon->reconstruct_3D_z(p, pgc, d, U, d->Ub, d->Ut);
    
    // reconstruct  V
    precon->reconstruct_3D_x(p, pgc, d, V, d->Vs, d->Vn);
    precon->reconstruct_3D_y(p, pgc, d, V, d->Ve, d->Vw);
    precon->reconstruct_3D_z(p, pgc, d, V, d->Vb, d->Vt);
    
    // reconstruct  W
    precon->reconstruct_3D_x(p, pgc, d, W, d->Ws, d->Wn);
    precon->reconstruct_3D_y(p, pgc, d, W, d->We, d->Ww);
    precon->reconstruct_3D_z(p, pgc, d, W, d->Wb, d->Wt);
    
    // reconstruct UH
    precon->reconstruct_3D_x(p, pgc, d, UH, d->UHs, d->UHn);
    precon->reconstruct_3D_y(p, pgc, d, UH, d->UHe, d->UHw);
    
    // reconstruct  VH
    precon->reconstruct_3D_x(p, pgc, d, VH, d->VHs, d->VHn);
    precon->reconstruct_3D_y(p, pgc, d, VH, d->VHe, d->VHw);
    
    // reconstruct  WH
    precon->reconstruct_3D_x(p, pgc, d, WH, d->WHs, d->WHn);
    precon->reconstruct_3D_y(p, pgc, d, WH, d->WHe, d->WHw);
    
    // wetdry
    pfsf->wetdry_fluxes(p,d,pgc,WL,U,V,W,UH,VH,WH);
    
    pss->signal_speed_update(p, pgc, d, d->Us, d->Un, d->Ve, d->Vw, d->Ds, d->Dn, d->De, d->Dw);
    
    p->recontime+=pgc->timer()-starttime;
}

void nhflow_momentum_RK2::velcalc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *UH, double *VH, double *WH, slice &WL)
{
    // Fr nuber limiter
    LOOP
    WETDRY
    {
    UH[IJK] = MIN(UH[IJK], p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    VH[IJK] = MIN(VH[IJK], p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    WH[IJK] = MIN(WH[IJK], p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));      
    
    UH[IJK] = MAX(UH[IJK], -p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    VH[IJK] = MAX(VH[IJK], -p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    WH[IJK] = MAX(WH[IJK], -p->A531*WL(i,j)*sqrt(9.81*WL(i,j))); 
    }
    
    
    LOOP
    {
    d->U[IJK] = UH[IJK]/WLVL;
    d->V[IJK] = VH[IJK]/WLVL;
    d->W[IJK] = WH[IJK]/WLVL;       
    }
    
    if(p->A520==0)
    LOOP
    {
    d->W[IJK] = 0.0;  
    //WH[IJK] = 0.0;
    }
    
    
    LOOP
    if(p->wet[IJ]==0)
    {
    d->U[IJK] = 0.0;
    d->V[IJK] = 0.0;
    d->W[IJK] = 0.0;
    }
    
    pgc->start4V(p,d->U,gcval_u);
    pgc->start4V(p,d->V,gcval_v);
    pgc->start4V(p,d->W,gcval_w);
}

void nhflow_momentum_RK2::irhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
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

void nhflow_momentum_RK2::jrhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
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

void nhflow_momentum_RK2::krhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
}

void nhflow_momentum_RK2::clearrhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
	n=0;
	LOOP
	{
	d->rhsvec.V[n]=0.0;
	++n;
	}
}

void nhflow_momentum_RK2::inidisc(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_fsf *pfsf)
{
    pfsf->wetdry(p,d,pgc,d->U,d->V,d->W,d->WL);
    sigma_update(p,d,pgc,d->WL);
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
    velcalc(p,d,pgc,d->UH,d->VH,d->WH,d->WL);
}
     



