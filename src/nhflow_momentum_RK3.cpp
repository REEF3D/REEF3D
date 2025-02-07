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

#include"nhflow_momentum_RK3.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_bcmom.h"
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
#include"net.h"
#include"6DOF.h"
#include"nhflow_forcing.h"
#include"wind_f.h"
#include"wind_v.h"

#define WLVL (fabs(WL(i,j))>(1.0*p->A544)?WL(i,j):1.0e20)

nhflow_momentum_RK3::nhflow_momentum_RK3(lexer *p, fdm_nhf *d, ghostcell *pgc, sixdof *pp6dof, vrans* ppvrans, 
                                                      vector<net*>& ppnet, nhflow_forcing *ppnhfdf)
                                                    : nhflow_momentum_func(p,d,pgc), WLRK1(p), WLRK2(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_uh=14;
	gcval_vh=15;
	gcval_wh=16;
    
    
    p->Darray(UHRK1,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VHRK1,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WHRK1,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(UHRK2,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VHRK2,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WHRK2,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(UHDIFF,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VHDIFF,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WHDIFF,p->imax*p->jmax*(p->kmax+2));
    
    p6dof=pp6dof;
    pnhfdf = ppnhfdf;
    pvrans = ppvrans;
    pnet = ppnet;
    
    if(p->A570==0)
    pwind = new wind_v(p);
    
    if(p->A570>0)
    pwind = new wind_f(p);
}

nhflow_momentum_RK3::~nhflow_momentum_RK3()
{
}

void nhflow_momentum_RK3::start(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, nhflow_signal_speed *pss,
                                     nhflow_reconstruct *precon, nhflow_convection *pconvec, nhflow_diffusion *pnhfdiff, 
                                     nhflow_pressure *ppress, solver *ppoissonsolv, solver *psolv, nhflow *pnhf, nhflow_fsf *pfsf, nhflow_turbulence *pnhfturb, vrans *pvrans)
{	
    pflow->discharge_nhflow(p,d,pgc);
    pflow->inflow_nhflow(p,d,pgc,d->U,d->V,d->W,d->UH,d->VH,d->WH);
    pflow->rkinflow_nhflow(p,d,pgc,d->U,d->V,d->W,UHRK1,VHRK1,WHRK1);
    pflow->rkinflow_nhflow(p,d,pgc,d->U,d->V,d->W,UHRK2,VHRK2,WHRK2);
		
//Step 1
//--------------------------------------------------------
    p->RK_alpha = 1.0;
    
    sigma_update(p,d,pgc,d->WL);
    reconstruct(p,d,pgc,pfsf,pss,precon,d->WL,d->U,d->V,d->W,d->UH,d->VH,d->WH);
    
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
    
    // FSF
    starttime=pgc->timer();
    
    pconvec->start(p,d,4,d->eta);
    pfsf->rk3_step1(p, d, pgc, pflow, d->UH, d->VH, d->WH, WLRK1, WLRK2, 1.0);
    omega_update(p,d,pgc,WLRK1,d->U,d->V,d->W);
    
    p->fsftime+=pgc->timer()-starttime;
    
	// U
	starttime=pgc->timer();

	pnhfturb->isource(p,d);
	pflow->isource_nhflow(p,d,pgc,pvrans); 
	ppress->upgrad(p,d,WLRK1);
    p6dof->isource(p,d,pgc,WLRK1);
	irhs(p,d,pgc);
    pwind->wind_forcing_nhf_x(p,d,pgc,d->U,d->V, d->F, WLRK1, d->eta);
    roughness_u(p,d,d->U,d->F,WLRK1);
	pconvec->start(p,d,1,WLRK1);
	pnhfdiff->diff_u(p,d,pgc,psolv,UHDIFF,d->UH,d->UH,d->VH,d->WH,WLRK1,1.0);

	LOOP
	UHRK1[IJK] = UHDIFF[IJK]
				+ p->dt*CPORNH*d->F[IJK];

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pnhfturb->jsource(p,d);
	pflow->jsource_nhflow(p,d,pgc,pvrans); 
    ppress->vpgrad(p,d,WLRK1);
    p6dof->jsource(p,d,pgc,WLRK1);
	jrhs(p,d,pgc);
    pwind->wind_forcing_nhf_y(p,d,pgc,d->U,d->V, d->G, WLRK1, d->eta);
    roughness_v(p,d,d->V,d->G,WLRK1);
    pconvec->start(p,d,2,WLRK1);
	pnhfdiff->diff_v(p,d,pgc,psolv,VHDIFF,d->VH,d->UH,d->VH,d->WH,WLRK1,1.0);

	LOOP
	VHRK1[IJK] = VHDIFF[IJK]
				+ p->dt*CPORNH*d->G[IJK];

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pnhfturb->ksource(p,d);
	//pflow->ksource_nhflow(p,d,pgc,pvrans); 
	ppress->wpgrad(p,d,WLRK1);
	krhs(p,d,pgc);
	pconvec->start(p,d,3,WLRK1);
	pnhfdiff->diff_w(p,d,pgc,psolv,WHDIFF,d->WH,d->UH,d->VH,d->WH,WLRK1,1.0);

	LOOP
	WHRK1[IJK] = WHDIFF[IJK]
				+ p->dt*CPORNH*d->H[IJK];
	
    p->wtime=pgc->timer()-starttime;
    
    velcalc(p,d,pgc,UHRK1,VHRK1,WHRK1,WLRK1);
    
    pnhfdf->forcing(p, d, pgc, p6dof, pvrans, pnet, 0, 1.0, UHRK1, VHRK1, WHRK1, WLRK1, 0);
    
	ppress->start(p,d,ppoissonsolv,pgc,pflow,WLRK1,UHRK1,VHRK1,WHRK1,1.0);
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
    p->RK_alpha = 0.25;
    
    sigma_update(p,d,pgc,WLRK1);
    reconstruct(p,d,pgc,pfsf,pss,precon,WLRK1,d->U,d->V,d->W,UHRK1,VHRK1,WHRK1);
	
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
    
    // FSF
    starttime=pgc->timer();
    
    pconvec->start(p,d,4,WLRK1);
    pfsf->rk3_step2(p, d, pgc, pflow, d->UH, d->VH, d->WH, WLRK1, WLRK2, 0.25);
    omega_update(p,d,pgc,WLRK2,d->U,d->V,d->W);
    
    p->fsftime+=pgc->timer()-starttime;
    
	// U
	starttime=pgc->timer();

	pnhfturb->isource(p,d);
	//pflow->isource(p,a,pgc,pvrans);
	ppress->upgrad(p,d,WLRK2);
    p6dof->isource(p,d,pgc,WLRK2);
	irhs(p,d,pgc);
    pwind->wind_forcing_nhf_x(p,d,pgc,d->U,d->V, d->F, WLRK2, d->eta);
    roughness_u(p,d,d->U,d->F,WLRK2);
    pconvec->start(p,d,1,WLRK2);
	pnhfdiff->diff_u(p,d,pgc,psolv,UHDIFF,UHRK1,UHRK1,VHRK1,WHRK1,WLRK2,0.25);

	LOOP
	UHRK2[IJK] = 0.75*d->UH[IJK] + 0.25*UHDIFF[IJK]
				+ 0.25*p->dt*CPORNH*d->F[IJK];
                
    p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	pnhfturb->jsource(p,d);
	//pflow->jsource(p,a,pgc,pvrans);
	ppress->vpgrad(p,d,WLRK2);
    p6dof->jsource(p,d,pgc,WLRK2);
	jrhs(p,d,pgc);
    pwind->wind_forcing_nhf_y(p,d,pgc,d->U,d->V, d->G, WLRK2, d->eta);
    roughness_v(p,d,d->V,d->G,WLRK2);
	pconvec->start(p,d,2,WLRK2);
	pnhfdiff->diff_v(p,d,pgc,psolv,VHDIFF,VHRK1,UHRK1,VHRK1,WHRK1,WLRK2,0.25);

	LOOP
	VHRK2[IJK] = 0.75*d->VH[IJK] + 0.25*VHDIFF[IJK]
				+ 0.25*p->dt*CPORNH*d->G[IJK];
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pnhfturb->ksource(p,d);
	//pflow->ksource(p,a,pgc,pvrans);
	ppress->wpgrad(p,d,WLRK2);
	krhs(p,d,pgc);
	pconvec->start(p,d,3,WLRK2);
	pnhfdiff->diff_w(p,d,pgc,psolv,WHDIFF,WHRK1,UHRK1,VHRK1,WHRK1,WLRK2,0.25);

	LOOP
	WHRK2[IJK] = 0.75*d->WH[IJK] + 0.25*WHDIFF[IJK]
				+ 0.25*p->dt*CPORNH*d->H[IJK];

    p->wtime+=pgc->timer()-starttime;
    
    velcalc(p,d,pgc,UHRK2,VHRK2,WHRK2,WLRK2);
    
    pnhfdf->forcing(p, d, pgc, p6dof, pvrans, pnet, 1, 0.25, UHRK2, VHRK2, WHRK2, WLRK2, 0);
    
	ppress->start(p,d,ppoissonsolv,pgc,pflow,WLRK2,UHRK2,VHRK2,WHRK2,0.25);
    velcalc(p,d,pgc,UHRK2,VHRK2,WHRK2,WLRK2);
    
	pflow->U_relax(p,pgc,d->U,UHRK2);
    pflow->V_relax(p,pgc,d->V,VHRK2);
    pflow->W_relax(p,pgc,d->W,WHRK2);

	pflow->P_relax(p,pgc,d->P);

	pgc->start4V(p,UHRK2,gcval_uh);
    pgc->start4V(p,VHRK2,gcval_vh);
    pgc->start4V(p,WHRK2,gcval_wh);
    clearrhs(p,d,pgc);

//Step 3
//--------------------------------------------------------
    p->RK_alpha = 2.0/3.0;
    
    sigma_update(p,d,pgc,WLRK2);
    reconstruct(p,d,pgc,pfsf,pss,precon,WLRK2,d->U,d->V,d->W,UHRK2,VHRK2,WHRK2);
    
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
    
    // FSF
    starttime=pgc->timer();
    
    pconvec->start(p,d,4,WLRK2);
    pfsf->rk3_step3(p, d, pgc, pflow, d->UH, d->VH, d->WH, WLRK1, WLRK2, 2.0/3.0);
    omega_update(p,d,pgc,d->WL,d->U,d->V,d->W);
    
    p->fsftime+=pgc->timer()-starttime;
    
	// U
	starttime=pgc->timer();

	pnhfturb->isource(p,d);
	//pflow->isource(p,a,pgc,pvrans);
	ppress->upgrad(p,d,d->WL);
    p6dof->isource(p,d,pgc,d->WL);
	irhs(p,d,pgc);
    pwind->wind_forcing_nhf_x(p,d,pgc,d->U,d->V, d->F, d->WL, d->eta);
    roughness_u(p,d,d->U,d->F,d->WL);
    pconvec->start(p,d,1,d->WL);
	pnhfdiff->diff_u(p,d,pgc,psolv,UHDIFF,UHRK2,UHRK2,VHRK2,WHRK2,d->WL,2.0/3.0);

	LOOP
	d->UH[IJK] = (1.0/3.0)*d->UH[IJK] + (2.0/3.0)*UHDIFF[IJK]
				+ (2.0/3.0)*p->dt*CPORNH*d->F[IJK];
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pnhfturb->jsource(p,d);
	//pflow->jsource(p,a,pgc,pvrans);
	ppress->vpgrad(p,d,d->WL);
    p6dof->jsource(p,d,pgc,d->WL);
	jrhs(p,d,pgc);
    pwind->wind_forcing_nhf_y(p,d,pgc,d->U,d->V, d->G, d->WL, d->eta);
    roughness_v(p,d,d->V,d->G,d->WL);
	pconvec->start(p,d,2,d->WL);
	pnhfdiff->diff_v(p,d,pgc,psolv,VHDIFF,VHRK2,UHRK2,VHRK2,WHRK2,d->WL,2.0/3.0);

	LOOP
	d->VH[IJK] = (1.0/3.0)*d->VH[IJK] + (2.0/3.0)*VHDIFF[IJK]
				+ (2.0/3.0)*p->dt*CPORNH*d->G[IJK];
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pnhfturb->ksource(p,d);
	//pflow->ksource(p,a,pgc,pvrans);
	ppress->wpgrad(p,d,d->WL);
	krhs(p,d,pgc);
	pconvec->start(p,d,3,d->WL);
	pnhfdiff->diff_u(p,d,pgc,psolv,WHDIFF,WHRK2,UHRK2,VHRK2,WHRK2,d->WL,2.0/3.0);

	LOOP
	d->WH[IJK] = (1.0/3.0)*d->WH[IJK] + (2.0/3.0)*WHDIFF[IJK]
				+ (2.0/3.0)*p->dt*CPORNH*d->H[IJK];
	
    p->wtime+=pgc->timer()-starttime;
    
    velcalc(p,d,pgc,d->UH,d->VH,d->WH,d->WL);
    
    pnhfdf->forcing(p, d, pgc, p6dof, pvrans, pnet, 2, 2.0/3.0, d->UH, d->VH, d->WH, d->WL, 1);
    
    ppress->start(p,d,ppoissonsolv,pgc,pflow,d->WL,d->UH,d->VH,d->WH,2.0/3.0);
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



