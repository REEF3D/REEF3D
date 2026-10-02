/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
#include"vrans_nhflow.h"
#include"6DOF.h"
#include"nhflow_forcing.h"
#include"wind_f.h"
#include"wind_v.h"

#define WLVL (fabs(WL(i,j))>(1.0*p->A544)?WL(i,j):1.0e20)

nhflow_momentum_RK3::nhflow_momentum_RK3(lexer *p, fdm_nhf *d, ghostcell *pgc, sixdof *pp6dof, vrans_nhflow* ppvrans, 
                                                      nhflow_forcing *ppnhfdf)
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
                                     nhflow_pressure *ppress, solver *ppoissonsolv, solver *psolv, nhflow *pnhf, nhflow_fsf *pfsf, 
                                     nhflow_turbulence *pnhfturb, vrans_nhflow *pvrans)
{	
    nhflow_stage_obj S = {pflow,pss,precon,pconvec,pnhfdiff,ppress,ppoissonsolv,psolv,pnhf,pfsf,pnhfturb,pvrans};
    
    if(prun!=nullptr)
    {
    prun->step(p,d,pgc,this,S);
    return;
    }
    
    step_begin(p,d,pgc,S);
    
    for(int s=0; s<3; ++s)
    {
    phase_F(p,d,pgc,S,s);
    phase_M(p,d,pgc,S,s);
    phase_P(p,d,pgc,S,s);
    phase_E(p,d,pgc,S,s);
    }
}

slice& nhflow_momentum_RK3::stage_WL(fdm_nhf *d, int s)
{
    if(s==0)
    return WLRK1;
    if(s==1)
    return WLRK2;
    return d->WL;
}

double* nhflow_momentum_RK3::stage_UH(fdm_nhf *d, int s, int m)
{
    if(s==0)
    return m==0?UHRK1:(m==1?VHRK1:WHRK1);
    if(s==1)
    return m==0?UHRK2:(m==1?VHRK2:WHRK2);
    return m==0?d->UH:(m==1?d->VH:d->WH);
}

void nhflow_momentum_RK3::step_begin(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S)
{
    S.pflow->discharge_nhflow(p,d,pgc);
    S.pflow->inflow_nhflow(p,d,pgc,d->U,d->V,d->W,d->UH,d->VH,d->WH,d->WL);
    S.pflow->rkinflow_nhflow(p,d,pgc,d->U,d->V,d->W,UHRK1,VHRK1,WHRK1,WLRK1);
    S.pflow->rkinflow_nhflow(p,d,pgc,d->U,d->V,d->W,UHRK2,VHRK2,WHRK2,WLRK1);
}

// stage s: sigma, reconstruction, continuity flux, water level, omega
void nhflow_momentum_RK3::phase_F(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S, int s)
{
    ioflow *pflow = S.pflow;
    nhflow_fsf *pfsf = S.pfsf;
    nhflow_convection *pconvec = S.pconvec;
    
    if(s==0)
    {
    sigma_update(p,d,pgc,d->WL);
    S.pvrans->update(p,d,pgc,2.0/3.0,0);
    reconstruct(p,d,pgc,pfsf,S.pss,S.precon,d->WL,d->U,d->V,d->W,d->UH,d->VH,d->WH);
    
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
    
    // FSF
    starttime=pgc->timer();
    pconvec->start(p,d,4,d->eta,UHRK1);
    
    pfsf->rk3_step1(p, d, pgc, pflow, d->UH, d->VH, d->WH, WLRK1, WLRK2, 1.0);
    omega_update(p,d,pgc,WLRK1,d->U,d->V,d->W);
    p->fsftime+=pgc->timer()-starttime;
    }
    
    if(s==1)
    {
    sigma_update(p,d,pgc,WLRK1);
    S.pvrans->update(p,d,pgc,1.0,1);
    reconstruct(p,d,pgc,pfsf,S.pss,S.precon,WLRK1,d->U,d->V,d->W,UHRK1,VHRK1,WHRK1);
	
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
    
    // FSF
    starttime=pgc->timer();
    
    pconvec->start(p,d,4,WLRK1,UHRK2);
    pfsf->rk3_step2(p, d, pgc, pflow, d->UH, d->VH, d->WH, WLRK1, WLRK2, 0.25);
    omega_update(p,d,pgc,WLRK2,d->U,d->V,d->W);
    
    p->fsftime+=pgc->timer()-starttime;
    }
    
    if(s==2)
    {
    sigma_update(p,d,pgc,WLRK2);
    S.pvrans->update(p,d,pgc,0.25,2);
    reconstruct(p,d,pgc,pfsf,S.pss,S.precon,WLRK2,d->U,d->V,d->W,UHRK2,VHRK2,WHRK2);
    
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
    
    // FSF
    starttime=pgc->timer();
    
    pconvec->start(p,d,4,WLRK2,d->UH);
    pfsf->rk3_step3(p, d, pgc, pflow, d->UH, d->VH, d->WH, WLRK1, WLRK2, 2.0/3.0);
    omega_update(p,d,pgc,d->WL,d->U,d->V,d->W);
    
    p->fsftime+=pgc->timer()-starttime;
    }
}

// stage s: momentum fluxes and RK update
void nhflow_momentum_RK3::phase_M(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S, int s)
{
    ioflow *pflow = S.pflow;
    nhflow_convection *pconvec = S.pconvec;
    nhflow_diffusion *pnhfdiff = S.pdiff;
    nhflow_pressure *ppress = S.ppress;
    nhflow_turbulence *pnhfturb = S.pturb;
    vrans_nhflow *pvrans = S.pvrans;
    solver *psolv = S.psolv;
    
    // stage input (diffusion) and output
    double *UHi = s==0?d->UH:(s==1?UHRK1:UHRK2);
    double *VHi = s==0?d->VH:(s==1?VHRK1:VHRK2);
    double *WHi = s==0?d->WH:(s==1?WHRK1:WHRK2);
    double *UHo = stage_UH(d,s,0);
    double *VHo = stage_UH(d,s,1);
    double *WHo = stage_UH(d,s,2);
    slice &WL = stage_WL(d,s);
    const double alpha = s==0?1.0:(s==1?0.25:2.0/3.0);
    
	// U
	starttime=pgc->timer();

	pnhfturb->isource(p,d);
	pflow->isource_nhflow(p,d,pgc,pvrans,WL); 
	ppress->upgrad(p,d,WL);
    p6dof->isource(p,d,pgc,WL);
    pwind->wind_forcing_nhf_x(p,d,pgc,d->U,d->V, d->F, WL, d->eta);
    roughness_u(p,d,d->U,d->F,WL);
    irhs(p,d,pgc);
    pconvec->start(p,d,1,WL,UHo);
    pnhfdiff->diff_u(p,d,pgc,pflow,psolv,UHDIFF,UHi,UHi,VHi,WHi,WL,alpha);

    if(s==0)
	LOOP
	UHRK1[IJK] = UHDIFF[IJK]
				+ p->dt*CPORNH*d->F[IJK];
    
    if(s==1)
    LOOP
	UHRK2[IJK] = 0.75*d->UH[IJK] + 0.25*UHDIFF[IJK]
				+ 0.25*p->dt*CPORNH*d->F[IJK];
    
    if(s==2)
    LOOP
	d->UH[IJK] = (1.0/3.0)*d->UH[IJK] + (2.0/3.0)*UHDIFF[IJK]
				+ (2.0/3.0)*p->dt*CPORNH*d->F[IJK];

    if(s==0)
    p->utime=pgc->timer()-starttime;
    else
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pnhfturb->jsource(p,d);
	pflow->jsource_nhflow(p,d,pgc,pvrans,WL); 
    ppress->vpgrad(p,d,WL);
    p6dof->jsource(p,d,pgc,WL);
    pwind->wind_forcing_nhf_y(p,d,pgc,d->U,d->V, d->G, WL, d->eta);
    roughness_v(p,d,d->V,d->G,WL);
    jrhs(p,d,pgc);
    pconvec->start(p,d,2,WL,VHo);
    pnhfdiff->diff_v(p,d,pgc,pflow,psolv,VHDIFF,VHi,UHi,VHi,WHi,WL,alpha);

    if(s==0)
	LOOP
	VHRK1[IJK] = VHDIFF[IJK]
				+ p->dt*CPORNH*d->G[IJK];
    
    if(s==1)
    LOOP
	VHRK2[IJK] = 0.75*d->VH[IJK] + 0.25*VHDIFF[IJK]
				+ 0.25*p->dt*CPORNH*d->G[IJK];
    
    if(s==2)
    LOOP
	d->VH[IJK] = (1.0/3.0)*d->VH[IJK] + (2.0/3.0)*VHDIFF[IJK]
				+ (2.0/3.0)*p->dt*CPORNH*d->G[IJK];

    if(s==0)
    p->vtime=pgc->timer()-starttime;
    else
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

    pnhfturb->ksource(p,d);
    pflow->ksource_nhflow(p,d,pgc,pvrans,WL); 
    ppress->wpgrad(p,d,WL);
    krhs(p,d,pgc);
    pconvec->start(p,d,3,WL,WHo);
    pnhfdiff->diff_w(p,d,pgc,pflow,psolv,WHDIFF,WHi,UHi,VHi,WHi,WL,alpha);

    if(s==0)
	LOOP
	WHRK1[IJK] = WHDIFF[IJK]
				+ p->dt*CPORNH*d->H[IJK];
    
    if(s==1)
    LOOP
	WHRK2[IJK] = 0.75*d->WH[IJK] + 0.25*WHDIFF[IJK]
				+ 0.25*p->dt*CPORNH*d->H[IJK];
    
    if(s==2)
    LOOP
	d->WH[IJK] = (1.0/3.0)*d->WH[IJK] + (2.0/3.0)*WHDIFF[IJK]
				+ (2.0/3.0)*p->dt*CPORNH*d->H[IJK];
	
    if(s==0)
    p->wtime=pgc->timer()-starttime;
    else
    p->wtime+=pgc->timer()-starttime;
}

// stage s: velocities, forcing (before the pressure projection)
void nhflow_momentum_RK3::phase_P1(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S, int s)
{
    double *UHo = stage_UH(d,s,0);
    double *VHo = stage_UH(d,s,1);
    double *WHo = stage_UH(d,s,2);
    slice &WL = stage_WL(d,s);
    const double alpha = stage_alpha(s);
    const int fin = s==2?1:0;
    
    velcalc(p,d,pgc,UHo,VHo,WHo,WL,alpha);
    
    pnhfdf->forcing(p, d, pgc, p6dof, s, alpha, UHo, VHo, WHo, WL, fin);
}

// stage s: velocities, reforcing (after the pressure projection)
void nhflow_momentum_RK3::phase_P2(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S, int s)
{
    double *UHo = stage_UH(d,s,0);
    double *VHo = stage_UH(d,s,1);
    double *WHo = stage_UH(d,s,2);
    slice &WL = stage_WL(d,s);
    const double alpha = stage_alpha(s);
    const int fin = s==2?1:0;
    
    velcalc(p,d,pgc,UHo,VHo,WHo,WL,alpha);
    
    pnhfdf->reforcing(p, d, pgc, p6dof, s, alpha, UHo, VHo, WHo, WL, fin);
}

// stage s: relaxation zones, ghost cells
void nhflow_momentum_RK3::phase_E(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S, int s)
{
    ioflow *pflow = S.pflow;
    double *UHo = stage_UH(d,s,0);
    double *VHo = stage_UH(d,s,1);
    double *WHo = stage_UH(d,s,2);
    
    pflow->U_relax(p,pgc,d->U,UHo);
    pflow->V_relax(p,pgc,d->V,VHo);
    pflow->W_relax(p,pgc,d->W,WHo);

	pflow->P_relax(p,pgc,d->P);

	pgc->start4V(p,UHo,gcval_uh);
    pgc->start4V(p,VHo,gcval_vh);
    pgc->start4V(p,WHo,gcval_wh);
    
    clearrhs(p,d,pgc);
}
