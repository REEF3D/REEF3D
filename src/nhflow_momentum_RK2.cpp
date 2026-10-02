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

#include"nhflow_momentum_RK2.h"
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
#include"sediment_f.h"

#define WLVL (fabs(WL(i,j))>(1.0*p->A544)?WL(i,j):1.0e20)

nhflow_momentum_RK2::nhflow_momentum_RK2(lexer *p, fdm_nhf *d, ghostcell *pgc, sixdof *pp6dof,vrans_nhflow* ppvrans, 
                                                      nhflow_forcing *ppnhfdf, sediment *ppsed)
                                                    : nhflow_momentum_func(p,d,pgc), nhflow_breaking(p,d,pgc), WLRK1(p)
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
    
    p->Darray(UHDIFF,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VHDIFF,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WHDIFF,p->imax*p->jmax*(p->kmax+2));
    
    sigma_ini(p,d,pgc,d->eta);
    
    p6dof = pp6dof;
    pnhfdf = ppnhfdf;
    pvrans = ppvrans;
    psed = ppsed;
    
    // wind forcing
    if(p->A570==0)
    pwind = new wind_v(p);
    
    if(p->A570>0)
    pwind = new wind_f(p);
}

nhflow_momentum_RK2::~nhflow_momentum_RK2()
{
}

void nhflow_momentum_RK2::start(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, nhflow_signal_speed *pss, 
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
    
    for(int s=0; s<2; ++s)
    {
    phase_F(p,d,pgc,S,s);
    phase_M(p,d,pgc,S,s);
    phase_P(p,d,pgc,S,s);
    phase_E(p,d,pgc,S,s);
    }
}

slice& nhflow_momentum_RK2::stage_WL(fdm_nhf *d, int s)
{
    if(s==0)
    return WLRK1;
    return d->WL;
}

double* nhflow_momentum_RK2::stage_UH(fdm_nhf *d, int s, int m)
{
    if(s==0)
    return m==0?UHRK1:(m==1?VHRK1:WHRK1);
    return m==0?d->UH:(m==1?d->VH:d->WH);
}

void nhflow_momentum_RK2::step_begin(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S)
{
    S.pflow->discharge_nhflow(p,d,pgc);
    S.pflow->inflow_nhflow(p,d,pgc,d->U,d->V,d->W,d->UH,d->VH,d->WH,d->WL);
    S.pflow->rkinflow_nhflow(p,d,pgc,d->U,d->V,d->W,UHRK1,VHRK1,WHRK1,WLRK1);
}

// stage s: sigma, reconstruction, continuity flux, water level, omega, breaking
void nhflow_momentum_RK2::phase_F(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S, int s)
{
    ioflow *pflow = S.pflow;
    nhflow_fsf *pfsf = S.pfsf;
    nhflow_convection *pconvec = S.pconvec;
    
    if(s==0)
    {
    sigma_update(p,d,pgc,d->WL);
    S.pvrans->update(p,d,pgc,0.5,0);
    reconstruct(p,d,pgc,pfsf,S.pss,S.precon,d->WL,d->U,d->V,d->W,d->UH,d->VH,d->WH);
    
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);   
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);

    // FSF
    starttime=pgc->timer();
    pconvec->start(p,d,4,d->WL,UHRK1);

    pfsf->rk2_step1(p, d, pgc, pflow, d->UH, d->VH, d->WH, WLRK1, WLRK1, 1.0);
    omega_update(p,d,pgc,WLRK1,d->U,d->V,d->W);
    breaking(p,d,pgc,d->eta,d->eta_n,WLRK1,1.0);
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
    
    pconvec->start(p,d,4,WLRK1,d->UH);
    pfsf->rk2_step2(p, d, pgc, pflow, UHRK1,VHRK1,WHRK1, WLRK1, WLRK1, 0.5);
    omega_update(p,d,pgc,d->WL,d->U,d->V,d->W);
    breaking(p,d,pgc,d->eta,d->eta_n,d->WL,0.5);
    p->fsftime+=pgc->timer()-starttime;
    }
}

// stage s: momentum fluxes and RK update
void nhflow_momentum_RK2::phase_M(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S, int s)
{
    ioflow *pflow = S.pflow;
    nhflow_convection *pconvec = S.pconvec;
    nhflow_diffusion *pnhfdiff = S.pdiff;
    nhflow_pressure *ppress = S.ppress;
    nhflow_turbulence *pnhfturb = S.pturb;
    vrans_nhflow *pvrans = S.pvrans;
    solver *psolv = S.psolv;
    
    // stage input (diffusion) and output
    double *UHi = s==0?d->UH:UHRK1;
    double *VHi = s==0?d->VH:VHRK1;
    double *WHi = s==0?d->WH:WHRK1;
    double *UHo = stage_UH(d,s,0);
    double *VHo = stage_UH(d,s,1);
    double *WHo = stage_UH(d,s,2);
    slice &WL = stage_WL(d,s);
    const double alpha = s==0?1.0:0.5;
    
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
	d->UH[IJK] = 0.5*d->UH[IJK] + 0.5*UHDIFF[IJK]
				+ 0.5*p->dt*CPORNH*d->F[IJK];

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
	d->VH[IJK] = 0.5*d->VH[IJK] + 0.5*VHDIFF[IJK]
                + 0.5*p->dt*CPORNH*d->G[IJK];

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
    
    if(p->A520!=3 && s==0)
	LOOP
	WHRK1[IJK] = WHDIFF[IJK]
				+ p->dt*CPORNH*d->H[IJK];
    
    if(p->A520!=3 && s==1)
	LOOP
	d->WH[IJK] = 0.5*d->WH[IJK] + 0.5*WHDIFF[IJK]
				+ 0.5*p->dt*CPORNH*d->H[IJK];
	
    if(s==0)
    p->wtime=pgc->timer()-starttime;
    else
    p->wtime+=pgc->timer()-starttime;
}

// stage s: velocities, forcing (before the pressure projection)
void nhflow_momentum_RK2::phase_P1(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S, int s)
{
    double *UHo = stage_UH(d,s,0);
    double *VHo = stage_UH(d,s,1);
    double *WHo = stage_UH(d,s,2);
    slice &WL = stage_WL(d,s);
    const double alpha = stage_alpha(s);
    const int fin = s;
    
    velcalc(p,d,pgc,UHo,VHo,WHo,WL,alpha);
    
    pnhfdf->forcing(p, d, pgc, p6dof, s, alpha, UHo, VHo, WHo, WL, fin);
}

// stage s: velocities, reforcing (after the pressure projection)
void nhflow_momentum_RK2::phase_P2(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S, int s)
{
    double *UHo = stage_UH(d,s,0);
    double *VHo = stage_UH(d,s,1);
    double *WHo = stage_UH(d,s,2);
    slice &WL = stage_WL(d,s);
    const double alpha = stage_alpha(s);
    const int fin = s;
    
    velcalc(p,d,pgc,UHo,VHo,WHo,WL,alpha);
    
    pnhfdf->reforcing(p, d, pgc, p6dof, s, alpha, UHo, VHo, WHo, WL, fin);
}

// stage s: relaxation zones, ghost cells, sediment and depth
void nhflow_momentum_RK2::phase_E(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S, int s)
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
    
    if(s==0)
    {
    psed->RK2_step1_nhflow(p,d,pgc,pflow);
    S.pfsf->depth_update(p,d,pgc,pflow);
    }
    
    if(s==1)
    {
    psed->RK2_step2_nhflow(p,d,pgc,pflow);
    S.pfsf->depth_update(p,d,pgc,pflow);
    bed_acceleration(p,d,pgc,d->WL,d->U,d->V,d->W);
    }
}
