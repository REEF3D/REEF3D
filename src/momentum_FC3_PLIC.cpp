/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MEFCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Fabian Knoblauch
--------------------------------------------------------------------*/

#include"momentum_FC3_PLIC.h"
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
#include"reini.h"
#include"picard.h"
#include"fluid_update_vof.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"nhflow.h"
#include"heat.h"
#include"concentration.h"
#include"density_vof.h"
#include"VOF_PLIC.h"

momentum_FC3_PLIC::momentum_FC3_PLIC(lexer *p, fdm *a, ghostcell *pgc, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, ioflow *pioflow,
                                                    heat *&pheat, concentration *&pconc, reini *ppreini,
                                                    fsi *ppfsi)
                                                    :momentum_forcing(p),bcmom(p),udiff(p),vdiff(p),wdiff(p),urk1(p),urk2(p),vrk1(p),vrk2(p),wrk1(p),wrk2(p),
                                                    VoF(p),vof_rk1(p),vof_rk2(p),fx(p),fy(p),fz(p)
                                                    
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    gcval_ro=1;
    gcval_visc=1;
    gcval_vof=1;
    
    if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;

	pconvec=pconvection;
	pdiff=pdiffusion;
	ppress=ppressure;
	ppois=ppoisson;
	pturb=pturbulence;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
    preini=ppreini;
    pfsi=ppfsi;
    pplic= new VOF_PLIC(p,a,pgc,pheat);
    pupdate = new fluid_update_vof(p,a,pgc);
	pd = new density_vof(p);
    
	if(p->F46==2)
	ppicard = new picard_f(p);

	if(p->F46==3)
	ppicard = new picard_lsm(p);

	if(p->F46!=2 && p->F46!=3)
	ppicard = new picard_void(p);
    
}

momentum_FC3_PLIC::~momentum_FC3_PLIC()
{
}

void momentum_FC3_PLIC::start(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof, vector<net*>& pnet)
{	
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
	pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);
    
//********************************************************
//Step 1
//********************************************************
    // get vectorized face density from density_f
    
    //-------------------------------------------
    // FSF
    LOOP
    {
	a->L(i,j,k)=0.0;
    VoF(i,j,k)=a->vof(i,j,k);
    }

	pplic->RKcalcL(a,p,pgc,a->u,a->v,a->w);
	
	LOOP
    {
        vof_rk1(i,j,k) = VoF(i,j,k) + a->L(i,j,k);
        
        if(vof_rk1(i,j,k)<0.0)
            vof_rk1(i,j,k)=0.0;
        if(vof_rk1(i,j,k)>1.0)
            vof_rk1(i,j,k)=1.0;
    }
    
     pgc->start4(p,vof_rk1,gcval_vof);
   // pplic->updatePlaneData(p,a,pgc,vof_rk1);
	pflow->vof_relax(p,a,pgc,vof_rk1);
	pgc->start4(p,vof_rk1,gcval_vof);
    
    LOOP
    {
     a->vof(i,j,k) = vof_rk1(i,j,k);
     a->L(i,j,k)=0.0;
    }
    pgc->start4(p,a->vof,gcval_vof); 
    
    //!no update yet -> update after diffusion!
    
    if(p->F92==1)
    {
        pplic->RK_redistance(a,p,pgc);
        pgc->start4(p,a->phi,gcval_phi);
        
        p->F44=3;
        preini->start(a,p,a->phi, pgc, pflow);
        ppicard->correct_ls(p,a,pgc,a->phi);
    }
    
    //-------------------------------------------
    
    pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
    
    // advect U
    pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
    pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
    pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
    
    ULOOP
	urk1(i,j,k) = a->u(i,j,k)
				+ p->dt*CPOR1*a->F(i,j,k);
                
    VLOOP
	vrk1(i,j,k) = a->v(i,j,k)
				+ p->dt*CPOR2*a->G(i,j,k);
                
    WLOOP
	wrk1(i,j,k) = a->w(i,j,k)
				+ p->dt*CPOR3*a->H(i,j,k);
                
    // clear_FGH
    clear_FGH(p,a);
    
    pgc->start1(p,urk1,gcval_u);
	pgc->start2(p,vrk1,gcval_v);
	pgc->start3(p,wrk1,gcval_w);
    
    //-------------------------------------------
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans); 
	bcmomPLIC_start(a,p,pgc,pturb,pplic,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,a->u,a->u,a->v,a->w,1.0);
	pdiff->diff_u(p,a,pgc,psolv,udiff,urk1,a->u,a->v,a->w,1.0);

	ULOOP
	urk1(i,j,k) = udiff(i,j,k)
				+ p->dt*CPOR1*a->F(i,j,k);

    p->utime=pgc->timer()-starttime;
    
    //-------------------------------------------
	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmomPLIC_start(a,p,pgc,pturb,pplic,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,a->v,a->u,a->v,a->w,1.0);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,vrk1,a->u,a->v,a->w,1.0);

	VLOOP
	vrk1(i,j,k) = vdiff(i,j,k)
				+ p->dt*CPOR2*a->G(i,j,k);

    p->vtime=pgc->timer()-starttime;
    
    //-------------------------------------------
	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmomPLIC_start(a,p,pgc,pturb,pplic,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,a->w,a->u,a->v,a->w,1.0);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,wrk1,a->u,a->v,a->w,1.0);

	WLOOP
	wrk1(i,j,k) = wdiff(i,j,k)
				+ p->dt*CPOR3*a->H(i,j,k);
	
    p->wtime=pgc->timer()-starttime;
    
    pgc->start1(p,urk1,gcval_u);
	pgc->start2(p,vrk1,gcval_v);
    pgc->start3(p,wrk1,gcval_w);
    clear_FGH(p,a);
    
    //update rho(vof) after diffusion but before pressure
    if(p->F92==3)
        pplic->calculateSubFractions(p,a,pgc,a->vof);
    pupdate->start(p,a,pgc);
    pgc->start4(p,a->ro,gcval_ro);
    pgc->start4(p,a->visc,gcval_visc);
    
    momentum_forcing_start(a, p, pgc, p6dof, pvrans, pnet, pfsi,
                           urk1, vrk1, wrk1, fx, fy, fz, 0, 1.0, false);
    
    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk1, vrk1, wrk1, 1.0);
	
	pflow->u_relax(p,a,pgc,urk1);
	pflow->v_relax(p,a,pgc,vrk1);
	pflow->w_relax(p,a,pgc,wrk1);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,urk1,gcval_u);
	pgc->start2(p,vrk1,gcval_v);
	pgc->start3(p,wrk1,gcval_w);
    
    clear_FGH(p,a);
    
    pplic->updatePhasemarkersCorrection(p,a,pgc,vof_rk1);
    LOOP
        a->vof(i,j,k)=vof_rk1(i,j,k);
    pgc->start4(p,a->vof,gcval_vof);
    pgc->start4(p,vof_rk1,gcval_vof);
    
    if(p->F92==3)
        pplic->calculateSubFractions(p,a,pgc,a->vof);
    pupdate->start(p,a,pgc);
    pgc->start4(p,a->ro,gcval_ro);
    pgc->start4(p,a->visc,gcval_visc);

	
//********************************************************
//Step 2
//********************************************************

    //-------------------------------------------
    // FSF
    
    pplic->RKcalcL(a,p,pgc,urk1,vrk1,wrk1);
	
	LOOP
    {
        vof_rk2(i,j,k) = 0.75*VoF(i,j,k) + 0.25*vof_rk1(i,j,k)+0.25*a->L(i,j,k);
        
        if(vof_rk2(i,j,k)<0.0)
            vof_rk2(i,j,k)=0.0;
        if(vof_rk2(i,j,k)>1.0)
            vof_rk2(i,j,k)=1.0;
    }
    
   // pplic->updatePlaneData(p,a,pgc,vof_rk2);
	pflow->vof_relax(p,a,pgc,vof_rk2);
    pgc->start4(p,vof_rk2,gcval_vof);
    
    LOOP
    {
        a->vof(i,j,k) = vof_rk2(i,j,k);
        a->L(i,j,k)=0.0;
    }
    pgc->start4(p,a->vof,gcval_vof);
    
     if(p->F92==1)
    {
        pplic->RK_redistance(a,p,pgc);
        pgc->start4(p,a->phi,gcval_phi);
        
        p->F44=3;
        preini->start(a,p,a->phi, pgc, pflow);
        ppicard->correct_ls(p,a,pgc,a->phi);
    }
    
    // advect U
    pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
    pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
    pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
    
    ULOOP
	urk2(i,j,k) = 0.75*a->u(i,j,k) + 0.25*urk1(i,j,k)
				+ 0.25*p->dt*CPOR1*a->F(i,j,k);
                
    VLOOP
	vrk2(i,j,k) = 0.75*a->v(i,j,k) + 0.25*vrk1(i,j,k)
				+ 0.25*p->dt*CPOR2*a->G(i,j,k);
                
    WLOOP
	wrk2(i,j,k) = 0.75*a->w(i,j,k) + 0.25*wrk1(i,j,k)
				+ 0.25*p->dt*CPOR3*a->H(i,j,k);
    
    pgc->start1(p,urk2,gcval_u);
	pgc->start2(p,vrk2,gcval_v);
	pgc->start3(p,wrk2,gcval_w);
    clear_FGH(p,a);
    
    //-------------------------------------------
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmomPLIC_start(a,p,pgc,pturb,pplic,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,urk1,urk1,vrk1,wrk1,0.25);
	pdiff->diff_u(p,a,pgc,psolv,udiff,urk2,urk1,vrk1,wrk1,1.0);

	ULOOP
	urk2(i,j,k) = udiff(i,j,k)
				+ 0.25*p->dt*CPOR1*a->F(i,j,k);
                
    p->utime+=pgc->timer()-starttime;
	
    //-------------------------------------------
	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmomPLIC_start(a,p,pgc,pturb,pplic,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,vrk1,urk1,vrk1,wrk1,0.25);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,vrk2,urk1,vrk1,wrk1,1.0);

	VLOOP
	vrk2(i,j,k) = vdiff(i,j,k)
				+ 0.25*p->dt*CPOR2*a->G(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;
    
    //-------------------------------------------
	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmomPLIC_start(a,p,pgc,pturb,pplic,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,wrk1,urk1,vrk1,wrk1,0.25);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,wrk2,urk1,vrk1,wrk1,1.0);

	WLOOP
	wrk2(i,j,k) = wdiff(i,j,k)
				+ 0.25*p->dt*CPOR3*a->H(i,j,k);

    p->wtime+=pgc->timer()-starttime;
    
    pgc->start1(p,urk2,gcval_u);
	pgc->start2(p,vrk2,gcval_v);
    pgc->start3(p,wrk2,gcval_w);
    clear_FGH(p,a);
    
    if(p->F92==3)
        pplic->calculateSubFractions(p,a,pgc,a->vof);
    pupdate->start(p,a,pgc);
    pgc->start4(p,a->ro,gcval_ro);
    pgc->start4(p,a->visc,gcval_visc);
    
    momentum_forcing_start(a, p, pgc, p6dof, pvrans, pnet, pfsi,
                           urk2, vrk2, wrk2, fx, fy, fz, 1, 0.25, false);

    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk2, vrk2, wrk2, 0.25);
	
	pflow->u_relax(p,a,pgc,urk2);
	pflow->v_relax(p,a,pgc,vrk2);
	pflow->w_relax(p,a,pgc,wrk2);
	pflow->p_relax(p,a,pgc,a->press);
	
	pgc->start1(p,urk2,gcval_u);
	pgc->start2(p,vrk2,gcval_v);
	pgc->start3(p,wrk2,gcval_w);
    clear_FGH(p,a);
    
    pplic->updatePhasemarkersCorrection(p,a,pgc,vof_rk2);
    LOOP
        a->vof(i,j,k)=vof_rk2(i,j,k);
    pgc->start4(p,a->vof,gcval_vof);
    pgc->start4(p,vof_rk1,gcval_vof);
    
    if(p->F92==3)
        pplic->calculateSubFractions(p,a,pgc,a->vof);
    pupdate->start(p,a,pgc);
    pgc->start4(p,a->ro,gcval_ro);
    pgc->start4(p,a->visc,gcval_visc);
    

//********************************************************
//Step 3
//********************************************************
    //-------------------------------------------
    // FSF
    
    pplic->RKcalcL(a,p,pgc,urk2,vrk2,wrk2);
	
	LOOP
    {
        a->vof(i,j,k) = (1.0/3.0)*VoF(i,j,k)
                            + (2.0/3.0)*vof_rk2(i,j,k)
                            + (2.0/3.0)*a->L(i,j,k);
        
        if(a->vof(i,j,k)<0.0)
            a->vof(i,j,k)=0.0;
        if(a->vof(i,j,k)>1.0)
            a->vof(i,j,k)=1.0;
    }
    
  //  pplic->updatePlaneData(p,a,pgc,a->vof);
    pflow->vof_relax(p,a,pgc,a->vof);
    pgc->start4(p,a->vof,gcval_vof);
    
    LOOP
        a->L(i,j,k)=0.0;
    
    if(p->F92==1)
    {
        pplic->RK_redistance(a,p,pgc);
        pgc->start4(p,a->phi,gcval_phi);
        
        p->F44=4;
        preini->start(a,p,a->phi, pgc, pflow);
        ppicard->correct_ls(p,a,pgc,a->phi);
    }
    
    //-------------------------------------------
    pconvec->start(p,a,urk2,1,urk2,vrk2,wrk2);
    pconvec->start(p,a,vrk2,2,urk2,vrk2,wrk2);
    pconvec->start(p,a,wrk2,3,urk2,vrk2,wrk2);
    
    ULOOP
	a->u(i,j,k) = (1.0/3.0)*a->u(i,j,k) + (2.0/3.0)*urk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR1*a->F(i,j,k);
                
    VLOOP
	a->v(i,j,k) = (1.0/3.0)*a->v(i,j,k) + (2.0/3.0)*vrk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR2*a->G(i,j,k);
                
    WLOOP
	a->w(i,j,k) = (1.0/3.0)*a->w(i,j,k) + (2.0/3.0)*wrk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR3*a->H(i,j,k);
    
    pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
    clear_FGH(p,a);
    
    //-------------------------------------------
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmomPLIC_start(a,p,pgc,pturb,pplic,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,urk2,urk2,vrk2,wrk2,2.0/3.0);
	pdiff->diff_u(p,a,pgc,psolv,udiff,a->u,urk2,vrk2,wrk2,1.0);

	ULOOP
	a->u(i,j,k) = udiff(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR1*a->F(i,j,k);
	
    p->utime+=pgc->timer()-starttime;
    
    //-------------------------------------------
	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmomPLIC_start(a,p,pgc,pturb,pplic,a->v,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
	jrhs(p,a,pgc,vrk2,urk2,vrk2,wrk2,2.0/3.0);
	pdiff->diff_v(p,a,pgc,psolv,vdiff,a->v,urk2,vrk2,wrk2,1.0);

	VLOOP
	a->v(i,j,k) = vdiff(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR2*a->G(i,j,k);
	
    p->vtime+=pgc->timer()-starttime;
    
    //-------------------------------------------
	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmomPLIC_start(a,p,pgc,pturb,pplic,a->w,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
	krhs(p,a,pgc,wrk2,urk2,vrk2,wrk2,2.0/3.0);
	pdiff->diff_w(p,a,pgc,psolv,wdiff,a->w,urk2,vrk2,wrk2,1.0);

	WLOOP
	a->w(i,j,k) = wdiff(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR3*a->H(i,j,k);
	
    p->wtime+=pgc->timer()-starttime;
    
    pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
    clear_FGH(p,a);
    
    if(p->F92==3)
        pplic->calculateSubFractions(p,a,pgc,a->vof);
    pupdate->start(p,a,pgc);
    pgc->start4(p,a->ro,gcval_ro);
    pgc->start4(p,a->visc,gcval_visc);
    
    momentum_forcing_start(a, p, pgc, p6dof, pvrans, pnet, pfsi,
                           a->u, a->v, a->w, fx, fy, fz, 2, 2.0/3.0, true);

	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, a->u, a->v,a->w,2.0/3.0);
	
	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
    clear_FGH(p,a);
    
    pplic->updatePhasemarkersCorrection(p,a,pgc,a->vof);
    pgc->start4(p,a->vof,gcval_vof);
    
    pupdate->start(p,a,pgc);
    pgc->start4(p,a->ro,gcval_ro);
    pgc->start4(p,a->visc,gcval_visc);
    
    LOOP
    {
        if(a->vof(i,j,k)>p->F94)
            a->phi(i,j,k)=1.0;
        else if(a->vof(i,j,k)<p->F93)
            a->phi(i,j,k)=-1.0;
        else
            a->phi(i,j,k)=(a->vof(i,j,k)-0.5)*p->DZN[KP];
    }
    pgc->start4(p,a->phi,1);
    
    
    double vofchecksum;
    vofchecksum=0.0;
    LOOP
        vofchecksum+=a->vof(i,j,k)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
    vofchecksum=pgc->globalsum(vofchecksum);
    cout<<"Total water volume:"<<vofchecksum<<endl;
        
    
}

void momentum_FC3_PLIC::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	ULOOP
	{
    a->maxF=MAX(fabs(a->rhsvec.V[n] + a->gi),a->maxF);
	a->F(i,j,k) += (a->rhsvec.V[n] + a->gi + p->W29_x)*PORVAL1;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void momentum_FC3_PLIC::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	VLOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n] + a->gj),a->maxG);
	a->G(i,j,k) += (a->rhsvec.V[n] + a->gj + p->W29_y)*PORVAL2;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void momentum_FC3_PLIC::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	WLOOP
	{
    a->maxH=MAX(fabs(a->rhsvec.V[n] + a->gk),a->maxH);
	a->H(i,j,k) += (a->rhsvec.V[n] + a->gk + p->W29_z)*PORVAL3;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}


void momentum_FC3_PLIC::utimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FC3_PLIC::vtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FC3_PLIC::wtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FC3_PLIC::clear_FGH(lexer *p, fdm *a)
{
    ULOOP
    a->F(i,j,k) = 0.0;
    
    VLOOP
    a->G(i,j,k) = 0.0;
    
    WLOOP
    a->H(i,j,k) = 0.0;
}

