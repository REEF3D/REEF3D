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

#include"momentum_FCC3_PLIC.h"
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

momentum_FCC3_PLIC::momentum_FCC3_PLIC(lexer *p, fdm *a, ghostcell *pgc, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, ioflow *pioflow,
                                                    heat *&pheat, concentration *&pconc, reini *ppreini,
                                                    fsi *ppfsi)
                                                    :momentum_forcing(p),bcmom(p),udiff(p),vdiff(p),wdiff(p),ur(p),vr(p),wr(p),urk1(p),urk2(p),vrk1(p),vrk2(p),wrk1(p),wrk2(p),ls(p),frk1(p),frk2(p),
                                                    Mx(p),rox(p),My(p),roy(p),Mz(p),roz(p),
                                                    Mx_rk1(p),Mx_rk2(p),My_rk1(p),My_rk2(p),Mz_rk1(p),Mz_rk2(p),
                                                    rox_rk1(p),rox_rk2(p),roy_rk1(p),roy_rk2(p),roz_rk1(p),roz_rk2(p),
                                                    fx(p),fy(p),fz(p),VoF(p),vof_rk1(p),vof_rk2(p),vof_rk3(p)
                                                    
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
    
   ro_threshold = p->F91*p->W1 + p->W3;
}

momentum_FCC3_PLIC::~momentum_FCC3_PLIC()
{
}

void momentum_FCC3_PLIC::start(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof)
{	
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
	pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);
    
//********************************************************
//Step 1
//********************************************************
    // get vectorized face density from density_f
    pgc->start4(p,a->ro,gcval_ro);
    if(p->F92==3)
        pplic->calculateSubFractions(p,a,pgc,a->vof);
    face_density(p,a,pgc,rox,roy,roz);
    
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
        
    // get M form M = rho * U
    ULOOP
    Mx(i,j,k) = rox(i,j,k)*a->u(i,j,k);
    
    VLOOP
    My(i,j,k) = roy(i,j,k)*a->v(i,j,k);
    
    WLOOP
    Mz(i,j,k) = roz(i,j,k)*a->w(i,j,k);
    
    pgc->start1(p,Mx,gcval_u);
	pgc->start2(p,My,gcval_v);
	pgc->start3(p,Mz,gcval_w);
    
    // advect M    
    pconvec->start(p,a,Mx,1,a->u,a->v,a->w);
    pconvec->start(p,a,My,2,a->u,a->v,a->w);
    pconvec->start(p,a,Mz,3,a->u,a->v,a->w);
    
    ULOOP
	Mx_rk1(i,j,k) = Mx(i,j,k)
				+ p->dt*CPOR1*a->F(i,j,k);
                
    VLOOP
	My_rk1(i,j,k) = My(i,j,k)
				+ p->dt*CPOR2*a->G(i,j,k);
                
    WLOOP
	Mz_rk1(i,j,k) = Mz(i,j,k)
				+ p->dt*CPOR3*a->H(i,j,k);
    
   
    clear_FGH(p,a);
    
    // advect rho
    pgc->start1(p,rox,gcval_u);
    pgc->start2(p,roy,gcval_v);
    pgc->start3(p,roz,gcval_w);
    
    pconvec->start(p,a,rox,1,a->u,a->v,a->w);
    pconvec->start(p,a,roy,2,a->u,a->v,a->w);
    pconvec->start(p,a,roz,3,a->u,a->v,a->w);
    
    ULOOP
	rox_rk1(i,j,k) = rox(i,j,k)
				+ p->dt*CPOR1*a->F(i,j,k);
                
    VLOOP
	roy_rk1(i,j,k) = roy(i,j,k)
				+ p->dt*CPOR2*a->G(i,j,k);
                
    WLOOP
	roz_rk1(i,j,k) = roz(i,j,k)
				+ p->dt*CPOR3*a->H(i,j,k);
    
    clear_FGH(p,a);
    
    // reconstruct U
    ULOOP
    ur(i,j,k) = vel_limiter(p,a,urk1,Mx_rk1,rox_rk1,rox);
    
    VLOOP
    vr(i,j,k) = vel_limiter(p,a,vrk1,My_rk1,roy_rk1,roy);
    
    WLOOP
    wr(i,j,k) = vel_limiter(p,a,wrk1,Mz_rk1,roz_rk1,roz);
    
    pgc->start1(p,ur,gcval_u);
	pgc->start2(p,vr,gcval_v);
	pgc->start3(p,wr,gcval_w);
    
    
    //-------------------------------------------
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans); 
	bcmomPLIC_start(a,p,pgc,pturb,pplic,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,a->u,a->u,a->v,a->w,1.0);
	pdiff->diff_u(p,a,pgc,psolv,udiff,ur,a->u,a->v,a->w,1.0);

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
	pdiff->diff_v(p,a,pgc,psolv,vdiff,vr,a->u,a->v,a->w,1.0);

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
	pdiff->diff_w(p,a,pgc,psolv,wdiff,wr,a->u,a->v,a->w,1.0);

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
    
    momentum_forcing_start(a, p, pgc, p6dof, pvrans, pfsi,
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
    
    //pplic->updatePhasemarkersCompression(p,a,pgc,vof_rk1);
    pplic->updatePhasemarkersCorrection(p,a,pgc,vof_rk1);
    pgc->start4(p,vof_rk1,gcval_vof);
    LOOP
        a->vof(i,j,k)=vof_rk1(i,j,k);
    pgc->start4(p,a->vof,gcval_vof);
    if(p->F92==3)
        pplic->calculateSubFractions(p,a,pgc,a->vof);
    pupdate->start(p,a,pgc);
    pgc->start4(p,a->ro,gcval_ro);
    pgc->start4(p,a->visc,gcval_visc);
	
//********************************************************
//Step 2
//********************************************************
   // face_density(p,a,pgc,rox_rk1,roy_rk1,roz_rk1); 
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
        a->L(i,j,k) = 0.0;
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
                
    clear_FGH(p,a);
        
    ULOOP
    Mx_rk1(i,j,k) = rox_rk1(i,j,k)*urk1(i,j,k);
    
    VLOOP
    My_rk1(i,j,k) = roy_rk1(i,j,k)*vrk1(i,j,k);
    
    WLOOP
    Mz_rk1(i,j,k) = roz_rk1(i,j,k)*wrk1(i,j,k);
    
    pgc->start1(p,Mx_rk1,gcval_u);
	pgc->start2(p,My_rk1,gcval_v);
	pgc->start3(p,Mz_rk1,gcval_w);
      
    pconvec->start(p,a,Mx_rk1,1,urk1,vrk1,wrk1);
    pconvec->start(p,a,My_rk1,2,urk1,vrk1,wrk1);
    pconvec->start(p,a,Mz_rk1,3,urk1,vrk1,wrk1);
    
    ULOOP
	Mx_rk2(i,j,k) = 0.75*Mx(i,j,k) + 0.25*Mx_rk1(i,j,k)
				+ 0.25*p->dt*CPOR1*a->F(i,j,k);
                
    VLOOP
	My_rk2(i,j,k) = 0.75*My(i,j,k) + 0.25*My_rk1(i,j,k)
				+ 0.25*p->dt*CPOR2*a->G(i,j,k);
                
    WLOOP
	Mz_rk2(i,j,k) = 0.75*Mz(i,j,k) + 0.25*Mz_rk1(i,j,k)
				+ 0.25*p->dt*CPOR3*a->H(i,j,k);
    
    clear_FGH(p,a);

    pgc->start1(p,rox_rk1,gcval_u);
	pgc->start2(p,roy_rk1,gcval_v);
	pgc->start3(p,roz_rk1,gcval_w);
    
    pconvec->start(p,a,rox_rk1,1,urk1,vrk1,wrk1);
    pconvec->start(p,a,roy_rk1,2,urk1,vrk1,wrk1);
    pconvec->start(p,a,roz_rk1,3,urk1,vrk1,wrk1);
    
    ULOOP
	rox_rk2(i,j,k) = 0.75*rox(i,j,k) + 0.25*rox_rk1(i,j,k)
				+ 0.25*p->dt*CPOR1*a->F(i,j,k);
                
    VLOOP
	roy_rk2(i,j,k) = 0.75*roy(i,j,k) + 0.25*roy_rk1(i,j,k)
				+ 0.25*p->dt*CPOR1*a->G(i,j,k);
                
    WLOOP
	roz_rk2(i,j,k) = 0.75*roz(i,j,k) + 0.25*roz_rk1(i,j,k)
				+ 0.25*p->dt*CPOR1*a->H(i,j,k);
    
    clear_FGH(p,a);
    
    ULOOP
    ur(i,j,k) = vel_limiter(p,a,urk2,Mx_rk2,rox_rk2,rox_rk1);
    
    VLOOP
    vr(i,j,k) = vel_limiter(p,a,vrk2,My_rk2,roy_rk2,roy_rk1);
    
    WLOOP
    wr(i,j,k) = vel_limiter(p,a,wrk2,Mz_rk2,roz_rk2,roz_rk1);
    
    pgc->start1(p,ur,gcval_u);
	pgc->start2(p,vr,gcval_v);
	pgc->start3(p,wr,gcval_w);
    
    //-------------------------------------------
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmomPLIC_start(a,p,pgc,pturb,pplic,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,urk1,urk1,vrk1,wrk1,0.25);
	pdiff->diff_u(p,a,pgc,psolv,udiff,ur,urk1,vrk1,wrk1,1.0);

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
	pdiff->diff_v(p,a,pgc,psolv,vdiff,vr,urk1,vrk1,wrk1,1.0);

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
	pdiff->diff_w(p,a,pgc,psolv,wdiff,wr,urk1,vrk1,wrk1,1.0);

	WLOOP
	wrk2(i,j,k) = wdiff(i,j,k)
				+ 0.25*p->dt*CPOR3*a->H(i,j,k);

    p->wtime+=pgc->timer()-starttime;
    
    pgc->start1(p,urk2,gcval_u);
	pgc->start2(p,vrk2,gcval_v);
    pgc->start3(p,wrk2,gcval_w);
    clear_FGH(p,a);
    
    if(p->F92==3)
    {
        pplic->calculateSubFractions(p,a,pgc,a->vof);
    }
    pupdate->start(p,a,pgc);
    pgc->start4(p,a->ro,gcval_ro);
    pgc->start4(p,a->visc,gcval_visc);
    
    momentum_forcing_start(a, p, pgc, p6dof, pvrans, pfsi,
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
    
   // pplic->updatePhasemarkersCompression(p,a,pgc,vof_rk2);
    pplic->updatePhasemarkersCorrection(p,a,pgc,vof_rk2);
    pgc->start4(p,vof_rk2,gcval_vof);
    LOOP
        a->vof(i,j,k)=vof_rk2(i,j,k);
    pgc->start4(p,a->vof,gcval_vof);
    if(p->F92==3)
    {
        pplic->calculateSubFractions(p,a,pgc,a->vof);
    }
    pupdate->start(p,a,pgc);
    pgc->start4(p,a->ro,gcval_ro);
    pgc->start4(p,a->visc,gcval_visc);
    

//********************************************************
//Step 3
//********************************************************
    
   // face_density(p,a,pgc,rox_rk2,roy_rk2,roz_rk2); 
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
    
   // pplic->updatePlaneData(p,a,pgc,a->vof);
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
                
    clear_FGH(p,a);

    ULOOP
    Mx_rk2(i,j,k) = rox_rk2(i,j,k)*urk2(i,j,k);
    
    VLOOP
    My_rk2(i,j,k) = roy_rk2(i,j,k)*vrk2(i,j,k);
    
    WLOOP
    Mz_rk2(i,j,k) = roz_rk2(i,j,k)*wrk2(i,j,k);
    
    pgc->start1(p,Mx_rk2,gcval_u);
	pgc->start2(p,My_rk2,gcval_v);
	pgc->start3(p,Mz_rk2,gcval_w);
    
    pconvec->start(p,a,Mx_rk2,1,urk2,vrk2,wrk2);
    pconvec->start(p,a,My_rk2,2,urk2,vrk2,wrk2);
    pconvec->start(p,a,Mz_rk2,3,urk1,vrk1,wrk1);
    
    ULOOP
    Mx(i,j,k) = (1.0/3.0)*Mx(i,j,k) + (2.0/3.0)*Mx_rk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR1*a->F(i,j,k);
    
    VLOOP
    My(i,j,k) = (1.0/3.0)*My(i,j,k) + (2.0/3.0)*My_rk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR3*a->G(i,j,k);
    
    WLOOP
    Mz(i,j,k) = (1.0/3.0)*Mz(i,j,k) + (2.0/3.0)*Mz_rk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR3*a->H(i,j,k);
    
    clear_FGH(p,a);
    
    pgc->start1(p,rox_rk2,gcval_u);
	pgc->start2(p,roy_rk2,gcval_v);
	pgc->start3(p,roz_rk2,gcval_w);
    
    pconvec->start(p,a,rox_rk2,1,urk2,vrk2,wrk2);
    pconvec->start(p,a,roy_rk2,2,urk2,vrk2,wrk2);
    pconvec->start(p,a,roz_rk2,3,urk2,vrk2,wrk2);
    
    ULOOP
	rox(i,j,k) = (1.0/3.0)*rox(i,j,k) + (2.0/3.0)*rox_rk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR1*a->F(i,j,k);
                
    VLOOP
	roy(i,j,k) = (1.0/3.0)*roy(i,j,k) + (2.0/3.0)*roy_rk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR2*a->G(i,j,k);
                
    WLOOP
	roz(i,j,k) = (1.0/3.0)*roz(i,j,k) + (2.0/3.0)*roz_rk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR3*a->H(i,j,k);

    clear_FGH(p,a);
    
    ULOOP
    ur(i,j,k) = vel_limiter(p,a,a->u,Mx,rox,rox_rk2);
    
    VLOOP
    vr(i,j,k) = vel_limiter(p,a,a->v,My,roy,roy_rk2);
    
    WLOOP
    wr(i,j,k) = vel_limiter(p,a,a->w,Mz,roz,roz_rk2);
    
    pgc->start1(p,ur,gcval_u);
	pgc->start2(p,vr,gcval_v);
	pgc->start3(p,wr,gcval_w);
    
    //-------------------------------------------
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmomPLIC_start(a,p,pgc,pturb,pplic,a->u,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
	irhs(p,a,pgc,urk2,urk2,vrk2,wrk2,2.0/3.0);
	pdiff->diff_u(p,a,pgc,psolv,udiff,ur,urk2,vrk2,wrk2,1.0);

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
	pdiff->diff_v(p,a,pgc,psolv,vdiff,vr,urk2,vrk2,wrk2,1.0);

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
	pdiff->diff_w(p,a,pgc,psolv,wdiff,wr,urk2,vrk2,wrk2,1.0);

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
    
    momentum_forcing_start(a, p, pgc, p6dof, pvrans, pfsi,
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
    
   // pplic->updatePhasemarkersCompression(p,a,pgc,a->vof);
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

void momentum_FCC3_PLIC::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
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

void momentum_FCC3_PLIC::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
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

void momentum_FCC3_PLIC::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
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


void momentum_FCC3_PLIC::utimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FCC3_PLIC::vtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FCC3_PLIC::wtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FCC3_PLIC::clear_FGH(lexer *p, fdm *a)
{
    ULOOP
    a->F(i,j,k) = 0.0;
    
    VLOOP
    a->G(i,j,k) = 0.0;
    
    WLOOP
    a->H(i,j,k) = 0.0;
}

void momentum_FCC3_PLIC::face_density(lexer *p, fdm *a, ghostcell *pgc, field &rox, field &roy, field &roz)
{
    ULOOP
    rox(i,j,k) = pd->roface(p,a,1,0,0);
    
    VLOOP
    roy(i,j,k) = pd->roface(p,a,0,1,0);
    
    WLOOP
    roz(i,j,k) = pd->roface(p,a,0,0,1);
    
    pgc->start1(p,rox,50);
	pgc->start2(p,roy,50);
	pgc->start3(p,roz,50);
}

double momentum_FCC3_PLIC::vel_limiter(lexer *p, fdm *a, field &vel, field &M, field &ro,field &ro_n)
{
    if(ro(i,j,k)>=ro_threshold)
    val = M(i,j,k)/ro(i,j,k);
    
    else
    if(ro(i,j,k)>p->W3 && ro(i,j,k)<ro_threshold && ro(i,j,k)<ro_n(i,j,k))
    val = (M(i,j,k)/ro_filter(p,a,ro))*(ro_filter(p,a,ro)/ro_threshold) + vel(i,j,k)*(ro_threshold-ro_filter(p,a,ro))/ro_threshold;
    
    else
    val = vel(i,j,k);  

    return val;
}

double momentum_FCC3_PLIC::ro_filter(lexer *p, fdm *a, field &ro)
{
    if(ro(i,j,k)<p->W3)
    val = p->W3;
    
    else
    if(ro(i,j,k)>=p->W3)
    val = ro(i,j,k);

    return val;
}
