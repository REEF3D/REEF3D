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
Author: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"momentum_FCLS3.h"
#include"momentum_FCLS3.h"
#include"vrans.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"convection.h"
#include"diffusion.h"
#include"density.h"
#include"ediff2.h"
#include"reini.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"pressure.h"
#include"poisson.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"6DOF.h"
#include"FSI.h"
#include"picard.h"
#include"fluid_update_fsf.h"
#include"fluid_update_fsf_heat.h"
#include"fluid_update_fsf_heat_Bouss.h"
#include"fluid_update_fsf_comp.h"
#include"fluid_update_void.h"
#include"fluid_update_fsf_concentration.h"
#include"fluid_update_rheology.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"heat.h"

momentum_FCLS3::momentum_FCLS3(lexer *p, fdm *a, ghostcell *pgc, convection *pconvection, convection *ppfsfdisc, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, ioflow *pioflow,
                                                    heat *&pheat, concentration *&pconc, reini *ppreini,
                                                    fsi *ppfsi) :
                                                    momentum_forcing(p),bcmom(p),urk(p),vrk(p),wrk(p),
                                                    Cu(p),Cv(p),Cw(p),Cf(p),fx(p),fy(p),fz(p)
{
    
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;
	
	pconvec=pconvection;
    pfsfdisc=ppfsfdisc;
	pdiff=pdiffusion;
	ppress=ppressure;
	ppois=ppoisson;
	pturb=pturbulence;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
    preini=ppreini;
    pfsi=ppfsi;
    
    if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90==0)
	pupdate = new fluid_update_fsf(p,a,pgc);
	
	if(p->F30>0 && p->H10==0 && p->W30==1 && p->F300==0 && p->W90==0)
	pupdate = new fluid_update_fsf_comp(p,a,pgc);
	
	if(p->F30>0 && p->H10>0 && p->W90==0 && p->F300==0 && p->H3==1)
	pupdate = new fluid_update_fsf_heat(p,a,pgc,pheat);
    
    if(p->F30>0 && p->H10>0 && p->W90==0 && p->F300==0 && p->H3==2)
    pupdate = new fluid_update_fsf_heat_Bouss(p,a,pgc,pheat);
    
    if(p->F30>0 && p->C10>0 && p->W90==0 && p->F300==0)
    pupdate = new fluid_update_fsf_concentration(p,a,pgc,pconc);
    
    if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90>0)
    pupdate = new fluid_update_rheology(p);
    
    if(p->F300>0)
	pupdate = new fluid_update_void();

	if(p->F46==2)
	ppicard = new picard_f(p);

	if(p->F46==3)
	ppicard = new picard_lsm(p);

	if(p->F46!=2 && p->F46!=3)
	ppicard = new picard_void(p);

    alpha << 4.0/15.0, 1.0/15.0, 1.0/6.0;
    gamma << 8.0/15.0, 5.0/12.0, 3.0/4.0;
    zeta << 0.0, -17.0/60.0, -5.0/12.0;
}

momentum_FCLS3::~momentum_FCLS3(){}


void momentum_FCLS3::start(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans, sixdof *p6dof)
{	

    // Set inflow 
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
		
    bool final = false;

    for (int loop=0; loop<3; loop++)
    {
        if (loop==2) final = true;
        
        pflow->rkinflow(p,a,pgc,urk,vrk,wrk);
        
        // FSF
        pfsfdisc->start(p,a,a->phi,4,a->u,a->v,a->w);
        
        LOOP
        a->phi(i,j,k) += gamma(loop)*p->dt*CPOR1*a->L(i,j,k) + zeta(loop)*p->dt*CPOR1*Cf(i,j,k);
                    
        LOOP
        Cf(i,j,k)=a->L(i,j,k);
        
        pflow->phi_relax(p,pgc,a->phi);
        
        pgc->start4(p,a->phi,gcval_phi);
        
        
        p->F44=2;
        preini->start(a,p,a->phi, pgc, pflow);
        ppicard->correct_ls(p,a,pgc,a->phi);
        
        pupdate->start(p,a,pgc,a->u,a->v,a->w);
        
        
    // -------------------
        // U
        starttime=pgc->timer();

        // Fill F
        pturb->isource(p,a);
        pflow->isource(p,a,pgc,pvrans); 
        bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
        ppress->upgrad(p,a,a->eta,a->eta_n);
        irhs(p,a,pgc,a->u,a->u,a->v,a->w,2.0*alpha(loop));
        pdiff->diff_u(p,a,pgc,psolv,urk,a->u,a->u,a->v,a->w,2.0*alpha(loop));

        ULOOP
        urk(i,j,k) += 2.0*alpha(loop)*p->dt*CPOR1*a->F(i,j,k);

        // Add convection
        ULOOP
        a->F(i,j,k)=0.0;

        pconvec->start(p,a,a->u,1,a->u,a->v,a->w);

        ULOOP
        urk(i,j,k) += gamma(loop)*p->dt*CPOR1*a->F(i,j,k) + zeta(loop)*p->dt*CPOR1*Cu(i,j,k);
        
        ULOOP
        Cu(i,j,k)=a->F(i,j,k);

        p->utime+=pgc->timer()-starttime;
        
    // -------------------
        // V
        starttime=pgc->timer();

        // Add source
        pturb->jsource(p,a);
        pflow->jsource(p,a,pgc,pvrans);
        bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
        ppress->vpgrad(p,a,a->eta,a->eta_n);
        jrhs(p,a,pgc,a->v,a->u,a->v,a->w,2.0*alpha(loop));
        pdiff->diff_v(p,a,pgc,psolv,vrk,a->v,a->u,a->v,a->w,2.0*alpha(loop));
        
        VLOOP
        vrk(i,j,k) += 2.0*alpha(loop)*p->dt*CPOR2*a->G(i,j,k);

        // Add convection
        VLOOP
        a->G(i,j,k)=0.0;

        pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
        
        VLOOP
        vrk(i,j,k) += gamma(loop)*p->dt*CPOR2*a->G(i,j,k) + zeta(loop)*p->dt*CPOR2*Cv(i,j,k);
        
        VLOOP
        Cv(i,j,k)=a->G(i,j,k);

        p->vtime+=pgc->timer()-starttime;

    // -------------------
        // W
        starttime=pgc->timer();

        // Add source
        pturb->ksource(p,a);
        pflow->ksource(p,a,pgc,pvrans);
        bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
        ppress->wpgrad(p,a,a->eta,a->eta_n);
        krhs(p,a,pgc,a->w,a->u,a->v,a->w,2.0*alpha(loop));
        pdiff->diff_w(p,a,pgc,psolv,wrk,a->w,a->u,a->v,a->w,2.0*alpha(loop));

        WLOOP
        wrk(i,j,k) += 2.0*alpha(loop)*p->dt*CPOR3*a->H(i,j,k);
        
        // Add convection
        WLOOP
        a->H(i,j,k)=0.0;

        pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
        
        WLOOP
        wrk(i,j,k) += gamma(loop)*p->dt*CPOR3*a->H(i,j,k) + zeta(loop)*p->dt*CPOR3*Cw(i,j,k);
        
        WLOOP
        Cw(i,j,k)=a->H(i,j,k);

        p->wtime+=pgc->timer()-starttime;

        // Forcing
        momentum_forcing_start(a, p, pgc, p6dof, pfsi,
                           urk,vrk,wrk, fx, fy, fz, loop, 2.0*alpha(loop), final);
                           
        // Direct Forcing
        ULOOP
        {
        a->u(i,j,k) = urk(i,j,k);
        
        if(p->count<10)
        a->maxF = MAX(fabs(2.0*alpha(loop)*CPOR1*fx(i,j,k)), a->maxF);
        
        p->sfmax = MAX(fabs(2.0*alpha(loop)*CPOR1*fx(i,j,k)), p->sfmax);
        }
        
        VLOOP
        {
        a->v(i,j,k) = vrk(i,j,k);
        
        if(p->count<10)
        a->maxG = MAX(fabs(2.0*alpha(loop)*CPOR2*fy(i,j,k)), a->maxG);
        
        p->sfmax = MAX(fabs(2.0*alpha(loop)*CPOR2*fy(i,j,k)), p->sfmax);
        }
        
        WLOOP
        {
        a->w(i,j,k) = wrk(i,j,k);
        
        if(p->count<10)
        a->maxH = MAX(fabs(2.0*alpha(loop)*CPOR3*fz(i,j,k)), a->maxH);
        
        p->sfmax = MAX(fabs(2.0*alpha(loop)*CPOR3*fz(i,j,k)), p->sfmax);
        }

        // Pressure
        pflow->pressure_io(p,a,pgc);
        ppress->start(a,p,ppois,ppoissonsolv,pgc, pflow, a->u, a->v, a->w, 2.0*alpha(loop));
        
        pflow->u_relax(p,a,pgc,a->u);
        pflow->v_relax(p,a,pgc,a->v);
        pflow->w_relax(p,a,pgc,a->w);
        pflow->p_relax(p,a,pgc,a->press);

        pgc->start1(p,a->u,gcval_u);
        pgc->start2(p,a->v,gcval_v);
        pgc->start3(p,a->w,gcval_w);
    }
}

void momentum_FCLS3::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;

        ULOOP
        {
            a->maxF = MAX(fabs(a->rhsvec.V[n] + a->gi), a->maxF);
            a->F(i,j,k) += (a->rhsvec.V[n] + a->gi + p->W29_x + a->Fext(i,j,k))*PORVAL1;
            
            a->rhsvec.V[n] = 0.0;
            a->Fext(i,j,k) = 0.0;
            ++n;
        }
        
}

void momentum_FCLS3::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
    
        VLOOP
        {
            a->maxG = MAX(fabs(a->rhsvec.V[n] + a->gj), a->maxG);
            a->G(i,j,k) += (a->rhsvec.V[n] + a->gj + p->W29_y + a->Gext(i,j,k))*PORVAL2;
            
            a->rhsvec.V[n] = 0.0;
            a->Gext(i,j,k) = 0.0;
            ++n;
        }
}

void momentum_FCLS3::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;

        WLOOP
        {
            a->maxH = MAX(fabs(a->rhsvec.V[n] + a->gk), a->maxH);
            a->H(i,j,k) += (a->rhsvec.V[n] + a->gk + p->W29_z + a->Hext(i,j,k))*PORVAL3;
            
            a->rhsvec.V[n] = 0.0;
            a->Hext(i,j,k) = 0.0;
            ++n;
        }
}
