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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#include"momentum_RKLS3_sf.h"
#include"vrans.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"convection.h"
#include"diffusion.h"
#include"density_df.h"
#include"density_sf.h"
#include"ediff2.h"
#include"pressure.h"
#include"poisson.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"

momentum_RKLS3_sf::momentum_RKLS3_sf
(
    lexer *p, 
    fdm *a, 
    ghostcell *pgc, 
    convection *pconvection, 
    diffusion *pdiffusion, 
    pressure* ppressure, 
    poisson* ppoisson,
    turbulence *pturbulence, 
    solver *psolver, 
    solver *ppoissonsolver, 
    ioflow *pioflow
):bcmom(p),urk(p),vrk(p),wrk(p),Cu(p),Cv(p),Cw(p),Du(p),Dv(p),Dw(p),fx(p),fy(p),fz(p)
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

	pdiff_e=new ediff2(p);
    
    pdensity = new density_sf(p);

    alpha << 4.0/15.0, 1.0/15.0, 1.0/6.0;
    gamma << 8.0/15.0, 5.0/12.0, 3.0/4.0;
    zeta << 0.0, -17.0/60.0, -5.0/12.0;
}

momentum_RKLS3_sf::~momentum_RKLS3_sf()
{
}

void momentum_RKLS3_sf::start(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans, sixdof *p6dof, vector<net*>& pnet)
{
}

void momentum_RKLS3_sf::starti(lexer* p, fdm* a, ghostcell* pgc, sixdof* p6dof, vrans* pvrans, vector<net*>& pnet, fsi* pfsi)
{	
    // Set inflow 
    double udisctime=0.0;
    double udiscstart=0.0;
    
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	//pflow->rkinflow(p,a,pgc,urk,vrk,wrk);
		
    bool final = false;

    for (int loop=0; loop<3; loop++)
    {
        if (loop==2) final = true;
        
        pflow->rkinflow(p,a,pgc,urk,vrk,wrk);
        pflow->pressure_io(p,a,pgc);
        
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

        udiscstart=pgc->timer();
        pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
        udisctime=pgc->timer()-udiscstart;

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
        pdiff->diff_v(p,a,pgc,psolv,vrk,a->v,a->v,a->v,a->w,2.0*alpha(loop));
        
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

        pgc->start1(p,urk,gcval_u);
        pgc->start2(p,vrk,gcval_v);
        pgc->start3(p,wrk,gcval_w);

        
    // -------------------
    // Forcing
        starttime=pgc->timer();
        
        ULOOP
        fx(i,j,k) = 0.0;
       
        VLOOP
        fy(i,j,k) = 0.0;
      
        WLOOP
        fz(i,j,k) = 0.0;
        
        pgc->start1(p,fx,10);
        pgc->start2(p,fy,11);
        pgc->start3(p,fz,12);           
        
        pgc->solid_forcing(p,a,2.0*alpha(loop),urk,vrk,wrk,fx,fy,fz);

        ULOOP
        {
        a->u(i,j,k) = urk(i,j,k) + 2.0*alpha(loop)*p->dt*CPOR1*fx(i,j,k);
        
        if(p->count<10)
        a->maxF = MAX(fabs(2.0*alpha(loop)*CPOR1*fx(i,j,k)), a->maxF);
        
        p->sfmax = MAX(fabs(2.0*alpha(loop)*CPOR1*fx(i,j,k)), p->sfmax);
        }
        
        VLOOP
        {
        a->v(i,j,k) = vrk(i,j,k) + 2.0*alpha(loop)*p->dt*CPOR2*fy(i,j,k);
        
        if(p->count<10)
        a->maxG = MAX(fabs(2.0*alpha(loop)*CPOR2*fy(i,j,k)), a->maxG);
        
        p->sfmax = MAX(fabs(2.0*alpha(loop)*CPOR2*fy(i,j,k)), p->sfmax);
        }
        
        WLOOP
        {
        a->w(i,j,k) = wrk(i,j,k) + 2.0*alpha(loop)*p->dt*CPOR3*fz(i,j,k);
        
        if(p->count<10)
        a->maxH = MAX(fabs(2.0*alpha(loop)*CPOR3*fz(i,j,k)), a->maxH);
        
        p->sfmax = MAX(fabs(2.0*alpha(loop)*CPOR3*fz(i,j,k)), p->sfmax);
        }

        pgc->start1(p,a->u,gcval_u);
        pgc->start2(p,a->v,gcval_v);
        pgc->start3(p,a->w,gcval_w);
        
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

void momentum_RKLS3_sf::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	if(p->D20<3)
    {
        ULOOP
        {
            a->maxF = MAX(fabs(a->rhsvec.V[n] + a->gi), a->maxF);
            
            a->F(i,j,k) += (a->rhsvec.V[n] + a->gi + p->W29_x)*PORVAL1;
            
            a->rhsvec.V[n] = 0.0;
            
            ++n;
        }
    }
	
	n=0;
	if(p->D20==3)
    {
        ULOOP
        {
            a->rhsvec.V[n] += a->gi;
            
            ++n;
        }
    }
}

void momentum_RKLS3_sf::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	if(p->D20<3)
    {
        VLOOP
        {
            a->maxG = MAX(fabs(a->rhsvec.V[n] + a->gj), a->maxG);
            
            a->G(i,j,k) += (a->rhsvec.V[n] + a->gj + p->W29_y)*PORVAL2;
            
            a->rhsvec.V[n]=0.0;
            
            ++n;
        }
    }
	
	n=0;
	if(p->D20==3)
    {
        VLOOP
        {
            a->rhsvec.V[n] += a->gj;
            
            ++n;
        }
    }
}

void momentum_RKLS3_sf::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	if(p->D20<3)
    {
        WLOOP
        {
            a->maxH = MAX(fabs(a->rhsvec.V[n] + a->gk), a->maxH);
            
            a->H(i,j,k) += (a->rhsvec.V[n] + a->gk + p->W29_z)*PORVAL3;
            
            a->rhsvec.V[n]=0.0;
            
            ++n;
        }
    }
	
	n=0;
	if(p->D20==3)
    {
        WLOOP
        {
            a->rhsvec.V[n] += a->gk;
            
            ++n;
        }
    }
}

void momentum_RKLS3_sf::utimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_RKLS3_sf::vtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_RKLS3_sf::wtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_RKLS3_sf::fillaij1(lexer *p, fdm *a, ghostcell* pgc, solver *psolv){}
void momentum_RKLS3_sf::fillaij2(lexer *p, fdm *a, ghostcell* pgc, solver *psolv){}
void momentum_RKLS3_sf::fillaij3(lexer *p, fdm *a, ghostcell* pgc, solver *psolv){}
