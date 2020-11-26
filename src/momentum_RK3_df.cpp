/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"momentum_RK3_df.h"
#include"vrans.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"convection.h"
#include"diffusion.h"
#include"density_f.h"
#include"ediff2.h"
#include"pressure.h"
#include"poisson.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"6DOF_df.h"
#include"net.h"

momentum_RK3_df::momentum_RK3_df
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
):bcmom(p),urk1(p),vrk1(p),wrk1(p),urk2(p),vrk2(p),wrk2(p),fx(p),fy(p),fz(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
	
	gcval_urk=20;
	gcval_vrk=21;
	gcval_wrk=22;

	pconvec=pconvection;
	pdiff=pdiffusion;
	ppress=ppressure;
	ppois=ppoisson;
	pturb=pturbulence;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;

    pdensity = new density_f(p);
    
    // Ini Forcing
    ULOOP
    {
        fx(i,j,k) = 0.0;
        a->fbh1(i,j,k) = 0.0;
    }
    VLOOP
    {
        fy(i,j,k) = 0.0;
        a->fbh2(i,j,k) = 0.0;
    }
    WLOOP
    {
        fz(i,j,k) = 0.0;
        a->fbh3(i,j,k) = 0.0;
    }
    LOOP
        a->fbh4(i,j,k) = 0.0;
    
    pgc->start1(p,fx,10);
    pgc->start2(p,fy,11);
    pgc->start3(p,fz,12);           
    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);
}

momentum_RK3_df::~momentum_RK3_df(){}


void momentum_RK3_df::start(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans){}

void momentum_RK3_df::starti(lexer* p, fdm* a, ghostcell* pgc, sixdof_df* p6dof_df, vrans* pvrans, vector<net*>& pnet)
{	
    // Set inflow 
    double udisctime=0.0;
    double udiscstart=0.0;
    
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
	pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);
		

//Step 1
//--------------------------------------------------------

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans); 
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a,pgc,a->u,a->u,a->v,a->w,1.0);
    udiscstart=pgc->timer();
	pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
    udisctime=pgc->timer()-udiscstart;
	pdiff->diff_u(p,a,pgc,psolv,a->u,a->v,a->w,1.0);

	ULOOP
	urk1(i,j,k) = a->u(i,j,k)
				+ p->dt*CPOR1*(a->F(i,j,k));

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a,pgc,a->v,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
	pdiff->diff_v(p,a,pgc,psolv,a->u,a->v,a->w,1.0);

	VLOOP
	vrk1(i,j,k) = a->v(i,j,k)
				+ p->dt*CPOR2*(a->G(i,j,k));

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a,pgc,a->w,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
	pdiff->diff_w(p,a,pgc,psolv,a->u,a->v,a->w,1.0);

	WLOOP
	wrk1(i,j,k) = a->w(i,j,k)
				+ p->dt*CPOR3*(a->H(i,j,k));
	
    p->wtime=pgc->timer()-starttime;
   
    // Forcing
    ULOOP
    {
        fx(i,j,k) = 0.0;
        a->fbh1(i,j,k) = 0.0;
    }
    VLOOP
    {
        fy(i,j,k) = 0.0;
        a->fbh2(i,j,k) = 0.0;
    }
    WLOOP
    {
        fz(i,j,k) = 0.0;
        a->fbh3(i,j,k) = 0.0;
    }
    LOOP
        a->fbh4(i,j,k) = 0.0;
    
    pgc->start1(p,fx,10);
    pgc->start2(p,fy,11);
    pgc->start3(p,fz,12);           
    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);
    
    p6dof_df->forcing(p,a,pgc,pvrans,pnet,1.0,urk1,vrk1,wrk1,fx,fy,fz);
	
    ULOOP
	urk1(i,j,k) += 1.0*p->dt*CPOR1*(fx(i,j,k));
	VLOOP
	vrk1(i,j,k) += 1.0*p->dt*CPOR2*(fy(i,j,k));
	WLOOP
	wrk1(i,j,k) += 1.0*p->dt*CPOR3*(fz(i,j,k));
    

	pgc->start1(p,urk1,gcval_urk);
	pgc->start2(p,vrk1,gcval_vrk);
	pgc->start3(p,wrk1,gcval_wrk);
	
	//urk1.ggcpol(p);
	//vrk1.ggcpol(p);
	//wrk1.ggcpol(p);
	
    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc, pflow, urk1, vrk1, wrk1, 1.0);
	
	pflow->u_relax(p,a,pgc,urk1);
	pflow->v_relax(p,a,pgc,vrk1);
	pflow->w_relax(p,a,pgc,wrk1);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,urk1,gcval_urk);
	pgc->start2(p,vrk1,gcval_vrk);
	pgc->start3(p,wrk1,gcval_wrk);
	
//Step 2
//--------------------------------------------------------
	
	
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a,pgc,urk1,urk1,vrk1,wrk1,0.25);
    udiscstart=pgc->timer();
	pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
    udisctime+=pgc->timer()-udiscstart;
	pdiff->diff_u(p,a,pgc,psolv,urk1,vrk1,wrk1,0.25);

	ULOOP
	urk2(i,j,k) = 0.75*a->u(i,j,k) + 0.25*urk1(i,j,k)
				+ 0.25*p->dt*CPOR1*(a->F(i,j,k));
                
    p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a,pgc,vrk1,urk1,vrk1,wrk1,0.25);
	pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
	pdiff->diff_v(p,a,pgc,psolv,urk1,vrk1,wrk1,0.25);

	VLOOP
	vrk2(i,j,k) = 0.75*a->v(i,j,k) + 0.25*vrk1(i,j,k)
				+ 0.25*p->dt*CPOR2*(a->G(i,j,k));
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a,pgc,wrk1,urk1,vrk1,wrk1,0.25);
	pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
	pdiff->diff_w(p,a,pgc,psolv,urk1,vrk1,wrk1,0.25);

	WLOOP
	wrk2(i,j,k) = 0.75*a->w(i,j,k) + 0.25*wrk1(i,j,k)
				+ 0.25*p->dt*CPOR3*(a->H(i,j,k));

    p->wtime+=pgc->timer()-starttime;

    // Forcing
    ULOOP
    fx(i,j,k) = 0.0; 
    VLOOP
    fy(i,j,k) = 0.0;
    WLOOP
    fz(i,j,k) = 0.0;
    
    pgc->start1(p,fx,10);
    pgc->start2(p,fy,11);
    pgc->start3(p,fz,12);           

    p6dof_df->forcing(p,a,pgc,pvrans,pnet,0.25,urk2,vrk2,wrk2,fx,fy,fz);

	ULOOP
	urk2(i,j,k) += 0.25*p->dt*CPOR1*(fx(i,j,k));
	VLOOP
	vrk2(i,j,k) += 0.25*p->dt*CPOR2*(fy(i,j,k));
	WLOOP
	wrk2(i,j,k) += 0.25*p->dt*CPOR3*(fz(i,j,k));
	
    pgc->start1(p,urk2,gcval_urk);
	pgc->start2(p,vrk2,gcval_vrk);
	pgc->start3(p,wrk2,gcval_wrk);
	
	//urk2.ggcpol(p);
	//vrk2.ggcpol(p);
	//wrk2.ggcpol(p);

    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc, pflow, urk2, vrk2, wrk2, 0.25);
	
	pflow->u_relax(p,a,pgc,urk2);
	pflow->v_relax(p,a,pgc,vrk2);
	pflow->w_relax(p,a,pgc,wrk2);
	pflow->p_relax(p,a,pgc,a->press);
	
	pgc->start1(p,urk2,gcval_urk);
	pgc->start2(p,vrk2,gcval_vrk);
	pgc->start3(p,wrk2,gcval_wrk);


//Step 3
//--------------------------------------------------------

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a,pgc,urk2,urk2,vrk2,wrk2,2.0/3.0);
    udiscstart=pgc->timer();
	pconvec->start(p,a,urk2,1,urk2,vrk2,wrk2);
    udisctime+=pgc->timer()-udiscstart;
	pdiff->diff_u(p,a,pgc,psolv,urk2,vrk2,wrk2,2.0/3.0);

	ULOOP
	a->u(i,j,k) = (1.0/3.0)*a->u(i,j,k) + (2.0/3.0)*urk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR1*(a->F(i,j,k));
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a,pgc,vrk2,urk2,vrk2,wrk2,2.0/3.0);
	pconvec->start(p,a,vrk2,2,urk2,vrk2,wrk2);
	pdiff->diff_v(p,a,pgc,psolv,urk2,vrk2,wrk2,2.0/3.0);

	VLOOP
	a->v(i,j,k) = (1.0/3.0)*a->v(i,j,k) + (2.0/3.0)*vrk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR2*(a->G(i,j,k));
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a,pgc,wrk2,urk2,vrk2,wrk2,2.0/3.0);
	pconvec->start(p,a,wrk2,3,urk2,vrk2,wrk2);
	pdiff->diff_w(p,a,pgc,psolv,urk2,vrk2,wrk2,2.0/3.0);

	WLOOP
	a->w(i,j,k) = (1.0/3.0)*a->w(i,j,k) + (2.0/3.0)*wrk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR3*(a->H(i,j,k));
	
    p->wtime+=pgc->timer()-starttime;

    // Forcing
    ULOOP
    fx(i,j,k) = 0.0; 
    VLOOP
    fy(i,j,k) = 0.0;
    WLOOP
    fz(i,j,k) = 0.0;
    
    pgc->start1(p,fx,10);
    pgc->start2(p,fy,11);
    pgc->start3(p,fz,12);   

    p6dof_df->forcing(p,a,pgc,pvrans,pnet,2.0/3.0,a->u,a->v,a->w,fx,fy,fz);

	ULOOP
	a->u(i,j,k) += 2.0/3.0*p->dt*CPOR1*(fx(i,j,k));
	VLOOP
	a->v(i,j,k) += 2.0/3.0*p->dt*CPOR2*(fy(i,j,k));
	WLOOP
	a->w(i,j,k) += 2.0/3.0*p->dt*CPOR3*(fz(i,j,k));
	

    pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
	
	a->u.ggcpol(p);
	a->v.ggcpol(p);
	a->w.ggcpol(p);

	//--------------------------------------------------------
	// pressure
	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc, pflow, a->u, a->v,a->w,2.0/3.0);

    pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
}

void momentum_RK3_df::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	if(p->D20 < 3)
    {
        ULOOP
        {
            a->maxF = MAX(fabs(a->rhsvec.V[n] + a->gi), a->maxF);
            
            a->F(i,j,k) += (a->rhsvec.V[n] + a->gi)*PORVAL1;
            
            a->rhsvec.V[n] = 0.0;
            
            ++n;
        }
    }
	
	n=0;
	if(p->D20 == 3)
    {
        ULOOP
        {
            a->rhsvec.V[n] += a->gi;
            
            ++n;
        }
    }
}


void momentum_RK3_df::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	if(p->D20 < 3)
    {
        VLOOP
        {
            a->maxG = MAX(fabs(a->rhsvec.V[n] + a->gj), a->maxG);
            
            a->G(i,j,k) += (a->rhsvec.V[n] + a->gj)*PORVAL2;
            
            a->rhsvec.V[n]=0.0;
            
            ++n;
        }
    }
	
	n=0;
	if(p->D20 == 3)
    {
        VLOOP
        {
            a->rhsvec.V[n] += a->gj;
            
            ++n;
        }
    }
}

void momentum_RK3_df::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	if(p->D20 < 3)
    {
        WLOOP
        {
            a->maxH = MAX(fabs(a->rhsvec.V[n] + a->gk), a->maxH);
            
            a->H(i,j,k) += (a->rhsvec.V[n] + a->gk)*PORVAL3;
            
            a->rhsvec.V[n]=0.0;
            
            ++n;
        }
    }
	
	n=0;
	if(p->D20 == 3)
    {
        WLOOP
        {
            a->rhsvec.V[n] += a->gk;
            
            ++n;
        }
    }
}


void momentum_RK3_df::utimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}


void momentum_RK3_df::vtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}


void momentum_RK3_df::wtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_RK3_df::fillaij1(lexer *p, fdm *a, ghostcell* pgc, solver *psolv){}
void momentum_RK3_df::fillaij2(lexer *p, fdm *a, ghostcell* pgc, solver *psolv){}
void momentum_RK3_df::fillaij3(lexer *p, fdm *a, ghostcell* pgc, solver *psolv){}
