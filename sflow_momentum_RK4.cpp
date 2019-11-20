/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"sflow_momentum_RK4.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sflow_convection.h"
#include"sflow_pressure.h"
#include"sflow_diffusion.h"
#include"sflow_fsf.h"
#include"sflow_boussinesq.h"
#include"sflow_rough_manning.h"
#include"sflow_rough_void.h"
#include"ioflow.h"
#include"solver2D.h"

sflow_momentum_RK4::sflow_momentum_RK4(lexer *p, fdm2D *b, sflow_convection *pconvection, sflow_diffusion *ppdiff, sflow_pressure* ppressure,
                                                    solver2D *psolver, solver2D *ppoissonsolver, ioflow *pioflow, sflow_fsf *pfreesurf,
                                                    sflow_boussinesq *ppbouss)
                                                    :Prk1(p),Prk2(p),Prk3(p),Qrk1(p),Qrk2(p),Qrk3(p),
														Prk(p),Qrk(p),wrk1(p),wrk2(p),wrk3(p),wrk(p),
                                                      etark1(p),etark2(p),etark3(p),etark(p),etan(p)
{
	gcval_u=10;
	gcval_v=11;
    gcval_w=12;

	gcval_urk=20;
	gcval_vrk=21;
    gcval_wrk=12;
    
    if(p->F50==1)
	gcval_eta = 51;
    
    if(p->F50==2)
	gcval_eta = 52;
    
    if(p->F50==3)
	gcval_eta = 53;
    
    if(p->F50==4)
	gcval_eta = 54;

	pconvec=pconvection;
	pdiff=ppdiff;
	ppress=ppressure;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
	pfsf=pfreesurf;
    pbouss=ppbouss;
    
    if(p->A218==0)
    prough = new sflow_rough_void(p);
    
    if(p->A218==1)
    prough = new sflow_rough_manning(p);
}

sflow_momentum_RK4::~sflow_momentum_RK4()
{
}

void sflow_momentum_RK4::start(lexer *p, fdm2D* b, ghostcell* pgc)
{	
    pflow->discharge2D(p,b,pgc);
    pflow->inflow2D(p,b,pgc,b->P,b->Q,b->bed,b->eta);
	pflow->rkinflow2D(p,b,pgc,Prk,Qrk,b->bed,b->eta);

//Step 1
//--------------------------------------------------------	
     // eta
    pfsf->wetdry(p,b,pgc,b->P,b->Q,b->ws);
                   
    SLICELOOP4
    {
    etark1(i,j) =      -p->dt*(b->P(i,j)*b->hx(i,j) - b->P(i-1,j)*b->hx(i-1,j)
                       +       b->Q(i,j)*b->hy(i,j) - b->Q(i,j-1)*b->hy(i,j-1))/p->dx;
                     
    etark(i,j)  = b->eta(i,j) + 0.5*etark1(i,j);
    }
                
    pgc->gcsl_start4(p,etark,gcval_eta);
    pfsf->depth_update(p,b,pgc,b->P,b->Q,b->ws,etark);
    pfsf->breaking(p,b,pgc,etark,b->eta,0.5);
    pflow->eta_relax(p,pgc,etark);
    pgc->gcsl_start4(p,etark,gcval_eta);
    
    // U
	starttime=pgc->timer();
	pflow->isource2D(p,b,pgc); 
    pbouss->psi1(p,b,pgc,b->P,b->Q,b->eta,0.5);
	ppress->upgrad(p,b,b->eta,b->eta);
	irhs(p,b,pgc,b->P,0.5);
	pconvec->start(p,b,b->P,1,b->P,b->Q);
	pdiff->diff_u(p,b,pgc,psolv,b->P,b->Q,0.5);

	SLICELOOP1
	{
	Prk1(i,j) = p->dt*b->F(i,j);
	Prk(i,j)  = b->P(i,j) + 0.5*Prk1(i,j);
	}
    
    pgc->gcsl_start1(p,Prk,gcval_urk);
    
    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pflow->jsource2D(p,b,pgc);
    pbouss->psi2(p,b,pgc,b->P,b->Q,b->eta,0.5);
	ppress->vpgrad(p,b,etark,b->eta);
	jrhs(p,b,pgc,b->Q,0.5);
	pconvec->start(p,b,b->Q,2,b->P,b->Q);
	pdiff->diff_v(p,b,pgc,psolv,b->P,b->Q,0.5);

	SLICELOOP2
	{
	Qrk1(i,j) = p->dt*b->G(i,j);
	Qrk(i,j)  = b->Q(i,j) + 0.5*Qrk1(i,j);
	}
	
	pgc->gcsl_start2(p,Qrk,gcval_vrk);
    
    p->vtime=pgc->timer()-starttime;
    
    // W
    SLICELOOP4
    b->L(i,j)=0.0;
    
    if(p->A214==1)
    pconvec->start(p,b,b->ws,4,b->P,b->Q);
    ppress->wpgrad(p,b,etark,b->eta);
    
    SLICELOOP4
    {
	wrk1(i,j) = p->dt*b->L(i,j);
    wrk(i,j)  = b->ws(i,j) + 0.5*wrk1(i,j);
    }
              
    pgc->gcsl_start4(p,wrk,12);
	
	//--------------------------------------------------------
	// pressure
    ppress->start(p,b,pgc,ppoissonsolv,pflow, Prk, Qrk, b->Pn, b->Qn, wrk, etark, 0.5);
    
	pflow->pm_relax(p,pgc,b->press);
	pflow->um_relax(p,pgc,Prk,b->bed,b->eta);
	pflow->vm_relax(p,pgc,Qrk,b->bed,b->eta);
    pflow->wm_relax(p,pgc,wrk,b->bed,b->eta);

	pgc->gcsl_start1(p,Prk,gcval_urk);
	pgc->gcsl_start2(p,Qrk,gcval_vrk);
    pgc->gcsl_start4(p,wrk,gcval_wrk);
	
//Step 2
//--------------------------------------------------------
	
    // eta
    pfsf->wetdry(p,b,pgc,Prk,Qrk,wrk);
    
    SLICELOOP4
    {
    etark2(i,j) =      -p->dt*(Prk(i,j)*b->hx(i,j) - Prk(i-1,j)*b->hx(i-1,j)
                       +       Qrk(i,j)*b->hy(i,j) - Qrk(i,j-1)*b->hy(i,j-1))/p->dx;
                    
    etan(i,j) = etark(i,j); 
    etark(i,j)  = b->eta(i,j) + 0.5*etark2(i,j);
    }
                
    pgc->gcsl_start4(p,etark,gcval_eta);
    pfsf->depth_update(p,b,pgc,Prk,Qrk,wrk,etark);
    pfsf->breaking(p,b,pgc,etan,b->eta,0.5);
    pflow->eta_relax(p,pgc,etark);
    pgc->gcsl_start4(p,etark,gcval_eta);
	
	// U
	starttime=pgc->timer();
	pflow->isource2D(p,b,pgc); 
    pbouss->psi1(p,b,pgc,Prk,Qrk,etark,0.5);
	ppress->upgrad(p,b,b->eta,etan);
	irhs(p,b,pgc,Prk,0.5);
	pconvec->start(p,b,Prk,1,Prk,Qrk);
	pdiff->diff_u(p,b,pgc,psolv,Prk,Qrk,0.5);

	SLICELOOP1
	{
	Prk2(i,j) = p->dt*b->F(i,j);
	Prk(i,j)  = b->P(i,j) + 0.5*Prk2(i,j);
	}
    
    pgc->gcsl_start1(p,Prk,gcval_urk);
                
    p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	pflow->jsource2D(p,b,pgc);
    pbouss->psi2(p,b,pgc,Prk,Qrk,etark,0.5);
	ppress->vpgrad(p,b,etark,etan);
	jrhs(p,b,pgc,Qrk,0.5);
	pconvec->start(p,b,Qrk,2,Prk,Qrk);
	pdiff->diff_v(p,b,pgc,psolv,Prk,Qrk,0.5);

	SLICELOOP2
	{
	Qrk2(i,j) = p->dt*b->G(i,j);
	Qrk(i,j)  = b->Q(i,j) + 0.5*Qrk2(i,j);
	}
	
	pgc->gcsl_start2(p,Qrk,gcval_vrk);

    p->vtime+=pgc->timer()-starttime;
	
	// W
    SLICELOOP4
    b->L(i,j)=0.0;
    
    if(p->A214==1)
    pconvec->start(p,b,b->ws,4,b->P,b->Q);
    ppress->wpgrad(p,b,etark,b->eta);
    
    SLICELOOP4
    {
	wrk2(i,j) = p->dt*b->L(i,j);
    wrk(i,j)  = b->ws(i,j) + 0.5*wrk2(i,j);
    }
              
    pgc->gcsl_start4(p,wrk,12);
	
	//--------------------------------------------------------
	// pressure
	ppress->start(p,b,pgc,ppoissonsolv,pflow, Prk, Qrk, b->Pn, b->Qn, wrk, etark, 0.5);
	
    pflow->pm_relax(p,pgc,b->press);
	pflow->um_relax(p,pgc,Prk,b->bed,b->eta);
	pflow->vm_relax(p,pgc,Qrk,b->bed,b->eta);
    pflow->wm_relax(p,pgc,wrk,b->bed,b->eta);

	pgc->gcsl_start1(p,Prk,gcval_urk);
	pgc->gcsl_start2(p,Qrk,gcval_vrk);
    pgc->gcsl_start4(p,wrk,gcval_wrk);

//Step 3
//--------------------------------------------------------
	
    // eta
    pfsf->wetdry(p,b,pgc,Prk,Qrk,wrk);
    
    SLICELOOP4
    {
    etark3(i,j) =      -p->dt*(Prk(i,j)*b->hx(i,j) - Prk(i-1,j)*b->hx(i-1,j)
                       +       Qrk(i,j)*b->hy(i,j) - Qrk(i,j-1)*b->hy(i,j-1))/p->dx;

    etan(i,j) = etark(i,j); 
    etark(i,j)  = b->eta(i,j) + etark3(i,j);
    }
                
    pgc->gcsl_start4(p,etark,gcval_eta);
    pfsf->depth_update(p,b,pgc,Prk,Qrk,wrk,etark);
    pfsf->breaking(p,b,pgc,etark,b->eta,1.0);
    pflow->eta_relax(p,pgc,etark);
    pgc->gcsl_start4(p,etark,gcval_eta);
    
	// U
	starttime=pgc->timer();
	pflow->isource2D(p,b,pgc); 
    pbouss->psi1(p,b,pgc,Prk,Qrk,etark,1.0);
	ppress->upgrad(p,b,b->eta,etan);
	irhs(p,b,pgc,Prk,1.0);
	pconvec->start(p,b,Prk,1,Prk,Qrk);
	pdiff->diff_u(p,b,pgc,psolv,Prk,Qrk,1.0);

	SLICELOOP1
	{
	Prk3(i,j) = p->dt*b->F(i,j);
	Prk(i,j)  = b->P(i,j) + Prk3(i,j);
	}
    
    pgc->gcsl_start1(p,Prk,gcval_urk);
    
    p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	pflow->jsource2D(p,b,pgc);
    pbouss->psi2(p,b,pgc,Prk,Qrk,etark,1.0);
	ppress->vpgrad(p,b,etark,etan);
	jrhs(p,b,pgc,Qrk,1.0);
	pconvec->start(p,b,Qrk,2,Prk,Qrk);
	pdiff->diff_v(p,b,pgc,psolv,Prk,Qrk,1.0);

	SLICELOOP2
	{
	Qrk3(i,j) = p->dt*b->G(i,j);
	Qrk(i,j)  = b->Q(i,j) + Qrk3(i,j);
	}
	
	pgc->gcsl_start2(p,Qrk,gcval_vrk);
    
    p->vtime+=pgc->timer()-starttime;
		
	// W
    SLICELOOP4
    b->L(i,j)=0.0;
    
    if(p->A214==1)
    pconvec->start(p,b,b->ws,4,b->P,b->Q);
    ppress->wpgrad(p,b,etark,b->eta);
    
    SLICELOOP4
    {
	wrk3(i,j) = p->dt*b->L(i,j);
    wrk(i,j)  = b->ws(i,j) + wrk3(i,j);
    }
              
    pgc->gcsl_start4(p,wrk,12);
	
	//--------------------------------------------------------
	// pressure
	ppress->start(p,b,pgc,ppoissonsolv,pflow, Prk, Qrk, b->Pn, b->Qn, wrk, etark, 1.0);
	
    pflow->pm_relax(p,pgc,b->press);
	pflow->um_relax(p,pgc,Prk,b->bed,b->eta);
	pflow->vm_relax(p,pgc,Qrk,b->bed,b->eta);
    pflow->wm_relax(p,pgc,wrk,b->bed,b->eta);

	pgc->gcsl_start1(p,Prk,gcval_urk);
	pgc->gcsl_start2(p,Qrk,gcval_vrk);
    pgc->gcsl_start4(p,wrk,gcval_wrk);


//Step 4
//--------------------------------------------------------
    
    // eta
    pfsf->wetdry(p,b,pgc,Prk,Qrk,wrk);
    
    SLICELOOP4
    {
    etan(i,j) = b->eta(i,j);
    b->eta(i,j) = b->eta(i,j) + (1.0/6.0)*(etark1(i,j) + 2.0*etark2(i,j) + 2.0*etark3(i,j)  
                       - p->dt*(Prk(i,j)*b->hx(i,j) - Prk(i-1,j)*b->hx(i-1,j)
                       +        Qrk(i,j)*b->hy(i,j) - Qrk(i,j-1)*b->hy(i,j-1))/p->dx);
    }
    
                
    pgc->gcsl_start4(p,b->eta,gcval_eta);
    pfsf->depth_update(p,b,pgc,Prk,Qrk,wrk,b->eta);
    pfsf->breaking(p,b,pgc,b->eta,etan,1.0);
    pflow->eta_relax(p,pgc,b->eta);
    pgc->gcsl_start4(p,b->eta,gcval_eta);
    
	// U
	starttime=pgc->timer();
	pflow->isource2D(p,b,pgc); 
    pbouss->psi1(p,b,pgc,Prk,Qrk,b->eta,1.0);
	ppress->upgrad(p,b,etan,etark);
	irhs(p,b,pgc,Prk,1.0);
	pconvec->start(p,b,Prk,1,Prk,Qrk);
	pdiff->diff_u(p,b,pgc,psolv,Prk,Qrk,1.0);

	SLICELOOP1
	b->P(i,j) = b->P(i,j) + (1.0/6.0)*(Prk1(i,j) + 2.0*Prk2(i,j) + 2.0*Prk3(i,j) + p->dt*b->F(i,j));
	
    pgc->gcsl_start1(p,b->P,gcval_u);
    
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pflow->jsource2D(p,b,pgc);
    pbouss->psi2(p,b,pgc,Prk,Qrk,b->eta,1.0);
	ppress->vpgrad(p,b,b->eta,etark);
	jrhs(p,b,pgc,Qrk,1.0);
	pconvec->start(p,b,Qrk,2,Prk,Qrk);
	pdiff->diff_v(p,b,pgc,psolv,Prk,Qrk,1.0);

	SLICELOOP2
	b->Q(i,j) = b->Q(i,j) + (1.0/6.0)*(Qrk1(i,j) + 2.0*Qrk2(i,j) + 2.0*Qrk3(i,j) + p->dt*b->G(i,j));
	
	pgc->gcsl_start2(p,b->Q,gcval_v);
	
    p->vtime+=pgc->timer()-starttime;
    
	// W
    SLICELOOP4
    b->L(i,j)=0.0;
    
    if(p->A214==1)
    pconvec->start(p,b,wrk,4,Prk,Qrk);
    ppress->wpgrad(p,b,etark,b->eta);
    
    SLICELOOP4
    b->ws(i,j) = b->ws(i,j) + (1.0/6.0)*(wrk1(i,j) + 2.0*wrk2(i,j) + 2.0*wrk3(i,j) + p->dt*b->L(i,j));
              
    pgc->gcsl_start4(p,b->ws,12);
	
	//--------------------------------------------------------
	// pressure
	ppress->start(p,b,pgc,ppoissonsolv,pflow, b->P, b->Q, b->Pn, b->Qn, b->ws, b->eta, 1.0);
	
    pflow->pm_relax(p,pgc,b->press);
	pflow->um_relax(p,pgc,b->P,b->bed,b->eta);
	pflow->vm_relax(p,pgc,b->Q,b->bed,b->eta);
    pflow->wm_relax(p,pgc,b->ws,b->bed,b->eta);

	pgc->gcsl_start1(p,b->P,gcval_urk);
	pgc->gcsl_start2(p,b->Q,gcval_vrk);
    pgc->gcsl_start4(p,b->ws,gcval_wrk);
    
    SLICELOOP4
    etan(i,j) = b->eta(i,j);
    
    pgc->gcsl_start4(p,etan,gcval_eta);
}

void sflow_momentum_RK4::irhs(lexer *p, fdm2D *b, ghostcell *pgc, slice &f, double alpha)
{
    n=0;
	if(p->D20<3)
	SLICELOOP1
	{
	b->maxF=MAX(fabs(b->F(i,j)),b->maxF);
	b->rhsvec.V[n]=0.0;
	++n;
	}
}

void sflow_momentum_RK4::jrhs(lexer *p, fdm2D *b, ghostcell *pgc, slice &f, double alpha)
{
	n=0;
	if(p->D20<3)
	SLICELOOP2
	{
	b->maxG=MAX(fabs(b->G(i,j)),b->maxG);
	b->rhsvec.V[n]=0.0;
	++n;
	}
}

