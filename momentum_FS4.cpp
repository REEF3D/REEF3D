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

#include"momentum_FS4.h"
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

momentum_FS4::momentum_FS4(lexer *p, fdm *a, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, ioflow *pioflow)
                                                    :bcmom(p),urk1(p),vrk1(p),wrk1(p),urk2(p),vrk2(p),wrk2(p),urk3(p),vrk3(p),wrk3(p),
                                                    urk(p),vrk(p),wrk(p)
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
}

momentum_FS4::~momentum_FS4()
{
}

void momentum_FS4::start(lexer *p, fdm* a, ghostcell* pgc, momentum *pmom)
{
	pflow->discharge(p,a,pgc);
	pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
	pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);
	pflow->rkinflow(p,a,pgc,urk3,vrk3,wrk3);
	pflow->rkinflow(p,a,pgc,urk,vrk,wrk);

//Step 1
//--------------------------------------------------------

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a);
	pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
    pdiff->diff_u(p,a,pgc,psolv,a->u,a->v,a->w,0.5);

	ULOOP
	{
	urk1(i,j,k) = a->F(i,j,k);
	urk(i,j,k)  = a->u(i,j,k) + 0.5*p->dt*CPOR1*urk1(i,j,k);
	}
	
    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a);
	pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
	pdiff->diff_v(p,a,pgc,psolv,a->u,a->v,a->w,0.5);

	VLOOP
	{
	vrk1(i,j,k) = a->G(i,j,k);
	vrk(i,j,k)  = a->v(i,j,k) +  0.5*p->dt*CPOR2*vrk1(i,j,k);
	}
	
    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a);
	pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
	pdiff->diff_w(p,a,pgc,psolv,a->u,a->v,a->w,0.5);

	WLOOP
	{
	wrk1(i,j,k) = a->H(i,j,k);
	wrk(i,j,k)  = a->w(i,j,k) + 0.5*p->dt*CPOR3*wrk1(i,j,k);
	}
	
    p->wtime=pgc->timer()-starttime;

    ucorr(p,a,0.5,urk);
	vcorr(p,a,0.5,vrk);
	wcorr(p,a,0.5,wrk);

	pflow->u_relax(p,a,pgc,urk);
	pflow->v_relax(p,a,pgc,vrk);
	pflow->w_relax(p,a,pgc,wrk);

	pgc->start1(p,urk,gcval_urk);
	pgc->start2(p,vrk,gcval_vrk);
	pgc->start3(p,wrk,gcval_wrk);


//Step 2
//--------------------------------------------------------

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,urk,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a);
	pconvec->start(p,a,urk,1,urk,vrk,wrk);
	pdiff->diff_u(p,a,pgc,psolv,urk,vrk,wrk,0.5);

	ULOOP
    {
	urk2(i,j,k) = a->F(i,j,k);
	urk(i,j,k)  = a->u(i,j,k) + 0.5*p->dt*CPOR1*urk2(i,j,k);
    }

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,vrk,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a);
	pconvec->start(p,a,vrk,2,urk,vrk,wrk);
	pdiff->diff_v(p,a,pgc,psolv,urk,vrk,wrk,0.5);

	VLOOP
    {
	vrk2(i,j,k) = a->G(i,j,k);
	vrk(i,j,k)  = a->v(i,j,k) +  0.5*p->dt*CPOR2*vrk2(i,j,k);
    }

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a);
	pconvec->start(p,a,wrk,3,urk,vrk,wrk);
	pdiff->diff_w(p,a,pgc,psolv,urk,vrk,wrk,0.5);

	WLOOP
	{
	wrk2(i,j,k) = a->H(i,j,k);
	wrk(i,j,k)  = a->w(i,j,k) + 0.5*p->dt*CPOR3*wrk2(i,j,k);
	}
	
    p->wtime=pgc->timer()-starttime;

	ucorr(p,a,0.5,urk);
	vcorr(p,a,0.5,vrk);
	wcorr(p,a,0.5,wrk);

	pflow->u_relax(p,a,pgc,urk);
	pflow->v_relax(p,a,pgc,vrk);
	pflow->w_relax(p,a,pgc,wrk);

	pgc->start1(p,urk,gcval_urk);
	pgc->start2(p,vrk,gcval_vrk);
	pgc->start3(p,wrk,gcval_wrk);


//Step 3
//--------------------------------------------------------

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,urk,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a);
	pconvec->start(p,a,urk,1,urk,vrk,wrk);
	pdiff->diff_u(p,a,pgc,psolv,urk,vrk,wrk,1.0);

	ULOOP
	{
	urk3(i,j,k) = a->F(i,j,k);
	urk(i,j,k)  = a->u(i,j,k) + p->dt*CPOR1*urk3(i,j,k);
    }
	
    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,vrk,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a);
	pconvec->start(p,a,vrk,2,urk,vrk,wrk);
	pdiff->diff_v(p,a,pgc,psolv,urk,vrk,wrk,1.0);

	VLOOP
	{
	vrk3(i,j,k) = a->G(i,j,k);
	vrk(i,j,k)  = a->v(i,j,k) +  p->dt*CPOR2*vrk3(i,j,k);
    }
	
    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,wrk,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a);
	pconvec->start(p,a,wrk,3,urk,vrk,wrk);
	pdiff->diff_w(p,a,pgc,psolv,urk,vrk,wrk,1.0);

	WLOOP
	{
	wrk3(i,j,k) = a->H(i,j,k);
	wrk(i,j,k)  = a->w(i,j,k) + p->dt*CPOR3*wrk3(i,j,k);
	}
	
    p->wtime=pgc->timer()-starttime;

    ucorr(p,a,1.0,urk);
	vcorr(p,a,1.0,vrk);
	wcorr(p,a,1.0,wrk);

	pflow->u_relax(p,a,pgc,urk);
	pflow->v_relax(p,a,pgc,vrk);
	pflow->w_relax(p,a,pgc,wrk);

	pgc->start1(p,urk,gcval_urk);
	pgc->start2(p,vrk,gcval_vrk);
	pgc->start3(p,wrk,gcval_wrk);


//Step 4
//--------------------------------------------------------

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,urk,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a);
	pconvec->start(p,a,urk,1,urk,vrk,wrk);
	pdiff->diff_u(p,a,pgc,psolv,urk,vrk,wrk,1.0);

	ULOOP
	a->u(i,j,k) = a->u(i,j,k) + (1.0/6.0)*p->dt*CPOR1*(urk1(i,j,k) + 2.0*urk2(i,j,k) + 2.0*urk3(i,j,k) + a->F(i,j,k));

    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,vrk,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a);
	pconvec->start(p,a,vrk,2,urk,vrk,wrk);
	pdiff->diff_v(p,a,pgc,psolv,urk,vrk,wrk,1.0);

    VLOOP
	a->v(i,j,k) = a->v(i,j,k) + (1.0/6.0)*p->dt*CPOR2*(vrk1(i,j,k) + 2.0*vrk2(i,j,k) + 2.0*vrk3(i,j,k) + a->G(i,j,k));
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,wrk,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a);
	pconvec->start(p,a,wrk,3,urk,vrk,wrk);
	pdiff->diff_w(p,a,pgc,psolv,urk,vrk,wrk,1.0);

	WLOOP
	a->w(i,j,k) = a->w(i,j,k) + (1.0/6.0)*p->dt*CPOR3*(wrk1(i,j,k) + 2.0*wrk2(i,j,k) + 2.0*wrk3(i,j,k) + a->H(i,j,k));
	
    p->wtime+=pgc->timer()-starttime;

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);

	//--------------------------------------------------------
	// pressure
    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pmom,pflow,a->u,a->v,a->w,6.0/6.0);

	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
}

void momentum_FS4::ucorr(lexer *p, fdm *a, double alpha, field& uvel)
{
    ULOOP
	uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*((a->press(i+1,j,k)-a->press(i,j,k))
	/(p->DXP[IP]*roface(p,a,1,0,0)));
}

void momentum_FS4::vcorr(lexer *p, fdm *a, double alpha, field& vvel)
{
    VLOOP
	vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*((a->press(i,j+1,k)-a->press(i,j,k))
	/(p->DYP[JP]*roface(p,a,0,1,0)));
}

void momentum_FS4::wcorr(lexer *p, fdm *a, double alpha, field& wvel)
{
    WLOOP
	wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((a->press(i,j,k+1)-a->press(i,j,k))
	/(p->DZP[KP]*roface(p,a,0,0,1)));
}

void momentum_FS4::irhs(lexer *p, fdm *a)
{
	n=0;
	if(p->D20<3)
	ULOOP
	{
    a->maxF=MAX(fabs(a->rhsvec.V[n]),a->maxF);
	a->F(i,j,k) += (a->rhsvec.V[n] + a->gi)*PORVAL1;
	a->rhsvec.V[n]=0.0;
	++n;
	}
	
	n=0;
	if(p->D20==3)
	ULOOP
	{
	a->rhsvec.V[n]+=a->gi;
	++n;
	}
}

void momentum_FS4::jrhs(lexer *p, fdm *a)
{
	n=0;
	if(p->D20<3)
	VLOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n]),a->maxG);
	a->G(i,j,k) += (a->rhsvec.V[n] + a->gj)*PORVAL2;
	a->rhsvec.V[n]=0.0;
	++n;
	}
	
	n=0;
	if(p->D20==3)
	VLOOP
	{
	a->rhsvec.V[n]+=a->gj;
	++n;
	}
}

void momentum_FS4::krhs(lexer *p, fdm *a)
{
	n=0;
	if(p->D20<3)
	WLOOP
	{
    a->maxH=MAX(fabs(a->rhsvec.V[n]),a->maxH);
	a->H(i,j,k) += (a->rhsvec.V[n] + a->gk)*PORVAL3;
	a->rhsvec.V[n]=0.0;
	++n;
	}
	
	n=0;
	if(p->D20==3)
	WLOOP
	{
	a->rhsvec.V[n]+=a->gk;
	++n;
	}
}

void momentum_FS4::utimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FS4::vtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FS4::wtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_FS4::fillaij1(lexer *p, fdm *a, ghostcell* pgc, solver *psolv)
{
}

void momentum_FS4::fillaij2(lexer *p, fdm *a, ghostcell* pgc, solver *psolv)
{
}

void momentum_FS4::fillaij3(lexer *p, fdm *a, ghostcell* pgc, solver *psolv)
{
}

