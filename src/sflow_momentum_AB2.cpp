/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"sflow_momentum_AB2.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sflow_convection.h"
#include"sflow_pressure.h"
#include"sflow_diffusion.h"
#include"sflow_fsf.h"
#include"ioflow.h"
#include"solver2D.h"
#include"6DOF.h"

sflow_momentum_AB2::sflow_momentum_AB2(lexer *p, fdm2D *b, sflow_convection *pconvection, sflow_diffusion *ppdiff, sflow_pressure* ppressure,
                                                    solver2D *psolver, solver2D *ppoissonsolver, ioflow *pioflow, sflow_fsf *pfreesurf, sixdof *pp6dof)
                                                    :Pab(p),Qab(p)
{
	gcval_u=10;
	gcval_v=11;

	pconvec=pconvection;
	pdiff=ppdiff;
	ppress=ppressure;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
	pfsf=pfreesurf;
}

sflow_momentum_AB2::~sflow_momentum_AB2()
{
}

void sflow_momentum_AB2::start(lexer *p, fdm2D* b, ghostcell* pgc)
{	/*
    pflow->discharge2D(p,b,pgc);
    pflow->inflow2D(p,b,pgc,b->P,b->Q,b->bed,b->eta);

//--------------------------------------------------------

	// U
	starttime=pgc->timer();

	pflow->isource2D(p,b,pgc);
	ppress->upgrad(p,b,b->eta);
	irhs(p,b,pgc,b->P,1.0);
	pconvec->start(p,b,b->P,1,b->P,b->Q);
	pdiff->diff_u(p,b,pgc,psolv,b->P,1.0);
	
    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pflow->jsource2D(p,b,pgc);
	ppress->vpgrad(p,b,b->eta);
	jrhs(p,b,pgc,b->Q,1.0);
	pconvec->start(p,b,b->Q,2,b->P,b->Q);
	pdiff->diff_v(p,b,pgc,psolv,b->Q,1.0);

    p->vtime=pgc->timer()-starttime;


	pgc->gcsl_start1(p,b->P,gcval_u);
	pgc->gcsl_start2(p,b->Q,gcval_v);
	
	
	pflow->pm_relax(p,pgc,b->press);
	
	//--------------------------------------------------------
	//U 2
    starttime=pgc->timer();

	if(p->count==1)
	SLICELOOP1
	Pab(i,j)=b->F(i,j);

	SLICELOOP1
	{
	b->P(i,j)+=p->dt*0.5*(((p->dt+2.0*p->dt_old)/p->dt_old)*b->F(i,j) - (p->dt/p->dt_old)*Pab(i,j));

	Pab(i,j)=b->F(i,j);
	}
	pgc->gcsl_start1(p,b->P,gcval_u);
	
	p->utime+=pgc->timer()-starttime;
	
	//--------------------------------------------------------
    //V 2
    starttime=pgc->timer();

	if(p->count==1)
	SLICELOOP2
	Qab(i,j)=b->G(i,j);

	SLICELOOP2
	{
	b->Q(i,j)+=p->dt*0.5*(((p->dt+2.0*p->dt_old)/p->dt_old)*b->G(i,j) - (p->dt/p->dt_old)*Qab(i,j));
	
	Qab(i,j)=b->G(i,j);
	}
	pgc->gcsl_start2(p,b->Q,gcval_v);
	
	p->vtime+=pgc->timer()-starttime;
	
	ppress->wcalc(p,b,1.0,b->P,b->Q,b->ws);
	pflow->pm_relax(p,pgc,b->press);
	//--------------------------------------------------------
	// pressure
	ppress->start(p,b,pgc,ppoissonsolv, b->P, b->Q, b->ws, b->eta, 1.0);
	
	pflow->um_relax(p,pgc,b->P,b->bed,b->eta);
	pflow->vm_relax(p,pgc,b->Q,b->bed,b->eta);

	pgc->gcsl_start1(p,b->P,gcval_u);
	pgc->gcsl_start2(p,b->Q,gcval_v);
	
	//pfsf->start(p,b,pgc,pflow,b->P,b->Q,1.0);*/
}

void sflow_momentum_AB2::irhs(lexer *p, fdm2D *b, ghostcell *pgc, slice &f, double alpha)
{

	n=0;
	if(p->D20<4)
	SLICELOOP1
	{
    b->maxF=MAX(fabs(b->rhsvec.V[n]),b->maxF);
	b->F(i,j) += (b->gi);
	
	b->maxF=MAX(fabs(b->F(i,j)),b->maxF);
	b->rhsvec.V[n]=0.0;
	++n;
	}
	
	n=0;
	if(p->D20==4)
	SLICELOOP1
	{
	b->rhsvec.V[n]+=b->gi;
	++n;
	}
}

void sflow_momentum_AB2::jrhs(lexer *p, fdm2D *b, ghostcell *pgc, slice &f, double alpha)
{
    
	n=0;
	if(p->D20<4)
	SLICELOOP2
	{
    b->maxG=MAX(fabs(b->rhsvec.V[n]),b->maxG);
	b->G(i,j) += (b->gj);
	
	b->maxG=MAX(fabs(b->G(i,j)),b->maxG);
	b->rhsvec.V[n]=0.0;
	++n;
	}
	
	n=0;
	if(p->D20==4)
	SLICELOOP2
	{
	b->rhsvec.V[n]+=b->gj;
	++n;
	}
}

