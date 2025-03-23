/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Elyas Larkermani
--------------------------------------------------------------------*/
#include"heat_RK3CN.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"fluid_update_fsf_heat.h"

heat_RK3CN::heat_RK3CN(lexer* p, fdm* a, ghostcell *pgc, heat *&pheat) : bcheat(p), heat_print(p,a), thermdiff(p),ark1(p),ark2(p),Tdiff(p)
{
	gcval_heat=80;
}

heat_RK3CN::~heat_RK3CN()
{
}

void heat_RK3CN::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow)
{

// Step 1
    starttime=pgc->timer();
    diff_update(p,a,pgc);
    
    clearrhs(p,a,pgc);
    pconvec->start(p,a,T,4,a->u,a->v,a->w);
    addrhs(p,a,pgc,1.0);
    pdiff->diff_scalar(p,a,pgc,psolv,ark1,T,thermdiff,a->eddyv,p->sigT, 8.0/15.0);
    
    bcheat_start(p,a,pgc,ark1); 
    pgc->start4(p,ark1,gcval_heat);

// Step 2
    clearrhs(p,a,pgc);
    pconvec->start(p,a,ark1,4,a->u,a->v,a->w);
    addrhs(p,a,pgc,25.0/8.0);
    pconvec->start(p,a,T,4,a->u,a->v,a->w);
    addrhs(p,a,pgc,-17.0/8.0);
    pdiff->diff_scalar(p,a,pgc,psolv,ark2,ark1,thermdiff,a->eddyv,p->sigT, 2.0/15.0);
	
    bcheat_start(p,a,pgc,ark2);
    pgc->start4(p,ark2,gcval_heat);

// Step 3
    clearrhs(p,a,pgc);
    pconvec->start(p,a,ark2,4,a->u,a->v,a->w);
    addrhs(p,a,pgc,9.0/4.0);
    pconvec->start(p,a,ark1,4,a->u,a->v,a->w);
    addrhs(p,a,pgc,-5.0/4.0);
    pdiff->diff_scalar(p,a,pgc,psolv,T,ark2,thermdiff,a->eddyv,p->sigT, 1.0/3.0);
    
	
    bcheat_start(p,a,pgc,T);
    pgc->start4(p,T,gcval_heat);

    pupdate->start(p,a,pgc);
    p->susptime=pgc->timer()-starttime;

}

void heat_RK3CN::ttimesave(lexer *p, fdm* a)
{
}

void heat_RK3CN::diff_update(lexer *p, fdm *a, ghostcell *pgc)
{
    double alpha_1;
	double alpha_2;
    double H;
    double epsi=p->F45*p->DXM;
    
    if(p->H9==1)
    {
    alpha_1 = p->H1;
	alpha_2 = p->H2;
    }
    
    if(p->H9==2)
    {
    alpha_1 = p->H2;
	alpha_2 = p->H1;
    }
    
    LOOP
	{
		if(a->phi(i,j,k)>epsi)
		H=1.0;

		if(a->phi(i,j,k)<-epsi)
		H=0.0;

		if(fabs(a->phi(i,j,k))<=epsi)
		H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));

		thermdiff(i,j,k) = alpha_1*H + alpha_2*(1.0-H);
	}
    
    pgc->start4(p,thermdiff,1);
}


void heat_RK3CN::clearrhs(lexer *p, fdm *a, ghostcell *pgc)
{
    int n=0;
	LOOP
	{
    a->L(i,j,k)=0.0;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void heat_RK3CN::addrhs(lexer *p, fdm *a, ghostcell *pgc,double alpha)
{
    int n=0;
    LOOP
    {
    a->rhsvec.V[n] += alpha*a->L(i,j,k);
    a->L(i,j,k) = 0.0;
    ++n;
    }
}




