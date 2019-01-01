/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIM2ILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"heat_IM2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"fluid_update_fsf_heat.h"

heat_IM2::heat_IM2(lexer* p, fdm* a, ghostcell *pgc, heat *&pheat) : ibcheat(p), heat_print(p,a), Tn(p), Tnn(p), thermdiff(p)
{
	gcval_heat=80;

	
}

heat_IM2::~heat_IM2()
{
}

void heat_IM2::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow)
{	
    starttime=pgc->timer();
    diff_update(p,a,pgc);
    clearrhs(p,a,pgc);
    pconvec->start(p,a,T,4,a->u,a->v,a->w);
	pdiff->idiff_scalar(p,a,pgc,psolv,T,thermdiff,p->sigT,1.0);
	timesource(p,a,T);
	psolv->start(p,a,pgc,T,a->xvec,a->rhsvec,4,gcval_heat,p->N43);
	ibcheat_start(p,a,pgc,T);
	pgc->start4(p,T,gcval_heat);
	pupdate->start(p,a,pgc);
	p->heattime=pgc->timer()-starttime;
	p->heatiter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"heatiter: "<<p->suspiter<<"  heattime: "<<setprecision(3)<<p->susptime<<endl;
}

void heat_IM2::timesource(lexer* p, fdm* a, field& fn)
{
    int count=0;
    int q;

    count=0;
    LOOP
    {
        a->M.p[count]+= 1.5/PDT;

        a->rhsvec.V[count] += a->L(i,j,k) + (2.0*Tn(i,j,k))/PDT - Tnn(i,j,k)/(2.0*PDT);

	++count;
    }
}

void heat_IM2::ttimesave(lexer *p, fdm* a)
{
    LOOP
    {
    Tnn(i,j,k)=Tn(i,j,k);
    Tn(i,j,k)=T(i,j,k);
    }
}

void heat_IM2::diff_update(lexer *p, fdm *a, ghostcell *pgc)
{
    double alpha_air = p->H2;
	double alpha_water = p->H1;
    double H;
    double epsi=p->F45*p->dx;
    
    LOOP
	{
		if(a->phi(i,j,k)>epsi)
		H=1.0;

		if(a->phi(i,j,k)<-epsi)
		H=0.0;

		if(fabs(a->phi(i,j,k))<=epsi)
		H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));


		thermdiff(i,j,k)= alpha_water*H + alpha_air*(1.0-H);
	}
    
    pgc->start4(p,thermdiff,1);
}

void heat_IM2::clearrhs(lexer *p, fdm *a, ghostcell *pgc)
{
    int n=0;
	LOOP
	{
    a->L(i,j,k)=0.0;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

