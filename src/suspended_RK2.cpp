/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"suspended_RK2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"solver.h"
#include"sediment_fdm.h"

suspended_RK2::suspended_RK2(lexer* p, fdm* a) : wvel(p)
{
	gcval_susp=60;
}

suspended_RK2::~suspended_RK2()
{
}


void suspended_RK2::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow, sediment_fdm *s)
{
    field4 ark1(p);
    fill_wvel(p,a,pgc,s);
    
// Step 1
    starttime=pgc->timer();

    suspsource(p,a,a->conc,s);
    pconvec->start(p,a,a->conc,4,a->u,a->v,wvel);
	pdiff->diff_scalar(p,a,pgc,psolv,a->conc,a->visc,a->eddyv,1.0,1.0);

	LOOP
	ark1(i,j,k) = a->conc(i,j,k)
				+ p->dt*a->L(i,j,k);
	
	pdiff->idiff_scalar(p,a,pgc,psolv,ark1,a->eddyv,1.0,1.0);
    bcsusp_start(p,a,pgc,s,ark1);
    sedfsf(p,a,ark1);
	pgc->start4(p,ark1,gcval_susp);



// Step 2
    suspsource(p,a,a->conc,s);
    pconvec->start(p,a,ark1,4,a->u,a->v,wvel);
	pdiff->diff_scalar(p,a,pgc,psolv,ark1,a->visc,a->eddyv,1.0,0.5);

	LOOP
	a->conc(i,j,k) = 0.5*a->conc(i,j,k)
				+ 0.5*ark1(i,j,k)
				+ 0.5*p->dt*a->L(i,j,k);
	
	pdiff->idiff_scalar(p,a,pgc,psolv,a->conc,a->eddyv,1.0,0.5);
    bcsusp_start(p,a,pgc,s,a->conc);
    sedfsf(p,a,a->conc);
	pgc->start4(p,a->conc,gcval_susp);
    fillconc(p,a,s);

	p->susptime=pgc->timer()-starttime;

}

void suspended_RK2::ctimesave(lexer *p, fdm* a)
{
}

void suspended_RK2::fill_wvel(lexer *p, fdm* a, ghostcell *pgc, sediment_fdm *s)
{
    WLOOP
    wvel(i,j,k) = a->w(i,j,k) - s->ws;
    
    pgc->start3(p,wvel,12);
}

void suspended_RK2::suspsource(lexer* p,fdm* a,field& conc, sediment_fdm *s)
{
    LOOP
    {
    a->L(i,j,k)=0.0;

    // if(a->phi(i,j,k)>0.0)    //a->L(i,j,k)=-s->ws*(conc(i,j,k+1)-conc(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]);
    }

}

void suspended_RK2::bcsusp_start(lexer* p, fdm* a,ghostcell *pgc, sediment_fdm *s, field& conc)
{
    GC4LOOP
    if(p->gcb4[n][4]==5)
    {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];
        
        conc(i,j,k) =  s->cb(i,j);
    }
}

void suspended_RK2::fillconc(lexer* p, fdm* a, sediment_fdm *s)
{
    GC4LOOP
    if(p->gcb4[n][4]==5)
    {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];
        
        s->conc(i,j) = a->conc(i,j,k+1);
    }
}

void suspended_RK2::sedfsf(lexer* p,fdm* a,field& conc)
{
    LOOP
    if(a->phi(i,j,k)<0.0)
    conc(i,j,k)=0.0;
}

void suspended_RK2::clearrhs(lexer* p, fdm* a)
{
    count=0;
    LOOP
    {
    a->rhsvec.V[count]=0.0;
	++count;
    }
}