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

#include"suspended_IM1.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"solver.h"
#include"sediment_fdm.h"

suspended_IM1::suspended_IM1(lexer* p, fdm* a) : concn(p),wvel(p)
{
	gcval_susp=60;
}

suspended_IM1::~suspended_IM1()
{
}

void suspended_IM1::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow, sediment_fdm *s)
{
    starttime=pgc->timer();
    clearrhs(p,a);
    fill_wvel(p,a,pgc,s);
    pconvec->start(p,a,a->conc,4,a->u,a->v,wvel);
	pdiff->idiff_scalar(p,a,pgc,psolv,a->conc,a->eddyv,1.0,1.0);
	suspsource(p,a,a->conc,s);
	timesource(p,a,a->conc);
	psolv->start(p,a,pgc,a->conc,a->rhsvec,4);
	bcsusp_start(p,a,pgc,s,a->conc);
	sedfsf(p,a,a->conc);
	pgc->start4(p,a->conc,gcval_susp);
	p->susptime=pgc->timer()-starttime;
	p->suspiter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"suspiter: "<<p->suspiter<<"  susptime: "<<setprecision(3)<<p->susptime<<endl;
}

void suspended_IM1::timesource(lexer* p, fdm* a, field& fn)
{
    int count=0;
    int q;

    LOOP
    {
        a->M.p[count]+= 1.0/DT;

        a->rhsvec.V[count] += a->L(i,j,k) + concn(i,j,k)/DT;

	++count;
    }
}

void suspended_IM1::ctimesave(lexer *p, fdm* a)
{
    LOOP
    concn(i,j,k)=a->conc(i,j,k);

}

void suspended_IM1::fill_wvel(lexer *p, fdm* a, ghostcell *pgc, sediment_fdm *s)
{
    WLOOP
    wvel(i,j,k) = a->w(i,j,k) + s->ws;
    
    pgc->start3(p,wvel,12);
}

void suspended_IM1::suspsource(lexer* p,fdm* a,field& conc, sediment_fdm *s)
{
    /*
    count=0;
    LOOP
    {
	if(a->phi(i,j,k)>0.0)
	a->rhsvec.V[count]  += -s->ws*(conc(i,j,k+1)-conc(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]);
	
	++count;
    }*/
}

void suspended_IM1::bcsusp_start(lexer* p, fdm* a,ghostcell *pgc, sediment_fdm *s, field& conc)
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

void suspended_IM1::sedfsf(lexer* p,fdm* a,field& conc)
{
    LOOP
    if(a->phi(i,j,k)<0.0)
    conc(i,j,k)=0.0;
}

void suspended_IM1::clearrhs(lexer* p, fdm* a)
{
    count=0;
    LOOP
    {
    a->rhsvec.V[count]=0.0;
	++count;
    }
}
