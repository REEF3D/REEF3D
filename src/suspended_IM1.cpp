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
    bcsusp_start(p,a,pgc,s,a->conc);
	psolv->start(p,a,pgc,a->conc,a->rhsvec,4);
	sedfsf(p,a,a->conc);
	pgc->start4(p,a->conc,gcval_susp);
    fillconc(p,a,s);
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
    wvel(i,j,k) = a->w(i,j,k) - s->ws;
    
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
        
        conc(i,j,k) =    s->cb(i,j);
        conc(i,j,k-1) =  s->cb(i,j);
        conc(i,j,k-2) =  s->cb(i,j);
        conc(i,j,k-3) =  s->cb(i,j);
    }
    
        n=0;
        LOOP
        {
            if(p->flag4[Im1JK]<0)
            {
            a->rhsvec.V[n] -= a->M.s[n]*conc(i-1,j,k);
            a->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0)
            {
            a->rhsvec.V[n] -= a->M.n[n]*conc(i+1,j,k);
            a->M.n[n] = 0.0;
            }
            
            if(p->flag4[IJm1K]<0 && p->j_dir==1)
            {
            a->rhsvec.V[n] -= a->M.e[n]*conc(i,j-1,k);
            a->M.e[n] = 0.0;
            }
            
            if(p->flag4[IJp1K]<0 && p->j_dir==1)
            {
            a->rhsvec.V[n] -= a->M.w[n]*conc(i,j+1,k);
            a->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0)
            {
            a->rhsvec.V[n] -= a->M.b[n]*conc(i,j,k-1);
            a->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0)
            {
            a->rhsvec.V[n] -= a->M.t[n]*conc(i,j,k+1);
            a->M.t[n] = 0.0;
            }

        ++n;
        }
}

void suspended_IM1::fillconc(lexer* p, fdm* a, sediment_fdm *s)
{
    double dist;
    double d50=p->S20;
    double adist=0.5*d50;
    double deltab=3.0*d50;

    double cx,cy;
    
    if(p->S34==1)
    GC4LOOP
    if(p->gcb4[n][4]==5)
    {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];
        
        //s->conc(i,j) = a->conc(i,j,k+1);
        
        //dist = p->ZP[KP1]-s->bedzh(i,j)-adist;
        
        //s->conc(i,j) = (s->cbe(i,j)*(dist-deltab+adist) + a->conc(i,j,k+1)*(deltab-adist))/(dist);

    }

    
    if(p->S34==2)
    ILOOP
    JLOOP
    {
        cx=0.0;
        cy=0.0;
    
        KLOOP
        PCHECK
        {
        cx += 0.5*(a->u(i,j,k) + a->u(i-1,j,k))*a->conc(i,j,k)*p->DZN[KP];
        cy += 0.5*(a->v(i,j,k) + a->v(i,j-1,k))*a->conc(i,j,k)*p->DZN[KP];
        }
    s->conc(i,j) = sqrt(cx*cx + cy*cy);
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
    a->L(i,j,k)=0.0;
	++count;
    }
}
