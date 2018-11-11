/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"sip.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"

sip::sip(lexer* p,fdm* a,ghostcell *pgc):le(p),ls(p),lb(p),lp(p),un(p),uw(p),ut(p),res(p),epsi(1e-19)
{	
	margin=3;
}

sip::~sip()
{
}

void sip::start(lexer* p,fdm* a, ghostcell* pgc, field &xfield, vec& xvec, vec& rhsvec, int var, int gcv, double stop_crit)
{

	if(var==1)
    {
	fillxvec1(p,a,xfield);
	pgc->gcparaxvec(p,xvec,var); 
	setup(p,a,pgc,var,a->C1);
	solve(p,a,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,a->C1);
    }
	
	if(var==2)
    {
	fillxvec2(p,a,xfield);
	pgc->gcparaxvec(p,xvec,var); 
	setup(p,a,pgc,var,a->C2);
	solve(p,a,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,a->C2);
    }
	
	if(var==3)
    {
	fillxvec3(p,a,xfield);
	pgc->gcparaxvec(p,xvec,var); 
	setup(p,a,pgc,var,a->C3);
	solve(p,a,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,a->C3);
    }
	
	if(var==4)
    {
	fillxvec4(p,a,xfield);
	pgc->gcparaxvec(p,xvec,var); 
	setup(p,a,pgc,var,a->C4);
	solve(p,a,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,a->C4);
    }	
	
	finalize(p,a,xfield,xvec,var);
}

void sip::startF(lexer* p, ghostcell* pgc, double *f, vec& rhsvec, matrix_diag &M, int var, int gcv, double stop_crit)
{
}

void sip::setup(lexer* p,fdm* a, ghostcell* pgc, int var, cpt &C)
{
	if(var==1)
	sizeM=p->sizeM1;
	
	if(var==2)
	sizeM=p->sizeM2;
	
	if(var==3)
	sizeM=p->sizeM3;
	
	if(var==4)
	sizeM=p->sizeM4;
	
	alpha=0.9;

    NLOOP
	{
	ijk_s = C.s[n];
	ijk_e = C.e[n];
	ijk_b = C.b[n];
	
	lb.V[n] = a->M.b[n]/(1.0 + alpha*(un.V[ijk_b]+uw.V[ijk_b]));
	le.V[n] = a->M.e[n]/(1.0 + alpha*(un.V[ijk_e]+ut.V[ijk_e]));
	ls.V[n] = a->M.s[n]/(1.0 + alpha*(uw.V[ijk_s]+ut.V[ijk_s]));

	p1		  = alpha*(lb.V[n]*un.V[ijk_b] + le.V[n]*un.V[ijk_e]);
	p2		  = alpha*(lb.V[n]*uw.V[ijk_b] + ls.V[n]*uw.V[ijk_s]);
	p3		  = alpha*(le.V[n]*ut.V[ijk_e] + ls.V[n]*ut.V[ijk_s]);

	lp.V[n] = 1.0/(a->M.p[n]+p1+p2+p3-lb.V[n]*ut.V[ijk_b]
				-le.V[n]*uw.V[ijk_e]-ls.V[n]*un.V[ijk_s] + 1.0e-20);

	un.V[n] = (a->M.n[n]-p1)*lp.V[n];
	uw.V[n] = (a->M.w[n]-p2)*lp.V[n];
	ut.V[n] = (a->M.t[n]-p3)*lp.V[n];
	}
	
	pgc->gcparaxvec(p,un,var);
	pgc->gcparaxvec(p,uw,var);
	pgc->gcparaxvec(p,ut,var);
	pgc->gcparaxvec(p,lb,var);
	pgc->gcparaxvec(p,le,var);
	pgc->gcparaxvec(p,ls,var);
}

void sip::solve(lexer* p,fdm* a, ghostcell* pgc, vec& xvec, vec& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt &C)
{
	
	numiter=0;
	p->solveriter=0;
	
	do
	{		
		residual=0.0;
		NLOOP
		{	
		ijk_s = C.s[n];
		ijk_e = C.e[n];
		ijk_b = C.b[n];
	
		res.V[n]=resi=rhsvec.V[n]

			-(a->M.p[n]*xvec.V[n]
			+ a->M.n[n]*xvec.V[Ip1_J_K] 
			+ a->M.s[n]*xvec.V[Im1_J_K]
			+ a->M.w[n]*xvec.V[I_Jp1_K]
			+ a->M.e[n]*xvec.V[I_Jm1_K]
			+ a->M.t[n]*xvec.V[I_J_Kp1]
			+ a->M.b[n]*xvec.V[I_J_Km1]);

		residual+=resi*resi;

		res.V[n] = (res.V[n]-lb.V[n]*res.V[ijk_b]

					-le.V[n]*res.V[ijk_e]-ls.V[n]*res.V[ijk_s])*lp.V[n];
		}
		residual = sqrt(pgc->globalsum(residual))/double(p->cellnumtot);

		NLOOP
		{
		ijk_n = C.n[n];
		ijk_w = C.w[n];
		ijk_t = C.t[n];
		
		res.V[n] = res.V[n]-un.V[n]*res.V[ijk_n]-uw.V[n]*res.V[ijk_w]-ut.V[n]*res.V[ijk_t];

		xvec.V[n] += res.V[n];
		}

		++numiter;
		pgc->gcparaxvec(p,xvec,var); 
		gcupdate(p,a,pgc,xvec,var,gcv);

	}while((residual>=stop_crit) && (numiter<p->N46));
	

	p->solveriter=numiter;	
	
}

void sip::fillxvec1(lexer* p, fdm* a, field& f)
{
	int count,q;
	
	count=0;
    ULOOP
    {
    a->xvec.V[count]=f(i,j,k);
    ++count;
    }

    GC1LOOP
    {
    i=p->gcb1[n][0];
    j=p->gcb1[n][1];
    k=p->gcb1[n][2];

        if(p->gcb1[n][3]==1)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i-1-q,j,k);
        ++count;
        }

        if(p->gcb1[n][3]==2)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j+1+q,k);
        ++count;
        }

        if(p->gcb1[n][3]==3)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j-1-q,k);
        ++count;
        }

        if(p->gcb1[n][3]==4)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i+1+q,j,k);
        ++count;
        }

        if(p->gcb1[n][3]==5)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k-1-q);
        ++count;
        }

        if(p->gcb1[n][3]==6)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k+1+q);
        ++count;
        }
    }
}

void sip::fillxvec2(lexer* p, fdm* a, field& f)
{
	int count,q;
	
	count=0;
    VLOOP
    {
    a->xvec.V[count]=f(i,j,k);
    ++count;
    }

    GC2LOOP
    {
    i=p->gcb2[n][0];
    j=p->gcb2[n][1];
    k=p->gcb2[n][2];

        if(p->gcb2[n][3]==1)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i-1-q,j,k);
        ++count;
        }

        if(p->gcb2[n][3]==2)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j+1+q,k);
        ++count;
        }

        if(p->gcb2[n][3]==3)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j-1-q,k);
        ++count;
        }

        if(p->gcb2[n][3]==4)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i+1+q,j,k);
        ++count;
        }

        if(p->gcb2[n][3]==5)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k-1-q);
        ++count;
        }

        if(p->gcb2[n][3]==6)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k+1+q);
        ++count;
        }
    }
}

void sip::fillxvec3(lexer* p, fdm* a, field& f)
{
	int count,q;
	
	count=0;
    WLOOP
    {
    a->xvec.V[count]=f(i,j,k);
    ++count;
    }

    GC3LOOP
    {
    i=p->gcb3[n][0];
    j=p->gcb3[n][1];
    k=p->gcb3[n][2];

        if(p->gcb3[n][3]==1)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i-1-q,j,k);
        ++count;
        }

        if(p->gcb3[n][3]==2)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j+1+q,k);
        ++count;
        }

        if(p->gcb3[n][3]==3)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j-1-q,k);
        ++count;
        }

        if(p->gcb3[n][3]==4)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i+1+q,j,k);
        ++count;
        }

        if(p->gcb3[n][3]==5)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k-1-q);
        ++count;
        }

        if(p->gcb3[n][3]==6)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k+1+q);
        ++count;
        }
    }
}


void sip::fillxvec4(lexer* p, fdm* a, field& f)
{
	count=0;
    LOOP
    {
    a->xvec.V[count]=f(i,j,k);
    ++count;
    }

    GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

        if(p->gcb4[n][3]==1)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i-1-q,j,k);
        ++count;
        }

        if(p->gcb4[n][3]==2)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j+1+q,k);
        ++count;
        }

        if(p->gcb4[n][3]==3)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j-1-q,k);
        ++count;
        }

        if(p->gcb4[n][3]==4)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i+1+q,j,k);
        ++count;
        }

        if(p->gcb4[n][3]==5)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k-1-q);
        ++count;
        }

        if(p->gcb4[n][3]==6)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k+1+q);
        ++count;
        }
    }
	
}

void sip::finalize(lexer *p, fdm *a, field &xfield, vec &xvec, int var)
{
	if(var==1)
    {
        count=0;
        ULOOP
        {
        xfield(i,j,k)=xvec.V[count];
        ++count;
        }
    }
	
	if(var==2)
    {
        count=0;
        VLOOP
        {
        xfield(i,j,k)=xvec.V[count];
        ++count;
        }
    }
	
	if(var==3)
    {
        count=0;
        WLOOP
        {
        xfield(i,j,k)=xvec.V[count];
        ++count;
        }
    }
	
	if(var==4)
    {
        count=0;
        LOOP
        {
        xfield(i,j,k)=xvec.V[count];
        ++count;
        }
    }
}

void sip::gcpara_update(lexer* p, vec& vec, ghostcell* pgc)
{
}

void sip::gcupdate(lexer *p, fdm *a, ghostcell *pgc, vec &x, int var, int gcv)
{
	if(p->N15==3 && gcv>=40 && gcv<=45)
	pgc->start4V(p,x,241);
	
	if(gcv>=40 && gcv<=45 && (p->N15==1||p->N15==2) && p->count<=p->N16)
	pgc->start4V(p,x,240);
	
	if(p->N15==4 && gcv>=40 && gcv<=45)
	pgc->start4V(p,x,240);
}


