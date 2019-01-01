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

#include"jacobi_block.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"

jacobi_block::jacobi_block(lexer* p,fdm* a,ghostcell *pgc):epsi(1e-19)
{
	
	margin=3;	
}


void jacobi_block::start(lexer* p,fdm* a, ghostcell* pgc, field& xfield, vec& xvec, vec& rhsvec, int var, int gcv, double stop_crit)
{
    
    if(var==1)
    {
	fillxvec1(p,a,xfield);
	pgc->gcparaxvec(p,xvec,var); 
    sizeM=p->sizeM1;
	solve(p,a,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,a->C1);
    }
	
	if(var==2)
    {
	fillxvec2(p,a,xfield);
	pgc->gcparaxvec(p,xvec,var); 
    sizeM=p->sizeM2;
	solve(p,a,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,a->C2);
    }
	
	if(var==3)
    {
	fillxvec3(p,a,xfield);
	pgc->gcparaxvec(p,xvec,var); 
    sizeM=p->sizeM3;
	solve(p,a,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,a->C3);
    }
	
	if(var==4)
    {
	fillxvec4(p,a,xfield);
	pgc->gcparaxvec(p,xvec,var); 
    sizeM=p->sizeM4;
	solve(p,a,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,a->C4);
    }
	
	finalize(p,a,xfield,xvec,var);
}

void jacobi_block::startF(lexer* p, ghostcell* pgc, double *f, vec& rhsvec, matrix_diag &M, int var, int gcv, double stop_crit)
{
}

jacobi_block::~jacobi_block()
{
}

void jacobi_block::setup(lexer* p,fdm* a, ghostcell* pgc, int var, cpt &C)
{
}

void jacobi_block::solve(lexer* p,fdm* a, ghostcell* pgc, vec& xvec, vec& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt &C)
{	
    solveriter=0;
    residual = 1.0e9;
   
    do{
   
		for(qn=0;qn<p->N13;++qn)
		{			
			NLOOP
			{
				a->xvec.V[n] = (1.0-p->N17)*a->xvec.V[n]
				+ p->N17*(a->rhsvec.V[n]-(
				+ a->M.n[n]*a->xvec.V[Ip1_J_K] 
				+ a->M.s[n]*a->xvec.V[Im1_J_K]
				+ a->M.w[n]*a->xvec.V[I_Jp1_K]
				+ a->M.e[n]*a->xvec.V[I_Jm1_K]
				+ a->M.t[n]*a->xvec.V[I_J_Kp1]
				+ a->M.b[n]*a->xvec.V[I_J_Km1]))/(a->M.p[n]+epsi);

			}
		++solveriter;
		}
		
		pgc->gcparaxvec(p,a->xvec,var); 
		gcupdate(p,a,pgc,a->xvec,var,gcv);
		
		residual=rescalc(p,a,pgc,xvec,rhsvec,var,C);
		
	
	}while((residual>=stop_crit) && (solveriter<maxiter));
}

double jacobi_block::rescalc(lexer* p, fdm *a, ghostcell* pgc, vec& xvec, vec& rhsvec, int var, cpt &C)
{

    resi=residual=0.0;

		NLOOP
		{
			resi += 
			(a->rhsvec.V[n]-(
			+ a->M.p[n]*xvec.V[n] 
			+ a->M.n[n]*xvec.V[Ip1_J_K] 
			+ a->M.s[n]*xvec.V[Im1_J_K]
			+ a->M.w[n]*xvec.V[I_Jp1_K]
			+ a->M.e[n]*xvec.V[I_Jm1_K]
			+ a->M.t[n]*xvec.V[I_J_Kp1]
			+ a->M.b[n]*xvec.V[I_J_Km1]))/a->M.p[n];
			
			residual = resi*resi;
		}
		
		residual=sqrt(pgc->globalsum(residual))/double(p->cellnumtot);
		
		return residual;
}


void jacobi_block::fillxvec1(lexer* p, fdm* a, field& f)
{
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

void jacobi_block::fillxvec2( lexer* p, fdm* a, field& f)
{
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

void jacobi_block::fillxvec3( lexer* p, fdm* a, field& f)
{
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

void jacobi_block::fillxvec4( lexer* p, fdm* a, field& f)
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


void jacobi_block::finalize(lexer *p, fdm *a, field &xfield, vec &xvec, int var)
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


void jacobi_block::gcupdate(lexer *p, fdm *a, ghostcell *pgc, vec &x, int var, int gcv)
{
	if(p->N15==3 && gcv>=40 && gcv<=45)
	pgc->start4V(p,x,241);
	
	if(gcv>=40 && gcv<=45 && (p->N15==1||p->N15==2) && p->count<=p->N16)
	pgc->start4V(p,x,240);
	
	if(p->N15==4 && gcv>=40 && gcv<=45)
	pgc->start4V(p,x,240);	
}

void jacobi_block::gcpara_update(lexer* p, vec& vec, ghostcell* pgc)
{
    //update_vec(p,vec,A,pgc);
}


