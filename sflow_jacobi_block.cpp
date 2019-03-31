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

#include"sflow_jacobi_block.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

sflow_jacobi_block::sflow_jacobi_block(lexer* p,fdm2D* b,ghostcell *pgc) : epsi(1e-19)
{	
	margin=3;
}

sflow_jacobi_block::~sflow_jacobi_block()
{
}

void sflow_jacobi_block::setup(lexer* p,fdm2D* b, ghostcell* pgc, int var, cpt2D &C)
{
}

void sflow_jacobi_block::start(lexer* p, fdm2D* b, ghostcell* pgc, slice &f, vec2D& xvec, vec2D& rhsvec, int var, int gcv, double stop_crit)
{
	p->preconiter=0;
    
	if(var==1)
    {
	fillxvec1(p,b,f);
    sizeS=p->sizeS1;
	solve(p,b,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,b->C1);
    }
	
	if(var==2)
    {
	fillxvec2(p,b,f);
    sizeS=p->sizeS2;
	solve(p,b,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,b->C2);
    }
	
	if(var==3||var==4)
    {
	fillxvec4(p,b,f);
    sizeS=p->sizeS4;
	solve(p,b,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,b->C4);
    }
	
	finalize(p,b,f,xvec,var);
}

void sflow_jacobi_block::solve(lexer* p,fdm2D* b, ghostcell* pgc, vec2D& xvec, vec2D& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt2D &C)
{
	solveriter=0;
    residual = 1.0e9;
   
    do{
   
		for(qn=0;qn<p->N13;++qn)
		{			
			NSLICELOOP
			{
				b->xvec.V[n] = (1.0-p->N17)*b->xvec.V[n]
				+ p->N17*(b->rhsvec.V[n]-(
				+ b->M.n[n]*b->xvec.V[Ip1_J] 
				+ b->M.s[n]*b->xvec.V[Im1_J]
				+ b->M.w[n]*b->xvec.V[I_Jp1]
				+ b->M.e[n]*b->xvec.V[I_Jm1]))/(b->M.p[n]+epsi);

			}
		++solveriter;
		}
		
         pgc->gcparaxvec2D(p,b,b->xvec,var);		
		
		residual=res_calc(p,b,xvec,pgc,C);
		
	
	}while((residual>=stop_crit) && (solveriter<maxiter));	
}

double sflow_jacobi_block::res_calc(lexer *p, fdm2D *b, vec2D &x, ghostcell *pgc, cpt2D &C)
{
	double y;
	double resi=0.0;

	NSLICELOOP
	{	
	y  = b->rhsvec.V[n]

		-(b->M.p[n]*x.V[n]
		+ b->M.n[n]*x.V[Ip1_J] 
		+ b->M.s[n]*x.V[Im1_J]
		+ b->M.w[n]*x.V[I_Jp1]
		+ b->M.e[n]*x.V[I_Jm1]);

	resi+=y*y;
	}

	resi=sqrt(pgc->globalsum(resi));

	return resi/double(p->cellnumtot2D);	
}

void sflow_jacobi_block::fillxvec1(lexer* p, fdm2D* b, slice& f)
{
	int count,q;
	
	count=0;
    SLICELOOP1
    {
    b->xvec.V[count]=f(i,j);
    ++count;
    }

    GCSL1LOOP
    {
    i=p->gcbsl1[n][0];
    j=p->gcbsl1[n][1];

        if(p->gcbsl1[n][3]==1)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i-1-q,j);
        ++count;
        }

        if(p->gcbsl1[n][3]==2)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i,j+1+q);
        ++count;
        }

        if(p->gcbsl1[n][3]==3)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i,j-1-q);
        ++count;
        }

        if(p->gcbsl1[n][3]==4)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i+1+q,j);
        ++count;
        }
    }
}

void sflow_jacobi_block::fillxvec2(lexer* p, fdm2D* b, slice& f)
{
	int count,q;
	
	count=0;
    SLICELOOP2
    {
    b->xvec.V[count]=f(i,j);
    ++count;
    }

    GCSL2LOOP
    {
    i=p->gcbsl2[n][0];
    j=p->gcbsl2[n][1];

        if(p->gcbsl2[n][3]==1)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i-1-q,j);
        ++count;
        }

        if(p->gcbsl2[n][3]==2)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i,j+1+q);
        ++count;
        }

        if(p->gcbsl2[n][3]==3)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i,j-1-q);
        ++count;
        }

        if(p->gcbsl2[n][3]==4)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i+1+q,j);
        ++count;
        }
    }
}

void sflow_jacobi_block::fillxvec4(lexer* p, fdm2D* b, slice& f)
{
	int count,q;
	
	count=0;
    SLICELOOP4
    {
    b->xvec.V[count]=f(i,j);
    ++count;
    }

    GCSL4LOOP
    {
    i=p->gcbsl4[n][0];
    j=p->gcbsl4[n][1];

        if(p->gcbsl4[n][3]==1)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i-1-q,j);
        ++count;
        }

        if(p->gcbsl4[n][3]==2)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i,j+1+q);
        ++count;
        }

        if(p->gcbsl4[n][3]==3)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i,j-1-q);
        ++count;
        }

        if(p->gcbsl4[n][3]==4)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i+1+q,j);
        ++count;
        }
    }
}

void sflow_jacobi_block::finalize(lexer *p, fdm2D *b, slice &f, vec2D &xvec, int var)
{
	if(var==1)
    {
        count=0;
        SLICELOOP1
        {
        f(i,j)=xvec.V[count];
        ++count;
        }
    }
	
	if(var==2)
    {
        count=0;
        SLICELOOP2
        {
        f(i,j)=xvec.V[count];
        ++count;
        }
    }
	
	if(var==3 || var==4)
    {
        count=0;
        SLICELOOP4
        {
        f(i,j)=xvec.V[count];
        ++count;
        }
    }
}


