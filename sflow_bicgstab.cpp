/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"sflow_bicgstab.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

sflow_bicgstab::sflow_bicgstab(lexer* p,ghostcell *pgc):sj(p),rj(p),r0(p),vj(p),tj(p),pj(p),precoeff(p),
												ph(p),sh(p),aii(p),epsi(1e-19)
{	
	margin=3;
}

sflow_bicgstab::~sflow_bicgstab()
{
}

void sflow_bicgstab::setup(lexer* p, ghostcell* pgc, int var, cpt2D &C)
{
}

void sflow_bicgstab::start(lexer* p, ghostcell* pgc, slice &f, matrix2D &M, vec2D &xvec, vec2D &rhsvec, int var, int gcv, double stop_crit, cpt2D &C)
{
	p->preconiter=0;
    
	if(var==1)
    {
	fillxvec1(p,f,xvec);
    sizeS=p->sizeS1;
	solve(p,pgc,M,xvec,rhsvec,1,gcv,p->solveriter,p->N46,stop_crit,C);
    }
	
	if(var==2)
    {
	fillxvec2(p,f,xvec);
    sizeS=p->sizeS2;
	solve(p,pgc,M,xvec,rhsvec,2,gcv,p->solveriter,p->N46,stop_crit,C);
    }
	
	if(var==3||var==4)
    {
	fillxvec4(p,f,xvec);
    sizeS=p->sizeS4;
	solve(p,pgc,M,xvec,rhsvec,4,gcv,p->solveriter,p->N46,stop_crit,C);
    }
	
	finalize(p,f,xvec,var);
}

	
void sflow_bicgstab::solve(lexer* p, ghostcell* pgc, matrix2D &M, vec2D &xvec, vec2D &rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt2D &C)
{
	solveriter=0;
	residual = 1.0e9;

	// -----------------
	precon_setup(p,pgc,M,var,C);
	// -----------------

 restart:
    r_j=norm_r0=0.0;	
	pgc->gcparaxvec2D(p,xvec,var,C);
	
	matvec_axb(p,M,xvec,rhsvec,rj,C);
	
	NSLICELOOP
	{
		r0.V[n]=pj.V[n]=rj.V[n];
		r_j += rj.V[n]*r0.V[n];
    }

    r_j=pgc->globalsum(r_j);
    norm_r0=sqrt(r_j);

    if((residual>=stop_crit) && (solveriter<maxiter))
	{

	do{
	    sigma=0.0;
	    norm_vj=0.0;
	    norm_rj=0.0;
		
		// -------------------------
		precon_solve(p,pgc,ph,pj,var,240,p->preconiter,p->N13,p->N18,C);
		pgc->gcparaxvec2D(p,ph,var,C);		
		// -------------------------
		
		matvec_std(p,M,ph,vj,C);
		
		NSLICELOOP
		{
			sigma   += vj.V[n]*r0.V[n];
			norm_vj += vj.V[n]*vj.V[n];
			norm_rj += rj.V[n]*rj.V[n];
	    }
		
        sigma = pgc->globalsum(sigma);
		norm_vj = sqrt(pgc->globalsum(norm_vj));
		norm_rj = sqrt(pgc->globalsum(norm_rj));

	    alpha=r_j/sigma;

	if(fabs(sigma) <= (1.0e-12*(norm_vj*norm_r0)))
	{	
		residual=res_calc(p,pgc,M,xvec,rhsvec,C);
		++solveriter;

		goto restart;
	}

    if((fabs(alpha)*norm_vj/(norm_rj==0?1.0e-15:norm_rj))<=0.08)
	{
		residual=res_calc(p,pgc,M,xvec,rhsvec,C);
		++solveriter;

		goto restart;
	}

		norm_sj=0.0;
		
		NSLICELOOP
		{
		sj.V[n] = rj.V[n] - alpha*vj.V[n];
		norm_sj += sj.V[n]*sj.V[n];
		}

	    norm_sj=sqrt(pgc->globalsum(norm_sj));

    if(norm_sj>stop_crit)
	{
		// -------------------------
		precon_solve(p,pgc,sh,sj,var,240,p->preconiter,p->N13,p->N18,C);
        pgc->gcparaxvec2D(p,sh,var,C);		
		// -------------------------

		matvec_std(p,M,sh,tj,C);
		
		w1=w2=0.0;
		
		NSLICELOOP
		{
		    w1 += tj.V[n]*sj.V[n];
		    w2 += tj.V[n]*tj.V[n];
		}

		w1=pgc->globalsum(w1);
		w2=pgc->globalsum(w2);

		w=w1/(w2==0?1.0e-15:w2);

		r_j1=0.0;
		
		NSLICELOOP
		{
		xvec.V[n] += alpha*ph.V[n] + w*sh.V[n];
		rj.V[n]  = sj.V[n]-w*tj.V[n];
		r_j1 += rj.V[n]*r0.V[n];
		}

		r_j1=pgc->globalsum(r_j1);

		beta=alpha*r_j1/(w*r_j==0?1.0e-15:(w*r_j));
		
		NSLICELOOP
		pj.V[n] = rj.V[n] + beta*(pj.V[n]-w*vj.V[n]);
	}


	if(norm_sj<=stop_crit)
	{
	r_j1=0.0;
		
		NSLICELOOP
		{
		xvec.V[n] += alpha*ph.V[n];
		rj.V[n]=sj.V[n];
		r_j1 += rj.V[n]*r0.V[n];
		}

    r_j1=pgc->globalsum(r_j1);
	}

	    r_j = r_j1 ;

	    residual=0.0;
		
		NSLICELOOP
		residual += rj.V[n]*rj.V[n];

	    residual = sqrt(pgc->globalsum(residual))/double(p->cellnumtot2D);
		
	    ++solveriter;

	}while((residual>=stop_crit) && (solveriter<maxiter));

    } 
		
	NSLICELOOP4
	{
	ph.V[n]=0.0;
	sh.V[n]=0.0;
	}
	pgc->gcsl_start4V(p,ph,1,C);
	pgc->gcsl_start4V(p,sh,1,C);
    
}

void sflow_bicgstab::matvec_axb(lexer *p, matrix2D &M, vec2D &x, vec2D &rhsvec, vec2D &y, cpt2D &C)
{
	NSLICELOOP
	{
	y.V[n]  = rhsvec.V[n]

			-(M.p[n]*x.V[I_J]
			+ M.n[n]*x.V[Ip1_J] 
			+ M.s[n]*x.V[Im1_J]
			+ M.w[n]*x.V[I_Jp1]
			+ M.e[n]*x.V[I_Jm1]);
	}
}

void sflow_bicgstab::matvec_std(lexer *p, matrix2D &M, vec2D &x, vec2D &y, cpt2D &C)
{
	NSLICELOOP
	{
	y.V[n]   = M.p[n]*x.V[I_J]
            + M.n[n]*x.V[Ip1_J] 
            + M.s[n]*x.V[Im1_J]
            + M.w[n]*x.V[I_Jp1]
            + M.e[n]*x.V[I_Jm1];
	}
}

double sflow_bicgstab::res_calc(lexer *p, ghostcell *pgc, matrix2D &M, vec2D &x, vec2D &rhsvec, cpt2D &C)
{
	double y;
	double resi=0.0;

	NSLICELOOP
	{	
	y  = rhsvec.V[n]

		-(M.p[n]*x.V[n]
		+ M.n[n]*x.V[Ip1_J] 
		+ M.s[n]*x.V[Im1_J]
		+ M.w[n]*x.V[I_Jp1]
		+ M.e[n]*x.V[I_Jm1]);

	resi+=y*y;
	}

	resi=sqrt(pgc->globalsum(resi));

	return resi/double(p->cellnumtot2D);	
}

void sflow_bicgstab::fillxvec1(lexer* p, slice& f, vec2D &xvec)
{
	int count,q;
	
	count=0;
    SLICELOOP1
    {
    xvec.V[count]=f(i,j);
    ++count;
    }

    GCSL1LOOP
    {
    i=p->gcbsl1[n][0];
    j=p->gcbsl1[n][1];

        if(p->gcbsl1[n][3]==1)
        for(q=0;q<margin;++q)
        {
        xvec.V[count]=f(i-1-q,j);
        ++count;
        }

        if(p->gcbsl1[n][3]==2)
        for(q=0;q<margin;++q)
        {
        xvec.V[count]=f(i,j+1+q);
        ++count;
        }

        if(p->gcbsl1[n][3]==3)
        for(q=0;q<margin;++q)
        {
        xvec.V[count]=f(i,j-1-q);
        ++count;
        }

        if(p->gcbsl1[n][3]==4)
        for(q=0;q<margin;++q)
        {
        xvec.V[count]=f(i+1+q,j);
        ++count;
        }
    }
}

void sflow_bicgstab::fillxvec2(lexer* p, slice& f, vec2D &xvec)
{
	int count,q;
	
	count=0;
    SLICELOOP2
    {
    xvec.V[count]=f(i,j);
    ++count;
    }

    GCSL2LOOP
    {
    i=p->gcbsl2[n][0];
    j=p->gcbsl2[n][1];

        if(p->gcbsl2[n][3]==1)
        for(q=0;q<margin;++q)
        {
        xvec.V[count]=f(i-1-q,j);
        ++count;
        }

        if(p->gcbsl2[n][3]==2)
        for(q=0;q<margin;++q)
        {
        xvec.V[count]=f(i,j+1+q);
        ++count;
        }

        if(p->gcbsl2[n][3]==3)
        for(q=0;q<margin;++q)
        {
        xvec.V[count]=f(i,j-1-q);
        ++count;
        }

        if(p->gcbsl2[n][3]==4)
        for(q=0;q<margin;++q)
        {
        xvec.V[count]=f(i+1+q,j);
        ++count;
        }
    }
}

void sflow_bicgstab::fillxvec4(lexer* p, slice& f, vec2D &xvec)
{
	int count,q;
	
	count=0;
    SLICELOOP4
    {
    xvec.V[count]=f(i,j);
    ++count;
    }

    GCSL4LOOP
    {
    i=p->gcbsl4[n][0];
    j=p->gcbsl4[n][1];

        if(p->gcbsl4[n][3]==1)
        for(q=0;q<margin;++q)
        {
        xvec.V[count]=f(i-1-q,j);
        ++count;
        }

        if(p->gcbsl4[n][3]==2)
        for(q=0;q<margin;++q)
        {
        xvec.V[count]=f(i,j+1+q);
        ++count;
        }

        if(p->gcbsl4[n][3]==3)
        for(q=0;q<margin;++q)
        {
        xvec.V[count]=f(i,j-1-q);
        ++count;
        }

        if(p->gcbsl4[n][3]==4)
        for(q=0;q<margin;++q)
        {
        xvec.V[count]=f(i+1+q,j);
        ++count;
        }
    }
}

void sflow_bicgstab::finalize(lexer *p, slice &f, vec2D &xvec, int var)
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

void sflow_bicgstab::precon_setup(lexer* p, ghostcell* pgc, matrix2D &M, int var, cpt2D &C)
{
    NSLICELOOP
	aii.V[n] = -1.0/(M.p[n]+epsi);	
}

void sflow_bicgstab::precon_solve(lexer* p, ghostcell* pgc, vec2D& xvec, vec2D& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt2D &C)
{
	NSLICELOOP
	xvec.V[n] = rhsvec.V[n]*aii.V[n];	
	
}
