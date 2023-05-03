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

#include"sflow_bicgstab.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

sflow_bicgstab::sflow_bicgstab(lexer* p,ghostcell *pgc):epsi(1e-19)
{	
	margin=3;
    
    p->Darray(sj,p->imax*p->jmax*p->kmax);
    p->Darray(rj,p->imax*p->jmax*p->kmax);
    p->Darray(r0,p->imax*p->jmax*p->kmax);
    p->Darray(vj,p->imax*p->jmax*p->kmax);
    p->Darray(tj,p->imax*p->jmax*p->kmax);
    p->Darray(pj,p->imax*p->jmax*p->kmax); 
    p->Darray(ph,p->imax*p->jmax*p->kmax);
    p->Darray(sh,p->imax*p->jmax*p->kmax);
    p->Darray(aii,p->imax*p->jmax*p->kmax);
    p->Darray(x,p->imax*p->jmax*p->kmax);
    p->Darray(rhs,p->imax*p->jmax*p->kmax);
}

sflow_bicgstab::~sflow_bicgstab()
{
}


void sflow_bicgstab::start(lexer* p, ghostcell* pgc, slice &f, matrix2D &M, vec2D &xvec, vec2D &rhsvec, int var)
{
	p->preconiter=0;
    
	
	if(var==1)
    {
    flagslice = p->flagslice1;
    ulast=p->ulast;
    vlast=0;
    }
	
	if(var==2)
    {
    flagslice = p->flagslice2;
    ulast=0;
    vlast=p->vlast;
    }
	
	if(var==4||var==5)
    {
    flagslice = p->flagslice4;
    ulast=0;
    vlast=0;
    }
    
    fillxvec(p,f,rhsvec);
	solve(p,pgc,M,xvec,rhsvec,var,p->solveriter);
	
	finalize(p,f);
}

	
void sflow_bicgstab::solve(lexer* p, ghostcell* pgc, matrix2D &M, vec2D &xvec, vec2D &rhsvec, int var, int &solveriter)
{
	solveriter=0;
	residual = 1.0e9;

	// -----------------
	precon_setup(p,M,pgc);
	// -----------------

 restart:
    r_j=norm_r0=0.0;	
	pgc->gcslparaxijk_single(p,x,var);
	
	matvec_axb(p,M,x,rj);
	
	SLICEFLEXLOOP
	{
		r0[IJ]=pj[IJ]=rj[IJ];
		r_j += rj[IJ]*r0[IJ];
    }

    r_j=pgc->globalsum(r_j);
    norm_r0=sqrt(r_j);

    if((residual>=p->N44) && (solveriter<p->N46))
	{

	do{
	    sigma=0.0;
	    norm_vj=0.0;
	    norm_rj=0.0;
		
		// -------------------------
		precon_solve(p,pgc,ph,pj);
		pgc->gcslparaxijk_single(p,ph,var);				
		// -------------------------
		
		matvec_std(p,M,ph,vj);
		
		SLICEFLEXLOOP
		{
			sigma   += vj[IJ]*r0[IJ];
			norm_vj += vj[IJ]*vj[IJ];
			norm_rj += rj[IJ]*rj[IJ];
	    }
		
        sigma = pgc->globalsum(sigma);
		norm_vj = sqrt(pgc->globalsum(norm_vj));
		norm_rj = sqrt(pgc->globalsum(norm_rj));

	    alpha=r_j/sigma;

	if(fabs(sigma) <= (1.0e-12*(norm_vj*norm_r0)))
	{	
		residual=res_calc(p,M,pgc,x);
		++solveriter;

		goto restart;
	}

    if((fabs(alpha)*norm_vj/(norm_rj==0?1.0e-15:norm_rj))<=0.08)
	{
		residual=res_calc(p,M,pgc,x);
		++solveriter;

		goto restart;
	}

		norm_sj=0.0;
		
		SLICEFLEXLOOP
		{
		sj[IJ] = rj[IJ] - alpha*vj[IJ];
		norm_sj += sj[IJ]*sj[IJ];
		}

	    norm_sj=sqrt(pgc->globalsum(norm_sj));

    if(norm_sj>p->N44)
	{
		// -------------------------
		precon_solve(p,pgc,sh,sj);
        pgc->gcslparaxijk_single(p,sh,var);		
		// -------------------------

		matvec_std(p,M,sh,tj);
		
		w1=w2=0.0;
		
		SLICEFLEXLOOP
		{
		    w1 += tj[IJ]*sj[IJ];
		    w2 += tj[IJ]*tj[IJ];
		}

		w1=pgc->globalsum(w1);
		w2=pgc->globalsum(w2);

		w=w1/(w2==0?1.0e-15:w2);

		r_j1=0.0;
		
		SLICEFLEXLOOP
		{
		x[IJ] += alpha*ph[IJ] + w*sh[IJ];
		rj[IJ]  = sj[IJ]-w*tj[IJ];
		r_j1 += rj[IJ]*r0[IJ];
		}

		r_j1=pgc->globalsum(r_j1);

		beta=alpha*r_j1/(w*r_j==0?1.0e-15:(w*r_j));
		
		SLICEFLEXLOOP
		pj[IJ] = rj[IJ] + beta*(pj[IJ]-w*vj[IJ]);
	}


	if(norm_sj<=p->N44)
	{
	r_j1=0.0;
		
		SLICEFLEXLOOP
		{
		x[IJ] += alpha*ph[IJ];
		rj[IJ]=sj[IJ];
		r_j1 += rj[IJ]*r0[IJ];
		}

    r_j1=pgc->globalsum(r_j1);
	}

	    r_j = r_j1 ;

	    residual=0.0;
		
		SLICEFLEXLOOP
		residual += rj[IJ]*rj[IJ];

	    residual = sqrt(pgc->globalsum(residual))/double(p->cellnumtot2D);
		
	    ++solveriter;
        
	}while((residual>=p->N44) && (solveriter<p->N46));

    } 
		
	
	SLICELOOP4
	{
	ph[IJ]=0.0;
	sh[IJ]=0.0;
	}
    
    pgc->gcslparaxijk_single(p,ph,var);	
    pgc->gcslparaxijk_single(p,sh,var);	
}

void sflow_bicgstab::matvec_axb(lexer *p, matrix2D &M, double *x, double *y)
{
    n=0;
	SLICEFLEXLOOP
	{
	y[IJ]  = rhs[IJ]

			-(M.p[n]*x[IJ]
			+ M.n[n]*x[Ip1J] 
			+ M.s[n]*x[Im1J]
			+ M.w[n]*x[IJp1]
			+ M.e[n]*x[IJm1]);
    ++n;
	}
}

void sflow_bicgstab::matvec_std(lexer *p, matrix2D &M, double *x, double *y)
{
    n=0;
	SLICEFLEXLOOP
	{
	y[IJ]      = M.p[n]*x[IJ]
				+ M.n[n]*x[Ip1J] 
				+ M.s[n]*x[Im1J]
				+ M.w[n]*x[IJp1]
				+ M.e[n]*x[IJm1];
    ++n;
	}
}

double sflow_bicgstab::res_calc(lexer *p, matrix2D &M, ghostcell *pgc, double *x)
{
	double y;
	double resi=0.0;
    
    n=0;
	SLICEFLEXLOOP
	{	
	y  = rhs[IJ]

		-(M.p[n]*x[IJ]
		+ M.n[n]*x[Ip1J] 
		+ M.s[n]*x[Im1J]
		+ M.w[n]*x[IJp1]
		+ M.e[n]*x[IJm1]);

	resi+=y*y;
    
    ++n;
	}

	resi=sqrt(pgc->globalsum(resi));

	return resi/double(p->cellnumtot2D);	
}

void sflow_bicgstab::precon_setup(lexer* p, matrix2D &M, ghostcell* pgc)
{
    n=0;
	SLICEFLEXLOOP
    {
	aii[IJ]=-1.0/(M.p[n]+epsi);
    ++n;
    }
}

void sflow_bicgstab::precon_solve(lexer* p, ghostcell* pgc, double *f, double *b)
{
	SLICEFLEXLOOP
	f[IJ]=b[IJ]*aii[IJ];
}

void sflow_bicgstab::fillxvec(lexer* p, slice& f, vec2D &rhsvec)
{
    n=0;
	SLICEFLEXLOOP
	{
	x[IJ] = f(i,j);
    
    rhs[IJ] = rhsvec.V[n];

    ++n;
    }
}


void sflow_bicgstab::finalize(lexer *p, slice &f)
{  
        SLICEFLEXLOOP
        f(i,j)=x[IJ];
}
