/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
#include"bicgstab_ijk.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

bicgstab_ijk::bicgstab_ijk(lexer* p, fdm *a, ghostcell *pgc):epsi(1e-19)
{
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

bicgstab_ijk::~bicgstab_ijk()
{
}

void bicgstab_ijk::setup(lexer* p, ghostcell* pgc, int var)
{
}

void bicgstab_ijk::start(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& rhsvec, int var)
{
	p->preconiter=0;
    
	if(var==1)
    {
    flag = p->flag1;
    ulast=p->ulast;
    vlast=0;
    wlast=0;
    stop_crit=p->N43;
    }
	
	if(var==2)
    {
    flag = p->flag2;
    ulast=0;
    vlast=p->vlast;
    wlast=0;
    stop_crit=p->N43;
    }
	
	if(var==3)
    {
    flag = p->flag3;
    ulast=0;
    vlast=0;
    wlast=p->wlast;
    stop_crit=p->N43;
    }
	
	if(var==4||var==5)
    {
    flag = p->flag4;
    ulast=0;
    vlast=0;
    wlast=0;
    stop_crit=p->N44;
    }
    
    fillxvec(p,a,f,rhsvec);
	solve(p,pgc,rhsvec,a->M,var,p->solveriter,p->N46,stop_crit);
	
	finalize(p,a,f);
}

void bicgstab_ijk::startV(lexer* p, ghostcell* pgc, double *f, vec& rhsvec, matrix_diag &M, int var)
{
    p->preconiter=0;
    
    flag = p->flag4;
    ulast=0;
    vlast=0;
    wlast=0;
    stop_crit=p->N43;
    
    fillxvecV(p,f,rhsvec);
	solve(p,pgc,rhsvec,M,var,p->solveriter,p->N46,stop_crit);
	
	finalizeV(p,f);
}

void bicgstab_ijk::startf(lexer* p, ghostcell* pgc, field &f, vec& rhs, matrix_diag &M, int var)
{
}

void bicgstab_ijk::startM(lexer* p, ghostcell* pgc, double *x, double *rhs, double *M, int var)
{
}

void bicgstab_ijk::startF(lexer* p, ghostcell* pgc, double *f, vec& rhsvec, matrix_diag &M, int var)
{
}
	
void bicgstab_ijk::solve(lexer* p, ghostcell* pgc, vec& rhsvec, matrix_diag &M, int var, int &solveriter, int maxiter, double stop_crit)
{
	solveriter=0;
	residual = 1.0e9;

	// -----------------
	precon_setup(p,pgc,M);
	// -----------------

 restart:
    r_j=norm_r0=0.0;	
	pgc->gcparaxijk_single(p,x,var);
	
	matvec_axb(p,x,rj,M);
	
	FLEXLOOP
	{
		r0[IJK]=pj[IJK]=rj[IJK];
		r_j += rj[IJK]*r0[IJK];
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
		precon_solve(p,pgc,ph,pj,M);
		pgc->gcparaxijk_single(p,ph,var);				
		// -------------------------
		
		matvec_std(p,ph,vj,M);
		
		FLEXLOOP
		{
			sigma   += vj[IJK]*r0[IJK];
			norm_vj += vj[IJK]*vj[IJK];
			norm_rj += rj[IJK]*rj[IJK];
	    }
		
        sigma = pgc->globalsum(sigma);
		norm_vj = sqrt(pgc->globalsum(norm_vj));
		norm_rj = sqrt(pgc->globalsum(norm_rj));

	    alpha=r_j/sigma;

	if(fabs(sigma) <= (1.0e-12*(norm_vj*norm_r0)))
	{	
		residual=res_calc(p,pgc,x,M);
		++solveriter;

		goto restart;
	}

    if((fabs(alpha)*norm_vj/(norm_rj==0?1.0e-15:norm_rj))<=0.08)
	{
		residual=res_calc(p,pgc,x,M);
		++solveriter;

		goto restart;
	}

		norm_sj=0.0;
		
		FLEXLOOP
		{
		sj[IJK] = rj[IJK] - alpha*vj[IJK];
		norm_sj += sj[IJK]*sj[IJK];
		}

	    norm_sj=sqrt(pgc->globalsum(norm_sj));

    if(norm_sj>stop_crit)
	{
        // -------------------------
        precon_solve(p,pgc,sh,sj,M);
        pgc->gcparaxijk_single(p,sh,var);		
        // -------------------------

		matvec_std(p,sh,tj,M);
		
		w1=w2=0.0;
		
		FLEXLOOP
		{
		    w1 += tj[IJK]*sj[IJK];
		    w2 += tj[IJK]*tj[IJK];
		}

		w1=pgc->globalsum(w1);
		w2=pgc->globalsum(w2);

		w=w1/(w2==0?1.0e-15:w2);

		r_j1=0.0;
		
		FLEXLOOP
		{
		x[IJK] += alpha*ph[IJK] + w*sh[IJK];
		rj[IJK]  = sj[IJK]-w*tj[IJK];
		r_j1 += rj[IJK]*r0[IJK];
		}

		r_j1=pgc->globalsum(r_j1);

		beta=alpha*r_j1/(w*r_j==0?1.0e-15:(w*r_j));
		
		FLEXLOOP
		pj[IJK] = rj[IJK] + beta*(pj[IJK]-w*vj[IJK]);
	}


	if(norm_sj<=stop_crit)
	{
	r_j1=0.0;
		
		FLEXLOOP
		{
		x[IJK] += alpha*ph[IJK];
		rj[IJK]=sj[IJK];
		r_j1 += rj[IJK]*r0[IJK];
		}

    r_j1=pgc->globalsum(r_j1);
	}

	    r_j = r_j1 ;

	    residual=0.0;
		
		FLEXLOOP
		residual += rj[IJK]*rj[IJK];

	    residual = sqrt(pgc->globalsum(residual))/double(p->cellnumtot);
		
	    ++solveriter;

	}while((residual>=stop_crit) && (solveriter<maxiter));

    } 
		
	
	LOOP
	{
	ph[IJK]=0.0;
	sh[IJK]=0.0;
	}
}

void bicgstab_ijk::matvec_axb(lexer *p, double *x, double *y, matrix_diag &M)
{
    n=0;
	FLEXLOOP
	{
	y[IJK]  = rhs[IJK]

			-(M.p[n]*x[IJK]
			+ M.n[n]*x[Ip1JK] 
			+ M.s[n]*x[Im1JK]
			+ M.w[n]*x[IJp1K]
			+ M.e[n]*x[IJm1K]
			+ M.t[n]*x[IJKp1]
			+ M.b[n]*x[IJKm1]);
    ++n;
	}
}

void bicgstab_ijk::matvec_std(lexer *p, double *x, double *y, matrix_diag &M)
{
    n=0;
	FLEXLOOP
	{
	y[IJK]      = M.p[n]*x[IJK]
				+ M.n[n]*x[Ip1JK] 
				+ M.s[n]*x[Im1JK]
				+ M.w[n]*x[IJp1K]
				+ M.e[n]*x[IJm1K]
				+ M.t[n]*x[IJKp1]
				+ M.b[n]*x[IJKm1];
    ++n;
	}
}

double bicgstab_ijk::res_calc(lexer *p, ghostcell *pgc, double *x, matrix_diag &M)
{
	double y;
	double resi=0.0;
    
    n=0;
	FLEXLOOP
	{	
	y  = rhs[IJK]

		-(M.p[n]*x[IJK]
		+ M.n[n]*x[Ip1JK] 
		+ M.s[n]*x[Im1JK]
		+ M.w[n]*x[IJp1K]
		+ M.e[n]*x[IJm1K]
		+ M.t[n]*x[IJKp1]
		+ M.b[n]*x[IJKm1]);

	resi+=y*y;
    
    ++n;
	}

	resi=sqrt(pgc->globalsum(resi));

	return resi/double(p->cellnumtot);	
}

void bicgstab_ijk::precon_setup(lexer* p, ghostcell* pgc, matrix_diag &M)
{
    n=0;
	FLEXLOOP
    {
	aii[IJK]=-1.0/(M.p[n]+epsi);
    ++n;
    }
}

void bicgstab_ijk::precon_solve(lexer* p, ghostcell* pgc, double *f, double *b, matrix_diag &M)
{
	FLEXLOOP
	f[IJK]=b[IJK]*aii[IJK];
}

void bicgstab_ijk::fillxvec(lexer* p, fdm* a, field& f, vec &rhsvec)
{
    n=0;
	FLEXLOOP
	{
	x[IJK] = f(i,j,k);
    
    rhs[IJK] = rhsvec.V[n];

    ++n;
    }
}

void bicgstab_ijk::finalize(lexer *p, fdm *a, field &f)
{  
    FLEXLOOP
    f(i,j,k)=x[IJK];
}

void bicgstab_ijk::fillxvecV(lexer* p, double *f, vec &rhsvec)
{
    n=0;
	FLEXLOOP
	{
	x[IJK] = f[IJK];
    
    rhs[IJK] = rhsvec.V[n];

    ++n;
    }
}

void bicgstab_ijk::finalizeV(lexer *p, double *f)
{  
    FLEXLOOP
    f[IJK]=x[IJK];
}

