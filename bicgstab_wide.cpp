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

#include"bicgstab_wide.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver_void.h"
#include"jacobi_scaling.h"
#include"jacobi_block.h"
#include"sip.h"

bicgstab_wide::bicgstab_wide(lexer* p, fdm *a, ghostcell *pgc, int precon_select):sj(p),rj(p),r0(p),vj(p),tj(p),pj(p),precoeff(p),
												ph(p),sh(p),epsi(1e-19)
{
	precon = new jacobi_scaling(p,a,pgc);

	margin=3;
}

bicgstab_wide::~bicgstab_wide()
{
}

void bicgstab_wide::setup(lexer* p,fdm* a, ghostcell* pgc, int var, cpt &C)
{
}

void bicgstab_wide::start(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& xvec, vec& rhsvec, int var, int gcv, double stop_crit)
{
	
}

void bicgstab_wide::startF(lexer* p, fdm_fnpf* c, ghostcell* pgc, double *f, vec& rhsvec, matrix_diag &M, int var, int gcv, double stop_crit)
{
   /* p->preconiter=0;
	
	fillxvec4(p,a,f);
    sizeM=p->sizeM4;
	solve(p,a,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,a->C4);
	
	finalize(p,a,f,xvec,var);*/
}

void bicgstab_wide::solveF(lexer* p, fdm_fnpf* c, ghostcell* pgc, double *f, vec& rhsvec, matrix_diag &M, int var, int gcv, double stop_crit)
{
   /* solveriter=0;
	residual = 1.0e9;

	// -----------------
	precon->setup(p,a,pgc,var,C);
	// -----------------

 restart:
    r_j=norm_r0=0.0;
	gcupdate(p,a,pgc,xvec,var,gcv,solveriter);		
	pgc->gcparaxvec(p,xvec,var);
	
	matvec_axb(p,a,xvec,rj,C);
	
	NLOOP
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
		precon->solve(p,a,pgc,ph,pj,var,240,p->preconiter,p->N13,p->N18,C);
		pgc->gcparaxvec(p,ph,var);		
		gcupdate(p,a,pgc,ph,var,240,solveriter);		
		// -------------------------
		
		matvec_std(p,a,ph,vj,C);
		
		NLOOP
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
		residual=res_calc(p,a,xvec,pgc,C);
		++solveriter;

		goto restart;
	}

    if((fabs(alpha)*norm_vj/(norm_rj==0?1.0e-15:norm_rj))<=0.08)
	{
		residual=res_calc(p,a,xvec,pgc,C);
		++solveriter;

		goto restart;
	}

		norm_sj=0.0;
		
		NLOOP
		{
		sj.V[n] = rj.V[n] - alpha*vj.V[n];
		norm_sj += sj.V[n]*sj.V[n];
		}

	    norm_sj=sqrt(pgc->globalsum(norm_sj));

    if(norm_sj>stop_crit)
	{
		// -------------------------
		precon->solve(p,a,pgc,sh,sj,var,240,p->preconiter,p->N13,p->N18,C);
        pgc->gcparaxvec(p,sh,var);		
		gcupdate(p,a,pgc,sh,var,240,solveriter);
		// -------------------------

		matvec_std(p,a,sh,tj,C);
		
		w1=w2=0.0;
		
		NLOOP
		{
		    w1 += tj.V[n]*sj.V[n];
		    w2 += tj.V[n]*tj.V[n];
		}

		w1=pgc->globalsum(w1);
		w2=pgc->globalsum(w2);

		w=w1/(w2==0?1.0e-15:w2);

		r_j1=0.0;
		
		NLOOP
		{
		xvec.V[n] += alpha*ph.V[n] + w*sh.V[n];
		rj.V[n]  = sj.V[n]-w*tj.V[n];
		r_j1 += rj.V[n]*r0.V[n];
		}

		r_j1=pgc->globalsum(r_j1);

		beta=alpha*r_j1/(w*r_j==0?1.0e-15:(w*r_j));
		
		NLOOP
		pj.V[n] = rj.V[n] + beta*(pj.V[n]-w*vj.V[n]);
	}


	if(norm_sj<=stop_crit)
	{
	r_j1=0.0;
		
		NLOOP
		{
		xvec.V[n] += alpha*ph.V[n];
		rj.V[n]=sj.V[n];
		r_j1 += rj.V[n]*r0.V[n];
		}

    r_j1=pgc->globalsum(r_j1);
	}

	    r_j = r_j1 ;

	    residual=0.0;
		
		NLOOP
		residual += rj.V[n]*rj.V[n];

	    residual = sqrt(pgc->globalsum(residual))/double(p->cellnumtot);
		
	    ++solveriter;

	}while((residual>=stop_crit) && (solveriter<maxiter));

    } 
		
	
	NLOOP4
	{
	ph.V[n]=0.0;
	sh.V[n]=0.0;
	}
	pgc->start4V(p,ph,240);
	pgc->start4V(p,sh,240);*/
    
}
	
void bicgstab_wide::solve(lexer* p,fdm* a, ghostcell* pgc, vec& xvec, vec& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt &C)
{
}

void bicgstab_wide::matvec_axb(lexer *p,fdm* a, matrix_diag &M, vec &x, vec &y, cpt &C)
{
	NLOOP
	{
	y.V[n]  = a->rhsvec.V[n]

			-(M.p[n]*x.V[I_J_K]
			+ M.n[n]*x.V[Ip1_J_K] 
			+ M.s[n]*x.V[Im1_J_K]
			+ M.w[n]*x.V[I_Jp1_K]
			+ M.e[n]*x.V[I_Jm1_K]
			+ M.t[n]*x.V[I_J_Kp1]
			+ M.b[n]*x.V[I_J_Km1]
            
            + M.nn[n]*x.V[Ip2_J_K] 
            + M.ss[n]*x.V[Im2_J_K]
            + M.ww[n]*x.V[I_Jp2_K]
            + M.ee[n]*x.V[I_Jm2_K]
            + M.tt[n]*x.V[I_J_Kp2]
            + M.bb[n]*x.V[I_J_Km2]);
	}
}

void bicgstab_wide::matvec_std(lexer *p,fdm* a, matrix_diag &M, vec &x, vec &y, cpt &C)
{
	NLOOP
	{
	y.V[n]      = M.p[n]*x.V[I_J_K]
				+ M.n[n]*x.V[Ip1_J_K] 
				+ M.s[n]*x.V[Im1_J_K]
				+ M.w[n]*x.V[I_Jp1_K]
				+ M.e[n]*x.V[I_Jm1_K]
				+ M.t[n]*x.V[I_J_Kp1]
				+ M.b[n]*x.V[I_J_Km1]
                
                + M.nn[n]*x.V[Ip2_J_K] 
                + M.ss[n]*x.V[Im2_J_K]
                + M.ww[n]*x.V[I_Jp2_K]
                + M.ee[n]*x.V[I_Jm2_K]
                + M.tt[n]*x.V[I_J_Kp2]
                + M.bb[n]*x.V[I_J_Km2];
	}
}

double bicgstab_wide::res_calc(lexer *p, fdm *a, matrix_diag &M, vec &x, ghostcell *pgc, cpt &C)
{
	double y;
	double resi=0.0;

	NLOOP
	{	
	y  = a->rhsvec.V[n]

		-(M.p[n]*x.V[n]
		+ M.n[n]*x.V[Ip1_J_K] 
		+ M.s[n]*x.V[Im1_J_K]
		+ M.w[n]*x.V[I_Jp1_K]
		+ M.e[n]*x.V[I_Jm1_K]
		+ M.t[n]*x.V[I_J_Kp1]
		+ M.b[n]*x.V[I_J_Km1]
        
        + M.nn[n]*x.V[Ip2_J_K] 
        + M.ss[n]*x.V[Im2_J_K]
        + M.ww[n]*x.V[I_Jp2_K]
        + M.ee[n]*x.V[I_Jm2_K]
        + M.tt[n]*x.V[I_J_Kp2]
        + M.bb[n]*x.V[I_J_Km2]);

	resi+=y*y;
	}

	resi=sqrt(pgc->globalsum(resi));

	return resi/double(p->cellnumtot);	
}


void bicgstab_wide::gcpara_update(lexer* p, vec &x, ghostcell* pgc)
{
}

void bicgstab_wide::gcupdate(lexer *p, fdm *a, ghostcell *pgc, vec &x, int var, int gcv, int solveriter)
{
	if(p->N15==3 && gcv>=40 && gcv<=45)
	pgc->start4V(p,x,241);
	
	if(gcv>=40 && gcv<=45 && (p->N15==1||p->N15==2) && p->count<=p->N16)
	pgc->start4V(p,x,240);
	
	if(p->N15==4 && gcv>=40 && gcv<=45)
	pgc->start4V(p,x,240);
	
	if(p->N15==5 && gcv>=40 && gcv<=45 && solveriter<5) 
	pgc->start4V(p,x,240);
	
	if(gcv==49) 
	pgc->start4V(p,x,49);
}

void bicgstab_wide::fillxvec1(lexer* p, fdm* a, field& f)
{
}

void bicgstab_wide::fillxvec2(lexer* p, fdm* a, field& f)
{
}

void bicgstab_wide::fillxvec3(lexer* p, fdm* a, field& f)
{
}

void bicgstab_wide::fillxvec4(lexer* p, fdm* a, field& f)
{
	int count,q;
	
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

void bicgstab_wide::finalize(lexer *p, fdm *a, field &f, vec &xvec, int var)
{
	if(var==1)
    {
        count=0;
        ULOOP
        {
        f(i,j,k)=xvec.V[count];
        ++count;
        }
    }
	
	if(var==2)
    {
        count=0;
        VLOOP
        {
        f(i,j,k)=xvec.V[count];
        ++count;
        }
    }
	
	if(var==3)
    {
        count=0;
        WLOOP
        {
        f(i,j,k)=xvec.V[count];
        ++count;
        }
    }
	
	if(var==4 || var==5)
    {
        count=0;
        LOOP
        {
        f(i,j,k)=xvec.V[count];
        ++count;
        }
    }
}

