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
#include"grid.h"

bicgstab_wide::bicgstab_wide(lexer* p):xvec(p),sj(p),rj(p),r0(p),vj(p),tj(p),pj(p),precoeff(p),
												ph(p),sh(p),aii(p),epsi(1e-19)
{
	margin=3;
    
    count=0;
	LOOP
	++count;
	
	p->sizeM4[0]=0;
	p->sizeM4[1]=count;
    
    pgrid = new grid(p);
    
    C4.allocate(p);
    
    pgrid->column_pt4_update(p,C4);
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
	
	fillxvec(p,f,xvec);
    solveF(p,c,pgc,xvec,rhsvec,M,var,gcv,stop_crit);
	finalize(p,f,xvec);
}

void bicgstab_wide::solveF(lexer* p, fdm_fnpf* c, ghostcell* pgc, vec& x, vec& rhsvec, matrix_diag &M, int var, int gcv, double stop_crit)
{
    var=4;
    gcv=250;
    maxiter=p->A322;
    p->solveriter=0;
	residual = 1.0e9;

	// -----------------
	precon_setup(p,M,C4);
	// -----------------

 restart:
    r_j=norm_r0=0.0;	
    pgc->gcparaxvec_sr(p, xvec, C4, 4);
	
	matvec_axb(p,M,xvec,rj,rhsvec,C4);
	
	NLOOP4
	{
		r0.V[n]=pj.V[n]=rj.V[n];
		r_j += rj.V[n]*r0.V[n];
    }

    r_j=pgc->globalsum(r_j);
    norm_r0=sqrt(r_j);

    if((residual>=stop_crit) && (p->solveriter<maxiter))
	{

	do{
	    sigma=0.0;
	    norm_vj=0.0;
	    norm_rj=0.0;
		
		// -------------------------
		precon(p,ph,pj,C4);	
        pgc->gcparaxvec_sr(p, ph, C4, 4);
		// -------------------------
		
		matvec_std(p,M,ph,vj,C4);
		
		NLOOP4
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
		residual=res_calc(p,pgc,M,xvec,rhsvec,C4);
		++p->solveriter;

		goto restart;
	}

    if((fabs(alpha)*norm_vj/(norm_rj==0?1.0e-15:norm_rj))<=0.08)
	{
		residual=res_calc(p,pgc,M,xvec,rhsvec,C4);
		++p->solveriter;

		goto restart;
	}

		norm_sj=0.0;
		
		NLOOP4
		{
		sj.V[n] = rj.V[n] - alpha*vj.V[n];
		norm_sj += sj.V[n]*sj.V[n];
		}

	    norm_sj=sqrt(pgc->globalsum(norm_sj));

    if(norm_sj>stop_crit)
	{
		// -------------------------
		precon(p,sh,sj,C4);
        pgc->gcparaxvec_sr(p, sh, C4, 4);	
		// -------------------------

		matvec_std(p,M,sh,tj,C4);
		
		w1=w2=0.0;
		
		NLOOP4
		{
		    w1 += tj.V[n]*sj.V[n];
		    w2 += tj.V[n]*tj.V[n];
		}

		w1=pgc->globalsum(w1);
		w2=pgc->globalsum(w2);

		w=w1/(w2==0?1.0e-15:w2);

		r_j1=0.0;
		
		NLOOP4
		{
		xvec.V[n] += alpha*ph.V[n] + w*sh.V[n];
		rj.V[n]  = sj.V[n]-w*tj.V[n];
		r_j1 += rj.V[n]*r0.V[n];
		}

		r_j1=pgc->globalsum(r_j1);

		beta=alpha*r_j1/(w*r_j==0?1.0e-15:(w*r_j));
		
		NLOOP4
		pj.V[n] = rj.V[n] + beta*(pj.V[n]-w*vj.V[n]);
	}


	if(norm_sj<=stop_crit)
	{
	r_j1=0.0;
		
		NLOOP4
		{
		xvec.V[n] += alpha*ph.V[n];
		rj.V[n]=sj.V[n];
		r_j1 += rj.V[n]*r0.V[n];
		}

    r_j1=pgc->globalsum(r_j1);
	}

	    r_j = r_j1 ;

	    residual=0.0;
		
		NLOOP4
		residual += rj.V[n]*rj.V[n];

	    residual = sqrt(pgc->globalsum(residual))/double(p->cellnumtot);
		
	    ++p->solveriter;

	}while((residual>=stop_crit) && (p->solveriter<maxiter));

    } 
    
    p->final_res = residual;
	
	VECLOOP
	{
	ph.V[n]=0.0;
	sh.V[n]=0.0;
	}
}
	
void bicgstab_wide::solve(lexer* p,fdm* a, ghostcell* pgc, vec& xvec, vec& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt &C)
{
}

void bicgstab_wide::matvec_axb(lexer *p, matrix_diag &M, vec &x, vec &y, vec &rhs, cpt &C)
{
	NLOOP4
	{
	y.V[n]  = rhs.V[n]

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

void bicgstab_wide::matvec_std(lexer *p, matrix_diag &M, vec &x, vec &y, cpt &C)
{
	NLOOP4
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

void bicgstab_wide::precon(lexer *p, vec &x, vec &y, cpt &C)
{
    NLOOP4
	x.V[n]=y.V[n]*aii.V[n];
    
}

void bicgstab_wide::precon_setup(lexer *p, matrix_diag &M, cpt &C)
{
    
    NLOOP4
	aii.V[n]=-1.0/(M.p[n]+epsi);
    
}

double bicgstab_wide::res_calc(lexer *p, ghostcell *pgc, matrix_diag &M, vec &x, vec &rhs, cpt &C)
{
	double y;
	double resi=0.0;

	NLOOP4
	{	
	y  = rhs.V[n]

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


void bicgstab_wide::fillxvec(lexer* p, double *f, vec &x)
{
	int count,q;
	
	count=0;
    LOOP
    {
    x.V[count]=f[FIJK];
    ++count;
    }

    GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

        if(p->gcb4[n][3]==1)
        {
        x.V[count]=f[FIm1JK];
        ++count;
        
        x.V[count]=f[FIm2JK];
        ++count;
        
        x.V[count]=f[FIm3JK];
        ++count;
        }

        if(p->gcb4[n][3]==2)
        {
        x.V[count]=f[FIJp1K];
        ++count;
        
        x.V[count]=f[FIJp2K];
        ++count;
        
        x.V[count]=f[FIJp3K];
        ++count;
        }

        if(p->gcb4[n][3]==3)
        {
        x.V[count]=f[FIJm1K];
        ++count;
        
        x.V[count]=f[FIJm2K];
        ++count;
        
        x.V[count]=f[FIJm3K];
        ++count;
        }

        if(p->gcb4[n][3]==4)
        {
        x.V[count]=f[FIp1JK];
        ++count;
        
        x.V[count]=f[FIp2JK];
        ++count;
        
        x.V[count]=f[FIp3JK];
        ++count;
        }

        if(p->gcb4[n][3]==5)
        {
        x.V[count]=f[FIJKm1];
        ++count;
        
        x.V[count]=f[FIJKm2];
        ++count;
        
        x.V[count]=f[FIJKm3];
        ++count;
        }

        if(p->gcb4[n][3]==6)
        {
        x.V[count]=f[FIJKp1];
        ++count;
        
        x.V[count]=f[FIJKp2];
        ++count;
        
        x.V[count]=f[FIJKp3];
        ++count;
        }
    }
}

void bicgstab_wide::finalize(lexer *p, double *f, vec &x)
{
    count=0;
    LOOP
    {
    f[FIJK]=x.V[count];
    ++count;
    }
}

