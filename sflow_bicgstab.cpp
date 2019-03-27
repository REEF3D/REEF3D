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

#include"sflow_bicgstab.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

sflow_bicgstab::sflow_bicgstab(lexer* p,fdm2D* b,ghostcell *pgc):sj(p),rj(p),r0(p),vj(p),tj(p),pj(p),precoeff(p),
												ph(p),sh(p),epsi(1e-19)
{	
	
	margin=3;
}

sflow_bicgstab::~sflow_bicgstab()
{
}

void sflow_bicgstab::setup(lexer* p,fdm2D* b, ghostcell* pgc, int var, cpt2D &C)
{
}

void sflow_bicgstab::start(lexer* p, fdm2D* b, ghostcell* pgc, slice &f, vec2D& xvec, vec2D& rhsvec, int var, int gcv, double stop_crit)
{
	p->preconiter=0;
	
	if(var==1)
    {
	fillxvec1(p,b,f);
    sizeS=p->sizeS1;
	precon=precon123;
	solve(p,b,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,b->C1);
    }
	
	if(var==2)
    {
	fillxvec2(p,b,f);
    sizeS=p->sizeS2;
	precon=precon123;
	solve(p,b,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,b->C2);
    }
	
	if(var==3||var==4||var==5)
    {
	fillxvec4(p,b,f);
    sizeS=p->sizeS4;
	precon=precon4;
	solve(p,b,pgc,xvec,rhsvec,var,gcv,p->solveriter,p->N46,stop_crit,b->C4);
    }
	
	finalize(p,b,f,xvec,var);
}

	
void sflow_bicgstab::solve(lexer* p,fdm2D* b, ghostcell* pgc, vec2D& xvec, vec2D& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt2D &C)
{
	solveriter=0;
	residual = 1.0e9;

	// -----------------
	precon->setup(p,b,pgc,var,C);
	// -----------------

 restart:
    r_j=norm_r0=0.0;	
	pgc->gcparaxvec(p,xvec,var);
	
	matvec_axb(p,b,xvec,rj,C);
	
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
		precon->solve(p,b,pgc,ph,pj,var,240,p->preconiter,p->N13,p->N18,C);
		pgc->gcparaxvec(p,ph,var);		
		// -------------------------
		
		matvec_std(p,b,ph,vj,C);
		
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
		residual=res_calc(p,b,xvec,pgc,C);
		++solveriter;

		goto restart;
	}

    if((fabs(alpha)*norm_vj/(norm_rj==0?1.0e-15:norm_rj))<=0.08)
	{
		residual=res_calc(p,b,xvec,pgc,C);
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
		precon->solve(p,b,pgc,sh,sj,var,240,p->preconiter,p->N13,p->N18,C);
        pgc->gcparaxvec(p,sh,var);		
		// -------------------------

		matvec_std(p,b,sh,tj,C);
		
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

	    residual = sqrt(pgc->globalsum(residual))/double(p->cellnumtot);
		
	    ++solveriter;

	}while((residual>=stop_crit) && (solveriter<maxiter));

    } 
		
	
	NSLICELOOP4
	{
	ph.V[n]=0.0;
	sh.V[n]=0.0;
	}
	pgc->start4V(p,ph,240);
	pgc->start4V(p,sh,240);
}

void sflow_bicgstab::matvec_axb(lexer *p,fdm2D* b, vec &x, vec &y, cpt &C)
{
	NSLICELOOP
	{
	y.V[n]  = b->rhsvec.V[n]

			-(b->M.p[n]*x.V[I_J_K]
			+ b->M.n[n]*x.V[Ip1_J_K] 
			+ b->M.s[n]*x.V[Im1_J_K]
			+ b->M.w[n]*x.V[I_Jp1_K]
			+ b->M.e[n]*x.V[I_Jm1_K]);
	}
}

void sflow_bicgstab::matvec_std(lexer *p,fdm2D* b, vec &x, vec &y, cpt &C)
{
	NSLICELOOP
	{
	y.V[n]      = b->M.p[n]*x.V[I_J_K]
				+ b->M.n[n]*x.V[Ip1_J_K] 
				+ b->M.s[n]*x.V[Im1_J_K]
				+ b->M.w[n]*x.V[I_Jp1_K]
				+ b->M.e[n]*x.V[I_Jm1_K];
	}
}

double sflow_bicgstab::res_calc(lexer *p, fdm *a, vec &x, ghostcell *pgc, cpt &C)
{
	double y;
	double resi=0.0;

	NSLICELOOP
	{	
	y  = b->rhsvec.V[n]

		-(b->M.p[n]*x.V[n]
		+ b->M.n[n]*x.V[Ip1_J_K] 
		+ b->M.s[n]*x.V[Im1_J_K]
		+ b->M.w[n]*x.V[I_Jp1_K]
		+ b->M.e[n]*x.V[I_Jm1_K]);

	resi+=y*y;
	}

	resi=sqrt(pgc->globalsum(resi));

	return resi/double(p->cellnumtot);	
}


void sflow_bicgstab::gcpara_update(lexer* p, vec &x, ghostcell* pgc)
{
}

void sflow_bicgstab::fillxvec1(lexer* p, fdm2D* b, slice& f)
{
	int count,q;
	
	count=0;
    SLICELOOP1
    {
    b->xvec.V[count]=f(i,j);
    ++count;
    }

    GC1LOOP
    {
    i=p->gcb1[n][0];
    j=p->gcb1[n][1];

        if(p->gcb1[n][3]==1)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i-1-q,j);
        ++count;
        }

        if(p->gcb1[n][3]==2)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i,j+1+q);
        ++count;
        }

        if(p->gcb1[n][3]==3)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i,j-1-q);
        ++count;
        }

        if(p->gcb1[n][3]==4)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i+1+q,j);
        ++count;
        }
    }
}

void sflow_bicgstab::fillxvec2(lexer* p, fdm2D* b, slice& f)
{
	int count,q;
	
	count=0;
    SLICELOOP2
    {
    b->xvec.V[count]=f(i,j);
    ++count;
    }

    GC2LOOP
    {
    i=p->gcb2[n][0];
    j=p->gcb2[n][1];

        if(p->gcb2[n][3]==1)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i-1-q,j);
        ++count;
        }

        if(p->gcb2[n][3]==2)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i,j+1+q);
        ++count;
        }

        if(p->gcb2[n][3]==3)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i,j-1-q);
        ++count;
        }

        if(p->gcb2[n][3]==4)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i+1+q,j);
        ++count;
        }
    }
}

void sflow_bicgstab::fillxvec4(lexer* p, fdm2D* b, slice& f)
{
	int count,q;
	
	count=0;
    SLICELOOP4
    {
    b->xvec.V[count]=f(i,j);
    ++count;
    }

    GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];

        if(p->gcb4[n][3]==1)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i-1-q,j);
        ++count;
        }

        if(p->gcb4[n][3]==2)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i,j+1+q);
        ++count;
        }

        if(p->gcb4[n][3]==3)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i,j-1-q);
        ++count;
        }

        if(p->gcb4[n][3]==4)
        for(q=0;q<margin;++q)
        {
        b->xvec.V[count]=f(i+1+q,j);
        ++count;
        }
    }
}


void sflow_bicgstab::finalize(lexer *p, fdm2D *b, slice &f, vec2D &xvec, int var)
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

void hypre_struct2D::precon_setup(lexer* p,fdm2D* b, ghostcell* pgc, int var, cpt2D &C)
{
}

void hypre_struct2D::precon_solve(lexer* p,fdm2D* b, ghostcell* pgc, vec2D& xvec, vec2D& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt2D &C)
{
}
