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

#include"hypre_aij.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"


hypre_aij::hypre_aij(lexer* p,fdm* a,ghostcell *pgc) : xvec(p)
{	
    int vecsize=p->knox*p->knoy*p->knoz;
    
    if(p->A10==3)
    vecsize=p->knox*p->knoy*(p->knoz+1); 
    
    p->Iarray(rows,vecsize);
}

hypre_aij::~hypre_aij()
{
}

void hypre_aij::start(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& rhsvec, int var)
{
    make_grid(p,pgc);
    create_solvers(p,pgc);
    
    if(var<=5)
	fill_matrix_7p(p,a,pgc,f);
    
    if(var==6)
	fill_matrix_13p(p,a,pgc,f);
    
    if(var==7)
	fill_matrix_19p(p,a,pgc,f);
  

    if(p->N10==21)
    {
	HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
	HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

	HYPRE_PCGGetNumIterations(solver, &num_iterations);
	HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==22)
    {
	HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_b, par_x);
	HYPRE_ParCSRGMRESSolve(solver, parcsr_A, par_b, par_x);

	HYPRE_GMRESGetNumIterations(solver, &num_iterations);
	HYPRE_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==23)
    {
	HYPRE_ParCSRLGMRESSetup(solver, parcsr_A, par_b, par_x);
	HYPRE_ParCSRLGMRESSolve(solver, parcsr_A, par_b, par_x);

	HYPRE_LGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_LGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==24)
    {
	HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_A, par_b, par_x);
	HYPRE_ParCSRBiCGSTABSolve(solver, parcsr_A, par_b, par_x);

	HYPRE_BiCGSTABGetNumIterations(solver, &num_iterations);
	HYPRE_BiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==25)
    {
    HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
    HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);
    
    HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
	} 
    
	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
	
	fillbackvec(p,a,f,xvec,var);
    
    delete_solvers(p,pgc);
    delete_grid(p,pgc);
}

void hypre_aij::startf(lexer* p, ghostcell* pgc, field &f, vec& rhs, matrix_diag &M, int var)
{
    
}

void hypre_aij::startM(lexer* p, ghostcell* pgc, double *x, double *rhs, double *M, int var)
{
}

void hypre_aij::startV(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
{
    
}

void hypre_aij::solve(lexer* p,fdm* a, ghostcell* pgc, vec& xvec, vec& rhsvec, int var, int gcv, int &solveriter)
{
	
	numiter=0;
	p->solveriter=0;

	p->solveriter=numiter;		
}

void hypre_aij::fillbackvec(lexer *p, fdm *a, field &f, vec &xvec, int var)
{
	HYPRE_IJVectorGetValues(x, p->N4_row, rows, xvec.V);
	
        n=0;
        FLUIDLOOP
        {
        f(i,j,k)=xvec.V[n];
        ++n;
        }
}

void hypre_aij::fillxvec1(lexer* p, fdm* a, field& f)
{
}
	
void hypre_aij::fillxvec2(lexer* p, fdm* a, field& f)
{
}
	
void hypre_aij::fillxvec3(lexer* p, fdm* a, field& f)
{
}
	
void hypre_aij::fillxvec4(lexer* p, fdm* a, field& f)
{
}

#endif
