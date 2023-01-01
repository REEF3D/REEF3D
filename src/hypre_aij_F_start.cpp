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

#include"hypre_aij.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"

void hypre_aij::startF(lexer* p, ghostcell* pgc, double *f, vec& rhsvec, matrix_diag &M, int var)
{     
    double *xvec;
    
    p->Darray(xvec,p->knox*p->knoy*(p->knoz+1));
    
    if(var==7||var==9||var==11)
    make_grid_F(p,pgc);
    
    if(var==8||var==10||var==12)
    make_grid(p,pgc);
    
    create_solvers(p,pgc);
    
    if(var==7)
	fill_matrix_F_7p(p,pgc,M,f,xvec,rhsvec);
    
    if(var==8)
	fill_matrix_F_7p_v2(p,pgc,M,f,xvec,rhsvec);
    
    if(var==9)
	fill_matrix_F_13p(p,pgc,M,f,xvec,rhsvec);
    
    if(var==10)
	fill_matrix_F_13p_v2(p,pgc,M,f,xvec,rhsvec);
    
    if(var==11)
	fill_matrix_F_19p(p,pgc,M,f,xvec,rhsvec);
    
    if(var==12)
	fill_matrix_F_19p_v2(p,pgc,M,f,xvec,rhsvec);
  

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
    
	if(var==7||var==9||var==11)
	fillbackvec_F(p,f,xvec,var);
    
    if(var==8||var==10||var==12)
    fillbackvec_F_v2(p,f,xvec,var);
    
    delete_solvers(p,pgc);
    delete_grid(p,pgc); 
    
    p->del_Darray(xvec,p->knox*p->knoy*(p->knoz+1));
}


void hypre_aij::fillbackvec_F(lexer *p, double *f, double *xvec, int var)
{
	HYPRE_IJVectorGetValues(x, p->N7_row, rows, xvec);
	
        n=0;
        FLOOP
        {
        f[FIJK]=xvec[n];
        ++n;
        }
}

void hypre_aij::fillbackvec_F_v2(lexer *p, double *f, double *xvec, int var)
{
	HYPRE_IJVectorGetValues(x, p->N4_row, rows, xvec);
	
        n=0;
        LOOP
        {
        f[FIJK]=xvec[n];
        ++n;
        }
}



#endif
