/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"hypre_struct.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"

hypre_struct::hypre_struct(lexer* p,fdm* a,ghostcell *pgc)
{	
    int vecsize=p->knox*p->knoy*p->knoz; 
    
    if(p->A10==3 && p->A300==1)
    vecsize=p->knox*p->knoy*(p->knoz+1); 
    
    p->Iarray(ilower,3);
    p->Iarray(iupper,3);
    p->Darray(values,vecsize*7);
    
    make_grid(p,a,pgc);	
}

hypre_struct::~hypre_struct()
{
}

void hypre_struct::start(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& xvec, vec& rhsvec, int var, int gcv, double stop_crit)
{
    if(var>=1 && var<=4)
    solve1234(p,a,pgc,f,xvec,rhsvec,var);
    
    if(var==5)
    solve5(p,a,pgc,f,xvec,rhsvec,var);
}

void hypre_struct::startF(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var, int gcv, double stop_crit)
{
    solve7(p,pgc,f,rhs,M,var);
}

void hypre_struct::solve1234(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& xvec, vec& rhsvec, int var)
{
    numiter=0;
	p->solveriter=0;
	
	create_solver1234(p,pgc);
    
    if(var==1)
    fill_matrix1(p,a,pgc,f);
    
    if(var==2)
    fill_matrix2(p,a,pgc,f);
    
    if(var==3)
    fill_matrix3(p,a,pgc,f);
    
    if(var==4)
    fill_matrix4(p,a,pgc,f);
    
    
    HYPRE_StructBiCGSTABSetup(solver, A, b, x);
    HYPRE_StructBiCGSTABSolve(solver, A, b, x);
    
    HYPRE_StructBiCGSTABGetNumIterations(solver, &num_iterations);
	HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
    
    
    p->solveriter=num_iterations;
        
    
    if(var==1)
    fillbackvec1(p,f,xvec,var);
    
    if(var==2)
    fillbackvec2(p,f,xvec,var);
    
    if(var==3)
    fillbackvec3(p,f,xvec,var);
    
    if(var==4)
    fillbackvec4(p,f,xvec,var);
	
	delete_solver1234(p,pgc);
    
}

void hypre_struct::solve5(lexer* p,fdm* a, ghostcell* pgc, field &f, vec& xvec, vec& rhsvec, int var)
{
	numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);

    fill_matrix4(p,a,pgc,f);


    if(p->N10==11)
    {
    HYPRE_StructPCGSetup(solver, A, b, x);
    HYPRE_StructPCGSolve(solver, A, b, x);
    
    HYPRE_StructPCGGetNumIterations(solver, &num_iterations);
	HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==12)
    {
    HYPRE_StructGMRESSetup(solver, A, b, x);
    HYPRE_StructGMRESSolve(solver, A, b, x);
    
    HYPRE_StructGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_StructGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==13)
    {
    HYPRE_StructLGMRESSetup(solver, A, b, x);
    HYPRE_StructLGMRESSolve(solver, A, b, x);
    
    HYPRE_StructLGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_StructLGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==14)
    {
    HYPRE_StructBiCGSTABSetup(solver, A, b, x);
    HYPRE_StructBiCGSTABSolve(solver, A, b, x);
    
    HYPRE_StructBiCGSTABGetNumIterations(solver, &num_iterations);
	HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
	
	if(p->N10==15 || p->N10==16 || p->N10==17)
    {
    HYPRE_StructHybridSetup(solver, A, b, x);
    HYPRE_StructHybridSolve(solver, A, b, x);
    
    HYPRE_StructHybridGetNumIterations(solver, &num_iterations);
	HYPRE_StructHybridGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==18)
    {
    HYPRE_StructPFMGSetup(solver, A, b, x);
    HYPRE_StructPFMGSolve(solver, A, b, x);
    
    HYPRE_StructPFMGGetNumIterations(solver, &num_iterations);
	HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==19)
    {
    HYPRE_StructSMGSetup(solver, A, b, x);
    HYPRE_StructSMGSolve(solver, A, b, x);
    
    HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
	HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }

	

	p->solveriter=num_iterations;
        
    fillbackvec4(p,f,xvec,var);
	
	delete_solver5(p,pgc);
}


void hypre_struct::solve7(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
{
    numiter=0;
	p->solveriter=0;
	
    create_solver5(p,pgc);

    fill_matrix7(p,pgc,f,rhs,M);


    if(p->N10==11)
    {
    HYPRE_StructPCGSetup(solver, A, b, x);
    HYPRE_StructPCGSolve(solver, A, b, x);
    
    HYPRE_StructPCGGetNumIterations(solver, &num_iterations);
	HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==12)
    {
    HYPRE_StructGMRESSetup(solver, A, b, x);
    HYPRE_StructGMRESSolve(solver, A, b, x);
    
    HYPRE_StructGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_StructGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==13)
    {
    HYPRE_StructLGMRESSetup(solver, A, b, x);
    HYPRE_StructLGMRESSolve(solver, A, b, x);
    
    HYPRE_StructLGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_StructLGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==14)
    {
    HYPRE_StructBiCGSTABSetup(solver, A, b, x);
    HYPRE_StructBiCGSTABSolve(solver, A, b, x);
    
    HYPRE_StructBiCGSTABGetNumIterations(solver, &num_iterations);
	HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
	
	if(p->N10==15 || p->N10==16 || p->N10==17)
    {
    HYPRE_StructHybridSetup(solver, A, b, x);
    HYPRE_StructHybridSolve(solver, A, b, x);
    
    HYPRE_StructHybridGetNumIterations(solver, &num_iterations);
	HYPRE_StructHybridGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==18)
    {
    HYPRE_StructPFMGSetup(solver, A, b, x);
    HYPRE_StructPFMGSolve(solver, A, b, x);
    
    HYPRE_StructPFMGGetNumIterations(solver, &num_iterations);
	HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==19)
    {
    HYPRE_StructSMGSetup(solver, A, b, x);
    HYPRE_StructSMGSolve(solver, A, b, x);
    
    HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
	HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }


	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
        
    fillbackvec7(p,f,var);
	
	delete_solver5(p,pgc);
}

void hypre_struct::setup(lexer* p,fdm* a, ghostcell* pgc, int var, cpt &C)
{
}

void hypre_struct::solve(lexer* p,fdm* a, ghostcell* pgc, vec& xvec, vec& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt &C)
{
}

void hypre_struct::fillxvec1(lexer* p, fdm* a, field& f)
{
}
	
void hypre_struct::fillxvec2(lexer* p, fdm* a, field& f)
{
}
	
void hypre_struct::fillxvec3(lexer* p, fdm* a, field& f)
{
}
	
void hypre_struct::fillxvec4(lexer* p, fdm* a, field& f)
{
}

#endif
