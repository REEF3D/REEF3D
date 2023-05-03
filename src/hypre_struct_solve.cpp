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

#include"hypre_struct.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"

void hypre_struct::solve1234(lexer* p)
{
    p->solver_status=0;
    
	p->solveriter=0;
	    
    HYPRE_StructBiCGSTABSetup(solver, A, b, x);
    p->solver_status = HYPRE_StructBiCGSTABSolve(solver, A, b, x);
    
    HYPRE_StructBiCGSTABGetNumIterations(solver, &num_iterations);
	HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
    
    
    p->solveriter=num_iterations;
    p->final_res = final_res_norm;
}

void hypre_struct::solve(lexer* p, ghostcell *pgc)
{
    p->solver_status=0;
    
	p->solveriter=0;
    
    if(p->N10==11)
    {
    HYPRE_StructPCGSetup(solver, A, b, x);
    p->solver_status = HYPRE_StructPCGSolve(solver, A, b, x);
    
    HYPRE_StructPCGGetNumIterations(solver, &num_iterations);
	HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==12)
    {
    HYPRE_StructGMRESSetup(solver, A, b, x);
    p->solver_status = HYPRE_StructGMRESSolve(solver, A, b, x);
    
    HYPRE_StructGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_StructGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==13)
    {
    HYPRE_StructLGMRESSetup(solver, A, b, x);
    p->solver_status = HYPRE_StructLGMRESSolve(solver, A, b, x);
    
    HYPRE_StructLGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_StructLGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==14)
    {
    HYPRE_StructBiCGSTABSetup(solver, A, b, x);
    p->solver_status = HYPRE_StructBiCGSTABSolve(solver, A, b, x);
    
    HYPRE_StructBiCGSTABGetNumIterations(solver, &num_iterations);
	HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
	
	if(p->N10==15 || p->N10==16 || p->N10==17)
    {
    HYPRE_StructHybridSetup(solver, A, b, x);
    p->solver_status = HYPRE_StructHybridSolve(solver, A, b, x);
    
    HYPRE_StructHybridGetNumIterations(solver, &num_iterations);
	HYPRE_StructHybridGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==18)
    {
    HYPRE_StructPFMGSetup(solver, A, b, x);
    p->solver_status = HYPRE_StructPFMGSolve(solver, A, b, x);
    
    HYPRE_StructPFMGGetNumIterations(solver, &num_iterations);
	HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==19)
    {
    HYPRE_StructSMGSetup(solver, A, b, x);
    p->solver_status = HYPRE_StructSMGSolve(solver, A, b, x);
    
    HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
	HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
    
    
}

#endif
