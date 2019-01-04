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

#include"hypre_sstruct.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"

void hypre_sstruct::solve1234(lexer* p)
{
	p->solveriter=0;
	    
    HYPRE_SStructBiCGSTABSetup(solver, A, b, x);
    HYPRE_SStructBiCGSTABSolve(solver, A, b, x);
    
    HYPRE_SStructBiCGSTABGetNumIterations(solver, &num_iterations);
	HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
    
    
    p->solveriter=num_iterations;
    p->final_res = final_res_norm;
}

void hypre_sstruct::solve(lexer* p)
{
	p->solveriter=0;

    if(p->N10==11)
    {
    HYPRE_SStructPCGSetup(solver, A, b, x);
    HYPRE_SStructPCGSolve(solver, A, b, x);
    
    HYPRE_SStructPCGGetNumIterations(solver, &num_iterations);
	HYPRE_SStructPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==12)
    {
    HYPRE_SStructGMRESSetup(solver, A, b, x);
    HYPRE_SStructGMRESSolve(solver, A, b, x);
    
    HYPRE_SStructGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_SStructGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==13)
    {
    HYPRE_SStructLGMRESSetup(solver, A, b, x);
    HYPRE_SStructLGMRESSolve(solver, A, b, x);
    
    HYPRE_SStructLGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_SStructLGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==14)
    {
    HYPRE_SStructBiCGSTABSetup(solver, A, b, x);
    HYPRE_SStructBiCGSTABSolve(solver, A, b, x);
    
    HYPRE_SStructBiCGSTABGetNumIterations(solver, &num_iterations);
	HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
	
	if(p->N10==15 || p->N10==16 || p->N10==17)
    {
    HYPRE_SStructHybridSetup(solver, A, b, x);
    HYPRE_SStructHybridSolve(solver, A, b, x);
    
    HYPRE_SStructHybridGetNumIterations(solver, &num_iterations);
	HYPRE_SStructHybridGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==18)
    {
    HYPRE_SStructPFMGSetup(solver, A, b, x);
    HYPRE_SStructPFMGSolve(solver, A, b, x);
    
    HYPRE_SStructPFMGGetNumIterations(solver, &num_iterations);
	HYPRE_SStructPFMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==19)
    {
    HYPRE_SStructSMGSetup(solver, A, b, x);
    HYPRE_SStructSMGSolve(solver, A, b, x);
    
    HYPRE_SStructSMGGetNumIterations(solver, &num_iterations);
	HYPRE_SStructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
}

#endif
