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

#include"hypre_sstruct_fnpf.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"

void hypre_sstruct_fnpf::solve(lexer* p, ghostcell *pgc)
{
	p->solveriter=0;
    
    if(p->N10==31)
    {
    HYPRE_SStructPCGSetup(solver, A, b, x);
    HYPRE_SStructPCGSolve(solver, A, b, x);
    
    HYPRE_SStructPCGGetNumIterations(solver, &num_iterations);
	HYPRE_SStructPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==32)
    {
    HYPRE_SStructGMRESSetup(solver, A, b, x);
    HYPRE_SStructGMRESSolve(solver, A, b, x);
    
    HYPRE_SStructGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_SStructGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==33)
    {
    HYPRE_SStructLGMRESSetup(solver, A, b, x);
    HYPRE_SStructLGMRESSolve(solver, A, b, x);
    
    HYPRE_SStructLGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_SStructLGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==34)
    {
    HYPRE_SStructBiCGSTABSetup(solver, A, b, x);
    HYPRE_SStructBiCGSTABSolve(solver, A, b, x);
    
    HYPRE_SStructBiCGSTABGetNumIterations(solver, &num_iterations);
	HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
	
    if(p->N10==38)
    {
    HYPRE_SStructSysPFMGSetup(solver, A, b, x);
    HYPRE_SStructSysPFMGSolve(solver, A, b, x);
    
    HYPRE_SStructSysPFMGGetNumIterations(solver, &num_iterations);
	HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
	p->solveriter=num_iterations;
    p->final_res = final_res_norm;
    
}

#endif
