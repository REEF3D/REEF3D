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

#include"hypre_struct2D.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void hypre_struct2D::solve(lexer* p, ghostcell *pgc)
{
    int feedback=0;
    numiter=0;
	p->solveriter=0;
    
    if(p->N10==11)
    {
    HYPRE_StructPCGSetup(solver, A, rhs, x);
    feedback = HYPRE_StructPCGSolve(solver, A, rhs, x);
    
    HYPRE_StructPCGGetNumIterations(solver, &num_iterations);
	HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==12)
    {
    HYPRE_StructGMRESSetup(solver, A, rhs, x);
    feedback = HYPRE_StructGMRESSolve(solver, A, rhs, x);
    
    HYPRE_StructGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_StructGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==13)
    {
    HYPRE_StructLGMRESSetup(solver, A, rhs, x);
    feedback = HYPRE_StructLGMRESSolve(solver, A, rhs, x);
    
    HYPRE_StructLGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_StructLGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==14)
    {
    HYPRE_StructBiCGSTABSetup(solver, A, rhs, x);
    feedback = HYPRE_StructBiCGSTABSolve(solver, A, rhs, x);
    
    HYPRE_StructBiCGSTABGetNumIterations(solver, &num_iterations);
	HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
	
	if(p->N10==15 || p->N10==16 || p->N10==17)
    {
    HYPRE_StructHybridSetup(solver, A, rhs, x);
    feedback = HYPRE_StructHybridSolve(solver, A, rhs, x);
    
    HYPRE_StructHybridGetNumIterations(solver, &num_iterations);
	HYPRE_StructHybridGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==18)
    {
    HYPRE_StructPFMGSetup(solver, A, rhs, x);
    feedback = HYPRE_StructPFMGSolve(solver, A, rhs, x);
    
    HYPRE_StructPFMGGetNumIterations(solver, &num_iterations);
	HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }
    
    if(p->N10==19)
    {
    HYPRE_StructSMGSetup(solver, A, rhs, x);
    feedback = HYPRE_StructSMGSolve(solver, A, rhs, x);
    
    HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
	HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    }

	
	p->solveriter=num_iterations;
    
    feedback = pgc->globalimax(feedback);
    /*
    if(feedback>=1)
    {
    if(p->mpirank==0)
    cout<<endl<<" HYPRE solver broke down! Emergency Stop! "<<feedback<<endl<<endl;
    
    exit(0);
    }*/
}

#endif