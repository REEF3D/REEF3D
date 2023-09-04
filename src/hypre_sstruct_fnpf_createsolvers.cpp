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

void hypre_sstruct_fnpf::create_solver5(lexer* p, ghostcell* pgc)
{
    // solver for pressure poisson and potential laplace equation
    
    if(p->N10==31)
    {
    HYPRE_SStructPCGCreate(pgc->mpi_comm, &solver);
    HYPRE_SStructPCGSetMaxIter(solver, p->N46 );
    HYPRE_SStructPCGSetTol(solver, p->N44 );
    HYPRE_SStructPCGSetTwoNorm(solver, 1 );
    HYPRE_SStructPCGSetRelChange(solver, 0 );
    HYPRE_SStructPCGSetPrintLevel(solver, 0 ); 
    HYPRE_SStructPCGSetLogging(solver, 1);
    }
    
    if(p->N10==32)
    {
    HYPRE_SStructGMRESCreate(pgc->mpi_comm, &solver);
    HYPRE_SStructGMRESSetMaxIter(solver, p->N46);
    HYPRE_SStructGMRESSetKDim(solver,30);
    HYPRE_SStructGMRESSetTol(solver, p->N44);
    HYPRE_SStructGMRESSetPrintLevel(solver, 0);
    HYPRE_SStructGMRESSetLogging(solver, 1);
    }
    
    if(p->N10==33)
    {
    HYPRE_SStructLGMRESCreate(pgc->mpi_comm, &solver);
    HYPRE_SStructLGMRESSetMaxIter(solver, p->N46);
    HYPRE_SStructLGMRESSetKDim(solver,30);
    HYPRE_SStructLGMRESSetTol(solver, p->N44);
    HYPRE_SStructLGMRESSetPrintLevel(solver, 0);
    HYPRE_SStructLGMRESSetLogging(solver, 1);
    }
    
    if(p->N10==34)
    {
    HYPRE_SStructBiCGSTABCreate(pgc->mpi_comm, &solver);
    HYPRE_SStructBiCGSTABSetMaxIter(solver, p->N46);
    HYPRE_SStructBiCGSTABSetTol(solver, p->N44);
    HYPRE_SStructBiCGSTABSetPrintLevel(solver, 0 ); 
    HYPRE_SStructBiCGSTABSetLogging(solver, 1);
    }
    
    if(p->N10==34)
    {
    HYPRE_SStructBiCGSTABCreate(pgc->mpi_comm, &solver);
    HYPRE_SStructBiCGSTABSetMaxIter(solver, p->N46);
    HYPRE_SStructBiCGSTABSetTol(solver, p->N44);
    HYPRE_SStructBiCGSTABSetPrintLevel(solver, 0 ); 
    HYPRE_SStructBiCGSTABSetLogging(solver, 1);
    }
    
    if(p->N10==38)
    {
    HYPRE_SStructSysPFMGCreate(pgc->mpi_comm, &solver);
	HYPRE_SStructSysPFMGSetMaxIter(solver, p->N46);
	HYPRE_SStructSysPFMGSetTol(solver, p->N44);
	HYPRE_SStructSysPFMGSetZeroGuess(solver);		
	//HYPRE_SStructSysPFMGSetRAPType(solver, 0);
	HYPRE_SStructSysPFMGSetRelaxType(solver, 1);
	HYPRE_SStructSysPFMGSetNumPreRelax(solver, 1);
	HYPRE_SStructSysPFMGSetNumPostRelax(solver, 1);
	HYPRE_SStructSysPFMGSetSkipRelax(solver, 0);
	HYPRE_SStructSysPFMGSetPrintLevel(solver, 0);
	HYPRE_SStructSysPFMGSetLogging(solver, 0);
    }
    
    if(p->N11==31)
    {
    HYPRE_SStructSysPFMGCreate(pgc->mpi_comm, &precond);
	HYPRE_SStructSysPFMGSetMaxIter(precond, 1);
	HYPRE_SStructSysPFMGSetTol(precond, 0.0);
	HYPRE_SStructSysPFMGSetZeroGuess(precond);		
	//HYPRE_SStructSysPFMGSetRAPType(precond, 0);
	HYPRE_SStructSysPFMGSetRelaxType(precond, 3);
	HYPRE_SStructSysPFMGSetNumPreRelax(precond, 1);
	HYPRE_SStructSysPFMGSetNumPostRelax(precond, 1);
	HYPRE_SStructSysPFMGSetSkipRelax(precond, 0);
	HYPRE_SStructSysPFMGSetPrintLevel(precond, 0);
	HYPRE_SStructSysPFMGSetLogging(precond, 0);
    }
    
    if(p->N11==32)
    {
    HYPRE_SStructSplitCreate(pgc->mpi_comm, &precond);
    HYPRE_SStructSplitSetMaxIter(precond, 1);
    HYPRE_SStructSplitSetTol(precond, 0.0);
    HYPRE_SStructSplitSetZeroGuess(precond);
    //HYPRE_SStructSplitSetRAPType(precond, 0);
    HYPRE_SStructSplitSetStructSolver(precond, HYPRE_PFMG);
    }
    
    if(p->N11==33)
    {
    HYPRE_SStructSplitCreate(pgc->mpi_comm, &precond);
    HYPRE_SStructSplitSetMaxIter(precond, 1);
    HYPRE_SStructSplitSetTol(precond, 0.0);
    HYPRE_SStructSplitSetZeroGuess(precond);
    HYPRE_SStructSplitSetStructSolver(precond, HYPRE_SMG);
    }

    if(p->N10==31 && p->N11==31)
    HYPRE_SStructPCGSetPrecond(solver, HYPRE_SStructSysPFMGSolve, HYPRE_SStructSysPFMGSetup, precond);
  
    if(p->N10==32 && p->N11==31)
    HYPRE_SStructGMRESSetPrecond(solver, HYPRE_SStructSysPFMGSolve, HYPRE_SStructSysPFMGSetup, precond);

    if(p->N10==33 && p->N11==31)
    HYPRE_SStructLGMRESSetPrecond(solver, HYPRE_SStructSysPFMGSolve, HYPRE_SStructSysPFMGSetup, precond);
    
    if(p->N10==34 && p->N11==31)
    HYPRE_SStructBiCGSTABSetPrecond(solver, HYPRE_SStructSysPFMGSolve, HYPRE_SStructSysPFMGSetup, precond);
    
    
    if(p->N10==31 && (p->N11==32 || p->N11==33))
    HYPRE_SStructPCGSetPrecond(solver,  HYPRE_SStructSplitSolve, HYPRE_SStructSplitSetup, precond);
  
    if(p->N10==32 && (p->N11==32 || p->N11==33))
    HYPRE_SStructGMRESSetPrecond(solver, HYPRE_SStructSplitSolve, HYPRE_SStructSplitSetup, precond);

    if(p->N10==33 && (p->N11==32 || p->N11==33))
    HYPRE_SStructLGMRESSetPrecond(solver, HYPRE_SStructSplitSolve, HYPRE_SStructSplitSetup, precond);
    
    if(p->N10==34 && (p->N11==32 || p->N11==33))
    HYPRE_SStructBiCGSTABSetPrecond(solver, HYPRE_SStructSplitSolve, HYPRE_SStructSplitSetup, precond);
}

void hypre_sstruct_fnpf::delete_solver5(lexer* p,ghostcell* pgc)
{
     if(p->N10==31)
    HYPRE_SStructPCGDestroy(solver);
    
    if(p->N10==32)
    HYPRE_SStructGMRESDestroy(solver);
    
    if(p->N10==33)
    HYPRE_SStructLGMRESDestroy(solver);
    
    if(p->N10==34)
    HYPRE_SStructBiCGSTABDestroy(solver);
	
    if(p->N10==38)
    HYPRE_SStructSysPFMGDestroy(solver);
    
    if(p->N11==31)
    HYPRE_SStructSysPFMGDestroy(precond);
    
    if(p->N11==32)
    HYPRE_SStructSplitDestroy(precond);
    
    if(p->N11==33)
    HYPRE_SStructSplitDestroy(precond);
    
}

#endif
