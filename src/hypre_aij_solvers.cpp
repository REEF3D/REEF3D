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


void hypre_aij::create_solvers(lexer* p, ghostcell* pgc)
{
    if(p->N10==21)
    {
    HYPRE_ParCSRPCGCreate(pgc->mpi_comm, &solver);
    HYPRE_PCGSetMaxIter(solver, p->N46 );
    HYPRE_PCGSetTol(solver, p->N44 );
    HYPRE_PCGSetTwoNorm(solver, 1 );
    HYPRE_PCGSetRelChange(solver, 0 );
    HYPRE_PCGSetPrintLevel(solver, 0 ); 
    HYPRE_PCGSetLogging(solver, 1);
    }
    
    if(p->N10==22)
    {
    HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_GMRESSetMaxIter(solver, p->N46);
    HYPRE_GMRESSetKDim(solver,30);
    HYPRE_GMRESSetTol(solver, p->N44);
    HYPRE_GMRESSetPrintLevel(solver, 0);
    HYPRE_GMRESSetLogging(solver, 1);
    }
    
    if(p->N10==23)
    {
    HYPRE_ParCSRLGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_LGMRESSetMaxIter(solver, p->N46);
    HYPRE_LGMRESSetKDim(solver,30);
    HYPRE_LGMRESSetTol(solver, p->N44);
    HYPRE_LGMRESSetPrintLevel(solver, 0);
    HYPRE_LGMRESSetLogging(solver, 1);
    }
    
    if(p->N10==24)
    {
    HYPRE_ParCSRBiCGSTABCreate(pgc->mpi_comm, &solver);
    HYPRE_BiCGSTABSetMaxIter(solver, p->N46);
    HYPRE_BiCGSTABSetTol(solver, p->N44);
    HYPRE_BiCGSTABSetPrintLevel(solver, 0 ); 
    HYPRE_BiCGSTABSetLogging(solver, 1);
    }
    
    if(p->N10==25)
    {
    HYPRE_BoomerAMGCreate(&solver);
    HYPRE_BoomerAMGSetPrintLevel(solver, 0); 
    HYPRE_BoomerAMGSetCoarsenType(solver, 22);
    HYPRE_BoomerAMGSetRelaxType(solver, 3); 
    HYPRE_BoomerAMGSetNumSweeps(solver, 1);
    HYPRE_BoomerAMGSetTol(solver, p->N44);
    HYPRE_BoomerAMGSetMaxIter(solver, p->N46); 
    }

    if(p->N11==21)
    {
    HYPRE_BoomerAMGCreate(&precond);
    HYPRE_BoomerAMGSetPrintLevel(precond, 0); 
    HYPRE_BoomerAMGSetCoarsenType(precond, 22);
    HYPRE_BoomerAMGSetRelaxType(precond, 3); 
    HYPRE_BoomerAMGSetNumSweeps(precond, 1);
    HYPRE_BoomerAMGSetTol(precond, 0.0);
    HYPRE_BoomerAMGSetMaxIter(precond, 1); 
    }
    
	  
    
    if(p->N10==21 && p->N11==21)
    HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
    
    if(p->N10==22 && p->N11==21)
    HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
    
    if(p->N10==23 && p->N11==21)
    HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
    
    if(p->N10==24 && p->N11==21)
    HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
    
}

void hypre_aij::delete_solvers(lexer* p, ghostcell* pgc)
{
    if(p->N10==21)
    HYPRE_ParCSRPCGDestroy(solver);
    
    if(p->N10==22)
    HYPRE_ParCSRGMRESDestroy(solver);
    
    if(p->N10==23)
    HYPRE_ParCSRLGMRESDestroy(solver);
    
    if(p->N10==24)
    HYPRE_ParCSRBiCGSTABDestroy(solver);
    
    if(p->N10==25)
    HYPRE_BoomerAMGDestroy(solver);
    
    if(p->N11==21)
    HYPRE_BoomerAMGDestroy(precond);
}

#endif
