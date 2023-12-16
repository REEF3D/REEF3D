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

#include"hypre_struct2D.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"ghostcell.h"


void hypre_struct2D::create_solvers(lexer* p, ghostcell* pgc)
{
    if(p->N10==11)
    {
    HYPRE_StructPCGCreate(pgc->mpi_comm, &solver);
    HYPRE_StructPCGSetMaxIter(solver, p->N46 );
    HYPRE_StructPCGSetTol(solver, p->N44 );
    HYPRE_StructPCGSetTwoNorm(solver, 1 );
    HYPRE_StructPCGSetRelChange(solver, 0 );
    HYPRE_StructPCGSetPrintLevel(solver, 0 ); 
    HYPRE_StructPCGSetLogging(solver, 1);
    }
    
    if(p->N10==12)
    {
    HYPRE_StructGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_StructGMRESSetMaxIter(solver, p->N46);
    HYPRE_StructGMRESSetKDim(solver,30);
    HYPRE_StructGMRESSetTol(solver, p->N44);
    HYPRE_StructGMRESSetPrintLevel(solver, 0);
    HYPRE_StructGMRESSetLogging(solver, 1);
    }
    
    if(p->N10==13)
    {
    HYPRE_StructLGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_StructLGMRESSetMaxIter(solver, p->N46);
    HYPRE_StructLGMRESSetKDim(solver,30);
    HYPRE_StructLGMRESSetTol(solver, p->N44);
    HYPRE_StructLGMRESSetPrintLevel(solver, 0);
    HYPRE_StructLGMRESSetLogging(solver, 1);
    }
    
    if(p->N10==14)
    {
    HYPRE_StructBiCGSTABCreate(pgc->mpi_comm, &solver);
    HYPRE_StructBiCGSTABSetMaxIter(solver, p->N46);
    HYPRE_StructBiCGSTABSetTol(solver, p->N44);
    HYPRE_StructBiCGSTABSetPrintLevel(solver, 0 ); 
    HYPRE_StructBiCGSTABSetLogging(solver, 1);
    }
	
	if(p->N10==15)
    {
    HYPRE_StructHybridCreate(pgc->mpi_comm, &solver);
	HYPRE_StructHybridSetSolverType(solver,0);
    HYPRE_StructHybridSetPCGMaxIter(solver, p->N46);
	HYPRE_StructHybridSetDSCGMaxIter(solver,1);
    HYPRE_StructHybridSetTol(solver, p->N44);
    HYPRE_StructHybridSetPrintLevel(solver, 0 ); 
    HYPRE_StructHybridSetLogging(solver, 1);
    }
	
	if(p->N10==16)
    {
    HYPRE_StructHybridCreate(pgc->mpi_comm, &solver);
	HYPRE_StructHybridSetSolverType(solver,1);
    HYPRE_StructHybridSetPCGMaxIter(solver, p->N46);
	HYPRE_StructHybridSetDSCGMaxIter(solver,1);
    HYPRE_StructHybridSetTol(solver, p->N44);
    HYPRE_StructHybridSetPrintLevel(solver, 0 ); 
    HYPRE_StructHybridSetLogging(solver, 1);
    }
	
	if(p->N10==17)
    {
    HYPRE_StructHybridCreate(pgc->mpi_comm, &solver);
	HYPRE_StructHybridSetSolverType(solver,2);
    HYPRE_StructHybridSetPCGMaxIter(solver, p->N46);
	HYPRE_StructHybridSetDSCGMaxIter(solver,1);
    HYPRE_StructHybridSetTol(solver, p->N44);
    HYPRE_StructHybridSetPrintLevel(solver, 0 ); 
    HYPRE_StructHybridSetLogging(solver, 1);
    }
    
    if(p->N10==18)
    {
    HYPRE_StructPFMGCreate(pgc->mpi_comm, &solver);
	HYPRE_StructPFMGSetMaxIter(solver, p->N46);
	HYPRE_StructPFMGSetTol(solver, p->N44);
	HYPRE_StructPFMGSetZeroGuess(solver);		
	HYPRE_StructPFMGSetRAPType(solver, 0);
	HYPRE_StructPFMGSetRelaxType(solver, 1);
	HYPRE_StructPFMGSetNumPreRelax(solver, 1);
	HYPRE_StructPFMGSetNumPostRelax(solver, 1);
	HYPRE_StructPFMGSetSkipRelax(solver, 0);
	HYPRE_StructPFMGSetPrintLevel(solver, 0);
	HYPRE_StructPFMGSetLogging(solver, 0);
    }
    
    if(p->N10==19)
    {
    HYPRE_StructSMGCreate(pgc->mpi_comm, &solver);
    HYPRE_StructSMGSetMemoryUse(solver,0);
    HYPRE_StructSMGSetMaxIter(solver,p->N46);
    HYPRE_StructSMGSetTol(solver, p->N44);
    HYPRE_StructSMGSetZeroGuess(solver);
    HYPRE_StructSMGSetNumPreRelax(solver,1);
    HYPRE_StructSMGSetNumPostRelax(solver,1);
    }

    if(p->N11==11)
    {
    HYPRE_StructPFMGCreate(pgc->mpi_comm, &precond);
	HYPRE_StructPFMGSetMaxIter(precond, 1);
	HYPRE_StructPFMGSetTol(precond, 0.0);
	HYPRE_StructPFMGSetZeroGuess(precond);		
	HYPRE_StructPFMGSetRAPType(precond, 0);
	HYPRE_StructPFMGSetRelaxType(precond, 1);
	HYPRE_StructPFMGSetNumPreRelax(precond, 1);
	HYPRE_StructPFMGSetNumPostRelax(precond, 3);
	HYPRE_StructPFMGSetSkipRelax(precond, 0);
	HYPRE_StructPFMGSetPrintLevel(precond, 0);
	HYPRE_StructPFMGSetLogging(precond, 0);
    }
    
    if(p->N11==12)
    {
    HYPRE_StructSMGCreate(pgc->mpi_comm, &precond);
    HYPRE_StructSMGSetMemoryUse(precond,0);
    HYPRE_StructSMGSetMaxIter(precond,1);
    HYPRE_StructSMGSetTol(precond, 0.0);
    HYPRE_StructSMGSetZeroGuess(precond);
    HYPRE_StructSMGSetNumPreRelax(precond,1);
    HYPRE_StructSMGSetNumPostRelax(precond,1);
    }
	  
    
    if(p->N10==11 && p->N11==11)
    HYPRE_StructPCGSetPrecond(solver, HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup, precond);
    
    if(p->N10==11 && p->N11==12)
    HYPRE_StructPCGSetPrecond(solver, HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, precond);
    
    
    if(p->N10==12 && p->N11==11)
    HYPRE_StructGMRESSetPrecond(solver, HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup, precond);
    
    if(p->N10==12 && p->N11==12)
    HYPRE_StructGMRESSetPrecond(solver, HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, precond);
    
    
    if(p->N10==13 && p->N11==11)
    HYPRE_StructLGMRESSetPrecond(solver, HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup, precond);
    
    if(p->N10==13 && p->N11==12)
    HYPRE_StructLGMRESSetPrecond(solver, HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, precond);
    
    
    if(p->N10==14 && p->N11==11)
    HYPRE_StructBiCGSTABSetPrecond(solver, HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup, precond);
    
    if(p->N10==14 && p->N11==12)
    HYPRE_StructBiCGSTABSetPrecond(solver, HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, precond);
	
	if((p->N10==15 || p->N10==16 || p->N10==17) && p->N11==11)
    HYPRE_StructHybridSetPrecond(solver, HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup, precond);
    
    if((p->N10==15 || p->N10==16 || p->N10==17) && p->N11==12)
    HYPRE_StructHybridSetPrecond(solver, HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, precond);
	
}

void hypre_struct2D::delete_solvers(lexer* p, ghostcell* pgc)
{
    if(p->N10==11)
    HYPRE_StructPCGDestroy(solver);
    
    if(p->N10==12)
    HYPRE_StructGMRESDestroy(solver);
    
    if(p->N10==13)
    HYPRE_StructLGMRESDestroy(solver);
    
    if(p->N10==14)
    HYPRE_StructBiCGSTABDestroy(solver);
	
	if(p->N10==15 || p->N10==16 || p->N10==17)
	HYPRE_StructHybridDestroy(solver);
    
    if(p->N10==18)
    HYPRE_StructPFMGDestroy(solver);
    
    if(p->N10==19)
    HYPRE_StructSMGDestroy(solver);
    
    if(p->N11==11)
    HYPRE_StructPFMGDestroy(precond);
    
    if(p->N11==12)
    HYPRE_StructSMGDestroy(precond);
    
}

#endif
