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

void hypre_sstruct::create_solver1234(lexer* p,ghostcell* pgc)
{
    // solver for velocities and other scalar variables, e.g. turbulence
    HYPRE_SStructBiCGSTABCreate(pgc->mpi_comm, &solver);
    HYPRE_SStructBiCGSTABSetMaxIter(solver, p->N46);
    HYPRE_SStructBiCGSTABSetTol(solver, p->D29);
    HYPRE_SStructBiCGSTABSetPrintLevel(solver, 0 ); 
    HYPRE_SStructBiCGSTABSetLogging(solver, 1);
    
    HYPRE_SStructJacobiCreate(pgc->mpi_comm, &precond);
    HYPRE_SStructJacobiSetMaxIter(precond,1);

    
    HYPRE_SStructBiCGSTABSetPrecond(solver, HYPRE_SStructJacobiSolve, HYPRE_SStructJacobiSetup, precond);
}

void hypre_sstruct::delete_solver1234(lexer* p,ghostcell* pgc)
{
    
    HYPRE_SStructJacobiDestroy(precond);  
    HYPRE_SStructBiCGSTABDestroy(solver);    
}


void hypre_sstruct::create_solver5(lexer* p, ghostcell* pgc)
{
    // solver for pressure poisson and potential laplace equation
    
    if(p->N10==11)
    {
    HYPRE_SStructPCGCreate(pgc->mpi_comm, &solver);
    HYPRE_SStructPCGSetMaxIter(solver, p->N46 );
    HYPRE_SStructPCGSetTol(solver, p->N44 );
    HYPRE_SStructPCGSetTwoNorm(solver, 1 );
    HYPRE_SStructPCGSetRelChange(solver, 0 );
    HYPRE_SStructPCGSetPrintLevel(solver, 0 ); 
    HYPRE_SStructPCGSetLogging(solver, 1);
    }
    
    if(p->N10==12)
    {
    HYPRE_SStructGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_SStructGMRESSetMaxIter(solver, p->N46);
    HYPRE_SStructGMRESSetKDim(solver,30);
    HYPRE_SStructGMRESSetTol(solver, p->N44);
    HYPRE_SStructGMRESSetPrintLevel(solver, 0);
    HYPRE_SStructGMRESSetLogging(solver, 1);
    }
    
    if(p->N10==13)
    {
    HYPRE_SStructLGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_SStructLGMRESSetMaxIter(solver, p->N46);
    HYPRE_SStructLGMRESSetKDim(solver,30);
    HYPRE_SStructLGMRESSetTol(solver, p->N44);
    HYPRE_SStructLGMRESSetPrintLevel(solver, 0);
    HYPRE_SStructLGMRESSetLogging(solver, 1);
    }
    
    if(p->N10==14)
    {
    HYPRE_SStructBiCGSTABCreate(pgc->mpi_comm, &solver);
    HYPRE_SStructBiCGSTABSetMaxIter(solver, p->N46);
    HYPRE_SStructBiCGSTABSetTol(solver, p->N44);
    HYPRE_SStructBiCGSTABSetPrintLevel(solver, 0 ); 
    HYPRE_SStructBiCGSTABSetLogging(solver, 1);
    }
	
	if(p->N10==15)
    {
    HYPRE_SStructHybridCreate(pgc->mpi_comm, &solver);
	HYPRE_SStructHybridSetSolverType(solver,0);
    HYPRE_SStructHybridSetPCGMaxIter(solver, p->N46);
	HYPRE_SStructHybridSetDSCGMaxIter(solver,1);
    HYPRE_SStructHybridSetTol(solver, p->N44);
    HYPRE_SStructHybridSetPrintLevel(solver, 0 ); 
    HYPRE_SStructHybridSetLogging(solver, 1);
    }
	
	if(p->N10==16)
    {
    HYPRE_SStructHybridCreate(pgc->mpi_comm, &solver);
	HYPRE_SStructHybridSetSolverType(solver,1);
    HYPRE_SStructHybridSetPCGMaxIter(solver, p->N46);
	HYPRE_SStructHybridSetDSCGMaxIter(solver,1);
    HYPRE_SStructHybridSetTol(solver, p->N44);
    HYPRE_SStructHybridSetPrintLevel(solver, 0 ); 
    HYPRE_SStructHybridSetLogging(solver, 1);
    }
	
	if(p->N10==17)
    {
    HYPRE_SStructHybridCreate(pgc->mpi_comm, &solver);
	HYPRE_SStructHybridSetSolverType(solver,2);
    HYPRE_SStructHybridSetPCGMaxIter(solver, p->N46);
	HYPRE_SStructHybridSetDSCGMaxIter(solver,1);
    HYPRE_SStructHybridSetTol(solver, p->N44);
    HYPRE_SStructHybridSetPrintLevel(solver, 0 ); 
    HYPRE_SStructHybridSetLogging(solver, 1);
    }
    
    if(p->N10==18)
    {
    HYPRE_SStructPFMGCreate(pgc->mpi_comm, &solver);
	HYPRE_SStructPFMGSetMaxIter(solver, p->N46);
	HYPRE_SStructPFMGSetTol(solver, p->N44);
	HYPRE_SStructPFMGSetZeroGuess(solver);		
	HYPRE_SStructPFMGSetRAPType(solver, 0);
	HYPRE_SStructPFMGSetRelaxType(solver, 1);
	HYPRE_SStructPFMGSetNumPreRelax(solver, 1);
	HYPRE_SStructPFMGSetNumPostRelax(solver, 1);
	HYPRE_SStructPFMGSetSkipRelax(solver, 0);
	HYPRE_SStructPFMGSetPrintLevel(solver, 0);
	HYPRE_SStructPFMGSetLogging(solver, 0);
    }
    
    if(p->N10==19)
    {
    HYPRE_SStructSMGCreate(pgc->mpi_comm, &solver);
    HYPRE_SStructSMGSetMemoryUse(solver,0);
    HYPRE_SStructSMGSetMaxIter(solver,p->N46);
    HYPRE_SStructSMGSetTol(solver, p->N44);
    HYPRE_SStructSMGSetZeroGuess(solver);
    HYPRE_SStructSMGSetNumPreRelax(solver,1);
    HYPRE_SStructSMGSetNumPostRelax(solver,1);
    }
    
    if(p->N11==10)
    {
    HYPRE_SStructJacobiCreate(pgc->mpi_comm, &precond);
    HYPRE_SStructJacobiSetMaxIter(precond,1);
    }
    
    if(p->N11==11)
    {
    HYPRE_SStructPFMGCreate(pgc->mpi_comm, &precond);
	HYPRE_SStructPFMGSetMaxIter(precond, 1);
	HYPRE_SStructPFMGSetTol(precond, 0.0);
	HYPRE_SStructPFMGSetZeroGuess(precond);		
	HYPRE_SStructPFMGSetRAPType(precond, 0);
	HYPRE_SStructPFMGSetRelaxType(precond, 1);
	HYPRE_SStructPFMGSetNumPreRelax(precond, 1);
	HYPRE_SStructPFMGSetNumPostRelax(precond, 1);
	HYPRE_SStructPFMGSetSkipRelax(precond, 0);
	HYPRE_SStructPFMGSetPrintLevel(precond, 0);
	HYPRE_SStructPFMGSetLogging(precond, 0);
    }

    if(p->N11==11)
    {
    HYPRE_SStructPFMGCreate(pgc->mpi_comm, &precond);
	HYPRE_SStructPFMGSetMaxIter(precond, 1);
	HYPRE_SStructPFMGSetTol(precond, 0.0);
	HYPRE_SStructPFMGSetZeroGuess(precond);		
	HYPRE_SStructPFMGSetRAPType(precond, 0);
	HYPRE_SStructPFMGSetRelaxType(precond, 1);
	HYPRE_SStructPFMGSetNumPreRelax(precond, 1);
	HYPRE_SStructPFMGSetNumPostRelax(precond, 1);
	HYPRE_SStructPFMGSetSkipRelax(precond, 0);
	HYPRE_SStructPFMGSetPrintLevel(precond, 0);
	HYPRE_SStructPFMGSetLogging(precond, 0);
    }
    
    if(p->N11==12)
    {
    HYPRE_SStructSMGCreate(pgc->mpi_comm, &precond);
    HYPRE_SStructSMGSetMemoryUse(precond,0);
    HYPRE_SStructSMGSetMaxIter(precond,1);
    HYPRE_SStructSMGSetTol(precond, 0.0);
    HYPRE_SStructSMGSetZeroGuess(precond);
    HYPRE_SStructSMGSetNumPreRelax(precond,1);
    HYPRE_SStructSMGSetNumPostRelax(precond,1);
    }
	  
    
    if(p->N10==11 && p->N11==10)
    HYPRE_SStructPCGSetPrecond(solver, HYPRE_SStructJacobiSolve, HYPRE_SStructJacobiSetup, precond);
    
    if(p->N10==11 && p->N11==11)
    HYPRE_SStructPCGSetPrecond(solver, HYPRE_SStructPFMGSolve, HYPRE_SStructPFMGSetup, precond);
    
    if(p->N10==11 && p->N11==12)
    HYPRE_SStructPCGSetPrecond(solver, HYPRE_SStructSMGSolve, HYPRE_SStructSMGSetup, precond);
    
    
    if(p->N10==12 && p->N11==10)
    HYPRE_SStructGMRESSetPrecond(solver, HYPRE_SStructJacobiSolve, HYPRE_SStructJacobiSetup, precond);
    
    if(p->N10==12 && p->N11==11)
    HYPRE_SStructGMRESSetPrecond(solver, HYPRE_SStructPFMGSolve, HYPRE_SStructPFMGSetup, precond);
    
    if(p->N10==12 && p->N11==12)
    HYPRE_SStructGMRESSetPrecond(solver, HYPRE_SStructSMGSolve, HYPRE_SStructSMGSetup, precond);
    
    
    if(p->N10==13 && p->N11==10)
    HYPRE_SStructLGMRESSetPrecond(solver, HYPRE_SStructJacobiSolve, HYPRE_SStructJacobiSetup, precond);
    
    if(p->N10==13 && p->N11==11)
    HYPRE_SStructLGMRESSetPrecond(solver, HYPRE_SStructPFMGSolve, HYPRE_SStructPFMGSetup, precond);
    
    if(p->N10==13 && p->N11==12)
    HYPRE_SStructLGMRESSetPrecond(solver, HYPRE_SStructSMGSolve, HYPRE_SStructSMGSetup, precond);
    
    
    if(p->N10==14 && p->N11==10)
    HYPRE_SStructBiCGSTABSetPrecond(solver, HYPRE_SStructJacobiSolve, HYPRE_SStructJacobiSetup, precond);
    
    if(p->N10==14 && p->N11==11)
    HYPRE_SStructBiCGSTABSetPrecond(solver, HYPRE_SStructPFMGSolve, HYPRE_SStructPFMGSetup, precond);
    
    if(p->N10==14 && p->N11==12)
    HYPRE_SStructBiCGSTABSetPrecond(solver, HYPRE_SStructSMGSolve, HYPRE_SStructSMGSetup, precond);
	
	if((p->N10==15 || p->N10==16 || p->N10==17) && p->N11==11)
    HYPRE_SStructHybridSetPrecond(solver, HYPRE_SStructPFMGSolve, HYPRE_SStructPFMGSetup, precond);
    
    if((p->N10==15 || p->N10==16 || p->N10==17) && p->N11==12)
    HYPRE_SStructHybridSetPrecond(solver, HYPRE_SStructSMGSolve, HYPRE_SStructSMGSetup, precond);
}

void hypre_sstruct::delete_solver5(lexer* p,ghostcell* pgc)
{
    if(p->N10==11)
    HYPRE_SStructPCGDestroy(solver);
    
    if(p->N10==12)
    HYPRE_SStructGMRESDestroy(solver);
    
    if(p->N10==13)
    HYPRE_SStructLGMRESDestroy(solver);
    
    if(p->N10==14)
    HYPRE_SStructBiCGSTABDestroy(solver);
	
	if(p->N10==15 || p->N10==16 || p->N10==17)
	HYPRE_SStructHybridDestroy(solver);
    
    if(p->N10==18)
    HYPRE_SStructPFMGDestroy(solver);
    
    if(p->N10==19)
    HYPRE_SStructSMGDestroy(solver);
    
    if(p->N11==11)
    HYPRE_SStructPFMGDestroy(precond);
    
    if(p->N11==12)
    HYPRE_SStructSMGDestroy(precond);
    
}

#endif