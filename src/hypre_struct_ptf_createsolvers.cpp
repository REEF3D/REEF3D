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

#include"hypre_struct_ptf.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"

void hypre_struct_ptf::create_solver1234(lexer* p,ghostcell* pgc)
{
    // solver for velocities and other scalar variables, e.g. turbulence
    HYPRE_struct_ptfBiCGSTABCreate(pgc->mpi_comm, &solver);
    HYPRE_struct_ptfBiCGSTABSetMaxIter(solver, p->N46);
    HYPRE_struct_ptfBiCGSTABSetTol(solver, p->N43);
    HYPRE_struct_ptfBiCGSTABSetPrintLevel(solver, 0 ); 
    HYPRE_struct_ptfBiCGSTABSetLogging(solver, 1);
    
    HYPRE_struct_ptfJacobiCreate(pgc->mpi_comm, &precond);
    HYPRE_struct_ptfJacobiSetMaxIter(precond,1);

    
    HYPRE_struct_ptfBiCGSTABSetPrecond(solver, HYPRE_struct_ptfJacobiSolve, HYPRE_struct_ptfJacobiSetup, precond);
}

void hypre_struct_ptf::delete_solver1234(lexer* p,ghostcell* pgc)
{
    
    HYPRE_struct_ptfJacobiDestroy(precond);  
    HYPRE_struct_ptfBiCGSTABDestroy(solver);    
}

void hypre_struct_ptf::create_solver5(lexer* p, ghostcell* pgc)
{
    // solver for pressure poisson and potential laplace equation
    
    if(solve_type==11)
    {
    HYPRE_struct_ptfPCGCreate(pgc->mpi_comm, &solver);
    HYPRE_struct_ptfPCGSetMaxIter(solver, p->N46 );
    HYPRE_struct_ptfPCGSetTol(solver, p->N44 );
    HYPRE_struct_ptfPCGSetTwoNorm(solver, 1 );
    HYPRE_struct_ptfPCGSetRelChange(solver, 0 );
    HYPRE_struct_ptfPCGSetPrintLevel(solver, 0 ); 
    HYPRE_struct_ptfPCGSetLogging(solver, 1);
    }
    
    if(solve_type==12)
    {
    HYPRE_struct_ptfGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_struct_ptfGMRESSetMaxIter(solver, p->N46);
    HYPRE_struct_ptfGMRESSetKDim(solver,30);
    HYPRE_struct_ptfGMRESSetTol(solver, p->N44);
    HYPRE_struct_ptfGMRESSetPrintLevel(solver, 0);
    HYPRE_struct_ptfGMRESSetLogging(solver, 1);
    }
    
    if(solve_type==13)
    {
    HYPRE_struct_ptfLGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_struct_ptfLGMRESSetMaxIter(solver, p->N46);
    HYPRE_struct_ptfLGMRESSetKDim(solver,30);
    HYPRE_struct_ptfLGMRESSetTol(solver, p->N44);
    HYPRE_struct_ptfLGMRESSetPrintLevel(solver, 0);
    HYPRE_struct_ptfLGMRESSetLogging(solver, 1);
    }
    
    if(solve_type==14)
    {
    HYPRE_struct_ptfBiCGSTABCreate(pgc->mpi_comm, &solver);
    HYPRE_struct_ptfBiCGSTABSetMaxIter(solver, p->N46);
    HYPRE_struct_ptfBiCGSTABSetTol(solver, p->N44);
    HYPRE_struct_ptfBiCGSTABSetPrintLevel(solver, 0 ); 
    HYPRE_struct_ptfBiCGSTABSetLogging(solver, 1);
    }
	
	if(solve_type==15)
    {
    HYPRE_struct_ptfHybridCreate(pgc->mpi_comm, &solver);
	HYPRE_struct_ptfHybridSetSolverType(solver,0);
    HYPRE_struct_ptfHybridSetPCGMaxIter(solver, p->N46);
	HYPRE_struct_ptfHybridSetDSCGMaxIter(solver,1);
    HYPRE_struct_ptfHybridSetTol(solver, p->N44);
    HYPRE_struct_ptfHybridSetPrintLevel(solver, 0 ); 
    HYPRE_struct_ptfHybridSetLogging(solver, 1);
    }
	
	if(solve_type==16)
    {
    HYPRE_struct_ptfHybridCreate(pgc->mpi_comm, &solver);
	HYPRE_struct_ptfHybridSetSolverType(solver,1);
    HYPRE_struct_ptfHybridSetPCGMaxIter(solver, p->N46);
	HYPRE_struct_ptfHybridSetDSCGMaxIter(solver,1);
    HYPRE_struct_ptfHybridSetTol(solver, p->N44);
    HYPRE_struct_ptfHybridSetPrintLevel(solver, 0 ); 
    HYPRE_struct_ptfHybridSetLogging(solver, 1);
    }
	
	if(solve_type==17)
    {
    HYPRE_struct_ptfHybridCreate(pgc->mpi_comm, &solver);
	HYPRE_struct_ptfHybridSetSolverType(solver,2);
    HYPRE_struct_ptfHybridSetPCGMaxIter(solver, p->N46);
	HYPRE_struct_ptfHybridSetDSCGMaxIter(solver,1);
    HYPRE_struct_ptfHybridSetTol(solver, p->N44);
    HYPRE_struct_ptfHybridSetPrintLevel(solver, 0 ); 
    HYPRE_struct_ptfHybridSetLogging(solver, 1);
    }
    
    if(solve_type==18)
    {
    HYPRE_struct_ptfPFMGCreate(pgc->mpi_comm, &solver);
	HYPRE_struct_ptfPFMGSetMaxIter(solver, p->N46);
	HYPRE_struct_ptfPFMGSetTol(solver, p->N44);
	HYPRE_struct_ptfPFMGSetZeroGuess(solver);		
	HYPRE_struct_ptfPFMGSetRAPType(solver, 0);
	HYPRE_struct_ptfPFMGSetRelaxType(solver, 1);
	HYPRE_struct_ptfPFMGSetNumPreRelax(solver, 1);
	HYPRE_struct_ptfPFMGSetNumPostRelax(solver, 1);
	HYPRE_struct_ptfPFMGSetSkipRelax(solver, 0);
	HYPRE_struct_ptfPFMGSetPrintLevel(solver, 0);
	HYPRE_struct_ptfPFMGSetLogging(solver, 0);
    }
    
    if(solve_type==19)
    {
    HYPRE_struct_ptfSMGCreate(pgc->mpi_comm, &solver);
    HYPRE_struct_ptfSMGSetMemoryUse(solver,0);
    HYPRE_struct_ptfSMGSetMaxIter(solver,p->N46);
    HYPRE_struct_ptfSMGSetTol(solver, p->N44);
    HYPRE_struct_ptfSMGSetZeroGuess(solver);
    HYPRE_struct_ptfSMGSetNumPreRelax(solver,1);
    HYPRE_struct_ptfSMGSetNumPostRelax(solver,1);
    }
    
    if(precon_type==10)
    {
    HYPRE_struct_ptfJacobiCreate(pgc->mpi_comm, &precond);
    HYPRE_struct_ptfJacobiSetMaxIter(precond,1);
    }
    
    if(precon_type==11)
    {
    HYPRE_struct_ptfPFMGCreate(pgc->mpi_comm, &precond);
	HYPRE_struct_ptfPFMGSetMaxIter(precond, 1);
	HYPRE_struct_ptfPFMGSetTol(precond, 0.0);
	HYPRE_struct_ptfPFMGSetZeroGuess(precond);		
	HYPRE_struct_ptfPFMGSetRAPType(precond, 0);    // now: 0; before: 0
	HYPRE_struct_ptfPFMGSetRelaxType(precond, 3);  // now: 3; before: 1
	HYPRE_struct_ptfPFMGSetNumPreRelax(precond, 1);
	HYPRE_struct_ptfPFMGSetNumPostRelax(precond, 1);
	HYPRE_struct_ptfPFMGSetSkipRelax(precond, 0);  // now: 0; before: 0
	HYPRE_struct_ptfPFMGSetPrintLevel(precond, 0);
	HYPRE_struct_ptfPFMGSetLogging(precond, 0);
    }
    
    if(precon_type==12)
    {
    HYPRE_struct_ptfSMGCreate(pgc->mpi_comm, &precond);
    HYPRE_struct_ptfSMGSetMemoryUse(precond,0);
    HYPRE_struct_ptfSMGSetMaxIter(precond,1);
    HYPRE_struct_ptfSMGSetTol(precond, 0.0);
    HYPRE_struct_ptfSMGSetZeroGuess(precond);
    HYPRE_struct_ptfSMGSetNumPreRelax(precond,1);
    HYPRE_struct_ptfSMGSetNumPostRelax(precond,1);
    }
	  
    
    if(solve_type==11 && precon_type==10)
    HYPRE_struct_ptfPCGSetPrecond(solver, HYPRE_struct_ptfJacobiSolve, HYPRE_struct_ptfJacobiSetup, precond);
    
    if(solve_type==11 && precon_type==11)
    HYPRE_struct_ptfPCGSetPrecond(solver, HYPRE_struct_ptfPFMGSolve, HYPRE_struct_ptfPFMGSetup, precond);
    
    if(solve_type==11 && precon_type==12)
    HYPRE_struct_ptfPCGSetPrecond(solver, HYPRE_struct_ptfSMGSolve, HYPRE_struct_ptfSMGSetup, precond);
    
    
    if(solve_type==12 && precon_type==10)
    HYPRE_struct_ptfGMRESSetPrecond(solver, HYPRE_struct_ptfJacobiSolve, HYPRE_struct_ptfJacobiSetup, precond);
    
    if(solve_type==12 && precon_type==11)
    HYPRE_struct_ptfGMRESSetPrecond(solver, HYPRE_struct_ptfPFMGSolve, HYPRE_struct_ptfPFMGSetup, precond);
    
    if(solve_type==12 && precon_type==12)
    HYPRE_struct_ptfGMRESSetPrecond(solver, HYPRE_struct_ptfSMGSolve, HYPRE_struct_ptfSMGSetup, precond);
    
    
    if(solve_type==13 && precon_type==10)
    HYPRE_struct_ptfLGMRESSetPrecond(solver, HYPRE_struct_ptfJacobiSolve, HYPRE_struct_ptfJacobiSetup, precond);
    
    if(solve_type==13 && precon_type==11)
    HYPRE_struct_ptfLGMRESSetPrecond(solver, HYPRE_struct_ptfPFMGSolve, HYPRE_struct_ptfPFMGSetup, precond);
    
    if(solve_type==13 && precon_type==12)
    HYPRE_struct_ptfLGMRESSetPrecond(solver, HYPRE_struct_ptfSMGSolve, HYPRE_struct_ptfSMGSetup, precond);
    
    
    if(solve_type==14 && precon_type==10)
    HYPRE_struct_ptfBiCGSTABSetPrecond(solver, HYPRE_struct_ptfJacobiSolve, HYPRE_struct_ptfJacobiSetup, precond);
    
    if(solve_type==14 && precon_type==11)
    HYPRE_struct_ptfBiCGSTABSetPrecond(solver, HYPRE_struct_ptfPFMGSolve, HYPRE_struct_ptfPFMGSetup, precond);
    
    if(solve_type==14 && precon_type==12)
    HYPRE_struct_ptfBiCGSTABSetPrecond(solver, HYPRE_struct_ptfSMGSolve, HYPRE_struct_ptfSMGSetup, precond);
	
	if((solve_type==15 || solve_type==16 || solve_type==17) && precon_type==11)
    HYPRE_struct_ptfHybridSetPrecond(solver, HYPRE_struct_ptfPFMGSolve, HYPRE_struct_ptfPFMGSetup, precond);
    
    if((solve_type==15 || solve_type==16 || solve_type==17) && precon_type==12)
    HYPRE_struct_ptfHybridSetPrecond(solver, HYPRE_struct_ptfSMGSolve, HYPRE_struct_ptfSMGSetup, precond);
}

void hypre_struct_ptf::delete_solver5(lexer* p,ghostcell* pgc)
{
    if(solve_type==11)
    HYPRE_struct_ptfPCGDestroy(solver);
    
    if(solve_type==12)
    HYPRE_struct_ptfGMRESDestroy(solver);
    
    if(solve_type==13)
    HYPRE_struct_ptfLGMRESDestroy(solver);
    
    if(solve_type==14)
    HYPRE_struct_ptfBiCGSTABDestroy(solver);
	
	if(solve_type==15 || solve_type==16 || solve_type==17)
	HYPRE_struct_ptfHybridDestroy(solver);
    
    if(solve_type==18)
    HYPRE_struct_ptfPFMGDestroy(solver);
    
    if(solve_type==19)
    HYPRE_struct_ptfSMGDestroy(solver);
    
    if(precon_type==11)
    HYPRE_struct_ptfPFMGDestroy(precond);
    
    if(precon_type==12)
    HYPRE_struct_ptfSMGDestroy(precond);
    
}

#endif
