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

#include"hypre_aij.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void hypre_aij::make_grid_F(lexer* p, ghostcell* pgc)
{
    
    int* rownum7;
    p->Iarray(rownum7,p->imax*p->jmax*(p->kmax+2));
    
    pgc->rownum7_update(p,rownum7);
    pgc->flagx7(p,rownum7);
    
        
     p->range_col7[0]=0;
    for(n=1;n<=p->M10;++n)
    p->range_col7[n]=p->range_row7[n]-1;
    
        
    HYPRE_IJMatrixCreate(pgc->mpi_comm, p->range_row7[p->mpirank], p->range_col7[p->mpirank+1], p->range_row7[p->mpirank], p->range_col7[p->mpirank+1], &A);
	HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(A);

	HYPRE_IJVectorCreate(pgc->mpi_comm, p->range_row7[p->mpirank], p->range_col7[p->mpirank+1], &b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);
	
	HYPRE_IJVectorCreate(pgc->mpi_comm, p->range_row7[p->mpirank], p->range_col7[p->mpirank+1], &x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);
    
    p->del_Iarray(rownum7,p->imax*p->jmax*(p->kmax+2));
}



#endif
