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

void hypre_aij::make_grid(lexer* p,ghostcell* pgc)
{
    fieldint4 rownum4(p);
    
    pgc->rownum4_update(p,rownum4);
    pgc->gcparaxint(p,rownum4,4);
        
     p->range_col4[0]=0;
    for(n=1;n<=p->M10;++n)
    p->range_col4[n]=p->range_row4[n]-1;
    
        
    HYPRE_IJMatrixCreate(pgc->mpi_comm, p->range_row4[p->mpirank], p->range_col4[p->mpirank+1], p->range_row4[p->mpirank], p->range_col4[p->mpirank+1], &A);
	HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(A);

	HYPRE_IJVectorCreate(pgc->mpi_comm, p->range_row4[p->mpirank], p->range_col4[p->mpirank+1], &b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);
	
	HYPRE_IJVectorCreate(pgc->mpi_comm, p->range_row4[p->mpirank], p->range_col4[p->mpirank+1], &x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);
}

void hypre_aij::delete_grid(lexer* p, ghostcell* pgc)
{
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(b);
    HYPRE_IJVectorDestroy(x);
}


#endif
