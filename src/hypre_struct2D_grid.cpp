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
#include"ghostcell.h"

void hypre_struct2D::make_grid(lexer* p, ghostcell* pgc)
{
    // grid
    ilower[0] = p->origin_i;
    ilower[1] = p->origin_j;
    
    iupper[0] = p->knox+p->origin_i-1;
    iupper[1] = p->knoy+p->origin_j-1;
    
    HYPRE_StructGridCreate(pgc->mpi_comm, 2, &grid);
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridAssemble(grid);
    
    
    // stencil
    HYPRE_StructStencilCreate(2, 5, &stencil);

    int entry;
    int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};

    for (entry=0; entry<5; ++entry)
    HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
    
    // matrix
    HYPRE_StructMatrixCreate(pgc->mpi_comm, grid, stencil, &A);
    HYPRE_StructMatrixInitialize(A);
    
    // vec
    HYPRE_StructVectorCreate(pgc->mpi_comm, grid, &rhs);
    HYPRE_StructVectorCreate(pgc->mpi_comm, grid, &x);

    HYPRE_StructVectorInitialize(rhs);
    HYPRE_StructVectorInitialize(x);
}



#endif
