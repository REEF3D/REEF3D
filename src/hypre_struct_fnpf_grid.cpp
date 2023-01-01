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

#include"hypre_struct_fnpf.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void hypre_struct_fnpf::make_grid(lexer* p, ghostcell* pgc)
{
    int kend=0;
    
    if(p->A10==3)
    kend=0;
    
    // grid
    ilower[0] = p->origin_i;
    ilower[1] = p->origin_j;
    ilower[2] = p->origin_k;
    
    iupper[0] = p->knox+p->origin_i-1;
    iupper[1] = p->knoy+p->origin_j-1;
    iupper[2] = p->knoz+p->origin_k-1+kend;
    
    HYPRE_StructGridCreate(pgc->mpi_comm, 3, &grid);
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridAssemble(grid);
    
    
    // stencil
    HYPRE_StructStencilCreate(3, 15, &stencil);

    int entry;
    int offsets[15][3] = {{0,0,0}, {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1},
                          {-1,0,-1},{-1,0,1},{1,0,-1},{1,0,1},{0,-1,-1},{0,-1,1},{0,1,-1},{0,1,1}};

    for (entry=0; entry<15; ++entry)
    HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
    
    // matrix
    HYPRE_StructMatrixCreate(pgc->mpi_comm, grid, stencil, &A);
    HYPRE_StructMatrixInitialize(A);
    
    // vec
    HYPRE_StructVectorCreate(pgc->mpi_comm, grid, &b);
    HYPRE_StructVectorCreate(pgc->mpi_comm, grid, &x);

    HYPRE_StructVectorInitialize(b);
    HYPRE_StructVectorInitialize(x);
}

void hypre_struct_fnpf::make_grid_2Dvert(lexer* p,ghostcell* pgc)
{
    int kend=0;
    
    if(p->A10==3)
    kend=0;
    
    // grid
    ilower[0] = p->origin_i;
    ilower[1] = p->origin_k;
    
    iupper[0] = p->knox+p->origin_i-1;
    iupper[1] = p->knoz+p->origin_k-1+kend;
    
    HYPRE_StructGridCreate(pgc->mpi_comm, 2, &grid);
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridAssemble(grid);
    
    
    // stencil
    HYPRE_StructStencilCreate(2, 9, &stencil);

    int entry;
    int offsets[9][2] = {{0,0}, {-1,0}, {1,0},  {0,-1}, {0,1}, {-1,-1},{-1,1},{1,-1},{1,1}};

    for (entry=0; entry<9; ++entry)
    HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
    
    // matrix
    HYPRE_StructMatrixCreate(pgc->mpi_comm, grid, stencil, &A);
    HYPRE_StructMatrixInitialize(A);
    
    // vec
    HYPRE_StructVectorCreate(pgc->mpi_comm, grid, &b);
    HYPRE_StructVectorCreate(pgc->mpi_comm, grid, &x);

    HYPRE_StructVectorInitialize(b);
    HYPRE_StructVectorInitialize(x);
}

#endif


