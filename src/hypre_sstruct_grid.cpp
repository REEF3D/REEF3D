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

#include"hypre_sstruct.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void hypre_sstruct::make_grid_7p(lexer* p,fdm* a, ghostcell* pgc)
{
    kend=0;
    numparts=1;
    part=0;
    dimensions = 3;
    variable = 0;
    numvar = 1;
    object_type = HYPRE_SSTRUCT;
    
    if(p->A10==3)
    kend=0;
    
    // grid
    ilower[0] = p->origin_i;
    ilower[1] = p->origin_j;
    ilower[2] = p->origin_k;
    
    iupper[0] = p->knox+p->origin_i-1;
    iupper[1] = p->knoy+p->origin_j-1;
    iupper[2] = p->knoz+p->origin_k-1+kend;
    
    vartypes[0] = HYPRE_SSTRUCT_VARIABLE_CELL;
    
    HYPRE_SStructGridCreate(pgc->mpi_comm, dimensions, numparts, &grid);
    HYPRE_SStructGridSetExtents(grid, part, ilower, iupper);
    HYPRE_SStructGridSetVariables(grid, part, numvar, vartypes);
    HYPRE_SStructGridAssemble(grid);
    
    // stencil
    HYPRE_SStructStencilCreate(3, 7, &stencil);

    int entry;
    int offsets[7][3] = {{0,0,0}, {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}};

    for (entry=0; entry<7; ++entry)
    HYPRE_SStructStencilSetEntry(stencil, entry, offsets[entry], variable);

    // graph
    HYPRE_SStructGraphCreate(pgc->mpi_comm, grid, &graph);
    HYPRE_SStructGraphSetStencil(graph, part, variable, stencil);
    HYPRE_SStructGraphAssemble(graph);

    // matrix
    HYPRE_SStructMatrixCreate(pgc->mpi_comm, graph, &A);
    HYPRE_SStructMatrixSetObjectType(A, object_type);
    HYPRE_SStructMatrixInitialize(A);
    
    // vec
    HYPRE_SStructVectorCreate(pgc->mpi_comm, grid, &b);
    HYPRE_SStructVectorCreate(pgc->mpi_comm, grid, &x);
    
    HYPRE_SStructVectorSetObjectType(b, object_type);
    HYPRE_SStructVectorSetObjectType(x, object_type);

    HYPRE_SStructVectorInitialize(b);
    HYPRE_SStructVectorInitialize(x);
}


void hypre_sstruct::make_grid_13p(lexer* p,fdm* a, ghostcell* pgc)
{
    kend=0;
    numparts=1;
    part=0;
    dimensions = 3;
    variable = 0;
    numvar = 1;
    object_type = HYPRE_PARCSR;
    
    if(p->A10==3)
    kend=0;
    
    // grid
    ilower[0] = p->origin_i;
    ilower[1] = p->origin_j;
    ilower[2] = p->origin_k;
    
    iupper[0] = p->knox+p->origin_i-1;
    iupper[1] = p->knoy+p->origin_j-1;
    iupper[2] = p->knoz+p->origin_k-1+kend;
    
    vartypes[0] = HYPRE_SSTRUCT_VARIABLE_CELL;
    
    HYPRE_SStructGridCreate(pgc->mpi_comm, dimensions, numparts, &grid);
    HYPRE_SStructGridSetExtents(grid, part, ilower, iupper);
    HYPRE_SStructGridSetVariables(grid, part, numvar, vartypes);
    HYPRE_SStructGridAssemble(grid);
    
    
    // stencil
    HYPRE_SStructStencilCreate(3, 13, &stencil);

    int entry;
    
    int offsets[13][3] = {{0,0,0}, {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1},
                                 {-2,0,0}, {2,0,0}, {0,-2,0}, {0,2,0}, {0,0,-2}, {0,0,2}};

    for (entry=0; entry<13; ++entry)
    HYPRE_SStructStencilSetEntry(stencil, entry, offsets[entry], variable);
    
    // graph
    HYPRE_SStructGraphCreate(pgc->mpi_comm, grid, &graph);
    HYPRE_SStructGraphSetStencil(graph, part, variable, stencil);
    HYPRE_SStructGraphAssemble(graph);

    // matrix
    HYPRE_SStructMatrixCreate(pgc->mpi_comm, graph, &A);
    HYPRE_SStructMatrixSetObjectType(A, object_type);
    HYPRE_SStructMatrixInitialize(A);
    
    // vec
    HYPRE_SStructVectorCreate(pgc->mpi_comm, grid, &b);
    HYPRE_SStructVectorCreate(pgc->mpi_comm, grid, &x);
    
    HYPRE_SStructVectorSetObjectType(b, object_type);
    HYPRE_SStructVectorSetObjectType(x, object_type);

    HYPRE_SStructVectorInitialize(b);
    HYPRE_SStructVectorInitialize(x);
}

void hypre_sstruct::make_grid_15p(lexer* p,fdm* a, ghostcell* pgc)
{
    kend=0;
    numparts=1;
    part=0;
    dimensions = 3;
    variable = 0;
    numvar = 1;
    object_type = HYPRE_SSTRUCT;
    
    // grid
    ilower[0] = p->origin_i;
    ilower[1] = p->origin_j;
    ilower[2] = p->origin_k;
    
    iupper[0] = p->knox+p->origin_i-1;
    iupper[1] = p->knoy+p->origin_j-1;
    iupper[2] = p->knoz+p->origin_k-1+kend;
    
    vartypes[0] = HYPRE_SSTRUCT_VARIABLE_CELL;
    
    HYPRE_SStructGridCreate(pgc->mpi_comm, dimensions, numparts, &grid);
    HYPRE_SStructGridSetExtents(grid, part, ilower, iupper);
    HYPRE_SStructGridSetVariables(grid, part, numvar, vartypes);
    HYPRE_SStructGridAssemble(grid);
    
    
    // stencil
    HYPRE_SStructStencilCreate(3, 15, &stencil);

    int entry;
    int offsets[15][3] = {{0,0,0}, {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1},
                          {-1,0,-1},{-1,0,1},{1,0,-1},{1,0,1},{0,-1,-1},{0,-1,1},{0,1,-1},{0,1,1}};

    for (entry=0; entry<15; ++entry)
    HYPRE_SStructStencilSetEntry(stencil, entry, offsets[entry], variable);
    
    // graph
    HYPRE_SStructGraphCreate(pgc->mpi_comm, grid, &graph);
    HYPRE_SStructGraphSetStencil(graph, part, variable, stencil);
    HYPRE_SStructGraphAssemble(graph);

    // matrix
    HYPRE_SStructMatrixCreate(pgc->mpi_comm, graph, &A);
    HYPRE_SStructMatrixSetObjectType(A, object_type);
    HYPRE_SStructMatrixInitialize(A);
    
    // vec
    HYPRE_SStructVectorCreate(pgc->mpi_comm, grid, &b);
    HYPRE_SStructVectorCreate(pgc->mpi_comm, grid, &x);
    
    HYPRE_SStructVectorSetObjectType(b, object_type);
    HYPRE_SStructVectorSetObjectType(x, object_type);

    HYPRE_SStructVectorInitialize(b);
    HYPRE_SStructVectorInitialize(x);
}

void hypre_sstruct::make_grid_2Dvert_9p(lexer* p,fdm* a, ghostcell* pgc)
{ 
    kend=0;
    numparts=1;
    part=0;
    dimensions = 2;
    variable = 0;
    numvar = 1;
    object_type = HYPRE_SSTRUCT;
     
    // grid
    ilower[0] = p->origin_i;
    ilower[1] = p->origin_k;
    
    iupper[0] = p->knox+p->origin_i-1;
    iupper[1] = p->knoz+p->origin_k-1+kend;
    
    vartypes[0] = HYPRE_SSTRUCT_VARIABLE_CELL;
    
    HYPRE_SStructGridCreate(pgc->mpi_comm, dimensions, numparts, &grid);
    HYPRE_SStructGridSetExtents(grid, part, ilower, iupper);
    HYPRE_SStructGridSetVariables(grid, part, numvar, vartypes);
    HYPRE_SStructGridAssemble(grid);
    
    // stencil
    HYPRE_SStructStencilCreate(2, 9, &stencil);

    int entry;
    int offsets[9][2] = {{0,0}, {-1,0},{1,0}, {0,-1},{0,1}, {-1,-1},{-1,1},{1,-1},{1,1}};

    for (entry=0; entry<9; ++entry)
    HYPRE_SStructStencilSetEntry(stencil, entry, offsets[entry], variable);
    
    // graph
    HYPRE_SStructGraphCreate(pgc->mpi_comm, grid, &graph);
    HYPRE_SStructGraphSetStencil(graph, part, variable, stencil);
    HYPRE_SStructGraphAssemble(graph);

    // matrix
    HYPRE_SStructMatrixCreate(pgc->mpi_comm, graph, &A);
    HYPRE_SStructMatrixSetObjectType(A, object_type);
    HYPRE_SStructMatrixInitialize(A);
    
    // vec
    HYPRE_SStructVectorCreate(pgc->mpi_comm, grid, &b);
    HYPRE_SStructVectorCreate(pgc->mpi_comm, grid, &x);
    
    HYPRE_SStructVectorSetObjectType(b, object_type);
    HYPRE_SStructVectorSetObjectType(x, object_type);

    HYPRE_SStructVectorInitialize(b);
    HYPRE_SStructVectorInitialize(x);
}

#endif


