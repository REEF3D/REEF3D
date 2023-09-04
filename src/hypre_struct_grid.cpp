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

#include"hypre_struct.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void hypre_struct::make_grid(lexer* p, ghostcell* pgc)
{

    // grid
    ilower[0] = p->origin_i;
    ilower[1] = p->origin_j;
    ilower[2] = p->origin_k;
    
    iupper[0] = p->knox+p->origin_i-1;
    iupper[1] = p->knoy+p->origin_j-1;
    iupper[2] = p->knoz+p->origin_k-1;
    
    // periodic BC
    periodic[0]=0;
    periodic[1]=0;
    periodic[2]=0;
    
    if(p->periodic1>0)
    periodic[0]=p->gknox;
    
    if(p->periodic2>0)
    periodic[1]=p->gknoy;
    
    if(p->periodic3>0)
    periodic[2]=p->gknoz;
    
    HYPRE_StructGridCreate(pgc->mpi_comm, 3, &grid);
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridSetPeriodic (grid, periodic);
    HYPRE_StructGridAssemble(grid);
    
    
    // stencil
    HYPRE_StructStencilCreate(3, 7, &stencil);

    int entry;
    int offsets[7][3] = {{0,0,0}, {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}};

    for (entry=0; entry<7; ++entry)
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

void hypre_struct::make_grid_2Dvert(lexer* p,ghostcell* pgc)
{    
    // grid
    ilower[0] = p->origin_i;
    ilower[1] = p->origin_k;
    
    iupper[0] = p->knox+p->origin_i-1;
    iupper[1] = p->knoz+p->origin_k-1;
    
    // periodic BC
    periodic[0]=0;
    periodic[1]=0;
    periodic[2]=0;
    
    if(p->periodic1>0)
    periodic[0]=p->gknox;
    
    if(p->periodic3>0)
    periodic[1]=p->gknoz;
    
    HYPRE_StructGridCreate(pgc->mpi_comm, 2, &grid);
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridSetPeriodic(grid, periodic);
    HYPRE_StructGridAssemble(grid);
    
    
    // stencil
    HYPRE_StructStencilCreate(2, 5, &stencil);

    int entry;
    int offsets[5][2] = {{0,0}, {-1,0}, {1,0},  {0,-1}, {0,1}};

    for (entry=0; entry<5; ++entry)
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

void hypre_struct::make_grid_15pt(lexer* p, ghostcell* pgc)
{    
    // grid
    ilower[0] = p->origin_i;
    ilower[1] = p->origin_j;
    ilower[2] = p->origin_k;
    
    iupper[0] = p->knox+p->origin_i-1;
    iupper[1] = p->knoy+p->origin_j-1;
    iupper[2] = p->knoz+p->origin_k-1;
    
    // periodic BC
    periodic[0]=0;
    periodic[1]=0;
    periodic[2]=0;
    
    if(p->periodic1>0)
    periodic[0]=p->gknox;
    
    if(p->periodic2>0)
    periodic[1]=p->gknoy;
    
    if(p->periodic3>0)
    periodic[2]=p->gknoz;
    
    HYPRE_StructGridCreate(pgc->mpi_comm, 3, &grid);
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridSetPeriodic (grid, periodic);
    HYPRE_StructGridAssemble(grid);
    
    
    // stencil
    HYPRE_StructStencilCreate(3, 15, &stencil);
    
    /*
sb 7
st 8
nb 9
nt 10
eb 11
et 12
wb 13
wt 14
*/

    int entry;
    int offsets[15][3] = {{0,0,0}, {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1},
                         {-1,0,-1}, {-1,0,1}, {1,0,-1}, {1,0,1}, {0,-1,-1}, {0,-1,1}, {0,1,-1}, {0,1,1}};

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


void hypre_struct::make_grid_2D_9pt(lexer* p,ghostcell* pgc)
{    
    // grid
    ilower[0] = p->origin_i;
    ilower[1] = p->origin_k;
    
    iupper[0] = p->knox+p->origin_i-1;
    iupper[1] = p->knoz+p->origin_k-1;
    
    // periodic BC
    periodic[0]=0;
    periodic[1]=0;
    periodic[2]=0;
    
    if(p->periodic1>0)
    periodic[0]=p->gknox;
    
    if(p->periodic3>0)
    periodic[1]=p->gknoz;
    
    HYPRE_StructGridCreate(pgc->mpi_comm, 2, &grid);
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridSetPeriodic(grid, periodic);
    HYPRE_StructGridAssemble(grid);
    
    
    // stencil
    HYPRE_StructStencilCreate(2, 9, &stencil);

    int entry;
    int offsets[9][2] = {{0,0}, {-1,0}, {1,0},  {0,-1}, {0,1}, {-1,-1}, {-1,1}, {1,-1}, {1,1}};

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


