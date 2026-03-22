/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_fnpf.h"
#include"fdm_nhf.h"
#include<sstream>

void ghostcell::gcx_cart_topology(lexer* p)
{
    // 1D
    if(p->mx>=1 && p->my==1 && p->mz==1)
    {
    ndims = 1;
    
    if(cart_comm != MPI_COMM_NULL)
    MPI_Comm_free(&cart_comm);

    int dims[1] = {p->mx};
    MPI_Dims_create(p->mpi_size, ndims, dims);

    int periods[1] = {p->periodic1};
    MPI_Cart_create(mpi_comm, ndims, dims, periods, false, &cart_comm);

    int cart_neg = MPI_PROC_NULL;
    int cart_pos = MPI_PROC_NULL;

    MPI_Cart_shift(cart_comm, 0, 1, &cart_neg, &cart_pos);
    neighbors[0] = cart_neg;
    neighbors[1] = cart_pos;

    neighbors[2] = MPI_PROC_NULL;
    neighbors[3] = MPI_PROC_NULL;

    neighbors[4] = MPI_PROC_NULL;
    neighbors[5] = MPI_PROC_NULL;
    }
    
    
    // 2D
    if(p->mx>1 && p->my>1 && p->mz==1)
    {
    ndims = 2;
    
    if(cart_comm != MPI_COMM_NULL)
    MPI_Comm_free(&cart_comm);

    int dims[2] = {p->mx, p->my};
    MPI_Dims_create(p->mpi_size, ndims, dims);

    int periods[2] = {p->periodic1,p->periodic2};
    MPI_Cart_create(mpi_comm, ndims, dims, periods, false, &cart_comm);

    int cart_neg = MPI_PROC_NULL;
    int cart_pos = MPI_PROC_NULL;

    MPI_Cart_shift(cart_comm, 0, 1, &cart_neg, &cart_pos);
    neighbors[0] = cart_neg;
    neighbors[1] = cart_pos;

    MPI_Cart_shift(cart_comm, 1, 1, &cart_neg, &cart_pos);
    neighbors[2] = cart_neg;
    neighbors[3] = cart_pos;

    neighbors[4] = MPI_PROC_NULL;
    neighbors[5] = MPI_PROC_NULL;
    }
    
    
    // 3D
    if(p->mx>1 && p->my>1 && p->mz>1)
    {
    ndims = 3;
    
    if(cart_comm != MPI_COMM_NULL)
    MPI_Comm_free(&cart_comm);

    int dims[3] = {p->mx, p->my, p->mz};
    MPI_Dims_create(p->mpi_size, ndims, dims);

    int periods[3] = {p->periodic1,p->periodic2,p->periodic3};
    MPI_Cart_create(mpi_comm, ndims, dims, periods, false, &cart_comm);

    int cart_neg = MPI_PROC_NULL;
    int cart_pos = MPI_PROC_NULL;

    MPI_Cart_shift(cart_comm, 0, 1, &cart_neg, &cart_pos);
    neighbors[0] = cart_neg;
    neighbors[1] = cart_pos;

    MPI_Cart_shift(cart_comm, 1, 1, &cart_neg, &cart_pos);
    neighbors[2] = cart_neg;
    neighbors[3] = cart_pos;

    MPI_Cart_shift(cart_comm, 2, 1, &cart_neg, &cart_pos);
    neighbors[4] = cart_neg;
    neighbors[5] = cart_pos;
    }

    
}
