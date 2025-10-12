/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"

void ghostcell::Sendrecv_double(int count1, int count2, int count3, int count4, int count5, int count6)
{
    const double* send_ptrs[6] = {send1, send4, send3, send2, send5, send6};
    double* recv_ptrs[6] = {recv1, recv4, recv3, recv2, recv5, recv6};

    int sendcounts[6] = {
        count1,
        count4,
        count3,
        count2,
        count5,
        count6
    };

    int recvcounts[6] = {
        count1,
        count4,
        count3,
        count2,
        count5,
        count6
    };

    if(cart_comm == MPI_COMM_NULL)
    {
        for(int qn=0; qn<6; ++qn)
            if(sendcounts[qn]>0)
            {
                MPI_Isend(send_ptrs[qn],sendcounts[qn],MPI_DOUBLE,nb[qn],stag[qn],mpi_comm,&sreq[qn]);
                MPI_Irecv(recv_ptrs[qn],recvcounts[qn],MPI_DOUBLE,nb[qn],rtag[qn],mpi_comm,&rreq[qn]);
            }

        for(int qn=0;qn<6;++qn)
            if(sendcounts[qn]>0)
            {
                MPI_Wait(&sreq[qn],&status);
                MPI_Wait(&rreq[qn],&status);
            }

        return;
    }

    int neighbors[6] = {MPI_PROC_NULL};
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

    MPI_Datatype sendtypes[6];
    MPI_Datatype recvtypes[6];
    MPI_Aint sdispls[6];
    MPI_Aint rdispls[6];

    for(int dir=0; dir<6; ++dir)
    {
        if(neighbors[dir] == MPI_PROC_NULL)
        {
            sendcounts[dir] = 0;
            recvcounts[dir] = 0;
        }

        sendtypes[dir] = MPI_DOUBLE;
        recvtypes[dir] = MPI_DOUBLE;
        MPI_Get_address(send_ptrs[dir], &sdispls[dir]);
        MPI_Get_address(recv_ptrs[dir], &rdispls[dir]);
    }

    MPI_Neighbor_alltoallw(MPI_BOTTOM,
                           sendcounts, sdispls, sendtypes,
                           MPI_BOTTOM,
                           recvcounts, rdispls, recvtypes,
                           cart_comm);
}

void ghostcell::Sendrecv_int(int count1, int count2, int count3, int count4, int count5, int count6)
{
    const int* send_ptrs[6] = {isend1, isend4, isend3, isend2, isend5, isend6};
    int* recv_ptrs[6] = {irecv1, irecv4, irecv3, irecv2, irecv5, irecv6};

    int sendcounts[6] = {
        count1,
        count4,
        count3,
        count2,
        count5,
        count6
    };

    int recvcounts[6] = {
        count1,
        count4,
        count3,
        count2,
        count5,
        count6
    };

    if(cart_comm == MPI_COMM_NULL)
    {
        for(int qn=0; qn<6; ++qn)
            if(sendcounts[qn]>0)
            {
                MPI_Isend(send_ptrs[qn],sendcounts[qn],MPI_INT,nb[qn],stag[qn],mpi_comm,&sreq[qn]);
                MPI_Irecv(recv_ptrs[qn],recvcounts[qn],MPI_INT,nb[qn],rtag[qn],mpi_comm,&rreq[qn]);
            }

        for(int qn=0;qn<6;++qn)
            if(sendcounts[qn]>0)
            {
                MPI_Wait(&sreq[qn],&status);
                MPI_Wait(&rreq[qn],&status);
            }

        return;
    }

    int neighbors[6] = {MPI_PROC_NULL};
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

    MPI_Datatype sendtypes[6];
    MPI_Datatype recvtypes[6];
    MPI_Aint sdispls[6];
    MPI_Aint rdispls[6];

    for(int dir=0; dir<6; ++dir)
    {
        if(neighbors[dir] == MPI_PROC_NULL)
        {
            sendcounts[dir] = 0;
            recvcounts[dir] = 0;
        }

        sendtypes[dir] = MPI_INT;
        recvtypes[dir] = MPI_INT;
        MPI_Get_address(send_ptrs[dir], &sdispls[dir]);
        MPI_Get_address(recv_ptrs[dir], &rdispls[dir]);
    }

    MPI_Neighbor_alltoallw(MPI_BOTTOM,
                           sendcounts, sdispls, sendtypes,
                           MPI_BOTTOM,
                           recvcounts, rdispls, recvtypes,
                           cart_comm);
}
