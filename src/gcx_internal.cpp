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
    const void* send_ptrs[6] = {send1, send4, send3, send2, send5, send6};
    void* recv_ptrs[6] = {recv1, recv4, recv3, recv2, recv5, recv6};

    int sendcounts[6] = {count1, count4, count3, count2, count5, count6};
    int recvcounts[6] = {count1, count4, count3, count2, count5, count6};

    Sendrecv(send_ptrs, sendcounts, recv_ptrs, recvcounts, MPI_DOUBLE);
}

void ghostcell::Sendrecv_int(int count1, int count2, int count3, int count4, int count5, int count6)
{
    const void* send_ptrs[6] = {isend1, isend4, isend3, isend2, isend5, isend6};
    void* recv_ptrs[6] = {irecv1, irecv4, irecv3, irecv2, irecv5, irecv6};

    int sendcounts[6] = {count1, count4, count3, count2, count5, count6};
    int recvcounts[6] = {count1, count4, count3, count2, count5, count6};

    Sendrecv(send_ptrs, sendcounts, recv_ptrs, recvcounts, MPI_INT);
}

void ghostcell::Sendrecv(const void* send_ptrs[6], int sendcounts[6], void* recv_ptrs[6], int recvcounts[6], MPI_Datatype datatype)
{
    MPI_Datatype sendtypes[6] = {datatype, datatype, datatype, datatype, datatype, datatype};
    MPI_Datatype recvtypes[6] = {datatype, datatype, datatype, datatype, datatype, datatype};
    MPI_Aint sdispls[6];
    MPI_Aint rdispls[6];

    for(int dir=0; dir<6; ++dir)
    {
        if(neighbors[dir] == MPI_PROC_NULL)
        {
            sendcounts[dir] = 0;
            recvcounts[dir] = 0;
        }

        MPI_Get_address(send_ptrs[dir], &sdispls[dir]);
        MPI_Get_address(recv_ptrs[dir], &rdispls[dir]);
    }

    MPI_Neighbor_alltoallw(MPI_BOTTOM,
                           sendcounts, sdispls, sendtypes,
                           MPI_BOTTOM,
                           recvcounts, rdispls, recvtypes,
                           cart_comm);
}
