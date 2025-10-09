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
    if(cart_comm == MPI_COMM_NULL)
    {
        // Fallback to point-to-point exchanges if the cartesian communicator is unavailable.
        if(count1>0)
        {
            MPI_Isend(send1,count1,MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
            MPI_Irecv(recv1,count1,MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
        }

        if(count4>0)
        {
            MPI_Isend(send4,count4,MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
            MPI_Irecv(recv4,count4,MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
        }

        if(count3>0)
        {
            MPI_Isend(send3,count3,MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
            MPI_Irecv(recv3,count3,MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
        }

        if(count2>0)
        {
            MPI_Isend(send2,count2,MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
            MPI_Irecv(recv2,count2,MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
        }

        if(count5>0)
        {
            MPI_Isend(send5,count5,MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
            MPI_Irecv(recv5,count5,MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
        }

        if(count6>0)
        {
            MPI_Isend(send6,count6,MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
            MPI_Irecv(recv6,count6,MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
        }

        if(count1>0)
        {
            MPI_Wait(&sreq1,&status);
            MPI_Wait(&rreq1,&status);
        }

        if(count4>0)
        {
            MPI_Wait(&sreq4,&status);
            MPI_Wait(&rreq4,&status);
        }

        if(count3>0)
        {
            MPI_Wait(&sreq3,&status);
            MPI_Wait(&rreq3,&status);
        }

        if(count2>0)
        {
            MPI_Wait(&sreq2,&status);
            MPI_Wait(&rreq2,&status);
        }

        if(count5>0)
        {
            MPI_Wait(&sreq5,&status);
            MPI_Wait(&rreq5,&status);
        }

        if(count6>0)
        {
            MPI_Wait(&sreq6,&status);
            MPI_Wait(&rreq6,&status);
        }

        return;
    }

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
