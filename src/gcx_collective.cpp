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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::gather_int(int *sendbuf, int sendcount, int *recvbuf, int recvcount)
{
    MPI_Gather(sendbuf,sendcount, MPI_INT, recvbuf, recvcount, MPI_INT, 0, mpi_comm);
}

void ghostcell::gather_double(double *sendbuf, int sendcount, double *recvbuf, int recvcount)
{
    MPI_Gather(sendbuf,sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, 0, mpi_comm);
}

void ghostcell::gatherv_int(int *sendbuf, int sendcount, int *recvbuf, int *recvcount, int *recvdispl)
{
    MPI_Gatherv(sendbuf,sendcount, MPI_INT, recvbuf, recvcount, recvdispl, MPI_INT, 0, mpi_comm);
}

void ghostcell::gatherv_double(double *sendbuf, int sendcount, double *recvbuf, int *recvcount, int *recvdispl)
{
    MPI_Gatherv(sendbuf,sendcount, MPI_DOUBLE, recvbuf, recvcount, recvdispl, MPI_DOUBLE, 0, mpi_comm);
}

void ghostcell::allgather_int(int *sendbuf, int sendcount, int *recvbuf, int recvcount)
{
    MPI_Allgather(sendbuf,sendcount, MPI_INT, recvbuf, recvcount, MPI_INT, mpi_comm);
}

void ghostcell::allgatherv_int(int *sendbuf, int sendcount, int *recvbuf, int *recvcount, int *recvdispl)
{
    MPI_Allgatherv(sendbuf,sendcount, MPI_INT, recvbuf, recvcount, recvdispl, MPI_INT, mpi_comm);
}

void ghostcell::bcast_int(int *sendbuf, int sendcount)
{
    MPI_Bcast(sendbuf,sendcount, MPI_INT, 0, mpi_comm);
}

void ghostcell::bcast_double(double *sendbuf, int sendcount)
{
    MPI_Bcast(sendbuf,sendcount, MPI_DOUBLE, 0, mpi_comm);
}
