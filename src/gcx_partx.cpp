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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"

void ghostcell::gcpartx(lexer *p, int *sendnum, int *recvnum, double **send, double **recv)
{
    //  SEND / RECEIVE
    for(int qn=0;qn<6;++qn)
	{
	MPI_Isend(send[qn],sendnum[qn],MPI_DOUBLE,nb[qn],stag[qn],mpi_comm,&sreq[qn]);
	MPI_Irecv(recv[qn],recvnum[qn],MPI_DOUBLE,nb[qn],rtag[qn],mpi_comm,&rreq[qn]);
    }

    //  WAIT
	for(int qn=0;qn<6;++qn)
	{
    MPI_Wait(&sreq[qn],&status);
	MPI_Wait(&rreq[qn],&status);
	}
}