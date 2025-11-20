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

void ghostcell::gcpartx(int sendnum[6], int recvnum[6], double *send[6], double *recv[6])
{
    const void* send_ptrs[6] = {send[0], send[1], send[2], send[3], send[4], send[5]};
    void* recv_ptrs[6] = {recv[0], recv[1], recv[2], recv[3], recv[4], recv[5]};

    Sendrecv(send_ptrs, sendnum, recv_ptrs, recvnum, MPI_DOUBLE);
}