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

void ghostcell::gcpartnum(int sendnum[6], int recvnum[6])
{
    const void* send_ptrs[6] = {&(sendnum[0]), &(sendnum[1]), &(sendnum[2]), &(sendnum[3]), &(sendnum[4]), &(sendnum[5])};
    void* recv_ptrs[6] = {&(recvnum[0]), &(recvnum[1]), &(recvnum[2]), &(recvnum[3]), &(recvnum[4]), &(recvnum[5])};

    int counts[6] = {1,1,1,1,1,1};

    Sendrecv(send_ptrs, counts, recv_ptrs, counts, MPI_INT);
}