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
Author: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_fnpf.h"
#include"fdm_nhf.h"
#include<sstream>

void ghostcell::mpi_check(lexer* p)
{
    int check=1;

    check=globalisum(check);

    if(p->mpirank==0 && check!=p->mpi_size)
        cout<<"mpi - checksum: "<<check<<" vs "<<p->mpi_size<<" ... mismatch"<<endl;
}

void ghostcell::gcx_ini(lexer* p)
{
    
    int gcx_count[6];
    gcx_count[0] = (p->gcpara1_count+p->flast)*paramargin + p->gcparaco1_count*paramargin;
    gcx_count[1] = (p->gcpara2_count+p->flast)*paramargin + p->gcparaco2_count*paramargin;
    gcx_count[2] = (p->gcpara3_count+p->flast)*paramargin + p->gcparaco3_count*paramargin;
    gcx_count[3] = (p->gcpara4_count+p->flast)*paramargin + p->gcparaco4_count*paramargin;
    gcx_count[4] = p->gcpara5_count*paramargin + p->gcparaco5_count*paramargin;
    gcx_count[5] = p->gcpara6_count*paramargin + p->gcparaco6_count*paramargin;

    p->Darray(send1,gcx_count[0]);
    p->Darray(send2,gcx_count[1]);
    p->Darray(send3,gcx_count[2]);
    p->Darray(send4,gcx_count[3]);
    p->Darray(send5,gcx_count[4]);
    p->Darray(send6,gcx_count[5]);

    p->Darray(recv1,gcx_count[0]);
    p->Darray(recv2,gcx_count[1]);
    p->Darray(recv3,gcx_count[2]);
    p->Darray(recv4,gcx_count[3]);
    p->Darray(recv5,gcx_count[4]);
    p->Darray(recv6,gcx_count[5]);

    p->Iarray(isend1,gcx_count[0]);
    p->Iarray(isend2,gcx_count[1]);
    p->Iarray(isend3,gcx_count[2]);
    p->Iarray(isend4,gcx_count[3]);
    p->Iarray(isend5,gcx_count[4]);
    p->Iarray(isend6,gcx_count[5]);

    p->Iarray(irecv1,gcx_count[0]);
    p->Iarray(irecv2,gcx_count[1]);
    p->Iarray(irecv3,gcx_count[2]);
    p->Iarray(irecv4,gcx_count[3]);
    p->Iarray(irecv5,gcx_count[4]);
    p->Iarray(irecv6,gcx_count[5]);

    if(cart_comm != MPI_COMM_NULL)
        MPI_Comm_free(&cart_comm);

    int dims[3] = {p->mx, p->my, p->mz};
    MPI_Dims_create(p->mpi_size, 3, dims);

    int periods[3] = {p->periodic1,p->periodic2,p->periodic3};
    MPI_Cart_create(mpi_comm, 3, dims, periods, false, &cart_comm);

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

    bool error = false;
    const int nb[6] = {p->nb1, p->nb4, p->nb3, p->nb2, p->nb5, p->nb6};

    for (int dir = 0; dir < 6; ++dir)
    {
        const int expected = (nb[dir] == -2) ? MPI_PROC_NULL : nb[dir];
        if (neighbors[dir] != expected) {
            error = true;
            std::cerr << "Rank " << p->mpirank << " mismatch dir " << dir
                      << " cart=" << neighbors[dir] << " nb=" << expected << std::endl;
        }
    }
    if(cart_comm == MPI_COMM_NULL)
        error = true;

    if(error)
    {
        std::cerr << "MPI Cartesian topology does not match user-specified neighbours or doesn't exist. Exiting." << std::endl;
        exit(1);
    }
    
    
    // nb
    nb0[0] = p->nb1;
	nb0[1] = p->nb2;
	nb0[2] = p->nb3;
	nb0[3] = p->nb4;
	nb0[4] = p->nb5;
	nb0[5] = p->nb6;
    
    stag[0] = 1;
	stag[1] = 2;
	stag[2] = 3;
	stag[3] = 4;
	stag[4] = 5;
	stag[5] = 6;
	
	rtag[0] = 4;
	rtag[1] = 3;
	rtag[2] = 2;
	rtag[3] = 1;
	rtag[4] = 6;
	rtag[5] = 5;
}

void ghostcell::final(bool error)
{
    if(cart_comm != MPI_COMM_NULL)
        MPI_Comm_free(&cart_comm);
    if(mpi_comm != MPI_COMM_NULL)
        MPI_Comm_free(&mpi_comm);
    MPI_Finalize();
    exit(error);
}
