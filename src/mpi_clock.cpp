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
Author: Hans Bihs
--------------------------------------------------------------------*/

//  The time of this rank in blocking MPI calls, for the load statistics of the mesh refinement
//  (amr_load.h, G 43).  The calls REEF3D makes are wrapped through the MPI profiling interface:
//  each wrapper adds its wall-clock time to amr_clock::mpi and calls the PMPI_ version.  The
//  nonblocking posts (Isend, Irecv) are not wrapped; their completion (Wait, Waitall, Test) is.
//  The cost is two MPI_Wtime calls per MPI call.  Build with -DREEF3D_NO_MPI_CLOCK to leave MPI
//  unwrapped (G 43 then reports no MPI time: busy = wall).

#ifndef REEF3D_NO_MPI_CLOCK

#include"amr_load.h"
#include<mpi.h>

namespace
{
struct mclk
{
    double t;
    mclk() : t(PMPI_Wtime()) {}
    ~mclk() { amr_clock::mpi += PMPI_Wtime()-t; ++amr_clock::mpi_calls; }
};
}

extern "C" {

int MPI_Allreduce(const void *sb, void *rb, int n, MPI_Datatype t, MPI_Op op, MPI_Comm cm)
{
    mclk c;
    return PMPI_Allreduce(sb,rb,n,t,op,cm);
}

int MPI_Reduce(const void *sb, void *rb, int n, MPI_Datatype t, MPI_Op op, int root, MPI_Comm cm)
{
    mclk c;
    return PMPI_Reduce(sb,rb,n,t,op,root,cm);
}

int MPI_Exscan(const void *sb, void *rb, int n, MPI_Datatype t, MPI_Op op, MPI_Comm cm)
{
    mclk c;
    return PMPI_Exscan(sb,rb,n,t,op,cm);
}

int MPI_Bcast(void *b, int n, MPI_Datatype t, int root, MPI_Comm cm)
{
    mclk c;
    return PMPI_Bcast(b,n,t,root,cm);
}

int MPI_Barrier(MPI_Comm cm)
{
    mclk c;
    return PMPI_Barrier(cm);
}

int MPI_Gather(const void *sb, int sn, MPI_Datatype st, void *rb, int rn, MPI_Datatype rt, int root, MPI_Comm cm)
{
    mclk c;
    return PMPI_Gather(sb,sn,st,rb,rn,rt,root,cm);
}

int MPI_Gatherv(const void *sb, int sn, MPI_Datatype st, void *rb, const int rn[], const int rd[], MPI_Datatype rt, int root, MPI_Comm cm)
{
    mclk c;
    return PMPI_Gatherv(sb,sn,st,rb,rn,rd,rt,root,cm);
}

int MPI_Allgather(const void *sb, int sn, MPI_Datatype st, void *rb, int rn, MPI_Datatype rt, MPI_Comm cm)
{
    mclk c;
    return PMPI_Allgather(sb,sn,st,rb,rn,rt,cm);
}

int MPI_Allgatherv(const void *sb, int sn, MPI_Datatype st, void *rb, const int rn[], const int rd[], MPI_Datatype rt, MPI_Comm cm)
{
    mclk c;
    return PMPI_Allgatherv(sb,sn,st,rb,rn,rd,rt,cm);
}

int MPI_Alltoall(const void *sb, int sn, MPI_Datatype st, void *rb, int rn, MPI_Datatype rt, MPI_Comm cm)
{
    mclk c;
    return PMPI_Alltoall(sb,sn,st,rb,rn,rt,cm);
}

int MPI_Alltoallv(const void *sb, const int sn[], const int sd[], MPI_Datatype st, void *rb, const int rn[], const int rd[], MPI_Datatype rt, MPI_Comm cm)
{
    mclk c;
    return PMPI_Alltoallv(sb,sn,sd,st,rb,rn,rd,rt,cm);
}

int MPI_Neighbor_alltoallw(const void *sb, const int sn[], const MPI_Aint sd[], const MPI_Datatype st[], void *rb, const int rn[], const MPI_Aint rd[], const MPI_Datatype rt[], MPI_Comm cm)
{
    mclk c;
    return PMPI_Neighbor_alltoallw(sb,sn,sd,st,rb,rn,rd,rt,cm);
}

int MPI_Sendrecv(const void *sb, int sn, MPI_Datatype st, int dest, int stag, void *rb, int rn, MPI_Datatype rt, int src, int rtag, MPI_Comm cm, MPI_Status *s)
{
    mclk c;
    return PMPI_Sendrecv(sb,sn,st,dest,stag,rb,rn,rt,src,rtag,cm,s);
}

int MPI_Send(const void *b, int n, MPI_Datatype t, int dest, int tg, MPI_Comm cm)
{
    mclk c;
    return PMPI_Send(b,n,t,dest,tg,cm);
}

int MPI_Recv(void *b, int n, MPI_Datatype t, int src, int tg, MPI_Comm cm, MPI_Status *s)
{
    mclk c;
    return PMPI_Recv(b,n,t,src,tg,cm,s);
}

int MPI_Wait(MPI_Request *r, MPI_Status *s)
{
    mclk c;
    return PMPI_Wait(r,s);
}

int MPI_Waitall(int n, MPI_Request r[], MPI_Status s[])
{
    mclk c;
    return PMPI_Waitall(n,r,s);
}

int MPI_Waitany(int n, MPI_Request r[], int *idx, MPI_Status *s)
{
    mclk c;
    return PMPI_Waitany(n,r,idx,s);
}

int MPI_Test(MPI_Request *r, int *flag, MPI_Status *s)
{
    mclk c;
    return PMPI_Test(r,flag,s);
}

int MPI_Testsome(int n, MPI_Request r[], int *outc, int idx[], MPI_Status s[])
{
    mclk c;
    return PMPI_Testsome(n,r,outc,idx,s);
}

}

#endif
