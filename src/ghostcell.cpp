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
Author: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_fnpf.h"
#include"fdm_nhf.h"
#include<sstream>

ghostcell::ghostcell(int& argc, char **argv, lexer *p)
{
    MPI_Init(&argc,&argv);
    MPI_Comm_dup(MPI_COMM_WORLD,&mpi_comm);

    MPI_Comm_rank(mpi_comm,&p->mpirank);
    MPI_Comm_size(mpi_comm,&p->mpi_size);

    ghostcell::p=p;
}

void ghostcell::mpi_check(lexer* p)
{
    int check=1;

    check=globalisum(check);

    if(p->mpirank==0 && check!=p->mpi_size)
        cout<<"mpi - checksum: "<<check<<" vs "<<p->mpi_size<<" ... mismatch"<<endl;
}

void ghostcell::gcini(lexer* p)
{
    margin=p->margin;
    paramargin=p->margin;
    gamma=p->B29;
    orderext=2;
    orderext2=2;
    orderdir2=2;

    if(p->B23==1)
        orderdir=2;
    else if(p->B23==2)
        orderdir=3;

    imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
    kmin=p->kmin;
    kmax=p->kmax;

    if(p->B20==1)
    {
        gclabel_u=4;
        gclabel_v=4;
        gclabel_w=4;

        gclabel_utopo=4;
        gclabel_vtopo=4;
        gclabel_wtopo=4;
    }
    else if(p->B20==2)
    {
        gclabel_u=5;
        gclabel_v=5;
        gclabel_w=5;

        gclabel_utopo=5;
        gclabel_vtopo=5;
        gclabel_wtopo=5;
    }
    else if(p->B20==3)
    {
        gclabel_u=2;
        gclabel_v=2;
        gclabel_w=2;

        gclabel_utopo=2;
        gclabel_vtopo=2;
        gclabel_wtopo=2;
    }
    else if(p->B20==4)
    {
        gclabel_u=5;
        gclabel_v=5;
        gclabel_w=5;

        gclabel_utopo=4;
        gclabel_vtopo=4;
        gclabel_wtopo=4;
    }

    gclabel_k=4;
    gclabel_e=4;

    gclabel_u_orth=1;
    gclabel_v_orth=1;
    gclabel_w_orth=1;
    gclabel_vel=5;    

    // for reflective BC
    if(p->B23==2)
    {
        gclabel_u=12;
        gclabel_v=12;
        gclabel_w=12;

        gclabel_u_orth=11;
        gclabel_v_orth=11;
        gclabel_w_orth=11;
    }

    gclabel_lsm=4;

    awa_lable=0;
    if(p->B99>=3)
        awa_lable=1;


    if(p->B75==1)
    {
        gclabel_u_out=4;
        gclabel_v_out=4;
        gclabel_w_out=4;
    }
    else if(p->B75==2)
    {
        gclabel_u_out=6;
        gclabel_v_out=6;
        gclabel_w_out=6;
    }
    else if(p->B75==3)
    {
        gclabel_u_out=0;
        gclabel_v_out=6;
        gclabel_w_out=6;
    }

    gclabel_outflow=1;
    if(p->B60==3||p->B60==4)
        gclabel_outflow=0;

    gclabel_u_in=1;
    gclabel_v_in=1;
    gclabel_w_in=1;
    gclabel_lsm_in=gclabel_lsm;

    if(p->I230>0 || p->B98>=3 || p->B60>0)
    {
        gclabel_u_in=0;
        gclabel_v_in=0;
        gclabel_w_in=0;
        gclabel_lsm_in=0;
    }

    if(p->B98>=3)
    {
        gclabel_u_in=0;
        gclabel_v_in=0;
        gclabel_w_in=0;
        gclabel_lsm_in=0;
    }

    // pressure bc labels
    gclabel_press=4;
    gclabel_press_in=gclabel_press;

    if(p->B76==2 || p->B76==3)
        gclabel_press_in=0;

    // pressure inflow
    pressin_lable=0;
    if(p->B76!=1)
        pressin_lable=1;

    // pressure outflow
    pressout_lable=0;
    if(p->B77==1 || p->B77==10)
        pressout_lable=1;


    // sflow slip/no-slip
    if(p->A217==1 && p->A10==2)
    {
        gclabel_u=4;
        gclabel_v=4;
    }
    else if(p->A217==2 && p->A10==2)
    {
        gclabel_u=5;
        gclabel_v=5;
    }

    for(m=0;m<15;m++)
    {
        y[m]=0.0;
        x[m]=0.0;
    }
    
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
}

void ghostcell::fdm2D_update(fdm2D *bb)
{
    b=bb;
}

void ghostcell::fdm_fnpf_update(fdm_fnpf *cc)
{
    c=cc;
}

void ghostcell::fdm_nhf_update(fdm_nhf *dd)
{
    d=dd;
}

void ghostcell::fdm_update(fdm *aa)
{
    a=aa;
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
