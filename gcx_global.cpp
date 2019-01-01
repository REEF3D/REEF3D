/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"time.h"

double ghostcell::globalsum(double sendsum)
{
    if(gcx==0)
    return sendsum;

    if(gcx==1)
    MPI_Allreduce(&sendsum,&recvsum,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    return recvsum;
}

int ghostcell::globalisum(int sendisum)
{
    if(gcx==0)
    return sendisum;

    if(gcx==1)
    MPI_Allreduce(&sendisum,&recvisum,1,MPI_INT,MPI_SUM,mpi_comm);
    return recvisum;
}

double ghostcell::globalmin(double sendmin)
{
    if(gcx==0)
    return sendmin;

    if(gcx==1)
    MPI_Allreduce(&sendmin,&recvmin,1,MPI_DOUBLE,MPI_MIN,mpi_comm);
    return recvmin;
}

double ghostcell::globalmax(double sendmax)
{
    if(gcx==0)
    return sendmax;

    if(gcx==1)
    MPI_Allreduce(&sendmax,&recvmax,1,MPI_DOUBLE,MPI_MAX,mpi_comm);
    return recvmax;
}

int ghostcell::globalimin(int sendimin)
{
    if(gcx==0)
    return sendimin;

    if(gcx==1)
    MPI_Allreduce(&sendimin,&recvimin,1,MPI_INT,MPI_MIN,mpi_comm);

    return recvimin;
}

int ghostcell::globalimax(int sendimax)
{
    if(gcx==0)
    return sendimax;

    if(gcx==1)
    MPI_Allreduce(&sendimax,&recvimax,1,MPI_INT,MPI_MAX,mpi_comm);


    return recvimax;
}

double ghostcell::timesync(double t)
{
    if(gcx==0)
    return t;

    if(gcx==1)
    MPI_Bcast(&t,1,MPI_DOUBLE,0,mpi_comm);


    return t;
}

void ghostcell::globalctrl(lexer* p)
{
	if(gcx==1)
	{
    MPI_Bcast(p->ictrl,p->ctrlsize,MPI_INT,0,mpi_comm);
    MPI_Bcast(p->dctrl,p->ctrlsize,MPI_DOUBLE,0,mpi_comm);
	}

}

double ghostcell::timer()
{
    if(gcx==0)
    {
    time_t seconds;
    seconds=time(0);
    return seconds;
    }

    double t=0.0;
    if(gcx==1)
    t=MPI_Wtime();

    return t;
}


