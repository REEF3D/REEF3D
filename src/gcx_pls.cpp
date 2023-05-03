/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

void ghostcell::parapls(lexer* p, double** s,double** r, int* psend, int* precv)
{
    if(p->nb1>=0)
    {
	MPI_Isend(&psend[0],1,MPI_INT,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(&precv[0],1,MPI_INT,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->nb2>=0)
    {
	MPI_Isend(&psend[1],1,MPI_INT,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(&precv[1],1,MPI_INT,p->nb2,tag3,mpi_comm,&rreq2);
    }

    if(p->nb3>=0)
    {
	MPI_Isend(&psend[2],1,MPI_INT,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(&precv[2],1,MPI_INT,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->nb4>=0)
    {
	MPI_Isend(&psend[3],1,MPI_INT,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(&precv[3],1,MPI_INT,p->nb4,tag1,mpi_comm,&rreq4);
	}

    if(p->nb5>=0)
    {
	MPI_Isend(&psend[4],1,MPI_INT,p->nb5,tag5,mpi_comm,&sreq5);
	MPI_Irecv(&precv[4],1,MPI_INT,p->nb5,tag6,mpi_comm,&rreq5);
    }

    if(p->nb6>=0)
    {
	MPI_Isend(&psend[5],1,MPI_INT,p->nb6,tag6,mpi_comm,&sreq6);
	MPI_Irecv(&precv[5],1,MPI_INT,p->nb6,tag5,mpi_comm,&rreq6);
    }

//  WAIT

    gcwait(p);

    if(p->nb1>=0)
    {
    MPI_Isend(s[0],psend[0],MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(r[0],precv[0],MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
	}

    if(p->nb2>=0)
    {
	MPI_Isend(s[1],psend[1],MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(r[1],precv[1],MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
	}

    if(p->nb3>=0)
    {
	MPI_Isend(s[2],psend[2],MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(r[2],precv[2],MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
	}

    if(p->nb4>=0)
    {
	MPI_Isend(s[3],psend[3],MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(r[3],precv[3],MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
	}

    if(p->nb5>=0)
    {
	MPI_Isend(s[4],psend[4],MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
	MPI_Irecv(r[4],precv[4],MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
	}

    if(p->nb6>=0)
    {
	MPI_Isend(s[5],psend[5],MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
	MPI_Irecv(r[5],precv[5],MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
	}

	gcwait(p);
}
