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
#include "particles_obj.h"

void ghostcell::para_tracersobj(lexer* p ,particles_obj* s ,particles_obj* r)
{
    // Setup amount information
    int send[6], recv[6];
    for (int n=0;n<6;n++)
    {
        send[n]=s[n].loopindex;
        recv[n]=0;
    }

    // Transfer amount information
    if(p->nb1>=0)
    {
        MPI_Isend(&send[0],1,MPI_INT,p->nb1,tag1,mpi_comm,&sreq1);
        MPI_Irecv(&recv[0],1,MPI_INT,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->nb2>=0)
    {
        MPI_Isend(&send[1],1,MPI_INT,p->nb2,tag2,mpi_comm,&sreq2);
        MPI_Irecv(&recv[1],1,MPI_INT,p->nb2,tag3,mpi_comm,&rreq2);
    }

    if(p->nb3>=0)
    {
        MPI_Isend(&send[2],1,MPI_INT,p->nb3,tag3,mpi_comm,&sreq3);
        MPI_Irecv(&recv[2],1,MPI_INT,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->nb4>=0)
    {
        MPI_Isend(&send[3],1,MPI_INT,p->nb4,tag4,mpi_comm,&sreq4);
        MPI_Irecv(&recv[3],1,MPI_INT,p->nb4,tag1,mpi_comm,&rreq4);
	}

    if(p->nb5>=0)
    {
        MPI_Isend(&send[4],1,MPI_INT,p->nb5,tag5,mpi_comm,&sreq5);
        MPI_Irecv(&recv[4],1,MPI_INT,p->nb5,tag6,mpi_comm,&rreq5);
    }

    if(p->nb6>=0)
    {
        MPI_Isend(&send[5],1,MPI_INT,p->nb6,tag6,mpi_comm,&sreq6);
        MPI_Irecv(&recv[5],1,MPI_INT,p->nb6,tag5,mpi_comm,&rreq6);
    }

    gcwait(p);

    for (int n=0;n<6;n++)
        r[n].fill(recv[n],false,1);

    // Tracer data

    // Transfer X
    if(p->nb1>=0)
    {
        MPI_Isend(s[0].X,send[0],MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
        MPI_Irecv(r[0].X,recv[0],MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
	}

    if(p->nb2>=0)
    {
        MPI_Isend(s[1].X,send[1],MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
        MPI_Irecv(r[1].X,recv[1],MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
	}

    if(p->nb3>=0)
    {
        MPI_Isend(s[2].X,send[2],MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
        MPI_Irecv(r[2].X,recv[2],MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
	}

    if(p->nb4>=0)
    {
        MPI_Isend(s[3].X,send[3],MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
        MPI_Irecv(r[3].X,recv[3],MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
	}

    if(p->nb5>=0)
    {
        MPI_Isend(s[4].X,send[4],MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
        MPI_Irecv(r[4].X,recv[4],MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
	}

    if(p->nb6>=0)
    {
        MPI_Isend(s[5].X,send[5],MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
        MPI_Irecv(r[5].X,recv[5],MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
	}

    gcwait(p);
    

    // Transfer Y
    if(p->nb1>=0)
    {
        MPI_Isend(s[0].Y,send[0],MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
        MPI_Irecv(r[0].Y,recv[0],MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
	}

    if(p->nb2>=0)
    {
        MPI_Isend(s[1].Y,send[1],MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
        MPI_Irecv(r[1].Y,recv[1],MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
	}

    if(p->nb3>=0)
    {
        MPI_Isend(s[2].Y,send[2],MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
        MPI_Irecv(r[2].Y,recv[2],MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
	}

    if(p->nb4>=0)
    {
        MPI_Isend(s[3].Y,send[3],MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
        MPI_Irecv(r[3].Y,recv[3],MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
	}

    if(p->nb5>=0)
    {
        MPI_Isend(s[4].Y,send[4],MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
        MPI_Irecv(r[4].Y,recv[4],MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
	}

    if(p->nb6>=0)
    {
        MPI_Isend(s[5].Y,send[5],MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
        MPI_Irecv(r[5].Y,recv[5],MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
	}

    gcwait(p);


    // Transfer Z
    if(p->nb1>=0)
    {
        MPI_Isend(s[0].Z,send[0],MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
        MPI_Irecv(r[0].Z,recv[0],MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
	}

    if(p->nb2>=0)
    {
        MPI_Isend(s[1].Z,send[1],MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
        MPI_Irecv(r[1].Z,recv[1],MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
	}

    if(p->nb3>=0)
    {
        MPI_Isend(s[2].Z,send[2],MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
        MPI_Irecv(r[2].Z,recv[2],MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
	}

    if(p->nb4>=0)
    {
        MPI_Isend(s[3].Z,send[3],MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
        MPI_Irecv(r[3].Z,recv[3],MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
	}

    if(p->nb5>=0)
    {
        MPI_Isend(s[4].Z,send[4],MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
        MPI_Irecv(r[4].Z,recv[4],MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
	}

    if(p->nb6>=0)
    {
        MPI_Isend(s[5].Z,send[5],MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
        MPI_Irecv(r[5].Z,recv[5],MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
	}

    gcwait(p);

    // Particle data
    if(s->entries>s->tracers_obj::entries)
    {
        // Transfer U
        if(p->nb1>=0)
        {
            MPI_Isend(s[0].U,send[0],MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
            MPI_Irecv(r[0].U,recv[0],MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
        }

        if(p->nb2>=0)
        {
            MPI_Isend(s[1].U,send[1],MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
            MPI_Irecv(r[1].U,recv[1],MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
        }

        if(p->nb3>=0)
        {
            MPI_Isend(s[2].U,send[2],MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
            MPI_Irecv(r[2].U,recv[2],MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
        }

        if(p->nb4>=0)
        {
            MPI_Isend(s[3].U,send[3],MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
            MPI_Irecv(r[3].U,recv[3],MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
        }

        if(p->nb5>=0)
        {
            MPI_Isend(s[4].U,send[4],MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
            MPI_Irecv(r[4].U,recv[4],MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
        }

        if(p->nb6>=0)
        {
            MPI_Isend(s[5].U,send[5],MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
            MPI_Irecv(r[5].U,recv[5],MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
        }

        gcwait(p);


        // Transfer V
        if(p->nb1>=0)
        {
            MPI_Isend(s[0].V,send[0],MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
            MPI_Irecv(r[0].V,recv[0],MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
        }

        if(p->nb2>=0)
        {
            MPI_Isend(s[1].V,send[1],MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
            MPI_Irecv(r[1].V,recv[1],MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
        }

        if(p->nb3>=0)
        {
            MPI_Isend(s[2].V,send[2],MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
            MPI_Irecv(r[2].V,recv[2],MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
        }

        if(p->nb4>=0)
        {
            MPI_Isend(s[3].V,send[3],MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
            MPI_Irecv(r[3].V,recv[3],MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
        }

        if(p->nb5>=0)
        {
            MPI_Isend(s[4].V,send[4],MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
            MPI_Irecv(r[4].V,recv[4],MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
        }

        if(p->nb6>=0)
        {
            MPI_Isend(s[5].V,send[5],MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
            MPI_Irecv(r[5].V,recv[5],MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
        }

        gcwait(p);


        // Transfer W
        if(p->nb1>=0)
        {
            MPI_Isend(s[0].W,send[0],MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
            MPI_Irecv(r[0].W,recv[0],MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
        }

        if(p->nb2>=0)
        {
            MPI_Isend(s[1].W,send[1],MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
            MPI_Irecv(r[1].W,recv[1],MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
        }

        if(p->nb3>=0)
        {
            MPI_Isend(s[2].W,send[2],MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
            MPI_Irecv(r[2].W,recv[2],MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
        }

        if(p->nb4>=0)
        {
            MPI_Isend(s[3].W,send[3],MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
            MPI_Irecv(r[3].W,recv[3],MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
        }

        if(p->nb5>=0)
        {
            MPI_Isend(s[4].W,send[4],MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
            MPI_Irecv(r[4].W,recv[4],MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
        }

        if(p->nb6>=0)
        {
            MPI_Isend(s[5].W,send[5],MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
            MPI_Irecv(r[5].W,recv[5],MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
        }

        gcwait(p);


        // Transfer U
        if(p->nb1>=0)
        {
            MPI_Isend(s[0].PackingFactor,send[0],MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
            MPI_Irecv(r[0].PackingFactor,recv[0],MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
        }

        if(p->nb2>=0)
        {
            MPI_Isend(s[1].PackingFactor,send[1],MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
            MPI_Irecv(r[1].PackingFactor,recv[1],MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
        }

        if(p->nb3>=0)
        {
            MPI_Isend(s[2].PackingFactor,send[2],MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
            MPI_Irecv(r[2].PackingFactor,recv[2],MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
        }

        if(p->nb4>=0)
        {
            MPI_Isend(s[3].PackingFactor,send[3],MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
            MPI_Irecv(r[3].PackingFactor,recv[3],MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
        }

        if(p->nb5>=0)
        {
            MPI_Isend(s[4].PackingFactor,send[4],MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
            MPI_Irecv(r[4].PackingFactor,recv[4],MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
        }

        if(p->nb6>=0)
        {
            MPI_Isend(s[5].PackingFactor,send[5],MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
            MPI_Irecv(r[5].PackingFactor,recv[5],MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
        }

        gcwait(p);
    }
}