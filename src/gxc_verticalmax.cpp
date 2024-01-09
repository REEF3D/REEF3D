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
#include"time.h"

void ghostcell::verticalmax(lexer *p, fdm* a, double **vmax)
{
		
//  FILL SEND
   
    count=0;
	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];


        if(p->gcpara5[q][5]==1)
		{
        send5[count]=vmax[i][j];
        ++count;
		}

	}

    count=0;
	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];

        if(p->gcpara6[q][5]==1)
		{
        send6[count]=vmax[i][j];
        ++count;
        }
	}


//  SEND / RECEIVE
    if(p->gcpara5_count>0)
    {
	MPI_Isend(send5,p->gcpara5_count,MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
	MPI_Irecv(recv5,p->gcpara5_count,MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
    }

    if(p->gcpara6_count>0)
    {
	MPI_Isend(send6,p->gcpara6_count,MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
	MPI_Irecv(recv6,p->gcpara6_count,MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
    }

//  WAIT

    if(p->gcpara5_count>0)
    {
    MPI_Wait(&sreq5,&status);
	MPI_Wait(&rreq5,&status);
    }

    if(p->gcpara6_count>0)
    {
    MPI_Wait(&sreq6,&status);
	MPI_Wait(&rreq6,&status);
    }

//  FILL RECEIVE

	count=0;
    for(q=0;q<p->gcpara5_count;++q)
    {
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];

        if(p->gcpara5[q][5]==1)
		{
        vmax[i][j]=MAX(vmax[i][j],recv5[count]);
        ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcpara6_count;++q)
	{
    i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];

        if(p->gcpara6[q][5]==1)
		{
        vmax[i][j]=MAX(vmax[i][j],recv6[count]);
        ++count;
        }
	}
}
