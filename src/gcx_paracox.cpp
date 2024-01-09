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
#include"fdm.h"

void ghostcell::gcparacox(lexer* p,field& f,int gcv)
{
    pip=4;
//  FILL SEND
    for(q=0;q<p->gcparaco1_count;++q)
    {
    i=p->gcparaco1[q][0];
    j=p->gcparaco1[q][1];
    k=p->gcparaco1[q][2];
	send1[q]=f(i,j,k);
    }

    for(q=0;q<p->gcparaco3_count;++q)
    {
    i=p->gcparaco3[q][0];
    j=p->gcparaco3[q][1];
    k=p->gcparaco3[q][2];
	send3[q]=f(i,j,k);
    }

	for(q=0;q<p->gcparaco5_count;++q)
	{
    i=p->gcparaco5[q][0];
    j=p->gcparaco5[q][1];
    k=p->gcparaco5[q][2];
	send5[q]=f(i,j,k);
	}

	for(q=0;q<p->gcparaco4_count;++q)
	{
    i=p->gcparaco4[q][0];
    j=p->gcparaco4[q][1];
    k=p->gcparaco4[q][2];
	send4[q]=f(i,j,k);
	}

	for(q=0;q<p->gcparaco2_count;++q)
	{
    i=p->gcparaco2[q][0];
    j=p->gcparaco2[q][1];
    k=p->gcparaco2[q][2];
	send2[q]=f(i,j,k);
	}

	for(q=0;q<p->gcparaco6_count;++q)
	{
    i=p->gcparaco6[q][0];
    j=p->gcparaco6[q][1];
    k=p->gcparaco6[q][2];
	send6[q]=f(i,j,k);
	}


//  SEND / RECEIVE

    if(p->gcparaco1_count>0)
    {
	MPI_Isend(send1,p->gcparaco1_count,MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(recv1,p->gcparaco1_count,MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcparaco4_count>0)
    {
	MPI_Isend(send4,p->gcparaco4_count,MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(recv4,p->gcparaco4_count,MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
    }

    if(p->gcparaco3_count>0)
    {
	MPI_Isend(send3,p->gcparaco3_count,MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(recv3,p->gcparaco3_count,MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcparaco2_count>0)
    {
	MPI_Isend(send2,p->gcparaco2_count,MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(recv2,p->gcparaco2_count,MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
    }

    if(p->gcparaco5_count>0)
    {
	MPI_Isend(send5,p->gcparaco5_count,MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
	MPI_Irecv(recv5,p->gcparaco5_count,MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
    }

    if(p->gcparaco6_count>0)
    {
	MPI_Isend(send6,p->gcparaco6_count,MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
	MPI_Irecv(recv6,p->gcparaco6_count,MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
    }

//  WAIT
    gcwait(p);

//  FILL RECEIVE
    //pip=4;
    for(q=0;q<p->gcparaco1_count;++q)
    {
    i=p->gcparaco1[q][0];
    j=p->gcparaco1[q][1];
    k=p->gcparaco1[q][2];
	f(i-1,j,k)=recv1[q];
    }

	for(q=0;q<p->gcparaco3_count;++q)
	{
    i=p->gcparaco3[q][0];
    j=p->gcparaco3[q][1];
    k=p->gcparaco3[q][2];
	f(i,j-1,k)=recv3[q];
	}

    for(q=0;q<p->gcparaco5_count;++q)
    {
    i=p->gcparaco5[q][0];
    j=p->gcparaco5[q][1];
    k=p->gcparaco5[q][2];
	f(i,j,k-1)=recv5[q];
    }

	for(q=0;q<p->gcparaco4_count;++q)
	{
    i=p->gcparaco4[q][0];
    j=p->gcparaco4[q][1];
    k=p->gcparaco4[q][2];
	f(i+1,j,k)=recv4[q];
	}

	for(q=0;q<p->gcparaco2_count;++q)
	{
    i=p->gcparaco2[q][0];
    j=p->gcparaco2[q][1];
    k=p->gcparaco2[q][2];
	f(i,j+1,k)=recv2[q];
	}

	for(q=0;q<p->gcparaco6_count;++q)
	{
    i=p->gcparaco6[q][0];
    j=p->gcparaco6[q][1];
    k=p->gcparaco6[q][2];
   	f(i,j,k+1)=recv6[q];
	}
    
    pip=0;
}

