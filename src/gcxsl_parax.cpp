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
#include"slice.h"

void ghostcell::gcslparax(lexer* p,slice& f,int gcv)
{
    paramargin=margin;

//  FILL SEND
    count=0;
    for(q=0;q<p->gcslpara1_count;++q)
    {
    i=p->gcslpara1[q][0];
    j=p->gcslpara1[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        send1[count]=f(i+n,j);
        ++count;
        }
    }

    count=0;
    for(q=0;q<p->gcslpara3_count;++q)
    {
    i=p->gcslpara3[q][0];
    j=p->gcslpara3[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        send3[count]=f(i,j+n);
        ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcslpara4_count;++q)
	{
    i=p->gcslpara4[q][0];
    j=p->gcslpara4[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        send4[count]=f(i-n,j);
        ++count;
        }
	}
    
    count=0;
	for(q=0;q<p->gcslpara2_count;++q)
	{
    i=p->gcslpara2[q][0];
    j=p->gcslpara2[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        send2[count]=f(i,j-n);
        ++count;
        }
	}

//  SEND / RECEIVE

    if(p->gcslpara1_count>0)
    {
	MPI_Isend(send1,p->gcslpara1_count*paramargin,MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(recv1,p->gcslpara1_count*paramargin,MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcslpara4_count>0)
    {
	MPI_Isend(send4,p->gcslpara4_count*paramargin,MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(recv4,p->gcslpara4_count*paramargin,MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
    }

    if(p->gcslpara3_count>0)
    {
	MPI_Isend(send3,p->gcslpara3_count*paramargin,MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(recv3,p->gcslpara3_count*paramargin,MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcslpara2_count>0)
    {
	MPI_Isend(send2,p->gcslpara2_count*paramargin,MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(recv2,p->gcslpara2_count*paramargin,MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
    }


//  WAIT

    gcslwait(p);

//  FILL RECEIVE

    count=0;
    for(q=0;q<p->gcslpara1_count;++q)
    {
    i=p->gcslpara1[q][0];
    j=p->gcslpara1[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        f(i-n-1,j)=recv1[count];
        ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcslpara3_count;++q)
	{
    i=p->gcslpara3[q][0];
    j=p->gcslpara3[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        f(i,j-n-1)=recv3[count];
        ++count;
        }
	}


    count=0;
	for(q=0;q<p->gcslpara4_count;++q)
	{
    i=p->gcslpara4[q][0];
    j=p->gcslpara4[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        f(i+n+1,j)=recv4[count];
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcslpara2_count;++q)
	{
    i=p->gcslpara2[q][0];
    j=p->gcslpara2[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        f(i,j+n+1)=recv2[count];
        ++count;
        }
	}
}

