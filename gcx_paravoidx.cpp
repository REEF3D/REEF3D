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
#include"fdm.h"

void ghostcell::gcparavoidx(lexer* p,field& f,int gcv)
{
    paramargin=margin;

//  FILL SEND
    count=0;
    for(q=0;q<p->gcparavoid1_count;++q)
    {
    i=p->gcparavoid1[q][0];
    j=p->gcparavoid1[q][1];
    k=p->gcparavoid1[q][2];

        for(n=0;n<paramargin;++n)
        {
        send1[count]=f(i+n,j,k);
        ++count;
        }
    }

    count=0;
    for(q=0;q<p->gcparavoid3_count;++q)
    {
    i=p->gcparavoid3[q][0];
    j=p->gcparavoid3[q][1];
    k=p->gcparavoid3[q][2];

        for(n=0;n<paramargin;++n)
        {
        send3[count]=f(i,j+n,k);
        ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcparavoid5_count;++q)
	{
    i=p->gcparavoid5[q][0];
    j=p->gcparavoid5[q][1];
    k=p->gcparavoid5[q][2];

        for(n=0;n<paramargin;++n)
        {
        send5[count]=f(i,j,k+n);
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcparavoid4_count;++q)
	{
    i=p->gcparavoid4[q][0];
    j=p->gcparavoid4[q][1];
    k=p->gcparavoid4[q][2];

        for(n=0;n<paramargin;++n)
        {
        send4[count]=f(i-n,j,k);
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcparavoid2_count;++q)
	{
    i=p->gcparavoid2[q][0];
    j=p->gcparavoid2[q][1];
    k=p->gcparavoid2[q][2];

        for(n=0;n<paramargin;++n)
        {
        send2[count]=f(i,j-n,k);
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcparavoid6_count;++q)
	{
	i=p->gcparavoid6[q][0];
    j=p->gcparavoid6[q][1];
    k=p->gcparavoid6[q][2];

        for(n=0;n<paramargin;++n)
        {
        send6[count]=f(i,j,k-n);
        ++count;
        }
	}


//  SEND / RECEIVE

    if(p->gcparavoid1_count>0)
    {
	MPI_Isend(send1,p->gcparavoid1_count*paramargin,MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(recv1,p->gcparavoid1_count*paramargin,MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcparavoid4_count>0)
    {
	MPI_Isend(send4,p->gcparavoid4_count*paramargin,MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(recv4,p->gcparavoid4_count*paramargin,MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
    }

    if(p->gcparavoid3_count>0)
    {
	MPI_Isend(send3,p->gcparavoid3_count*paramargin,MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(recv3,p->gcparavoid3_count*paramargin,MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcparavoid2_count>0)
    {
	MPI_Isend(send2,p->gcparavoid2_count*paramargin,MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(recv2,p->gcparavoid2_count*paramargin,MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
    }

    if(p->gcparavoid5_count>0)
    {
	MPI_Isend(send5,p->gcparavoid5_count*paramargin,MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
	MPI_Irecv(recv5,p->gcparavoid5_count*paramargin,MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
    }

    if(p->gcparavoid6_count>0)
    {
	MPI_Isend(send6,p->gcparavoid6_count*paramargin,MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
	MPI_Irecv(recv6,p->gcparavoid6_count*paramargin,MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
    }

//  WAIT

    gcwait(p);

//  FILL RECEIVE

    count=0;
    for(q=0;q<p->gcparavoid1_count;++q)
    {
    i=p->gcparavoid1[q][0];
    j=p->gcparavoid1[q][1];
    k=p->gcparavoid1[q][2];

        for(n=0;n<paramargin;++n)
        {
        if(p->gcparavoid1[q][3]==0)
        f(i-n-1,j,k)=recv1[count];
        ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcparavoid3_count;++q)
	{
    i=p->gcparavoid3[q][0];
    j=p->gcparavoid3[q][1];
    k=p->gcparavoid3[q][2];

        for(n=0;n<paramargin;++n)
        {
        if(p->gcparavoid3[q][3]==0)
        f(i,j-n-1,k)=recv3[count];
        ++count;
        }
	}

	count=0;
    for(q=0;q<p->gcparavoid5_count;++q)
    {
    i=p->gcparavoid5[q][0];
    j=p->gcparavoid5[q][1];
    k=p->gcparavoid5[q][2];

        for(n=0;n<paramargin;++n)
        {
        if(p->gcparavoid5[q][3]==0)
        f(i,j,k-n-1)=recv5[count];
        ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcparavoid4_count;++q)
	{
    i=p->gcparavoid4[q][0];
    j=p->gcparavoid4[q][1];
    k=p->gcparavoid4[q][2];

        for(n=0;n<paramargin;++n)
        {
        if(p->gcparavoid4[q][3]==0)
        f(i+n+1,j,k)=recv4[count];
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcparavoid2_count;++q)
	{
    i=p->gcparavoid2[q][0];
    j=p->gcparavoid2[q][1];
    k=p->gcparavoid2[q][2];

        for(n=0;n<paramargin;++n)
        {
        if(p->gcparavoid2[q][3]==0)
        f(i,j+n+1,k)=recv2[count];
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcparavoid6_count;++q)
	{
    i=p->gcparavoid6[q][0];
    j=p->gcparavoid6[q][1];
    k=p->gcparavoid6[q][2];

        for(n=0;n<paramargin;++n)
        {
        if(p->gcparavoid6[q][3]==0)
        f(i,j,k+n+1)=recv6[count];
        ++count;
        }
	}

}

