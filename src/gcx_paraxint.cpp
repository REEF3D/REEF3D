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

void ghostcell::gcparaxint(lexer* p,fieldint& f,int gcv)
{
    paramargin=margin;

//  FILL SEND
    count=0;
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->gcpara1[q][2+gcv]==1)
        for(n=0;n<paramargin;++n)
        {
        isend1[count]=f(i+n,j,k);
        ++count;
        }
    }

    count=0;
    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->gcpara3[q][2+gcv]==1)
        for(n=0;n<paramargin;++n)
        {
        isend3[count]=f(i,j+n,k);
        ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->gcpara5[q][2+gcv]==1)
        for(n=0;n<paramargin;++n)
        {
        isend5[count]=f(i,j,k+n);
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->gcpara4[q][2+gcv]==1)
        for(n=0;n<paramargin;++n)
        {
        isend4[count]=f(i-n,j,k);
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
        if(p->gcpara2[q][2+gcv]==1)
        for(n=0;n<paramargin;++n)
        {
        isend2[count]=f(i,j-n,k);
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
        if(p->gcpara6[q][2+gcv]==1)
        for(n=0;n<paramargin;++n)
        {
        isend6[count]=f(i,j,k-n);
        ++count;
        }
	}


//  SEND / RECEIVE

    if(p->gcpara1_count>0)
    {
	MPI_Isend(isend1,p->gcpara1_count*paramargin,MPI_INT,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(irecv1,p->gcpara1_count*paramargin,MPI_INT,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcpara4_count>0)
    {
	MPI_Isend(isend4,p->gcpara4_count*paramargin,MPI_INT,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(irecv4,p->gcpara4_count*paramargin,MPI_INT,p->nb4,tag1,mpi_comm,&rreq4);
    }

    if(p->gcpara3_count>0)
    {
	MPI_Isend(isend3,p->gcpara3_count*paramargin,MPI_INT,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(irecv3,p->gcpara3_count*paramargin,MPI_INT,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcpara2_count>0)
    {
	MPI_Isend(isend2,p->gcpara2_count*paramargin,MPI_INT,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(irecv2,p->gcpara2_count*paramargin,MPI_INT,p->nb2,tag3,mpi_comm,&rreq2);
    }

    if(p->gcpara5_count>0)
    {
	MPI_Isend(isend5,p->gcpara5_count*paramargin,MPI_INT,p->nb5,tag5,mpi_comm,&sreq5);
	MPI_Irecv(irecv5,p->gcpara5_count*paramargin,MPI_INT,p->nb5,tag6,mpi_comm,&rreq5);
    }

    if(p->gcpara6_count>0)
    {
	MPI_Isend(isend6,p->gcpara6_count*paramargin,MPI_INT,p->nb6,tag6,mpi_comm,&sreq6);
	MPI_Irecv(irecv6,p->gcpara6_count*paramargin,MPI_INT,p->nb6,tag5,mpi_comm,&rreq6);
    }

//  WAIT

    gcwait(p);

//  FILL RECEIVE

    count=0;
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->gcpara1[q][2+gcv]==1)
        for(n=0;n<paramargin;++n)
        {
        f(i-n-1,j,k)=irecv1[count];
        ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcpara3_count;++q)
	{
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->gcpara3[q][2+gcv]==1)
        for(n=0;n<paramargin;++n)
        {
        f(i,j-n-1,k)=irecv3[count];
        ++count;
        }
	}

	count=0;
    for(q=0;q<p->gcpara5_count;++q)
    {
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->gcpara5[q][2+gcv]==1)
        for(n=0;n<paramargin;++n)
        {
        f(i,j,k-n-1)=irecv5[count];
        ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->gcpara4[q][2+gcv]==1)
        for(n=0;n<paramargin;++n)
        {
        f(i+n+1,j,k)=irecv4[count];
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
        if(p->gcpara2[q][2+gcv]==1)
        for(n=0;n<paramargin;++n)
        {
        f(i,j+n+1,k)=irecv2[count];
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara6_count;++q)
	{
    i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
        if(p->gcpara6[q][2+gcv]==1)
        for(n=0;n<paramargin;++n)
        {
        f(i,j,k+n+1)=irecv6[count];
        ++count;
        }
	}

}

