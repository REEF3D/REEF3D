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

void ghostcell::flagx7(lexer* p,int *f)
{
    paramargin=margin;

//  FILL SEND
    count=0;
    for(q=0;q<p->gcx7_count[0];++q)
    {
    i=p->gcx7[0][q][0];
    j=p->gcx7[0][q][1];
    k=p->gcx7[0][q][2];
        

        isend1[count] = f[FIJK];  
        ++count;
        isend1[count] = f[FIp1JK];  
        ++count;
        isend1[count] = f[FIp2JK];  
        ++count;

    }

    count=0;
    for(q=0;q<p->gcx7_count[2];++q)
    {
    i=p->gcx7[2][q][0];
    j=p->gcx7[2][q][1];
    k=p->gcx7[2][q][2];
        
        isend3[count] = f[FIJK];  
        ++count;
        isend3[count] = f[FIJp1K];  
        ++count;
        isend3[count] = f[FIJp2K];  
        ++count;
    }

    count=0;
	for(q=0;q<p->gcx7_count[3];++q)
	{
    i=p->gcx7[3][q][0];
    j=p->gcx7[3][q][1];
    k=p->gcx7[3][q][2];
        
        isend4[count] = f[FIJK];  
        ++count;
        isend4[count] = f[FIm1JK];  
        ++count;
        isend4[count] = f[FIm2JK];  
        ++count;
	}

    count=0;
	for(q=0;q<p->gcx7_count[1];++q)
	{
    i=p->gcx7[1][q][0];
    j=p->gcx7[1][q][1];
    k=p->gcx7[1][q][2];
        
        isend2[count] = f[FIJK];  
        ++count;
        isend2[count] = f[FIJm1K];  
        ++count;
        isend2[count] = f[FIJm2K];  
        ++count;
	}


//  SEND / RECEIVE

    if(p->gcx7_count[0]>0)
    {
	MPI_Isend(isend1,p->gcx7_count[0]*paramargin,MPI_INT,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(irecv1,p->gcx7_count[0]*paramargin,MPI_INT,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcx7_count[3]>0)
    {
	MPI_Isend(isend4,p->gcx7_count[3]*paramargin,MPI_INT,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(irecv4,p->gcx7_count[3]*paramargin,MPI_INT,p->nb4,tag1,mpi_comm,&rreq4);
    }

    if(p->gcx7_count[2]>0)
    {
	MPI_Isend(isend3,p->gcx7_count[2]*paramargin,MPI_INT,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(irecv3,p->gcx7_count[2]*paramargin,MPI_INT,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcx7_count[1]>0)
    {
	MPI_Isend(isend2,p->gcx7_count[1]*paramargin,MPI_INT,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(irecv2,p->gcx7_count[1]*paramargin,MPI_INT,p->nb2,tag3,mpi_comm,&rreq2);
    }



//  WAIT

    gcwait7(p);

//  FILL RECEIVE

    count=0;
    for(q=0;q<p->gcx7_count[0];++q)
    {
    i=p->gcx7[0][q][0];
    j=p->gcx7[0][q][1];
    k=p->gcx7[0][q][2];
        

        f[FIm1JK] = irecv1[count];
        ++count;
        f[FIm2JK] = irecv1[count];
        ++count;
        f[FIm3JK] = irecv1[count];
        ++count;
        
    }

    count=0;
	for(q=0;q<p->gcx7_count[2];++q)
	{
    i=p->gcx7[2][q][0];
    j=p->gcx7[2][q][1];
    k=p->gcx7[2][q][2];
        
        f[FIJm1K] = irecv3[count];
        ++count;
        f[FIJm2K] = irecv3[count];
        ++count;
        f[FIJm3K] = irecv3[count];
        ++count;
	}

    count=0;
	for(q=0;q<p->gcx7_count[3];++q)
	{
    i=p->gcx7[3][q][0];
    j=p->gcx7[3][q][1];
    k=p->gcx7[3][q][2];
        
        f[FIp1JK] = irecv4[count];
        ++count;
        f[FIp2JK] = irecv4[count];
        ++count;
        f[FIp3JK] = irecv4[count];
        ++count;
	}

    count=0;
	for(q=0;q<p->gcx7_count[1];++q)
	{
    i=p->gcx7[1][q][0];
    j=p->gcx7[1][q][1];
    k=p->gcx7[1][q][2];
        
        f[FIJp1K] = irecv2[count];
        ++count;
        f[FIJp2K] = irecv2[count];
        ++count;
        f[FIJp3K] = irecv2[count];
        ++count;
	}

}

