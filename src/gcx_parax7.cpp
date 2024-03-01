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

void ghostcell::gcparax7(lexer* p,double *&f,int gcv)
{
    paramargin=2;//margin;

//  FILL SEND
    count=0;
    for(q=0;q<p->gcx7_count[0];++q)
    {
    i=p->gcx7[0][q][0];
    j=p->gcx7[0][q][1];
    k=p->gcx7[0][q][2];
    
    //cout<<p->mpirank<<" GCX7: "<<i<<" "<<j<<" "<<k<<" | "<<q<<endl;

        send1[count] = f[FIJK];  
        ++count;
        send1[count] = f[FIp1JK];  
        ++count;
        //send1[count] = f[FIp2JK];  
        //++count;

    }

    count=0;
    for(q=0;q<p->gcx7_count[2];++q)
    {
    i=p->gcx7[2][q][0];
    j=p->gcx7[2][q][1];
    k=p->gcx7[2][q][2];
        
        send3[count] = f[FIJK];  
        ++count;
        send3[count] = f[FIJp1K];  
        ++count;
        //send3[count] = f[FIJp2K];  
        //++count;
    }

    count=0;
	for(q=0;q<p->gcx7_count[3];++q)
	{
    i=p->gcx7[3][q][0];
    j=p->gcx7[3][q][1];
    k=p->gcx7[3][q][2];
        
        send4[count] = f[FIJK];  
        ++count;
        send4[count] = f[FIm1JK];  
        ++count;
        //send4[count] = f[FIm2JK];  
        //++count;
	}

    count=0;
	for(q=0;q<p->gcx7_count[1];++q)
	{
    i=p->gcx7[1][q][0];
    j=p->gcx7[1][q][1];
    k=p->gcx7[1][q][2];
        
        send2[count] = f[FIJK];  
        ++count;
        send2[count] = f[FIJm1K];  
        ++count;
        //send2[count] = f[FIJm2K];  
        //++count;
	}


//  SEND / RECEIVE

    if(p->gcx7_count[0]>0)
    {
	MPI_Isend(send1,p->gcx7_count[0]*paramargin,MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(recv1,p->gcx7_count[0]*paramargin,MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcx7_count[3]>0)
    {
	MPI_Isend(send4,p->gcx7_count[3]*paramargin,MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(recv4,p->gcx7_count[3]*paramargin,MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
    }

    if(p->gcx7_count[2]>0)
    {
	MPI_Isend(send3,p->gcx7_count[2]*paramargin,MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(recv3,p->gcx7_count[2]*paramargin,MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcx7_count[1]>0)
    {
	MPI_Isend(send2,p->gcx7_count[1]*paramargin,MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(recv2,p->gcx7_count[1]*paramargin,MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
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
        
        f[FIm1JK] = recv1[count];
        ++count;
        f[FIm2JK] = recv1[count];
        ++count;
        //f[FIm3JK] = recv1[count];
        //++count;
    }

    count=0;
	for(q=0;q<p->gcx7_count[2];++q)
	{
    i=p->gcx7[2][q][0];
    j=p->gcx7[2][q][1];
    k=p->gcx7[2][q][2];
        
        f[FIJm1K] = recv3[count];
        ++count;
        f[FIJm2K] = recv3[count];
        ++count;
        //f[FIJm3K] = recv3[count];
        //++count;
	}

    count=0;
	for(q=0;q<p->gcx7_count[3];++q)
	{
    i=p->gcx7[3][q][0];
    j=p->gcx7[3][q][1];
    k=p->gcx7[3][q][2];
        
        f[FIp1JK] = recv4[count];
        ++count;
        f[FIp2JK] = recv4[count];
        ++count;
        //f[FIp3JK] = recv4[count];
        //++count;
	}

    count=0;
	for(q=0;q<p->gcx7_count[1];++q)
	{
    i=p->gcx7[1][q][0];
    j=p->gcx7[1][q][1];
    k=p->gcx7[1][q][2];
        
        f[FIJp1K] = recv2[count];
        ++count;
        f[FIJp2K] = recv2[count];
        ++count;
        //f[FIJp3K] = recv2[count];
        //++count;
	}

}

