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
#include"fieldint4.h"

void ghostcell::gcparax_test(lexer* p,int gcv)
{
    int testmargin=2;
    
    fieldint4 f(p);

//  FILL SEND
    count=0;
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];


        if(p->gcpara1[q][2+gcv]==1)
        {
        isend1[count]=j;
        ++count;
        isend1[count]=k;
        ++count;
        }
    }

    count=0;
    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];

        for(n=0;n<testmargin;++n)
        {
        if(p->gcpara3[q][2+gcv]==1)
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

        for(n=0;n<testmargin;++n)
        {
        if(p->gcpara5[q][2+gcv]==1)
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
        {
        isend4[count]=j;
        ++count;
        
        isend4[count]=k;
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];

        for(n=0;n<testmargin;++n)
        {
        if(p->gcpara2[q][2+gcv]==1)
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

        for(n=0;n<testmargin;++n)
        {
        if(p->gcpara6[q][2+gcv]==1)
        isend6[count]=f(i,j,k-n);
        ++count;
        }
	}


//  SEND / RECEIVE

    if(p->gcpara1_count>0)
    {
	MPI_Isend(isend1,p->gcpara1_count*testmargin,MPI_INT,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(irecv1,p->gcpara1_count*testmargin,MPI_INT,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcpara4_count>0)
    {
	MPI_Isend(isend4,p->gcpara4_count*testmargin,MPI_INT,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(irecv4,p->gcpara4_count*testmargin,MPI_INT,p->nb4,tag1,mpi_comm,&rreq4);
    }

    if(p->gcpara3_count>0)
    {
	MPI_Isend(isend3,p->gcpara3_count*testmargin,MPI_INT,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(irecv3,p->gcpara3_count*testmargin,MPI_INT,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcpara2_count>0)
    {
	MPI_Isend(isend2,p->gcpara2_count*testmargin,MPI_INT,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(irecv2,p->gcpara2_count*testmargin,MPI_INT,p->nb2,tag3,mpi_comm,&rreq2);
    }

    if(p->gcpara5_count>0)
    {
	MPI_Isend(isend5,p->gcpara5_count*testmargin,MPI_INT,p->nb5,tag5,mpi_comm,&sreq5);
	MPI_Irecv(irecv5,p->gcpara5_count*testmargin,MPI_INT,p->nb5,tag6,mpi_comm,&rreq5);
    }

    if(p->gcpara6_count>0)
    {
	MPI_Isend(isend6,p->gcpara6_count*testmargin,MPI_INT,p->nb6,tag6,mpi_comm,&sreq6);
	MPI_Irecv(irecv6,p->gcpara6_count*testmargin,MPI_INT,p->nb6,tag5,mpi_comm,&rreq6);
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

        for(n=0;n<testmargin;++n)
        {
        if(p->gcpara1[q][2+gcv]==1)
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

        for(n=0;n<testmargin;++n)
        {
        if(p->gcpara3[q][2+gcv]==1)
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

        for(n=0;n<testmargin;++n)
        {
        if(p->gcpara5[q][2+gcv]==1)
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
        for(n=0;n<testmargin;++n)
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

        for(n=0;n<testmargin;++n)
        {
        if(p->gcpara2[q][2+gcv]==1)
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

        for(n=0;n<testmargin;++n)
        {
        if(p->gcpara6[q][2+gcv]==1)
        f(i,j,k+n+1)=irecv6[count];
        ++count;
        }
	}
    
    int coin=0;
    count=0;
    
    for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
    
         if(p->gcpara4[q][2+gcv]==1)
         {
         if(isend4[count]!=f(i+1,j,k) || isend4[count+1]!=f(i+2,j,k)) 
         cout<<p->mpirank<<" PARAX: send: "<<isend4[count]<<" "<<isend4[count+1]<<" recv: "<<irecv4[count]<<" "<<irecv4[count+1]<<" f: "<<f(i+1,j,k)<<" "<<f(i+2,j,k)<<endl;
        
        if(isend4[count]==irecv4[count] && isend4[count+1]==irecv4[count+1])
        ++coin;
        
         count+=2;
         }
    }
    

    //cout<<p->mpirank<<" COIN: "<<coin<<endl;
    
    
    
    /*
    int count1,count2,count3,count4;
    
    count1=count2=count3=count4=0;
    if(p->mpirank==2)
    {
        for(n=0;n<p->gcpara1_count;++n)
        {
        if(p->gcpara1[n][3]==1)
        ++count1;
        
        if(p->gcpara1[n][4]==1)
        ++count2;
        
        if(p->gcpara1[n][5]==1)
        ++count3;
        
        if(p->gcpara1[n][6]==1)
        ++count4;
        }
     cout<<p->mpirank<<" GCX_COUNT: "<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<endl;           
    }
    
    count1=count2=count3=count4=0;
    if(p->mpirank==1)
    {
        for(n=0;n<p->gcpara4_count;++n)
        {
        if(p->gcpara4[n][3]==1)
        ++count1;
        
        if(p->gcpara4[n][4]==1)
        ++count2;
        
        if(p->gcpara4[n][5]==1)
        ++count3;
        
        if(p->gcpara4[n][6]==1)
        ++count4;
        }
    cout<<p->mpirank<<" GCX_COUNT: "<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<endl;
    }*/

}

