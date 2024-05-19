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

void ghostcell::gcperiodicx(lexer* p,field& f,int gcv)
{
    paramargin=margin;
    
    int aa,bb,cc;
    aa=bb=cc=0;
    /*
    if(gcv==1)
    aa=1;
    
    if(gcv==2)
    bb=1;
    
    if(gcv==3)
    cc=1;*/
    
    //cout<<p->mpirank<<"  p->periodicX1: "<<p->periodicX1<<" "<<p->gcpara1_count<<endl;
    //cout<<p->mpirank<<"  p->periodicX4: "<<p->periodicX4<<" "<<p->gcpara4_count<<endl;
    
    //cout<<p->mpirank<<" nb1: "<<p->nb1<<" nb4: "<<p->nb4<<endl;

//  FILL SEND
    count=0;
    for(q=p->periodicX1;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0]-1;
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->gcpara1[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
        send1[count]=f(i+n+1,j,k);
        ++count;
        }
    }

    count=0;
    for(q=p->periodicX3;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1]-1;
    k=p->gcpara3[q][2];
        
        if(p->gcpara3[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
        send3[count]=f(i,j+n+1,k);
        ++count;
        }
    }

    count=0;
	for(q=p->periodicX5;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2]-1;
        
        if(p->gcpara5[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
        send5[count]=f(i,j,k+n+1);
        ++count;
        }
	}

    count=0;
	for(q=p->periodicX4;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0]+1-aa;
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->gcpara4[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
        send4[count]=f(i-n-1,j,k);
        ++count;
        }
	}

    count=0;
	for(q=p->periodicX2;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1]+1-bb;
    k=p->gcpara2[q][2];
        
        if(p->gcpara2[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
        send2[count]=f(i,j-n-1,k);
        ++count;
        }
	}

    count=0;
	for(q=p->periodicX6;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2]+1-cc;
        
        if(p->gcpara6[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
        send6[count]=f(i,j,k-n-1);
        ++count;
        }
	}

//  SEND / RECEIVE

    if(p->gcpara1_count>0)
    {
	MPI_Isend(send1,(p->gcpara1_count-p->periodicX1)*paramargin,MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(recv1,(p->gcpara1_count-p->periodicX1)*paramargin,MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcpara4_count>0)
    {
	MPI_Isend(send4,(p->gcpara4_count-p->periodicX4)*paramargin,MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(recv4,(p->gcpara4_count-p->periodicX4)*paramargin,MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
    }

    if(p->gcpara3_count>0)
    {
	MPI_Isend(send3,(p->gcpara3_count-p->periodicX3)*paramargin,MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(recv3,(p->gcpara3_count-p->periodicX3)*paramargin,MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcpara2_count>0)
    {
	MPI_Isend(send2,(p->gcpara2_count-p->periodicX2)*paramargin,MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(recv2,(p->gcpara2_count-p->periodicX2)*paramargin,MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
    }

    if(p->gcpara5_count>0)
    {
	MPI_Isend(send5,(p->gcpara5_count-p->periodicX5)*paramargin,MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
	MPI_Irecv(recv5,(p->gcpara5_count-p->periodicX5)*paramargin,MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
    }

    if(p->gcpara6_count>0)
    {
	MPI_Isend(send6,(p->gcpara6_count-p->periodicX6)*paramargin,MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
	MPI_Irecv(recv6,(p->gcpara6_count-p->periodicX6)*paramargin,MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
    }

//  WAIT

    gcwait(p);

//  FILL RECEIVE

    count=0;
    for(q=p->periodicX1;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->gcpara1[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
            if(p->gcpara1[q][2+gcv]==1)
            f(i-n-1,j,k)=recv1[count];
            ++count;
        }
    }

    count=0;
	for(q=p->periodicX3;q<p->gcpara3_count;++q)
	{
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->gcpara3[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
            if(p->gcpara3[q][2+gcv]==1)
            f(i,j-n-1,k)=recv3[count];
            ++count;
        }
	}

	count=0;
    for(q=p->periodicX5;q<p->gcpara5_count;++q)
    {
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->gcpara5[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
            if(p->gcpara5[q][2+gcv]==1)
            f(i,j,k-n-1)=recv5[count];
            ++count;
        }
    }

    count=0;
	for(q=p->periodicX4;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0]-aa;
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->gcpara4[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
            if(p->gcpara4[q][2+gcv]==1)
            f(i+n+1,j,k)=recv4[count];
            ++count;
        }
	}

    count=0;
	for(q=p->periodicX2;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1]-bb;
    k=p->gcpara2[q][2];
        
        if(p->gcpara2[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
            if(p->gcpara2[q][2+gcv]==1)
            f(i,j+n+1,k)=recv2[count];
            ++count;
        }
	}

    count=0;
	for(q=p->periodicX6;q<p->gcpara6_count;++q)
	{
    i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2]-cc;
        
        if(p->gcpara6[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {  
            if(p->gcpara6[q][2+gcv]==1)
            f(i,j,k+n+1)=recv6[count];
            ++count;
        }
	}
    
    
    /*
    if(p->mpirank==5)
    {
    i=p->knox-1;
    j=0;
    k=6;
    
    GC4LOOP
    {
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];
    if(p->gcb4[n][3]==4)
    cout<<i<<" "<<p->gcb4[n][4]<<" PRESS PERIODX 4: "<<a->press(i-2,j,k)<<" "<<a->press(i-1,j,k)<<" "<<a->press(i,j,k)<<" . "<<a->press(i+1,j,k)<<" "<<a->press(i+2,j,k)<<" "<<a->press(i+3,j,k)<<" "<<endl;
    }
    }
    
    
    if(p->mpirank==0)
    {
    i=0;
    j=0;
    k=6;
    
    GC4LOOP
    {
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];
    if(p->gcb4[n][3]==1)
    cout<<i<<" "<<p->gcb4[n][4]<<" PRESS PERIODX 1: "<<a->press(i-3,j,k)<<" "<<a->press(i-2,j,k)<<" "<<a->press(i-1,j,k)<<" . "<<a->press(i,j,k)<<" "<<a->press(i+1,j,k)<<" "<<a->press(i+2,j,k)<<" "<<endl;
    }
    }*/

}

