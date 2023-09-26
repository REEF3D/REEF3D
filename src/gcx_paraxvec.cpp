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
#include"fdm.h"

void ghostcell::gcparaxvec(lexer* p, vec &x, int gcv)
{
		
	if(gcv==4)
	gcparaxvec_sr(p,x,a->C4,4);
	
	if(gcv==5)
	gcparaxvec_sr(p,x,a->C4a,5);
    
    if(gcv==6)
	gcparaxvec_sr(p,x,a->C6,6);
}
	
void ghostcell::gcparaxvec_sr(lexer* p, vec &x, cpt &C, int gcv)
{
	starttime=timer();
	
    paramargin=3;
	
//  FILL SEND
    count=0;
    for(q=0;q<p->gcpara1_count;++q)
    {
    n=p->gcpara1[q][8+gcv];
    
        
        if(p->gcpara1[q][2+gcv]==1)
        {
        send1[count]=x.V[I_J_K];
        ++count;

        send1[count]=x.V[Ip1_J_K];
        ++count;
        
        send1[count]=x.V[Ip2_J_K];
        ++count;
        }
        
    }
	
	count=0;
	for(q=0;q<p->gcpara2_count;++q)
	{
    n=p->gcpara2[q][8+gcv];

        
        if(p->gcpara2[q][2+gcv]==1)
        {
        send2[count]=x.V[I_J_K];
        ++count;

        send2[count]=x.V[I_Jm1_K];
        ++count;
  
        send2[count]=x.V[I_Jm2_K];
        ++count;
        }
        
	}

    count=0;
    for(q=0;q<p->gcpara3_count;++q)
    {
    n=p->gcpara3[q][8+gcv];

        
        if(p->gcpara3[q][2+gcv]==1)
        {
        send3[count]=x.V[I_J_K];
        ++count;
        
        send3[count]=x.V[I_Jp1_K];
        ++count;
     
        send3[count]=x.V[I_Jp2_K];
        ++count;
        }
        
    }
	
	count=0;
	for(q=0;q<p->gcpara4_count;++q)
	{
    n=p->gcpara4[q][8+gcv];

        if(p->gcpara4[q][2+gcv]==1)
        {
        send4[count]=x.V[I_J_K];
        ++count;

        send4[count]=x.V[Im1_J_K];
        ++count;

        send4[count]=x.V[Im2_J_K];
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara5_count;++q)
	{
    n=p->gcpara5[q][8+gcv];
        
        if(p->gcpara5[q][2+gcv]==1)
        {
        send5[count]=x.V[I_J_K];
        ++count;

        send5[count]=x.V[I_J_Kp1];
        ++count;

        send5[count]=x.V[I_J_Kp2];
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara6_count;++q)
	{
	n=p->gcpara6[q][8+gcv];
        
        if(p->gcpara6[q][2+gcv]==1)
        {
        send6[count]=x.V[I_J_K];
        ++count;

        send6[count]=x.V[I_J_Km1];
        ++count;

        send6[count]=x.V[I_J_Km2];
        ++count;
        }
	}


//  SEND / RECEIVE

    if(p->gcpara1_count>0)
    {
	MPI_Isend(send1,p->gcpara1_count*paramargin,MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(recv1,p->gcpara1_count*paramargin,MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcpara4_count>0)
    {
	MPI_Isend(send4,p->gcpara4_count*paramargin,MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(recv4,p->gcpara4_count*paramargin,MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
    }

    if(p->gcpara3_count>0)
    {
	MPI_Isend(send3,p->gcpara3_count*paramargin,MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(recv3,p->gcpara3_count*paramargin,MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcpara2_count>0)
    {
	MPI_Isend(send2,p->gcpara2_count*paramargin,MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(recv2,p->gcpara2_count*paramargin,MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
    }

    if(p->gcpara5_count>0)
    {
	MPI_Isend(send5,p->gcpara5_count*paramargin,MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
	MPI_Irecv(recv5,p->gcpara5_count*paramargin,MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
    }

    if(p->gcpara6_count>0)
    {
	MPI_Isend(send6,p->gcpara6_count*paramargin,MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
	MPI_Irecv(recv6,p->gcpara6_count*paramargin,MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
    }

//  WAIT
    gcwait(p);

//  FILL RECEIVE

    count=0;
    for(q=0;q<p->gcpara1_count;++q)
    {
    n=p->gcpara1[q][8+gcv];
        
        if(p->gcpara1[q][2+gcv]>=1)
        {
        if(p->gcpara1[q][2+gcv]==1)
        x.V[Im1_J_K]=recv1[count];
        ++count;

        if(p->gcpara1[q][2+gcv]==1)
        x.V[Im2_J_K]=recv1[count];
        ++count;
        if(p->gcpara1[q][2+gcv]==1)
        x.V[Im3_J_K]=recv1[count];
        ++count;
        }
        
    }

    count=0;
	for(q=0;q<p->gcpara2_count;++q)
	{
    n=p->gcpara2[q][8+gcv];

        if(p->gcpara2[q][2+gcv]>=1)
        {
        if(p->gcpara2[q][2+gcv]==1)
        x.V[I_Jp1_K]=recv2[count];
        ++count;

        if(p->gcpara2[q][2+gcv]==1)
        x.V[I_Jp2_K]=recv2[count];
        ++count;
        
        if(p->gcpara2[q][2+gcv]==1)
        x.V[I_Jp3_K]=recv2[count];
        ++count;
        }
	}	
	
	count=0;
	for(q=0;q<p->gcpara3_count;++q)
	{
    n=p->gcpara3[q][8+gcv];
        
        if(p->gcpara3[q][2+gcv]>=1)
        {
        if(p->gcpara3[q][2+gcv]==1)
        x.V[I_Jm1_K]=recv3[count];
        ++count;

        if(p->gcpara3[q][2+gcv]==1)
        x.V[I_Jm2_K]=recv3[count];
        ++count;

        if(p->gcpara3[q][2+gcv]==1)
        x.V[I_Jm3_K]=recv3[count];
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara4_count;++q)
	{
    n=p->gcpara4[q][8+gcv];
        
        if(p->gcpara4[q][2+gcv]>=1)
        {
        if(p->gcpara4[q][2+gcv]==1)
        x.V[Ip1_J_K]=recv4[count];
        ++count;

        if(p->gcpara4[q][2+gcv]==1)
        x.V[Ip2_J_K]=recv4[count];
        ++count;

        if(p->gcpara4[q][2+gcv]==1)        
        x.V[Ip3_J_K]=recv4[count];
        ++count;
        }
	}
	
	count=0;
    for(q=0;q<p->gcpara5_count;++q)
    {
    n=p->gcpara5[q][8+gcv];
        
        if(p->gcpara5[q][2+gcv]>=1)
        {
        if(p->gcpara5[q][2+gcv]==1)
        x.V[I_J_Km1]=recv5[count];
        ++count;
    
        if(p->gcpara5[q][2+gcv]==1)
        x.V[I_J_Km2]=recv5[count];
        ++count;

        if(p->gcpara5[q][2+gcv]==1)
        x.V[I_J_Km3]=recv5[count];
        ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcpara6_count;++q)
	{
    n=p->gcpara6[q][8+gcv];
        
        if(p->gcpara6[q][2+gcv]>=1)
        {
        if(p->gcpara6[q][2+gcv]==1)
        x.V[I_J_Kp1]=recv6[count];
        ++count;

        if(p->gcpara6[q][2+gcv]==1)
        x.V[I_J_Kp2]=recv6[count];
        ++count;

        if(p->gcpara6[q][2+gcv]==1)
        x.V[I_J_Kp3]=recv6[count];
        ++count;
        }
	}
	endtime=timer();
	p->xtime+=endtime-starttime;
}

