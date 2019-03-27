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
#include"fdm2D.h"

void ghostcell::gcparaxvec2D(lexer* p, fdm2D *b, vec2D &x, int gcv)
{
	
	if(gcv==1)
	gcparaxvec_slr(p,x,b->C1,1);
	
	if(gcv==2)
	gcparaxvec_slr(p,x,b->C2,2);
	
	if(gcv==3 || gcv==4)
	gcparaxvec_slr(p,x,b->C4,4);
}
	
void ghostcell::gcparaxvec_slr(lexer* p, vec2D &x, cpt2D &C, int gcv)
{
	starttime=timer();
	
    paramargin=3;
	
//  FILL SEND
    count=0;
    for(q=0;q<p->gcpara1_count;++q)
    {
    n=p->gcpara1[q][8+gcv];
    
        
        if(p->gcpara1[q][2+gcv]==1 || gcv==6)
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

        
        if(p->gcpara2[q][2+gcv]==1 || gcv==6)
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

        
        if(p->gcpara3[q][2+gcv]==1 || gcv==6)
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

        if(p->gcpara4[q][2+gcv]==1 || gcv==6)
        {
        send4[count]=x.V[I_J_K];
        ++count;

        send4[count]=x.V[Im1_J_K];
        ++count;

        send4[count]=x.V[Im2_J_K];
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

//  WAIT

    gcwait(p);

//  FILL RECEIVE

    count=0;
    for(q=0;q<p->gcpara1_count;++q)
    {
    n=p->gcpara1[q][8+gcv];
        
        if(p->gcpara1[q][2+gcv]==1 || gcv==6)
        {
        x.V[Im1_J_K]=recv1[count];
        ++count;

        x.V[Im2_J_K]=recv1[count];
        ++count;

        x.V[Im3_J_K]=recv1[count];
        ++count;
        }
        
    }

    count=0;
	for(q=0;q<p->gcpara2_count;++q)
	{
    n=p->gcpara2[q][8+gcv];

        if(p->gcpara2[q][2+gcv]==1 || gcv==6)
        {
        x.V[I_Jp1_K]=recv2[count];
        ++count;

        x.V[I_Jp2_K]=recv2[count];
        ++count;
        
        x.V[I_Jp3_K]=recv2[count];
        ++count;
        }
	}	
	
	count=0;
	for(q=0;q<p->gcpara3_count;++q)
	{
    n=p->gcpara3[q][8+gcv];
        
        if(p->gcpara3[q][2+gcv]==1 || gcv==6)
        {
        x.V[I_Jm1_K]=recv3[count];
        ++count;

        x.V[I_Jm2_K]=recv3[count];
        ++count;
        
        x.V[I_Jm3_K]=recv3[count];
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara4_count;++q)
	{
    n=p->gcpara4[q][8+gcv];
        
        if(p->gcpara4[q][2+gcv]==1 || gcv==6)
        {
        x.V[Ip1_J_K]=recv4[count];
        ++count;

        x.V[Ip2_J_K]=recv4[count];
        ++count;
        
        x.V[Ip3_J_K]=recv4[count];
        ++count;
        }
	}

	endtime=timer();
	p->xtime+=endtime-starttime;
}

