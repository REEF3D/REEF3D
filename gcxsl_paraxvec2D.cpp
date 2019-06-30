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

void ghostcell::gcparaxvec2D(lexer* p, vec2D &x, int gcv, cpt2D &C)
{
	if(gcv==1)
	gcslparaxvec_slr(p,x,C,1);

	if(gcv==2)
	gcslparaxvec_slr(p,x,C,2);

	if(gcv==3 || gcv==4)
	gcslparaxvec_slr(p,x,C,4);
}

void ghostcell::gcslparaxvec_slr(lexer* p, vec2D &x, cpt2D &C, int gcv)
{
	starttime=timer();

    paramargin=3;

//  FILL SEND
    count=0;
    for(q=0;q<p->gcslpara1_count;++q)
    {
    n=p->gcslpara1[q][8+gcv];


        if(p->gcslpara1[q][2+gcv]==1)
        {
        send1[count]=x.V[I_J];
        ++count;

        send1[count]=x.V[Ip1_J];
        ++count;

        send1[count]=x.V[Ip2_J];
        ++count;
        }
    }

	count=0;
	for(q=0;q<p->gcslpara2_count;++q)
	{
    n=p->gcslpara2[q][8+gcv];


        if(p->gcslpara2[q][2+gcv]==1)
        {
        send2[count]=x.V[I_J];
        ++count;

        send2[count]=x.V[I_Jm1];
        ++count;

        send2[count]=x.V[I_Jm2];
        ++count;
        }
	}

    count=0;
    for(q=0;q<p->gcslpara3_count;++q)
    {
    n=p->gcslpara3[q][8+gcv];


        if(p->gcslpara3[q][2+gcv]==1)
        {
        send3[count]=x.V[I_J];
        ++count;

        send3[count]=x.V[I_Jp1];
        ++count;

        send3[count]=x.V[I_Jp2];
        ++count;
        }
    }

	count=0;
	for(q=0;q<p->gcslpara4_count;++q)
	{
    n=p->gcslpara4[q][8+gcv];

        if(p->gcslpara4[q][2+gcv]==1)
        {
        send4[count]=x.V[I_J];
        ++count;

        send4[count]=x.V[Im1_J];
        ++count;

        send4[count]=x.V[Im2_J];
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

    gcwait(p);

//  FILL RECEIVE

    count=0;
    for(q=0;q<p->gcslpara1_count;++q)
    {
    n=p->gcslpara1[q][8+gcv];

        if(p->gcslpara1[q][2+gcv]==1)
        {
        x.V[Im1_J]=recv1[count];
        ++count;

        x.V[Im2_J]=recv1[count];
        ++count;

        x.V[Im3_J]=recv1[count];
        ++count;
        }

    }

    count=0;
	for(q=0;q<p->gcslpara2_count;++q)
	{
    n=p->gcslpara2[q][8+gcv];

        if(p->gcslpara2[q][2+gcv]==1)
        {
        x.V[I_Jp1]=recv2[count];
        ++count;

        x.V[I_Jp2]=recv2[count];
        ++count;

        x.V[I_Jp3]=recv2[count];
        ++count;
        }
	}

	count=0;
	for(q=0;q<p->gcslpara3_count;++q)
	{
    n=p->gcslpara3[q][8+gcv];

        if(p->gcslpara3[q][2+gcv]==1)
        {
        x.V[I_Jm1]=recv3[count];
        ++count;

        x.V[I_Jm2]=recv3[count];
        ++count;

        x.V[I_Jm3]=recv3[count];
        ++count;
        
        }
	}

    count=0;
	for(q=0;q<p->gcslpara4_count;++q)
	{
    n=p->gcslpara4[q][8+gcv];

        if(p->gcslpara4[q][2+gcv]==1)
        {
        x.V[Ip1_J]=recv4[count];
        ++count;

        x.V[Ip2_J]=recv4[count];
        ++count;

        x.V[Ip3_J]=recv4[count];
        ++count;
        
        //if(p->mpirank==0)
        //cout<<n<<" . "<<Ip1_J<<" "<<Ip2_J<<" "<<Ip3_J<<" "<<endl;
        }
	}
    
    //if(p->mpirank==0)
    //cout<<endl<<endl;

	endtime=timer();
	p->xtime+=endtime-starttime;
}
