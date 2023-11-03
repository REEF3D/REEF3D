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

void ghostcell::flagx(lexer* p, int *flag)
{

    count=0;
    for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];

        isend1[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
        ++count;
        isend1[count]=flag[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
        ++count;
        isend1[count]=flag[(i-p->imin+2)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];

        isend2[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
        ++count;
        isend2[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin];
        ++count;
        isend2[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-2)*p->kmax + k-p->kmin];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];

        isend3[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
        ++count;
        isend3[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin];
        ++count;
        isend3[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+2)*p->kmax + k-p->kmin];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];

        isend4[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
        ++count;
        isend4[count]=flag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
        ++count;
        isend4[count]=flag[(i-p->imin-2)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];

        isend5[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
        ++count;
        isend5[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1];
        ++count;
        isend5[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+2];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];

        isend6[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
        ++count;
        isend6[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1];
        ++count;
        isend6[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-2];
        ++count;
    }


//  Communication

    if(p->gcpara1_count>0)
    {
	MPI_Isend(isend1,p->gcpara1_count*paramargin,MPI_INT,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(irecv1,p->gcpara1_count*paramargin,MPI_INT,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcpara2_count>0)
    {
	MPI_Isend(isend2,p->gcpara2_count*paramargin,MPI_INT,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(irecv2,p->gcpara2_count*paramargin,MPI_INT,p->nb2,tag3,mpi_comm,&rreq2);
    }

    if(p->gcpara3_count>0)
    {
	MPI_Isend(isend3,p->gcpara3_count*paramargin,MPI_INT,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(irecv3,p->gcpara3_count*paramargin,MPI_INT,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcpara4_count>0)
    {
	MPI_Isend(isend4,p->gcpara4_count*paramargin,MPI_INT,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(irecv4,p->gcpara4_count*paramargin,MPI_INT,p->nb4,tag1,mpi_comm,&rreq4);
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

    gcwait(p);

//  Unpack

    count=0;
    for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];

        flag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=irecv1[count];
        ++count;
        flag[(i-p->imin-2)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=irecv1[count];
        ++count;
        flag[(i-p->imin-3)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=irecv1[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];

        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]=irecv2[count];
        ++count;
        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+2)*p->kmax + k-p->kmin]=irecv2[count];
        ++count;
        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+3)*p->kmax + k-p->kmin]=irecv2[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];

        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]=irecv3[count];
        ++count;
        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-2)*p->kmax + k-p->kmin]=irecv3[count];
        ++count;
        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-3)*p->kmax + k-p->kmin]=irecv3[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];

        flag[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=irecv4[count];
        ++count;
        flag[(i-p->imin+2)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=irecv4[count];
        ++count;
        flag[(i-p->imin+3)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=irecv4[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];

        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]=irecv5[count];
        ++count;
        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-2]=irecv5[count];
        ++count;
        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-3]=irecv5[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];

        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]=irecv6[count];
        ++count;
        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+2]=irecv6[count];
        ++count;
        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+3]=irecv6[count];
        ++count;
    }
    
    
// -- Paraco


//  FILL SEND
    for(q=0;q<p->gcparaco1_count;++q)
    {
    i=p->gcparaco1[q][0];
    j=p->gcparaco1[q][1];
    k=p->gcparaco1[q][2];
	isend1[q]=flag[IJK];
    }

    for(q=0;q<p->gcparaco3_count;++q)
    {
    i=p->gcparaco3[q][0];
    j=p->gcparaco3[q][1];
    k=p->gcparaco3[q][2];
	isend3[q]=flag[IJK];
    }

	for(q=0;q<p->gcparaco5_count;++q)
	{
    i=p->gcparaco5[q][0];
    j=p->gcparaco5[q][1];
    k=p->gcparaco5[q][2];
	isend5[q]=flag[IJK];
	}

	for(q=0;q<p->gcparaco4_count;++q)
	{
    i=p->gcparaco4[q][0];
    j=p->gcparaco4[q][1];
    k=p->gcparaco4[q][2];
	isend4[q]=flag[IJK];
	}

	for(q=0;q<p->gcparaco2_count;++q)
	{
    i=p->gcparaco2[q][0];
    j=p->gcparaco2[q][1];
    k=p->gcparaco2[q][2];
	isend2[q]=flag[IJK];
	}

	for(q=0;q<p->gcparaco6_count;++q)
	{
    i=p->gcparaco6[q][0];
    j=p->gcparaco6[q][1];
    k=p->gcparaco6[q][2];
	isend6[q]=flag[IJK];
	}


//  SEND / RECEIVE

    if(p->gcparaco1_count>0)
    {
	MPI_Isend(isend1,p->gcparaco1_count,MPI_INT,p->nb1,tag,mpi_comm,&sreq1);
	MPI_Irecv(irecv1,p->gcparaco1_count,MPI_INT,p->nb1,tag,mpi_comm,&rreq1);
    }

    if(p->gcparaco4_count>0)
    {
	MPI_Isend(isend4,p->gcparaco4_count,MPI_INT,p->nb4,tag,mpi_comm,&sreq4);
	MPI_Irecv(irecv4,p->gcparaco4_count,MPI_INT,p->nb4,tag,mpi_comm,&rreq4);
    }

    if(p->gcparaco3_count>0)
    {
	MPI_Isend(isend3,p->gcparaco3_count,MPI_INT,p->nb3,tag,mpi_comm,&sreq3);
	MPI_Irecv(irecv3,p->gcparaco3_count,MPI_INT,p->nb3,tag,mpi_comm,&rreq3);
    }

    if(p->gcparaco2_count>0)
    {
	MPI_Isend(isend2,p->gcparaco2_count,MPI_INT,p->nb2,tag,mpi_comm,&sreq2);
	MPI_Irecv(irecv2,p->gcparaco2_count,MPI_INT,p->nb2,tag,mpi_comm,&rreq2);
    }

    if(p->gcparaco5_count>0)
    {
	MPI_Isend(isend5,p->gcparaco5_count,MPI_INT,p->nb5,tag,mpi_comm,&sreq5);
	MPI_Irecv(irecv5,p->gcparaco5_count,MPI_INT,p->nb5,tag,mpi_comm,&rreq5);
    }

    if(p->gcparaco6_count>0)
    {
	MPI_Isend(isend6,p->gcparaco6_count,MPI_INT,p->nb6,tag,mpi_comm,&sreq6);
	MPI_Irecv(irecv6,p->gcparaco6_count,MPI_INT,p->nb6,tag,mpi_comm,&rreq6);
    }

//  WAIT
    gcwait(p);

//  FILL RECEIVE
    //pip=4;
    for(q=0;q<p->gcparaco1_count;++q)
    {
    i=p->gcparaco1[q][0];
    j=p->gcparaco1[q][1];
    k=p->gcparaco1[q][2];
	flag[Im1JK]=irecv1[q];
    }

	for(q=0;q<p->gcparaco3_count;++q)
	{
    i=p->gcparaco3[q][0];
    j=p->gcparaco3[q][1];
    k=p->gcparaco3[q][2];
	flag[IJm1K]=irecv3[q];
	}

    for(q=0;q<p->gcparaco5_count;++q)
    {
    i=p->gcparaco5[q][0];
    j=p->gcparaco5[q][1];
    k=p->gcparaco5[q][2];
	flag[IJKm1]=irecv5[q];
    }

	for(q=0;q<p->gcparaco4_count;++q)
	{
    i=p->gcparaco4[q][0];
    j=p->gcparaco4[q][1];
    k=p->gcparaco4[q][2];
	flag[Ip1JK]=irecv4[q];
	}

	for(q=0;q<p->gcparaco2_count;++q)
	{
    i=p->gcparaco2[q][0];
    j=p->gcparaco2[q][1];
    k=p->gcparaco2[q][2];
	flag[IJp1K]=irecv2[q];
	}

	for(q=0;q<p->gcparaco6_count;++q)
	{
    i=p->gcparaco6[q][0];
    j=p->gcparaco6[q][1];
    k=p->gcparaco6[q][2];
   	flag[IJKp1]=irecv6[q];
	}
    

}
