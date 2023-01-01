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
#include"sliceint.h"

void ghostcell::gcslparax_int(lexer* p,sliceint& f,int gcv)
{
    paramargin=margin;

//  FILL SEND
    count=0;
    for(q=0;q<p->gcslpara1_count;++q)
    {
    i=p->gcslpara1[q][0];
    j=p->gcslpara1[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        isend1[count]=f(i+n,j);
        ++count;
        }
    }

    count=0;
    for(q=0;q<p->gcslpara3_count;++q)
    {
    i=p->gcslpara3[q][0];
    j=p->gcslpara3[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        isend3[count]=f(i,j+n);
        ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcslpara4_count;++q)
	{
    i=p->gcslpara4[q][0];
    j=p->gcslpara4[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        isend4[count]=f(i-n,j);
        ++count;
        }
	}
    
    count=0;
	for(q=0;q<p->gcslpara2_count;++q)
	{
    i=p->gcslpara2[q][0];
    j=p->gcslpara2[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        isend2[count]=f(i,j-n);
        ++count;
        }
	}

//  SEND / RECEIVE

    if(p->gcslpara1_count>0)
    {
	MPI_Isend(isend1,p->gcslpara1_count*paramargin,MPI_INT,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(irecv1,p->gcslpara1_count*paramargin,MPI_INT,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcslpara4_count>0)
    {
	MPI_Isend(isend4,p->gcslpara4_count*paramargin,MPI_INT,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(irecv4,p->gcslpara4_count*paramargin,MPI_INT,p->nb4,tag1,mpi_comm,&rreq4);
    }

    if(p->gcslpara3_count>0)
    {
	MPI_Isend(isend3,p->gcslpara3_count*paramargin,MPI_INT,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(irecv3,p->gcslpara3_count*paramargin,MPI_INT,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcslpara2_count>0)
    {
	MPI_Isend(isend2,p->gcslpara2_count*paramargin,MPI_INT,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(irecv2,p->gcslpara2_count*paramargin,MPI_INT,p->nb2,tag3,mpi_comm,&rreq2);
    }


//  WAIT

    gcslwait(p);

//  FILL RECEIVE

    count=0;
    for(q=0;q<p->gcslpara1_count;++q)
    {
    i=p->gcslpara1[q][0];
    j=p->gcslpara1[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        f(i-n-1,j)=irecv1[count];
        ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcslpara3_count;++q)
	{
    i=p->gcslpara3[q][0];
    j=p->gcslpara3[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        f(i,j-n-1)=irecv3[count];
        ++count;
        }
	}


    count=0;
	for(q=0;q<p->gcslpara4_count;++q)
	{
    i=p->gcslpara4[q][0];
    j=p->gcslpara4[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        f(i+n+1,j)=irecv4[count];
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcslpara2_count;++q)
	{
    i=p->gcslpara2[q][0];
    j=p->gcslpara2[q][1];
        
        for(n=0;n<paramargin;++n)
        {
        f(i,j+n+1)=irecv2[count];
        ++count;
        }
	}
}

void ghostcell::gcslparaxV_int(lexer* p, int *f,int gcv)
{
    paramargin=margin;

//  FILL SEND
    count=0;
    for(q=0;q<p->gcslpara1_count;++q)
    {
    i=p->gcslpara1[q][0];
    j=p->gcslpara1[q][1];
        
        isend1[count]=f[IJ];
        ++count;
        isend1[count]=f[Ip1J];
        ++count;
        isend1[count]=f[Ip2J];
        ++count;
    }

    count=0;
    for(q=0;q<p->gcslpara3_count;++q)
    {
    i=p->gcslpara3[q][0];
    j=p->gcslpara3[q][1];
        
        isend3[count]=f[IJ];
        ++count;
        isend3[count]=f[IJp1];
        ++count;
        isend3[count]=f[IJp2];
        ++count;
    }

    count=0;
	for(q=0;q<p->gcslpara4_count;++q)
	{
    i=p->gcslpara4[q][0];
    j=p->gcslpara4[q][1];
        
        isend4[count]=f[IJ];
        ++count;
        isend4[count]=f[Im1J];
        ++count;
        isend4[count]=f[Im2J];
        ++count;
	}
    
    count=0;
	for(q=0;q<p->gcslpara2_count;++q)
	{
    i=p->gcslpara2[q][0];
    j=p->gcslpara2[q][1];
        
        isend2[count]=f[IJ];
        ++count;
        isend2[count]=f[IJm1];
        ++count;
        isend2[count]=f[IJm2];
        ++count;
	}

//  SEND / RECEIVE

    if(p->gcslpara1_count>0)
    {
	MPI_Isend(isend1,p->gcslpara1_count*paramargin,MPI_INT,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(irecv1,p->gcslpara1_count*paramargin,MPI_INT,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcslpara4_count>0)
    {
	MPI_Isend(isend4,p->gcslpara4_count*paramargin,MPI_INT,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(irecv4,p->gcslpara4_count*paramargin,MPI_INT,p->nb4,tag1,mpi_comm,&rreq4);
    }

    if(p->gcslpara3_count>0)
    {
	MPI_Isend(isend3,p->gcslpara3_count*paramargin,MPI_INT,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(irecv3,p->gcslpara3_count*paramargin,MPI_INT,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcslpara2_count>0)
    {
	MPI_Isend(isend2,p->gcslpara2_count*paramargin,MPI_INT,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(irecv2,p->gcslpara2_count*paramargin,MPI_INT,p->nb2,tag3,mpi_comm,&rreq2);
    }


//  WAIT

    gcslwait(p);

//  FILL RECEIVE

    count=0;
    for(q=0;q<p->gcslpara1_count;++q)
    {
    i=p->gcslpara1[q][0];
    j=p->gcslpara1[q][1];
        
        f[Im1J]=irecv1[count];
        ++count;
        f[Im2J]=irecv1[count];
        ++count;
        f[Im3J]=irecv1[count];
        ++count;
    }

    count=0;
	for(q=0;q<p->gcslpara3_count;++q)
	{
    i=p->gcslpara3[q][0];
    j=p->gcslpara3[q][1];
        
        f[IJm1]=irecv3[count];
        ++count;
        f[IJm2]=irecv3[count];
        ++count;
        f[IJm3]=irecv3[count];
        ++count;
	}


    count=0;
	for(q=0;q<p->gcslpara4_count;++q)
	{
    i=p->gcslpara4[q][0];
    j=p->gcslpara4[q][1];
        
        f[Ip1J]=irecv4[count];
        ++count;
        f[Ip2J]=irecv4[count];
        ++count;
        f[Ip3J]=irecv4[count];
        ++count;
	}

    count=0;
	for(q=0;q<p->gcslpara2_count;++q)
	{
    i=p->gcslpara2[q][0];
    j=p->gcslpara2[q][1];
        
        f[IJp1]=irecv2[count];
        ++count;
        f[IJp2]=irecv2[count];
        ++count;
        f[IJp3]=irecv2[count];
        ++count;
	}
}