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
#include"field.h"

void ghostcell::gcparacox_generic(lexer* p,field& f,int *gcx_count, int ***gcx)
{    
    int aa,bb,cc;
	int r;
	int count[6];
    
    paramargin=3;
	
	for(qn=0;qn<6;++qn)
	count[qn]=0;

//  FILL SEND

    for(n=0;n<6;++n)
    for(q=0;q<gcx_count[n];++q)
    {
    i=gcx[n][q][0];
    j=gcx[n][q][1];
    k=gcx[n][q][2];
    
    aa=bb=cc=0;
    
		for(r=1;r<=paramargin;++r)
		{
		if(n==0)
		aa=r;
		
		if(n==1)
		bb=r;
		
		if(n==2)
		bb=r;
		
		if(n==3)
		aa=r;
		
		if(n==4)
		cc=r;
		
		if(n==5)
		cc=-r;
			
        //cout<<"Xs: "<<i<<" "<<k<<" "<<f(i+aa,j+bb,k+cc)<<"  "<<n<<endl;
		send[n][count[n]]=f(i+aa,j+bb,k+cc);
		++count[n];
		}
    }
	

//  SEND / RECEIVE

    for(qn=0;qn<6;++qn)
    if(count[qn]>0)
	{
	MPI_Isend(send[qn],count[qn],MPI_DOUBLE,nb[qn],stag[qn],mpi_comm,&sreq[qn]);
	MPI_Irecv(recv[qn],count[qn],MPI_DOUBLE,nb[qn],rtag[qn],mpi_comm,&rreq[qn]);
    }


//  WAIT

	for(qn=0;qn<6;++qn)
    if(count[qn]>0)
	{
    MPI_Wait(&sreq[qn],&status);
	MPI_Wait(&rreq[qn],&status);
	}

//  FILL RECEIVE

	for(qn=0;qn<6;++qn)
	count[qn]=0;

	for(n=0;n<6;++n)
    for(q=0;q<gcx_count[n];++q)
    {
    i=gcx[n][q][0];
    j=gcx[n][q][1];
    k=gcx[n][q][2];
    
    aa=bb=cc=0;
    
		for(r=1;r<=paramargin;++r)
		{
		if(n==0)
		aa=-r;
		
		if(n==1)
		bb=r;
		
		if(n==2)
		bb=-r;
		
		if(n==3)
		aa=r;
		
		if(n==4)
		cc=-r;
		
		if(n==5)
		cc=r;
			

		f(i+aa,j+bb,k+cc)=recv[n][count[n]];
		++count[n];
		}
    }						
}
