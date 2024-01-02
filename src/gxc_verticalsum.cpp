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
#include"time.h"

void ghostcell::verticalsum(lexer *p, fdm* a, double **sum)
{
    double ***vsum;
    int **marker;
    double *s5,*s6,*r5,*r6;
    
    p->Darray(vsum,p->knox,p->knoy,p->mz);
    p->Iarray(marker,p->mz,2);
    p->Darray(s5,p->knox*p->knoy*(p->mz+1));
    p->Darray(s6,p->knox*p->knoy*(p->mz+1));
    p->Darray(r5,p->knox*p->knoy*(p->mz+1));
    p->Darray(r6,p->knox*p->knoy*(p->mz+1));
    
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    for(n=0;n<p->mz;++n)
    vsum[i][j][n]=0;
    
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    vsum[i][j][p->mk]=sum[i][j];
    
    //  FILL SEND MARKER
    marker[p->mk][0]=1;
    marker[p->mk][1]=1;
    
    for(int qn=0;qn<p->mz;++qn)
    {   
        if(p->gcpara5_count>0)
        for(n=0;n<p->mz;++n)
        s5[n]=marker[n][1];
        
        if(p->gcpara6_count>0)
        for(n=0;n<p->mz;++n)
        s6[n]=marker[n][0];
        
    //  SEND / RECEIVE MARKER
        if(p->gcpara5_count>0)
        {
        MPI_Isend(s5,p->mz,MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
        MPI_Irecv(r5,p->mz,MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
        }

        if(p->gcpara6_count>0)
        {
        MPI_Isend(s6,p->mz,MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
        MPI_Irecv(r6,p->mz,MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
        }
        
    // FILL RECEIVE MARKER
        if(p->gcpara5_count>0)
        for(n=0;n<p->mz;++n)
        if(r5[n]==1)
        marker[n][0]=1;
        
        if(p->gcpara6_count>0)
        for(n=0;n<p->mz;++n)
        if(r6[n]==1)
        marker[n][1]=1;
    
            
    //  FILL SEND
        count=0;
        for(q=0;q<p->gcpara5_count;++q)
        {
        i=p->gcpara5[q][0];
        j=p->gcpara5[q][1];

            for(n=0;n<p->mz;++n)
            {
            s5[count]=vsum[i][j][n];
            ++count;
            }
        }

        count=0;
        for(q=0;q<p->gcpara6_count;++q)
        {
        i=p->gcpara6[q][0];
        j=p->gcpara6[q][1];

            for(n=0;n<p->mz;++n)
            {
            s6[count]=vsum[i][j][n];
            ++count;
            }
        }


    //  SEND / RECEIVE
        if(p->gcpara5_count>0)
        {
        MPI_Isend(s5,p->gcpara5_count,MPI_DOUBLE,p->nb5,tag5,mpi_comm,&sreq5);
        MPI_Irecv(r5,p->gcpara5_count,MPI_DOUBLE,p->nb5,tag6,mpi_comm,&rreq5);
        }

        if(p->gcpara6_count>0)
        {
        MPI_Isend(s6,p->gcpara6_count,MPI_DOUBLE,p->nb6,tag6,mpi_comm,&sreq6);
        MPI_Irecv(r6,p->gcpara6_count,MPI_DOUBLE,p->nb6,tag5,mpi_comm,&rreq6);
        }

    //  WAIT

        if(p->gcpara5_count>0)
        {
        MPI_Wait(&sreq5,&status);
        MPI_Wait(&rreq5,&status);
        }

        if(p->gcpara6_count>0)
        {
        MPI_Wait(&sreq6,&status);
        MPI_Wait(&rreq6,&status);
        }

    //  FILL RECEIVE

        count=0;
        for(q=0;q<p->gcpara5_count;++q)
        {
        i=p->gcpara5[q][0];
        j=p->gcpara5[q][1];
            
            for(n=0;n<p->mz;++n)
            {
            if(marker[n][0]==1)
            vsum[i][j][k]=r5[count];
            ++count;
            }
        }

        count=0;
        for(q=0;q<p->gcpara6_count;++q)
        {
        i=p->gcpara6[q][0];
        j=p->gcpara6[q][1]; 


            for(n=0;n<p->mz;++n)
            {
            if(marker[n][1]==1)
            vsum[i][j][k]=r6[count];
            ++count;
            }

        }
    }
    
    // sum up and fill back
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    sum[i][j]=0.0;
    
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    for(n=0;n<p->mz;++n)
    sum[i][j] += vsum[i][j][n];
    
    p->del_Darray(vsum,p->knox,p->knoy,p->mz);
    p->del_Iarray(marker,p->mz,2);
    p->del_Darray(s5,p->knox*p->knoy*(p->mz+1));
    p->del_Darray(s6,p->knox*p->knoy*(p->mz+1));
    p->del_Darray(r5,p->knox*p->knoy*(p->mz+1));
    p->del_Darray(r6,p->knox*p->knoy*(p->mz+1));
}
