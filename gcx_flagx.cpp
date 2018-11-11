/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
    }

    count=0;
    for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];

        isend2[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
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
    }

    count=0;
    for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];

        isend4[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
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
    }

    count=0;
    for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];

        isend6[count]=flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
        ++count;
    }


//  Communication

    if(p->gcpara1_count>0)
    {
	MPI_Isend(isend1,p->gcpara1_count,MPI_INT,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(irecv1,p->gcpara1_count,MPI_INT,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcpara2_count>0)
    {
	MPI_Isend(isend2,p->gcpara2_count,MPI_INT,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(irecv2,p->gcpara2_count,MPI_INT,p->nb2,tag3,mpi_comm,&rreq2);
    }

    if(p->gcpara3_count>0)
    {
	MPI_Isend(isend3,p->gcpara3_count,MPI_INT,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(irecv3,p->gcpara3_count,MPI_INT,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcpara4_count>0)
    {
	MPI_Isend(isend4,p->gcpara4_count,MPI_INT,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(irecv4,p->gcpara4_count,MPI_INT,p->nb4,tag1,mpi_comm,&rreq4);
    }

    if(p->gcpara5_count>0)
    {
	MPI_Isend(isend5,p->gcpara5_count,MPI_INT,p->nb5,tag5,mpi_comm,&sreq5);
	MPI_Irecv(irecv5,p->gcpara5_count,MPI_INT,p->nb5,tag6,mpi_comm,&rreq5);
    }

    if(p->gcpara6_count>0)
    {
	MPI_Isend(isend6,p->gcpara6_count,MPI_INT,p->nb6,tag6,mpi_comm,&sreq6);
	MPI_Irecv(irecv6,p->gcpara6_count,MPI_INT,p->nb6,tag5,mpi_comm,&rreq6);
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
    }

    count=0;
    for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];

        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]=irecv2[count];
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
    }

    count=0;
    for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];

        flag[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=irecv4[count];
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
    }

    count=0;
    for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];

        flag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]=irecv6[count];
        ++count;
    }

}
