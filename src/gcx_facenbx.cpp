/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
#include"fieldint.h"

void ghostcell::facenbx(lexer* p, fieldint& cv, int* range)
{
    count=0;
    for(n=0;n<p->gcpara1_count;++n)
    {
        i=p->gcpara1[n][0];
        j=p->gcpara1[n][1];
        k=p->gcpara1[n][2];

        isend1[count]=cv(i,j,k);
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara2_count;++n)
    {
        i=p->gcpara2[n][0];
        j=p->gcpara2[n][1];
        k=p->gcpara2[n][2];

        isend2[count]=cv(i,j,k);
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara3_count;++n)
    {
        i=p->gcpara3[n][0];
        j=p->gcpara3[n][1];
        k=p->gcpara3[n][2];

        isend3[count]=cv(i,j,k);
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara4_count;++n)
    {
        i=p->gcpara4[n][0];
        j=p->gcpara4[n][1];
        k=p->gcpara4[n][2];

        isend4[count]=cv(i,j,k);
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara5_count;++n)
    {
        i=p->gcpara5[n][0];
        j=p->gcpara5[n][1];
        k=p->gcpara5[n][2];

        isend5[count]=cv(i,j,k);
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara6_count;++n)
    {
        i=p->gcpara6[n][0];
        j=p->gcpara6[n][1];
        k=p->gcpara6[n][2];

        isend6[count]=cv(i,j,k);
        ++count;
    }

    Sendrecv_int(p->gcpara1_count,p->gcpara2_count,p->gcpara3_count,p->gcpara4_count,p->gcpara5_count,p->gcpara6_count);

    //  Unpack

    count=0;
    for(n=0;n<p->gcpara1_count;++n)
    {
        i=p->gcpara1[n][0];
        j=p->gcpara1[n][1];
        k=p->gcpara1[n][2];

        cv(i-1,j,k)=irecv1[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara2_count;++n)
    {
        i=p->gcpara2[n][0];
        j=p->gcpara2[n][1];
        k=p->gcpara2[n][2];

        cv(i,j+1,k)=irecv2[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara3_count;++n)
    {
        i=p->gcpara3[n][0];
        j=p->gcpara3[n][1];
        k=p->gcpara3[n][2];

        cv(i,j-1,k)=irecv3[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara4_count;++n)
    {
        i=p->gcpara4[n][0];
        j=p->gcpara4[n][1];
        k=p->gcpara4[n][2];

        cv(i+1,j,k)=irecv4[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara5_count;++n)
    {
        i=p->gcpara5[n][0];
        j=p->gcpara5[n][1];
        k=p->gcpara5[n][2];

        cv(i,j,k-1)=irecv5[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcpara6_count;++n)
    {
        i=p->gcpara6[n][0];
        j=p->gcpara6[n][1];
        k=p->gcpara6[n][2];

        cv(i,j,k+1)=irecv6[count];
        ++count;
    }
}
