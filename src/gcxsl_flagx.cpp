/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

void ghostcell::gcslflagx(lexer* p, int *flag)
{
    // a patch of the mesh refinement (set_comms off) has no partition neighbours
    if(!do_comms)
    return;

    count=0;
    for(n=0;n<p->gcslpara1_count;++n)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];

        isend1[count]=flag[IJ];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcslpara2_count;++n)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];

        isend2[count]=flag[IJ];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcslpara3_count;++n)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];

        isend3[count]=flag[IJ];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcslpara4_count;++n)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];

        isend4[count]=flag[IJ];
        ++count;
    }

    Sendrecv_int(p->gcslpara1_count, p->gcslpara2_count, p->gcslpara3_count, p->gcslpara4_count, 0, 0);

//  Unpack

    count=0;
    for(n=0;n<p->gcslpara1_count;++n)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];

        flag[Im1J]=irecv1[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcslpara2_count;++n)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];

        flag[IJp1]=irecv2[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcslpara3_count;++n)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];

        flag[IJm1]=irecv3[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcslpara4_count;++n)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];

        flag[Ip1J]=irecv4[count];
        ++count;
    }

    //  Corners: second hop through the face ghosts filled above, as gcslparacox does for slices.
    //  Without it the diagonal ghost at an interior partition corner keeps its -10 from read_grid.

    for(n=0;n<p->gcslparaco1_count;++n)
    {
        i=p->gcslparaco1[n][0];
        j=p->gcslparaco1[n][1];
        isend1[n]=flag[IJ];
    }

    for(n=0;n<p->gcslparaco2_count;++n)
    {
        i=p->gcslparaco2[n][0];
        j=p->gcslparaco2[n][1];
        isend2[n]=flag[IJ];
    }

    for(n=0;n<p->gcslparaco3_count;++n)
    {
        i=p->gcslparaco3[n][0];
        j=p->gcslparaco3[n][1];
        isend3[n]=flag[IJ];
    }

    for(n=0;n<p->gcslparaco4_count;++n)
    {
        i=p->gcslparaco4[n][0];
        j=p->gcslparaco4[n][1];
        isend4[n]=flag[IJ];
    }

    Sendrecv_int(p->gcslparaco1_count, p->gcslparaco2_count, p->gcslparaco3_count, p->gcslparaco4_count, 0, 0);

    for(n=0;n<p->gcslparaco1_count;++n)
    {
        i=p->gcslparaco1[n][0];
        j=p->gcslparaco1[n][1];
        flag[Im1J]=irecv1[n];
    }

    for(n=0;n<p->gcslparaco2_count;++n)
    {
        i=p->gcslparaco2[n][0];
        j=p->gcslparaco2[n][1];
        flag[IJp1]=irecv2[n];
    }

    for(n=0;n<p->gcslparaco3_count;++n)
    {
        i=p->gcslparaco3[n][0];
        j=p->gcslparaco3[n][1];
        flag[IJm1]=irecv3[n];
    }

    for(n=0;n<p->gcslparaco4_count;++n)
    {
        i=p->gcslparaco4[n][0];
        j=p->gcslparaco4[n][1];
        flag[Ip1J]=irecv4[n];
    }
}
