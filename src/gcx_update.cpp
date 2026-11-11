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

void ghostcell::gcxupdate(lexer* p)
{
    for(n=0;n<p->gcpara1_count;++n)
    {
        i=p->gcpara1[n][0];
        j=p->gcpara1[n][1];
        k=p->gcpara1[n][2];

        // 1
        if(p->flag1(i-1,j,k)<AIR_FLAG || p->flag1(i,j,k)<AIR_FLAG)
            p->gcpara1[n][3]=p->Y71;

        if(p->flag1(i-1,j,k)>=AIR_FLAG && p->flag1(i,j,k)>=AIR_FLAG)
            p->gcpara1[n][3]=1;

        if(p->flag1(i,j,k)<AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG && p->flag1(i-1,j,k)>AIR_FLAG && p->flag4(i-1,j,k)>=AIR_FLAG)
            p->gcpara1[n][3]=1;

        if(p->flag1(i,j,k)>=AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG && p->flag1(i-1,j,k)<AIR_FLAG && p->flag4(i-1,j,k)>=AIR_FLAG)
            p->gcpara1[n][3]=2;

        // 2
        if(p->flag2(i-1,j,k)<AIR_FLAG || p->flag2(i,j,k)<AIR_FLAG)
            p->gcpara1[n][4]=p->Y72;

        if(p->flag2(i-1,j,k)>=AIR_FLAG && p->flag2(i,j,k)>=AIR_FLAG)
            p->gcpara1[n][4]=1;

        if(j + p->origin_j >= p->gknoy-1 && p->flag2(i-1,j-1,k)>=AIR_FLAG && p->flag2(i,j-1,k)>=AIR_FLAG)
            p->gcpara1[n][4]=1;

        // 3
        if(p->flag3(i-1,j,k)<AIR_FLAG || p->flag3(i,j,k)<AIR_FLAG)
            p->gcpara1[n][5]=p->Y73;

        if(p->flag3(i-1,j,k)>=AIR_FLAG && p->flag3(i,j,k)>=AIR_FLAG)
            p->gcpara1[n][5]=1;

        if(k + p->origin_k >= p->gknoz-1 && p->flag3(i-1,j,k-1)>=AIR_FLAG && p->flag3(i,j,k-1)>=AIR_FLAG)
            p->gcpara1[n][5]=1;

        // 4
        if(p->flag4(i-1,j,k)<AIR_FLAG || p->flag4(i,j,k)<AIR_FLAG)
            p->gcpara1[n][6]=p->Y74;

        if(p->flag4(i-1,j,k)>=AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG)
            p->gcpara1[n][6]=1;

        // 4a
        if(p->flag4(i-1,j,k)<=SOLID_FLAG || p->flag4(i,j,k)<=SOLID_FLAG)
            p->gcpara1[n][7]=p->Y74;

        if(p->flag4(i-1,j,k)>SOLID_FLAG && p->flag4(i,j,k)>SOLID_FLAG)
            p->gcpara1[n][7]=1;

        // 6
        p->gcpara1[n][8]=1;

    }

    for(n=0;n<p->gcpara2_count;++n)
    {
        i=p->gcpara2[n][0];
        j=p->gcpara2[n][1];
        k=p->gcpara2[n][2];

        // 1
        if(p->flag1(i,j+1,k)<AIR_FLAG || p->flag1(i,j,k)<AIR_FLAG)
            p->gcpara2[n][3]=p->Y71;

        if(p->flag1(i,j+1,k)>=AIR_FLAG && p->flag1(i,j,k)>=AIR_FLAG)
            p->gcpara2[n][3]=1;

        if(i + p->origin_i >= p->gknox-1 && p->flag1(i-1,j+1,k)>=AIR_FLAG && p->flag1(i-1,j,k)>=AIR_FLAG)
            p->gcpara2[n][3]=1;

        // 2
        if(p->flag2(i,j+1,k)<AIR_FLAG || p->flag2(i,j,k)<AIR_FLAG)
            p->gcpara2[n][4]=p->Y72;

        if(p->flag2(i,j+1,k)>=AIR_FLAG  && p->flag2(i,j,k)>=AIR_FLAG)
            p->gcpara2[n][4]=1;

        if(p->flag2(i,j,k)<AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG && p->flag2(i,j+1,k)>AIR_FLAG && p->flag4(i,j+1,k)>=AIR_FLAG)
            p->gcpara2[n][4]=1;

        if(p->flag2(i,j,k)>=AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG && p->flag2(i,j+1,k)<AIR_FLAG && p->flag4(i,j+1,k)>=AIR_FLAG)
            p->gcpara2[n][4]=2;

        // 3
        if(p->flag3(i,j+1,k)<AIR_FLAG || p->flag3(i,j,k)<AIR_FLAG)
            p->gcpara2[n][5]=p->Y73;

        if(p->flag3(i,j+1,k)>=AIR_FLAG && p->flag3(i,j,k)>=AIR_FLAG)
            p->gcpara2[n][5]=1;

        if(k + p->origin_k >= p->gknoz-1 && p->flag3(i,j+1,k-1)>=AIR_FLAG && p->flag3(i,j,k-1)>=AIR_FLAG)
            p->gcpara2[n][5]=1;

        // 4
        if(p->flag4(i,j+1,k)<AIR_FLAG || p->flag4(i,j,k)<AIR_FLAG)
            p->gcpara2[n][6]=p->Y74;

        if(p->flag4(i,j+1,k)>=AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG)
            p->gcpara2[n][6]=1;

        // 4a
        if(p->flag4(i,j+1,k)<=SOLID_FLAG || p->flag4(i,j,k)<=SOLID_FLAG)
            p->gcpara2[n][7]=p->Y74;

        if(p->flag4(i,j+1,k)>SOLID_FLAG && p->flag4(i,j,k)>SOLID_FLAG)
            p->gcpara2[n][7]=1;

        // 6
        p->gcpara2[n][8]=1;
    }

    for(n=0;n<p->gcpara3_count;++n)
    {
        i=p->gcpara3[n][0];
        j=p->gcpara3[n][1];
        k=p->gcpara3[n][2];

        // 1
        if(p->flag1(i,j-1,k)<AIR_FLAG || p->flag1(i,j,k)<AIR_FLAG)
            p->gcpara3[n][3]=p->Y71;

        if(p->flag1(i,j-1,k)>=AIR_FLAG && p->flag1(i,j,k)>=AIR_FLAG)
            p->gcpara3[n][3]=1;

        if(i + p->origin_i >= p->gknox-1 && p->flag1(i-1,j-1,k)>=AIR_FLAG && p->flag1(i-1,j,k)>=AIR_FLAG)
            p->gcpara3[n][3]=1;

        // 2
        if(p->flag2(i,j-1,k)<AIR_FLAG || p->flag2(i,j,k)<AIR_FLAG)
            p->gcpara3[n][4]=p->Y72;

        if(p->flag2(i,j-1,k)>=AIR_FLAG && p->flag2(i,j,k)>=AIR_FLAG)
            p->gcpara3[n][4]=1;

        if(p->flag2(i,j,k)<AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG && p->flag2(i,j-1,k)>AIR_FLAG && p->flag4(i,j-1,k)>=AIR_FLAG)
            p->gcpara3[n][4]=1;

        if(p->flag2(i,j,k)>=AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG && p->flag2(i,j-1,k)<AIR_FLAG && p->flag4(i,j-1,k)>=AIR_FLAG)
            p->gcpara3[n][4]=2;

        // 3
        if(p->flag3(i,j-1,k)<AIR_FLAG || p->flag3(i,j,k)<AIR_FLAG)
            p->gcpara3[n][5]=p->Y73;

        if(p->flag3(i,j-1,k)>=AIR_FLAG && p->flag3(i,j,k)>=AIR_FLAG)
            p->gcpara3[n][5]=1;

        if(k + p->origin_k >= p->gknoz-1 && p->flag3(i,j-1,k-1)>=AIR_FLAG && p->flag3(i,j,k-1)>=AIR_FLAG)
            p->gcpara3[n][5]=1;

        // 4
        if(p->flag4(i,j-1,k)<AIR_FLAG || p->flag4(i,j,k)<AIR_FLAG)
            p->gcpara3[n][6]=p->Y74;

        if(p->flag4(i,j-1,k)>=AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG)
            p->gcpara3[n][6]=1;

        // 4a
        if(p->flag4(i,j-1,k)<=SOLID_FLAG || p->flag4(i,j,k)<=SOLID_FLAG)
            p->gcpara3[n][7]=p->Y74;

        if(p->flag4(i,j-1,k)>SOLID_FLAG && p->flag4(i,j,k)>SOLID_FLAG)
            p->gcpara3[n][7]=1;

        // 6
        p->gcpara3[n][8]=1;
    }

    for(n=0;n<p->gcpara4_count;++n)
    {
        i=p->gcpara4[n][0];
        j=p->gcpara4[n][1];
        k=p->gcpara4[n][2];

        // 1
        if(p->flag1(i+1,j,k)<AIR_FLAG || p->flag1(i,j,k)<AIR_FLAG)
            p->gcpara4[n][3]=p->Y71;

        if(p->flag1(i+1,j,k)>=AIR_FLAG && p->flag1(i,j,k)>=AIR_FLAG)
            p->gcpara4[n][3]=1;

        if(p->flag1(i,j,k)<AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG && p->flag1(i+1,j,k)>AIR_FLAG && p->flag4(i+1,j,k)>=AIR_FLAG)
            p->gcpara4[n][3]=1;

        if(p->flag1(i,j,k)>=AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG && p->flag1(i+1,j,k)<AIR_FLAG && p->flag4(i+1,j,k)>=AIR_FLAG)
            p->gcpara4[n][3]=2;

        // 2
        if(p->flag2(i+1,j,k)<AIR_FLAG || p->flag2(i,j,k)<AIR_FLAG)
            p->gcpara4[n][4]=p->Y72;

        if(p->flag2(i+1,j,k)>=AIR_FLAG && p->flag2(i,j,k)>=AIR_FLAG)
            p->gcpara4[n][4]=1;

        if(j + p->origin_j >= p->gknoy-1 && p->flag2(i+1,j-1,k)>=AIR_FLAG && p->flag2(i,j-1,k)>=AIR_FLAG)
            p->gcpara4[n][4]=1;

        // 3
        if(p->flag3(i+1,j,k)<AIR_FLAG || p->flag3(i,j,k)<AIR_FLAG)
            p->gcpara4[n][5]=p->Y73;

        if(p->flag3(i+1,j,k)>=AIR_FLAG && p->flag3(i,j,k)>=AIR_FLAG)
            p->gcpara4[n][5]=1;

        if(k + p->origin_k >= p->gknoz-1 && p->flag3(i+1,j,k-1)>=AIR_FLAG && p->flag3(i,j,k-1)>=AIR_FLAG)
            p->gcpara4[n][5]=1;

        // 4
        if(p->flag4(i+1,j,k)<AIR_FLAG || p->flag4(i,j,k)<AIR_FLAG)
            p->gcpara4[n][6]=p->Y74;

        if(p->flag4(i+1,j,k)>=AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG)
            p->gcpara4[n][6]=1;

        // 4a
        if(p->flag4(i+1,j,k)<=SOLID_FLAG || p->flag4(i,j,k)<=SOLID_FLAG)
            p->gcpara4[n][7]=p->Y74;

        if(p->flag4(i+1,j,k)>SOLID_FLAG && p->flag4(i,j,k)>SOLID_FLAG)
            p->gcpara4[n][7]=1;

        // 6
        p->gcpara4[n][8]=1;
    }

    for(n=0;n<p->gcpara5_count;++n)
    {
        i=p->gcpara5[n][0];
        j=p->gcpara5[n][1];
        k=p->gcpara5[n][2];

        // 1
        if(p->flag1(i,j,k-1)<AIR_FLAG || p->flag1(i,j,k)<AIR_FLAG)
            p->gcpara5[n][3]=p->Y71;

        if(p->flag1(i,j,k-1)>=AIR_FLAG && p->flag1(i,j,k)>=AIR_FLAG)
            p->gcpara5[n][3]=1;

        if(i + p->origin_i >= p->gknox-1 && p->flag1(i-1,j,k-1)>=AIR_FLAG && p->flag1(i-1,j,k)>=AIR_FLAG)
            p->gcpara5[n][3]=1;

        // 2
        if(p->flag2(i,j,k-1)<AIR_FLAG || p->flag2(i,j,k)<AIR_FLAG)
            p->gcpara5[n][4]=p->Y72;

        if(p->flag2(i,j,k-1)>=AIR_FLAG && p->flag2(i,j,k)>=AIR_FLAG)
            p->gcpara5[n][4]=1;

        if(j + p->origin_j >= p->gknoy-1 && p->flag2(i,j-1,k-1)>=AIR_FLAG && p->flag2(i,j-1,k)>=AIR_FLAG)
            p->gcpara5[n][4]=1;

        // 3
        if(p->flag3(i,j,k-1)<AIR_FLAG || p->flag3(i,j,k)<AIR_FLAG)
            p->gcpara5[n][5]=p->Y73;

        if(p->flag3(i,j,k-1)>=AIR_FLAG && p->flag3(i,j,k)>=AIR_FLAG)
            p->gcpara5[n][5]=1;

        if(p->flag3(i,j,k)<AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG && p->flag3(i,j,k-1)>AIR_FLAG && p->flag4(i,j,k-1)>=AIR_FLAG)
            p->gcpara5[n][5]=1;

        if(p->flag3(i,j,k)>=AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG && p->flag3(i,j,k-1)<AIR_FLAG && p->flag4(i,j,k-1)>=AIR_FLAG)
            p->gcpara5[n][5]=2;

        // 4
        if(p->flag4(i,j,k-1)<AIR_FLAG || p->flag4(i,j,k)<AIR_FLAG)
            p->gcpara5[n][6]=p->Y74;

        if(p->flag4(i,j,k-1)>=AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG)
            p->gcpara5[n][6]=1;

        // 4a
        if(p->flag4(i,j,k-1)<=SOLID_FLAG || p->flag4(i,j,k)<=SOLID_FLAG)
            p->gcpara5[n][7]=p->Y74;

        if(p->flag4(i,j,k-1)>SOLID_FLAG && p->flag4(i,j,k)>SOLID_FLAG)
            p->gcpara5[n][7]=1;

        // 6
        p->gcpara5[n][8]=1;
    }

    for(n=0;n<p->gcpara6_count;++n)
    {
        i=p->gcpara6[n][0];
        j=p->gcpara6[n][1];
        k=p->gcpara6[n][2];

        // 1
        if(p->flag1(i,j,k+1)<AIR_FLAG || p->flag1(i,j,k)<AIR_FLAG)
            p->gcpara6[n][3]=p->Y71;

        if(p->flag1(i,j,k+1)>=AIR_FLAG && p->flag1(i,j,k)>=AIR_FLAG)
            p->gcpara6[n][3]=1;

        if(i + p->origin_i >= p->gknox-1 && p->flag1(i-1,j,k+1)>=AIR_FLAG && p->flag1(i-1,j,k)>=AIR_FLAG)
            p->gcpara6[n][3]=1;

        // 2
        if(p->flag2(i,j,k+1)<AIR_FLAG || p->flag2(i,j,k)<AIR_FLAG)
            p->gcpara6[n][4]=p->Y72;

        if(p->flag2(i,j,k+1)>=AIR_FLAG && p->flag2(i,j,k)>=AIR_FLAG)
            p->gcpara6[n][4]=1;

        if(j + p->origin_j >= p->gknoy-1 && p->flag2(i,j-1,k+1)>=AIR_FLAG && p->flag2(i,j-1,k)>=AIR_FLAG)
            p->gcpara6[n][4]=1;

        // 3
        if(p->flag3(i,j,k+1)<AIR_FLAG || p->flag3(i,j,k)<AIR_FLAG)
            p->gcpara6[n][5]=p->Y73;

        if(p->flag3(i,j,k+1)>=AIR_FLAG && p->flag3(i,j,k)>=AIR_FLAG)
            p->gcpara6[n][5]=1;

        if(p->flag3(i,j,k)<AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG && p->flag3(i,j,k+1)>AIR_FLAG && p->flag4(i,j,k+1)>=AIR_FLAG)
            p->gcpara6[n][5]=1;

        if(p->flag3(i,j,k)>=AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG && p->flag3(i,j,k+1)<AIR_FLAG && p->flag4(i,j,k+1)>=AIR_FLAG)
            p->gcpara6[n][5]=2;

        // 4
        if(p->flag4(i,j,k+1)<AIR_FLAG || p->flag4(i,j,k)<AIR_FLAG)
            p->gcpara6[n][6]=p->Y74;

        if(p->flag4(i,j,k+1)>=AIR_FLAG && p->flag4(i,j,k)>=AIR_FLAG)
            p->gcpara6[n][6]=1;

        // 4a
        if(p->flag4(i,j,k+1)<=SOLID_FLAG || p->flag4(i,j,k)<=SOLID_FLAG)
            p->gcpara6[n][7]=p->Y74;

        if(p->flag4(i,j,k+1)>SOLID_FLAG && p->flag4(i,j,k)>SOLID_FLAG)
            p->gcpara6[n][7]=1;

        // 6
        p->gcpara6[n][8]=1;
    }
}
