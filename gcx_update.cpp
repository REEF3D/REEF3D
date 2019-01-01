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

void ghostcell::gcxupdate(lexer* p)
{

    for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];

        // 1
        if(p->flag1[UIm1JK]<AIR || p->flag1[UIJK]<AIR)
        p->gcpara1[n][3]=p->Y71;

        if(p->flag1[UIm1JK]>=AIR && p->flag1[UIJK]>=AIR)
        p->gcpara1[n][3]=1;
        
        // 2
        if(p->flag2[VIm1JK]<AIR || p->flag2[VIJK]<AIR)
        p->gcpara1[n][4]=p->Y72;

        if(p->flag2[VIm1JK]>=AIR && p->flag2[VIJK]>=AIR)
        p->gcpara1[n][4]=1;
        
        // 3
        if(p->flag3[WIm1JK]<AIR || p->flag3[WIJK]<AIR)
        p->gcpara1[n][5]=p->Y73;

        if(p->flag3[WIm1JK]>=AIR && p->flag3[WIJK]>=AIR)
        p->gcpara1[n][5]=1;
        
        // 4
        if(p->flag4[Im1JK]<AIR || p->flag4[IJK]<AIR)
        p->gcpara1[n][6]=p->Y74;

        if(p->flag4[Im1JK]>=AIR && p->flag4[IJK]>=AIR)
        p->gcpara1[n][6]=1;
        
        // 4a
        if(p->flag4[Im1JK]<=SOLID || p->flag4[IJK]<=SOLID)
        p->gcpara1[n][7]=p->Y74;

        if(p->flag4[Im1JK]>SOLID && p->flag4[IJK]>SOLID)
        p->gcpara1[n][7]=1;

    }

    for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];
        
        // 1
        if(p->flag1[UIJp1K]<AIR || p->flag1[UIJK]<AIR)
        p->gcpara2[n][3]=p->Y71;

        if(p->flag1[UIJp1K]>=AIR && p->flag1[UIJK]>=AIR)
        p->gcpara2[n][3]=1;
        
        // 2
        if(p->flag2[VIJp1K]<AIR || p->flag2[VIJK]<AIR)
        p->gcpara2[n][4]=p->Y72;

        if(p->flag2[VIJp1K]>=AIR && p->flag2[VIJK]>=AIR)
        p->gcpara2[n][4]=1;
        
        // 3
        if(p->flag3[WIJp1K]<AIR || p->flag3[WIJK]<AIR)
        p->gcpara2[n][5]=p->Y73;

        if(p->flag3[WIJp1K]>=AIR && p->flag3[WIJK]>=AIR)
        p->gcpara2[n][5]=1;
        
        // 4
        if(p->flag4[IJp1K]<AIR || p->flag4[IJK]<AIR)
        p->gcpara2[n][6]=p->Y74;

        if(p->flag4[IJp1K]>=AIR && p->flag4[IJK]>=AIR)
        p->gcpara2[n][6]=1;
        
        // 4a
        if(p->flag4[IJp1K]<=SOLID || p->flag4[IJK]<=SOLID)
        p->gcpara2[n][7]=p->Y74;

        if(p->flag4[IJp1K]>SOLID && p->flag4[IJK]>SOLID)
        p->gcpara2[n][7]=1;

    }

    for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];

        // 1
        if(p->flag1[UIJm1K]<AIR || p->flag1[UIJK]<AIR)
        p->gcpara3[n][3]=p->Y71;

        if(p->flag1[UIJm1K]>=AIR && p->flag1[UIJK]>=AIR)
        p->gcpara3[n][3]=1;
        
        // 2
        if(p->flag2[VIJm1K]<AIR || p->flag2[VIJK]<AIR)
        p->gcpara3[n][4]=p->Y72;

        if(p->flag2[VIJm1K]>=AIR && p->flag2[VIJK]>=AIR)
        p->gcpara3[n][4]=1;
        
        // 3
        if(p->flag3[WIJm1K]<AIR || p->flag3[WIJK]<AIR)
        p->gcpara3[n][5]=p->Y73;

        if(p->flag3[WIJm1K]>=AIR && p->flag3[WIJK]>=AIR)
        p->gcpara3[n][5]=1;
        
        // 4
        if(p->flag4[IJm1K]<AIR || p->flag4[IJK]<AIR)
        p->gcpara3[n][6]=p->Y74;

        if(p->flag4[IJm1K]>=AIR && p->flag4[IJK]>=AIR)
        p->gcpara3[n][6]=1;
        
        // 4a
        if(p->flag4[IJm1K]<=SOLID || p->flag4[IJK]<=SOLID)
        p->gcpara3[n][7]=p->Y74;

        if(p->flag4[IJm1K]>SOLID && p->flag4[IJK]>SOLID)
        p->gcpara3[n][7]=1;

    }

    for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];

        // 1
        if(p->flag1[UIp1JK]<AIR || p->flag1[UIJK]<AIR)
        p->gcpara4[n][3]=p->Y71;

        if(p->flag1[UIp1JK]>=AIR && p->flag1[UIJK]>=AIR)
        p->gcpara4[n][3]=1;
        
        // 2
        if(p->flag2[VIp1JK]<AIR || p->flag2[VIJK]<AIR)
        p->gcpara4[n][4]=p->Y72;

        if(p->flag2[VIp1JK]>=AIR && p->flag2[VIJK]>=AIR)
        p->gcpara4[n][4]=1;
        
        // 3
        if(p->flag3[WIp1JK]<AIR || p->flag3[WIJK]<AIR)
        p->gcpara4[n][5]=p->Y73;

        if(p->flag3[WIp1JK]>=AIR && p->flag3[WIJK]>=AIR)
        p->gcpara4[n][5]=1;
        
        // 4
        if(p->flag4[Ip1JK]<AIR || p->flag4[IJK]<AIR)
        p->gcpara4[n][6]=p->Y74;

        if(p->flag4[Ip1JK]>=AIR && p->flag4[IJK]>=AIR)
        p->gcpara4[n][6]=1;
        
        // 4a
        if(p->flag4[Ip1JK]<=SOLID || p->flag4[IJK]<=SOLID)
        p->gcpara4[n][7]=p->Y74;

        if(p->flag4[Ip1JK]>SOLID && p->flag4[IJK]>SOLID)
        p->gcpara4[n][7]=1;
        
    }

    for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];

        // 1
        if(p->flag1[UIJKm1]<AIR || p->flag1[UIJK]<AIR)
        p->gcpara5[n][3]=p->Y71;

        if(p->flag1[UIJKm1]>=AIR && p->flag1[UIJK]>=AIR)
        p->gcpara5[n][3]=1;
        
        // 2
        if(p->flag2[VIJKm1]<AIR || p->flag2[VIJK]<AIR)
        p->gcpara5[n][4]=p->Y72;

        if(p->flag2[VIJKm1]>=AIR && p->flag2[VIJK]>=AIR)
        p->gcpara5[n][4]=1;
        
        // 3
        if(p->flag3[WIJKm1]<AIR || p->flag3[WIJK]<AIR)
        p->gcpara5[n][5]=p->Y73;

        if(p->flag3[WIJKm1]>=AIR && p->flag3[WIJK]>=AIR)
        p->gcpara5[n][5]=1;
        
        // 4
        if(p->flag4[IJKm1]<AIR || p->flag4[IJK]<AIR)
        p->gcpara5[n][6]=p->Y74;

        if(p->flag4[IJKm1]>=AIR && p->flag4[IJK]>=AIR)
        p->gcpara5[n][6]=1;
        
        // 4a
        if(p->flag4[IJKm1]<=SOLID || p->flag4[IJK]<=SOLID)
        p->gcpara5[n][7]=p->Y74;

        if(p->flag4[IJKm1]>SOLID && p->flag4[IJK]>SOLID)
        p->gcpara5[n][7]=1;
		
    }

    for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];

        // 1
        if(p->flag1[UIJKp1]<AIR || p->flag1[UIJK]<AIR)
        p->gcpara6[n][3]=p->Y71;

        if(p->flag1[UIJKp1]>=AIR && p->flag1[UIJK]>=AIR)
        p->gcpara6[n][3]=1;
        
        // 2
        if(p->flag2[VIJKp1]<AIR || p->flag2[VIJK]<AIR)
        p->gcpara6[n][4]=p->Y72;

        if(p->flag2[VIJKp1]>=AIR && p->flag2[VIJK]>=AIR)
        p->gcpara6[n][4]=1;
        
        // 3
        if(p->flag3[WIJKp1]<AIR || p->flag3[WIJK]<AIR)
        p->gcpara6[n][5]=p->Y73;

        if(p->flag3[WIJKp1]>=AIR && p->flag3[WIJK]>=AIR)
        p->gcpara6[n][5]=1;
        
        // 4
        if(p->flag4[IJKp1]<AIR || p->flag4[IJK]<AIR)
        p->gcpara6[n][6]=p->Y74;

        if(p->flag4[IJKp1]>=AIR && p->flag4[IJK]>=AIR)
        p->gcpara6[n][6]=1;
        
        // 4a
        if(p->flag4[IJKp1]<=SOLID || p->flag4[IJK]<=SOLID)
        p->gcpara6[n][7]=p->Y74;

        if(p->flag4[IJKp1]>SOLID && p->flag4[IJK]>SOLID)
        p->gcpara6[n][7]=1;
    }
}

