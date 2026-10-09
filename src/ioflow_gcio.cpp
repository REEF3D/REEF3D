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

#include"ioflow_gcio.h"
#include"lexer.h"
#include"patchBC_interface.h"
#include"patch_obj.h"
#include"looping.h"

namespace
{
inline bool gcio_inflow(int bc)
{
    return bc==1 || bc==6;
}

inline bool gcio_outflow(int bc)
{
    return bc==2 || bc==7 || bc==8;
}

// IO mark of the ghost cell next to cell (i,j,k) on side cs
inline void gcio_mark(lexer *p, int i, int j, int k, int cs, int val)
{
    if(cs==1)
    p->IO[Im1JK] = val;

    if(cs==4)
    p->IO[Ip1JK] = val;

    if(cs==3)
    p->IO[IJm1K] = val;

    if(cs==2)
    p->IO[IJp1K] = val;

    if(cs==5)
    p->IO[IJKm1] = val;

    if(cs==6)
    p->IO[IJKp1] = val;
}
}

void ioflow_gcio_lists(lexer *p, const int *flag, const int *wet)
{
    int i,j,k,n;
    int count1=0;
    int count2=0;

    GC4LOOP
    {
    i = p->gcb4[n][0];
    j = p->gcb4[n][1];
    k = p->gcb4[n][2];

        if(gcio_inflow(p->gcb4[n][4]) && flag[IJK]>0 && (wet==nullptr || wet[IJ]==1))
        ++count1;

        if(gcio_outflow(p->gcb4[n][4]) && flag[IJK]>0)
        ++count2;
    }

	p->Iresize(p->gcin,p->gcin_count, count1, 6, 6);
	p->Iresize(p->gcout,p->gcout_count, count2, 6, 6);

    count1=0;
    count2=0;
    GC4LOOP
    {
    i = p->gcb4[n][0];
    j = p->gcb4[n][1];
    k = p->gcb4[n][2];

        if(gcio_inflow(p->gcb4[n][4]) && flag[IJK]>0 && (wet==nullptr || wet[IJ]==1))
        {
        p->gcin[count1][0]=p->gcb4[n][0];
        p->gcin[count1][1]=p->gcb4[n][1];
        p->gcin[count1][2]=p->gcb4[n][2];
        p->gcin[count1][3]=p->gcb4[n][3];
        p->gcin[count1][5]=p->gcb4[n][5];
        ++count1;
        }

        if(gcio_outflow(p->gcb4[n][4]) && flag[IJK]>0)
        {
        p->gcout[count2][0]=p->gcb4[n][0];
        p->gcout[count2][1]=p->gcb4[n][1];
        p->gcout[count2][2]=p->gcb4[n][2];
        p->gcout[count2][3]=p->gcb4[n][3];
        p->gcout[count2][5]=p->gcb4[n][5];
        ++count2;
        }
    }

    p->gcin_count=count1;
    p->gcout_count=count2;
}

void ioflow_gcio_marks(lexer *p, const int *flag, patchBC_interface *pBC)
{
    int i,j,k,n;

    // reset first: the boundary cells change with sediment / topography updates
    MALOOP
    p->IO[IJK] = 0;

    GC4LOOP
    {
    i = p->gcb4[n][0];
    j = p->gcb4[n][1];
    k = p->gcb4[n][2];

        if(gcio_inflow(p->gcb4[n][4]) && flag[IJK]>0)
        gcio_mark(p,i,j,k,p->gcb4[n][3],1);

        if(gcio_outflow(p->gcb4[n][4]) && flag[IJK]>0)
        gcio_mark(p,i,j,k,p->gcb4[n][3],2);
    }

    // patch faces (B 4xx): IO 3. IO 1 would make the turbulence ghosts at the patch Dirichlet
    // (komega_bc, kepsilon_bc, NHFLOW start20V/30V/24V) without a profile written there.
    for(int qq=0;qq<pBC->obj_count;++qq)
    for(n=0;n<pBC->patch[qq]->gcb_count;++n)
    {
    i=pBC->patch[qq]->gcb[n][0];
    j=pBC->patch[qq]->gcb[n][1];
    k=pBC->patch[qq]->gcb[n][2];

    gcio_mark(p,i,j,k,pBC->patch[qq]->gcb[n][3],3);
    }
}
