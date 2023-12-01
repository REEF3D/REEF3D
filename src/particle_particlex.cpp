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

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void particle_f::particlex(lexer* p, fdm* a, ghostcell* pgc)
{
    xchange=0;

    for(n=0;n<6;++n)
    {
    pxs[n]=0;
    pxr[n]=0;
    }

    for(n=0;n<posactive;++n)
    if(posflag[n]==1)
    {
        // POS
        i = p->posc_i(pos[n][0]);
        j = p->posc_j(pos[n][1]);
        k = p->posc_k(pos[n][2]);


        if(p->flag5[IJK]<0 && p->flag5[IJK]>-10)
        {
        ++xchange;
        posflag[n]=2;
        //posmem[pcount]=n;
        //++pcount;

            if(p->flag5[IJK]==-1)
            {
            posxs[0][pxs[0]+0]=pos[n][0];
            posxs[0][pxs[0]+1]=pos[n][1];
            posxs[0][pxs[0]+2]=pos[n][2];
            posxs[0][pxs[0]+3]=pos[n][3];
            posxs[0][pxs[0]+4]=pos[n][4];
            pxs[0]+=5;
            }

            if(p->flag5[IJK]==-2)
            {
            posxs[1][pxs[1]+0]=pos[n][0];
            posxs[1][pxs[1]+1]=pos[n][1];
            posxs[1][pxs[1]+2]=pos[n][2];
            posxs[1][pxs[1]+3]=pos[n][3];
            posxs[1][pxs[1]+4]=pos[n][4];
            pxs[1]+=5;
            }

            if(p->flag5[IJK]==-3)
            {
            posxs[2][pxs[2]+0]=pos[n][0];
            posxs[2][pxs[2]+1]=pos[n][1];
            posxs[2][pxs[2]+2]=pos[n][2];
            posxs[2][pxs[2]+3]=pos[n][3];
            posxs[2][pxs[2]+4]=pos[n][4];
            pxs[2]+=5;
            }

            if(p->flag5[IJK]==-4)
            {
            posxs[3][pxs[3]+0]=pos[n][0];
            posxs[3][pxs[3]+1]=pos[n][1];
            posxs[3][pxs[3]+2]=pos[n][2];
            posxs[3][pxs[3]+3]=pos[n][3];
            posxs[3][pxs[3]+4]=pos[n][4];
            pxs[3]+=5;
            }

            if(p->flag5[IJK]==-5)
            {
            posxs[4][pxs[4]+0]=pos[n][0];
            posxs[4][pxs[4]+1]=pos[n][1];
            posxs[4][pxs[4]+2]=pos[n][2];
            posxs[4][pxs[4]+3]=pos[n][3];
            posxs[4][pxs[4]+4]=pos[n][4];
            pxs[4]+=5;
            }

            if(p->flag5[IJK]==-6)
            {
            posxs[5][pxs[5]+0]=pos[n][0];
            posxs[5][pxs[5]+1]=pos[n][1];
            posxs[5][pxs[5]+2]=pos[n][2];
            posxs[5][pxs[5]+3]=pos[n][3];
            posxs[5][pxs[5]+4]=pos[n][4];
            pxs[5]+=5;
            }
        }
    }

    pgc->parapls(p,posxs,posxr,pxs,pxr);

// ---------------------------------------------------
// FILL
// ---------------------------------------------------

	//pos
    // q=0;
    // while(pcount>0)
    // {
    //     if(pxr[q]>0)
    //     {
    //     pos[posmem[pcount]][0]=posxr[q][pxr[q]-5];
    //     pos[posmem[pcount]][1]=posxr[q][pxr[q]-4];
    //     pos[posmem[pcount]][2]=posxr[q][pxr[q]-3];
    //     pos[posmem[pcount]][3]=posxr[q][pxr[q]-2];
    //     pos[posmem[pcount]][4]=posxr[q][pxr[q]-1];
	// 	--pcount;
    //     posflag[posmem[pcount]]=3;
    //     pxr[q]-=5;
    //     }

    //     if(pxr[q]<=0)
    //     ++q;

    //     if(q>=6)
    //     break;
    // }



    q=0;
    while(posactive<maxparticle)
    {
        if(pxr[q]>0)
        {
        pos[posactive][0]=posxr[q][pxr[q]-5];
        pos[posactive][1]=posxr[q][pxr[q]-4];
        pos[posactive][2]=posxr[q][pxr[q]-3];
        pos[posactive][3]=posxr[q][pxr[q]-2];
        pos[posactive][4]=posxr[q][pxr[q]-1];
        posflag[posactive]=3;
        pxr[q]-=5;
        ++posactive;
        }

        if(pxr[q]<=0)
        ++q;

        if(q>=6)
        break;
    }
}


