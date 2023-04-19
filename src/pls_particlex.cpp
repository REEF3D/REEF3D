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

#include"particle_pls.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void particle_pls::particlex(lexer* p, fdm* a, ghostcell* pgc)
{
    xchange=0;

    for(n=0;n<6;++n)
    {
    pxs[n]=0;
    nxs[n]=0;
    pxr[n]=0;
    nxr[n]=0;
    }

    for(n=0;n<posactive;++n)
    if(posflag[n]==1)
    {
        // POS
        i=int(pos[n][0]/dx);
        j=int(pos[n][1]/dx);
        k=int(pos[n][2]/dx);


        if(p->flag5[IJK]<0 && p->flag5[IJK]>-10)
        {
        ++xchange;
        posflag[n]=2;
        //posmem[pcount]=n;
        //++pcount;

            if(p->flag5[IJK]==-1)
            {
            posxs[0][pxs[0]+0]=pos[n][0]+p->originx;
            posxs[0][pxs[0]+1]=pos[n][1]+p->originy;
            posxs[0][pxs[0]+2]=pos[n][2]+p->originz;
            posxs[0][pxs[0]+3]=pos[n][3];
            posxs[0][pxs[0]+4]=pos[n][4];
            pxs[0]+=5;
            }

            if(p->flag5[IJK]==-2)
            {
            posxs[1][pxs[1]+0]=pos[n][0]+p->originx;
            posxs[1][pxs[1]+1]=pos[n][1]+p->originy;
            posxs[1][pxs[1]+2]=pos[n][2]+p->originz;
            posxs[1][pxs[1]+3]=pos[n][3];
            posxs[1][pxs[1]+4]=pos[n][4];
            pxs[1]+=5;
            }

            if(p->flag5[IJK]==-3)
            {
            posxs[2][pxs[2]+0]=pos[n][0]+p->originx;
            posxs[2][pxs[2]+1]=pos[n][1]+p->originy;
            posxs[2][pxs[2]+2]=pos[n][2]+p->originz;
            posxs[2][pxs[2]+3]=pos[n][3];
            posxs[2][pxs[2]+4]=pos[n][4];
            pxs[2]+=5;
            }

            if(p->flag5[IJK]==-4)
            {
            posxs[3][pxs[3]+0]=pos[n][0]+p->originx;
            posxs[3][pxs[3]+1]=pos[n][1]+p->originy;
            posxs[3][pxs[3]+2]=pos[n][2]+p->originz;
            posxs[3][pxs[3]+3]=pos[n][3];
            posxs[3][pxs[3]+4]=pos[n][4];
            pxs[3]+=5;
            }

            if(p->flag5[IJK]==-5)
            {
            posxs[4][pxs[4]+0]=pos[n][0]+p->originx;
            posxs[4][pxs[4]+1]=pos[n][1]+p->originy;
            posxs[4][pxs[4]+2]=pos[n][2]+p->originz;
            posxs[4][pxs[4]+3]=pos[n][3];
            posxs[4][pxs[4]+4]=pos[n][4];
            pxs[4]+=5;
            }

            if(p->flag5[IJK]==-6)
            {
            posxs[5][pxs[5]+0]=pos[n][0]+p->originx;
            posxs[5][pxs[5]+1]=pos[n][1]+p->originy;
            posxs[5][pxs[5]+2]=pos[n][2]+p->originz;
            posxs[5][pxs[5]+3]=pos[n][3];
            posxs[5][pxs[5]+4]=pos[n][4];
            pxs[5]+=5;
            }
        }
    }

    for(n=0;n<negactive;++n)
    if(negflag[n]==1)
    {

        //NEG
        i=int(neg[n][0]/dx);
        j=int(neg[n][1]/dx);
        k=int(neg[n][2]/dx);

        if(p->flag5[IJK]<0 && p->flag5[IJK]>-10)
        {
        ++xchange;
        negflag[n]=2;
        //negmem[pcount]=n;
        //ncount++;

            if(p->flag5[IJK]==-1)
            {
            negxs[0][nxs[0]+0]=neg[n][0]+p->originx;
            negxs[0][nxs[0]+1]=neg[n][1]+p->originy;
            negxs[0][nxs[0]+2]=neg[n][2]+p->originz;
            negxs[0][nxs[0]+3]=neg[n][3];
            negxs[0][nxs[0]+4]=neg[n][4];
            nxs[0]+=5;
            }

            if(p->flag5[IJK]==-2)
            {
            negxs[1][nxs[1]+0]=neg[n][0]+p->originx;
            negxs[1][nxs[1]+1]=neg[n][1]+p->originy;
            negxs[1][nxs[1]+2]=neg[n][2]+p->originz;
            negxs[1][nxs[1]+3]=neg[n][3];
            negxs[1][nxs[1]+4]=neg[n][4];
            nxs[1]+=5;
            }

            if(p->flag5[IJK]==-3)
            {
            negxs[2][nxs[2]+0]=neg[n][0]+p->originx;
            negxs[2][nxs[2]+1]=neg[n][1]+p->originy;
            negxs[2][nxs[2]+2]=neg[n][2]+p->originz;
            negxs[2][nxs[2]+3]=neg[n][3];
            negxs[2][nxs[2]+4]=neg[n][4];
            nxs[2]+=5;
            }

            if(p->flag5[IJK]==-4)
            {
            negxs[3][nxs[3]+0]=neg[n][0]+p->originx;
            negxs[3][nxs[3]+1]=neg[n][1]+p->originy;
            negxs[3][nxs[3]+2]=neg[n][2]+p->originz;
            negxs[3][nxs[3]+3]=neg[n][3];
            negxs[3][nxs[3]+4]=neg[n][4];
            nxs[3]+=5;
            }

            if(p->flag5[IJK]==-5)
            {
            negxs[4][nxs[4]+0]=neg[n][0]+p->originx;
            negxs[4][nxs[4]+1]=neg[n][1]+p->originy;
            negxs[4][nxs[4]+2]=neg[n][2]+p->originz;
            negxs[4][nxs[4]+3]=neg[n][3];
            negxs[4][nxs[4]+4]=neg[n][4];
            nxs[4]+=5;
            }

            if(p->flag5[IJK]==-6)
            {
            negxs[5][nxs[5]+0]=neg[n][0]+p->originx;
            negxs[5][nxs[5]+1]=neg[n][1]+p->originy;
            negxs[5][nxs[5]+2]=neg[n][2]+p->originz;
            negxs[5][nxs[5]+3]=neg[n][3];
            negxs[5][nxs[5]+4]=neg[n][4];
            nxs[5]+=5;
            }
        }
    }

    pgc->parapls(p,posxs,posxr,pxs,pxr);
    pgc->parapls(p,negxs,negxr,nxs,nxr);

// ---------------------------------------------------
// FILL
// ---------------------------------------------------

	//pos
    q=0;
    while(pcount>0)
    {
        if(pxr[q]>0)
        {
        pos[posmem[pcount]][0]=posxr[q][pxr[q]-5]-p->originx;
        pos[posmem[pcount]][1]=posxr[q][pxr[q]-4]-p->originy;
        pos[posmem[pcount]][2]=posxr[q][pxr[q]-3]-p->originz;
        pos[posmem[pcount]][3]=posxr[q][pxr[q]-2];
        pos[posmem[pcount]][4]=posxr[q][pxr[q]-1];
		--pcount;
        posflag[posmem[pcount]]=3;
        pxr[q]-=5;
        }

        if(pxr[q]<=0)
        ++q;

        if(q>=6)
        break;
    }



    q=0;
    while(posactive<maxparticle && pcount<=0)
    {
        if(pxr[q]>0)
        {
        pos[posactive][0]=posxr[q][pxr[q]-5]-p->originx;
        pos[posactive][1]=posxr[q][pxr[q]-4]-p->originy;
        pos[posactive][2]=posxr[q][pxr[q]-3]-p->originz;
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

    // neg
    q=0;
    while(ncount>0)
    {
        if(nxr[q]>0)
        {
        neg[negmem[ncount]][0]=negxr[q][nxr[q]-5]-p->originx;
        neg[negmem[ncount]][1]=negxr[q][nxr[q]-4]-p->originy;
        neg[negmem[ncount]][2]=negxr[q][nxr[q]-3]-p->originz;
        neg[negmem[ncount]][3]=negxr[q][nxr[q]-2];
        neg[negmem[ncount]][4]=negxr[q][nxr[q]-1];
        negflag[negmem[ncount]]=3;
		--ncount;
        nxr[q]-=5;
        }

        if(nxr[q]<=0)
        ++q;

        if(q>=6)
        break;
    }


    q=0;
    while(negactive<maxparticle && ncount<=0)
    {
        if(nxr[q]>0)
        {
        neg[negactive][0]=negxr[q][nxr[q]-5]-p->originx;
        neg[negactive][1]=negxr[q][nxr[q]-4]-p->originy;
        neg[negactive][2]=negxr[q][nxr[q]-3]-p->originz;
        neg[negactive][3]=negxr[q][nxr[q]-2];
        neg[negactive][4]=negxr[q][nxr[q]-1];
        negflag[negactive]=3;
        nxr[q]-=5;
        ++negactive;
        }

        if(nxr[q]<=0)
        ++q;

        if(q>=6)
        break;
    }
}


