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

#include"force.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void force::cellnodes(lexer* p, fdm *a, ghostcell *pgc)
{
    const double eps = 1.0e-5*p->DXM;
    int cn[4][4];
    int r,s;
    int count, fcount;
	double px1,px2,px3;
	double py1,py2,py3;
	double pz1,pz2,pz3;
	double ntx,nty,ntz;
	
    fcount=0;

    // identify active and unactive cell nodes
    for(n=0;n<p->facetnum;++n)
    {
    i=p->facet[n][0];
    j=p->facet[n][1];
    k=p->facet[n][2];

    if(i>=is && i<=ie && j>=js && j<=je && k>=ks && k<=ke)
	if(fcheck[fcount]==1)
    {
        for(r=0;r<4;++r)
        for(s=0;s<4;++s)
        cn[r][s]=-1;
		
		for(q=0;q<p->facet[n][6];++q)
		{

        cx = p->XP[p->facet[n][0]+marge];
        cy = p->YP[p->facet[n][1]+marge];
        cz = p->ZP[p->facet[n][2]+marge];
		
		px=p->ccpoint[p->facet[n][7+q]][0];
        py=p->ccpoint[p->facet[n][7+q]][1];
        pz=p->ccpoint[p->facet[n][7+q]][2];
	

    // Side 0
        // 0
        if(px > cx-0.5*p->dx-eps && px < cx-0.5*p->dx+eps
        && py > cy-0.5*p->dx-eps && py < cy-0.5*p->dx+eps
        && pz > cz-0.5*p->dx+eps && pz < cz+0.5*p->dx+eps)
        cn[0][0]=q;

        // 1
        if(px > cx-0.5*p->dx-eps && px < cx+0.5*p->dx-eps
        && py > cy-0.5*p->dx-eps && py < cy-0.5*p->dx+eps
        && pz > cz-0.5*p->dx-eps && pz < cz-0.5*p->dx+eps)
        cn[0][1]=q;

        // 2
        if(px > cx-0.5*p->dx+eps && px < cx+0.5*p->dx-eps
        && py > cy-0.5*p->dx-eps && py < cy-0.5*p->dx+eps
        && pz > cz+0.5*p->dx-eps && pz < cz+0.5*p->dx+eps)
        cn[0][2]=q;

    // Side 1
        // 0
        if(px > cx+0.5*p->dx-eps && px < cx+0.5*p->dx+eps
        && py > cy-0.5*p->dx-eps && py < cy-0.5*p->dx+eps
        && pz > cz-0.5*p->dx+eps && pz < cz+0.5*p->dx+eps)
        cn[1][0]=q;

        // 1
        if(px > cx+0.5*p->dx-eps && px < cx+0.5*p->dx+eps
        && py > cy-0.5*p->dx-eps && py < cy+0.5*p->dx-eps
        && pz > cz-0.5*p->dx-eps && pz < cz-0.5*p->dx+eps)
        cn[1][1]=q;

        // 2
        if(px > cx+0.5*p->dx-eps && px < cx+0.5*p->dx+eps
        && py > cy-0.5*p->dx+eps && py < cy+0.5*p->dx-eps
        && pz > cz+0.5*p->dx-eps && pz < cz+0.5*p->dx+eps)
        cn[1][2]=q;

    // Side 2
        // 0
        if(px > cx+0.5*p->dx-eps && px < cx+0.5*p->dx+eps
        && py > cy+0.5*p->dx-eps && py < cy+0.5*p->dx+eps
        && pz > cz-0.5*p->dx+eps && pz < cz+0.5*p->dx+eps)
        cn[2][0]=q;

        // 1
        if(px > cx-0.5*p->dx+eps && px < cx+0.5*p->dx+eps
        && py > cy+0.5*p->dx-eps && py < cy+0.5*p->dx+eps
        && pz > cz-0.5*p->dx-eps && pz < cz-0.5*p->dx+eps)
        cn[2][1]=q;

        // 2
        if(px > cx-0.5*p->dx+eps && px < cx+0.5*p->dx-eps
        && py > cy+0.5*p->dx-eps && py < cy+0.5*p->dx+eps
        && pz > cz+0.5*p->dx-eps && pz < cz+0.5*p->dx+eps)
        cn[2][2]=q;

    // Side 3
        // 0
        if(px > cx-0.5*p->dx-eps && px < cx-0.5*p->dx+eps
        && py > cy+0.5*p->dx-eps && py < cy+0.5*p->dx+eps
        && pz > cz-0.5*p->dx+eps && pz < cz+0.5*p->dx+eps)
        cn[3][0]=q;

        // 1
        if(px > cx-0.5*p->dx-eps && px < cx-0.5*p->dx+eps
        && py > cy-0.5*p->dx+eps && py < cy+0.5*p->dx+eps
        && pz > cz-0.5*p->dx-eps && pz < cz-0.5*p->dx+eps)
        cn[3][1]=q;

        // 2
        if(px > cx-0.5*p->dx-eps && px < cx-0.5*p->dx+eps
        && py > cy-0.5*p->dx+eps && py < cy+0.5*p->dx-eps
        && pz > cz+0.5*p->dx-eps && pz < cz+0.5*p->dx+eps)
        cn[3][2]=q;
        }


//_____________________________________________________
    // order
    count=0;

    // Side 0
        if(cn[0][0]>=0)
        {
        ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[0][0]]][0];
        ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[0][0]]][1];
        ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[0][0]]][2];
        ++count;
        }
    //cout<<p->mpirank<<" fcount: "<<fcount<<"  count: "<<count<<endl;

        if(cn[0][1]>=0 && cn[0][2]<0)
        {
        ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[0][1]]][0];
        ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[0][1]]][1];
        ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[0][1]]][2];
        ++count;
        }

        if(cn[0][2]>=0 && cn[0][1]<0)
        {
        ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[0][2]]][0];
        ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[0][2]]][1];
        ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[0][2]]][2];
        ++count;
        }

        if(cn[0][1]>=0 && cn[0][2]>=0)
        {
            if(p->ccpoint[p->facet[n][7+cn[0][1]]][0] < p->ccpoint[p->facet[n][7+cn[0][2]]][0])
            {
            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[0][1]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[0][1]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[0][1]]][2];
            ++count;

            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[0][2]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[0][2]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[0][2]]][2];
            ++count;
            }

            if(p->ccpoint[p->facet[n][7+cn[0][1]]][0] >= p->ccpoint[p->facet[n][7+cn[0][2]]][0])
            {
            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[0][2]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[0][2]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[0][2]]][2];
            ++count;

            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[0][1]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[0][1]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[0][1]]][2];
            ++count;
            }
        }

    // Side 1
        if(cn[1][0]>=0)
        {
        ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[1][0]]][0];
        ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[1][0]]][1];
        ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[1][0]]][2];
        ++count;
        }

        if(cn[1][1]>=0 && cn[1][2]<0)
        {
        ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[1][1]]][0];
        ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[1][1]]][1];
        ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[1][1]]][2];
        ++count;
        }

        if(cn[1][2]>=0 && cn[1][1]<0)
        {
        ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[1][2]]][0];
        ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[1][2]]][1];
        ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[1][2]]][2];
        ++count;
        }

        if(cn[1][1]>=0 && cn[1][2]>=0)
        {
            if(p->ccpoint[p->facet[n][7+cn[1][1]]][1] < p->ccpoint[p->facet[n][7+cn[1][2]]][1])
            {
            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[1][1]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[1][1]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[1][1]]][2];
            ++count;

            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[1][2]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[1][2]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[1][2]]][2];
            ++count;
            }

            if(p->ccpoint[p->facet[n][7+cn[1][1]]][1] >= p->ccpoint[p->facet[n][7+cn[1][2]]][1])
            {
            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[1][2]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[1][2]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[1][2]]][2];
            ++count;

            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[1][1]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[1][1]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[1][1]]][2];
            ++count;
            }
        }

    // Side 2
        if(cn[2][0]>=0)
        {
        ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[2][0]]][0];
        ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[2][0]]][1];
        ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[2][0]]][2];
        ++count;
        }

        if(cn[2][1]>=0 && cn[2][2]<0)
        {
        ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[2][1]]][0];
        ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[2][1]]][1];
        ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[2][1]]][2];
        ++count;
        }

        if(cn[2][2]>=0 && cn[2][1]<0)
        {
        ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[2][2]]][0];
        ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[2][2]]][1];
        ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[2][2]]][2];
        ++count;
        }

        if(cn[2][1]>=0 && cn[2][2]>=0)
        {
            if(p->ccpoint[p->facet[n][7+cn[2][1]]][0] >= p->ccpoint[p->facet[n][7+cn[2][2]]][0])
            {
            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[2][1]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[2][1]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[2][1]]][2];
            ++count;

            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[2][2]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[2][2]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[2][2]]][2];
            ++count;
            }

            if(p->ccpoint[p->facet[n][7+cn[2][1]]][0] < p->ccpoint[p->facet[n][7+cn[2][2]]][0])
            {
            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[2][2]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[2][2]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[2][2]]][2];
            ++count;

            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[2][1]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[2][1]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[2][1]]][2];
            ++count;
            }
        }

    // Side 3
        if(cn[3][0]>=0)
        {
        ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[3][0]]][0];
        ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[3][0]]][1];
        ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[3][0]]][2];
        ++count;
        }

        if(cn[3][1]>=0 && cn[3][2]<0)
        {
        ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[3][1]]][0];
        ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[3][1]]][1];
        ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[3][1]]][2];
        ++count;
        }

        if(cn[3][2]>=0 && cn[3][1]<0)
        {
        ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[3][2]]][0];
        ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[3][2]]][1];
        ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[3][2]]][2];
        ++count;
        }

        if(cn[3][1]>=0 && cn[3][2]>=0)
        {
            if(p->ccpoint[p->facet[n][7+cn[3][1]]][1] < p->ccpoint[p->facet[n][7+cn[3][2]]][1])
            {
            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[3][1]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[3][1]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[3][1]]][2];
            ++count;

            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[3][2]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[3][2]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[3][2]]][2];
            ++count;
            }

            if(p->ccpoint[p->facet[n][7+cn[3][1]]][1] >= p->ccpoint[p->facet[n][7+cn[3][2]]][1])
            {
            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[3][2]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[3][2]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[3][2]]][2];
            ++count;

            ccpt[fcount][count][0] = p->ccpoint[p->facet[n][7+cn[3][1]]][0];
            ccpt[fcount][count][1] = p->ccpoint[p->facet[n][7+cn[3][1]]][1];
            ccpt[fcount][count][2] = p->ccpoint[p->facet[n][7+cn[3][1]]][2];
            ++count;
            }
        }

    fid[fcount][0] = i;
    fid[fcount][1] = j;
    fid[fcount][2] = k;
	
	fx[fcount] = p->XP[IP];
	fy[fcount] = p->YP[JP];
	fz[fcount] = p->ZP[KP];

    surfnum[fcount] = p->facet[n][6];

    fn[fcount][0] = p->gcn[n][0];
    fn[fcount][1] = p->gcn[n][1];
    fn[fcount][2] = p->gcn[n][2];
    ++fcount;
    }
    }
}


