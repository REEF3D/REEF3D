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
for more da->solidils.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void sixdof_obj::triangulation(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
	int negcount, poscount;
    double eps;
    
    NDBASELOOP
    fsf[IJK] = 0.125*(d->FB[IJK] + d->FB[Ip1JK] + d->FB[IJp1K] + d->FB[Ip1Jp1K]
                    + d->FB[IJKp1] + d->FB[Ip1JKp1] + d->FB[IJp1Kp1] + d->FB[Ip1Jp1Kp1]);
	
    NDBASELOOP
    vert[IJK]=-1;

    NDBASELOOP
    nflag[IJK]=0;
	

    BASELOOP
    {
        eps = interfac*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]*d->WL(i,j));
        
        if(fabs(d->FB[IJK])<eps)
        {
            check=1;

            if(fsf[IJK]<zero && fsf[Im1JK]<zero && fsf[Im1Jm1K]<zero && fsf[IJm1K]<zero &&
               fsf[IJKm1]<zero && fsf[Im1JKm1]<zero && fsf[Im1Jm1Km1]<zero && fsf[IJm1Km1]<zero)
            check=0;
            
            if(fsf[IJK]>zero && fsf[Im1JK]>zero && fsf[Im1Jm1K]>zero && fsf[IJm1K]>zero &&
               fsf[IJKm1]>zero && fsf[Im1JKm1]>zero && fsf[Im1Jm1Km1]>zero && fsf[IJm1Km1]>zero)
            check=0;

            if(check==1)
            {
            nflag[IJK]=1;
            nflag[Im1JK]=1;
            nflag[Im1Jm1K]=1;
            nflag[IJm1K]=1;
            nflag[IJKm1]=1;
            nflag[Im1JKm1]=1;
            nflag[Im1Jm1Km1]=1;
            nflag[IJm1Km1]=1;
            }
        }
    }

	//--------------------
    countM=0;
    NDBASELOOP
    if(nflag[IJK]==1)
    ++countM;

    numtri = 6*countM;
    numvert = countM;

    numtri_mem = numtri;
    numvert_mem = numvert;

    // allocate -------
    allocate(p,d,pgc);

    countM=0;
    NDBASELOOP
    if(nflag[IJK]==1)
    {
    pt[countM][0] = p->posnode_x();
    pt[countM][1] = p->posnode_y();
    pt[countM][2] = p->posnode_z();

    ls[countM] = fsf[IJK];

    vert[IJK] = countM;

    ++countM;
    }
    
	// p. 725, 956
    count=0;
    BASELOOP
    if(nflag[IJK]==1)
    if(nflag[Im1JK]==1)
    if(nflag[Im1Jm1K]==1)
    if(nflag[IJm1K]==1)
    if(nflag[IJKm1]==1)
    if(nflag[Im1JKm1]==1)
    if(nflag[Im1Jm1Km1]==1)
    if(nflag[IJm1Km1]==1)
    {
    // 1
    tri[count][0] = vert[Im1Jm1Km1];
    tri[count][1] = vert[Im1JKm1];
    tri[count][2] = vert[Im1Jm1K];
    tri[count][3] = vert[IJm1K];
    ++count;

    // 2
    tri[count][0] = vert[Im1Jm1Km1];
    tri[count][1] = vert[IJm1Km1];
    tri[count][2] = vert[Im1JKm1];
    tri[count][3] = vert[IJm1K];
    ++count;

    // 3
    tri[count][0] = vert[Im1JKm1];
    tri[count][1] = vert[IJKm1];
    tri[count][2] = vert[IJm1Km1];
    tri[count][3] = vert[IJm1K];
    ++count;

    // 4
    tri[count][0] = vert[IJKm1];
    tri[count][1] = vert[Im1JKm1];
    tri[count][2] = vert[IJm1K];
    tri[count][3] = vert[IJK];
    ++count;

    // 5
	tri[count][0] = vert[Im1JKm1];
    tri[count][1] = vert[Im1JK];
    tri[count][2] = vert[IJK];
    tri[count][3] = vert[IJm1K];
    ++count;

    // 6
    tri[count][0] = vert[Im1JKm1];
    tri[count][1] = vert[Im1Jm1K];
    tri[count][2] = vert[IJm1K];
    tri[count][3] = vert[Im1JK];
    ++count;
    }
	
    numtri=count;
}
