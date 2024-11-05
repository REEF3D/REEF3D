/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"nhflow_force.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"


void nhflow_force::triangulation(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
	int negcount, poscount;
    
    NDBASELOOP
    eta[IJK] = 0.125*(d->SOLID[IJK] + d->SOLID[Ip1JK] + d->SOLID[IJp1K] + d->SOLID[Ip1Jp1K]
                      + d->SOLID[IJKp1] + d->SOLID[Ip1JKp1] + d->SOLID[IJp1Kp1] + d->SOLID[Ip1Jp1Kp1]);
	
    NDBASELOOP
    vertice[IJK]=-1;

    NDBASELOOP
    nodeflag[IJK]=0;
	

    BASELOOP
    if(i>=is && i<=ie && j>=js && j<=je && k>=ks && k<=ke)
    {
        epsi = interfac*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]*d->WL(i,j));
        
        if(fabs(d->SOLID[IJK])<epsi)
        {
            check=1;

            if(eta[IJK]<zero && eta[Im1JK]<zero && eta[Im1Jm1K]<zero && eta[IJm1K]<zero &&
               eta[IJKm1]<zero && eta[Im1JKm1]<zero && eta[Im1Jm1Km1]<zero && eta[IJm1Km1]<zero)
            check=0;
            
            if(eta[IJK]>zero && eta[Im1JK]>zero && eta[Im1Jm1K]>zero && eta[IJm1K]>zero &&
               eta[IJKm1]>zero && eta[Im1JKm1]>zero && eta[Im1Jm1Km1]>zero && eta[IJm1Km1]>zero)
            check=0;

            if(check==1)
            {
            nodeflag[IJK]=1;
            nodeflag[Im1JK]=1;
            nodeflag[Im1Jm1K]=1;
            nodeflag[IJm1K]=1;
            nodeflag[IJKm1]=1;
            nodeflag[Im1JKm1]=1;
            nodeflag[Im1Jm1Km1]=1;
            nodeflag[IJm1Km1]=1;
            }
        }
    }

	//--------------------
    countM=0;
    NDBASELOOP
    if(nodeflag[IJK]==1)
    ++countM;

    numtri = 6*countM;
    numvert = countM;

    numtri_mem = numtri;
    numvert_mem = numvert;

    // allocate -------
    allocate(p,d,pgc);

    countM=0;
    NDBASELOOP
    if(nodeflag[IJK]==1)
    {
    pt[countM][0] = p->posnode_x();
    pt[countM][1] = p->posnode_y();
    pt[countM][2] = p->posnode_z();

    ls[countM] = eta[IJK];

    vertice[IJK] = countM;

    ++countM;
    }
    
	// p. 725, 956
    count=0;
    BASELOOP
    if(nodeflag[IJK]==1)
    if(nodeflag[Im1JK]==1)
    if(nodeflag[Im1Jm1K]==1)
    if(nodeflag[IJm1K]==1)
    if(nodeflag[IJKm1]==1)
    if(nodeflag[Im1JKm1]==1)
    if(nodeflag[Im1Jm1Km1]==1)
    if(nodeflag[IJm1Km1]==1)
    {
    // 1
    tri[count][0] = vertice[Im1Jm1Km1];
    tri[count][1] = vertice[Im1JKm1];
    tri[count][2] = vertice[Im1Jm1K];
    tri[count][3] = vertice[IJm1K];
    ++count;

    // 2
    tri[count][0] = vertice[Im1Jm1Km1];
    tri[count][1] = vertice[IJm1Km1];
    tri[count][2] = vertice[Im1JKm1];
    tri[count][3] = vertice[IJm1K];
    ++count;

    // 3
    tri[count][0] = vertice[Im1JKm1];
    tri[count][1] = vertice[IJKm1];
    tri[count][2] = vertice[IJm1Km1];
    tri[count][3] = vertice[IJm1K];
    ++count;

    // 4
    tri[count][0] = vertice[IJKm1];
    tri[count][1] = vertice[Im1JKm1];
    tri[count][2] = vertice[IJm1K];
    tri[count][3] = vertice[IJK];
    ++count;

    // 5
	tri[count][0] = vertice[Im1JKm1];
    tri[count][1] = vertice[Im1JK];
    tri[count][2] = vertice[IJK];
    tri[count][3] = vertice[IJm1K];
    ++count;

    // 6
    tri[count][0] = vertice[Im1JKm1];
    tri[count][1] = vertice[Im1Jm1K];
    tri[count][2] = vertice[IJm1K];
    tri[count][3] = vertice[Im1JK];
    ++count;
    }
	
    numtri=count;
}
