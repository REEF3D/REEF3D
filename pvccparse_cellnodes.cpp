/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"pvccparse.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void pvccparse::cellnodes(lexer* p, fdm* a, ghostcell *pgc)
{

    // identify active and unactive cell nodes
    for(q=0;q<p->facet[n][6];q++)
    {
        cx = p->XP[p->facet[n][0]+marge]-p->originx;
        cy = p->YP[p->facet[n][1]+marge]-p->originy;
        cz = p->ZP[p->facet[n][2]+marge]-p->originz;

        px=p->ccpoint[p->facet[n][7+q]][0] - p->originx;
        py=p->ccpoint[p->facet[n][7+q]][1] - p->originy;
        pz=p->ccpoint[p->facet[n][7+q]][2] - p->originz;

        i=p->facet[n][0];
        j=p->facet[n][1];
        k=p->facet[n][2];

        // 0
        check1=0;
        if(px > cx-0.5*p->DXN[IP]-eps && px < cx-0.5*p->DXN[IP]+eps
        && py > cy-0.5*p->DYN[JP]-eps && py < cy-0.5*p->DYN[JP]+eps
        && pz > cz-0.5*p->DZN[KP]-eps && pz < cz-0.5*p->DZN[KP]+eps)
        {
			pt[0]=1;
			
			if(p->facet[n][3]==1)
            pt[1]=1;

            if(p->facet[n][3]==4)
            pt[1]=0;
			
            if(p->facet[n][4]==3)
            pt[3]=1;

            if(p->facet[n][4]==2)
            pt[3]=0;

			if(p->facet[n][5]==5)
            pt[4]=1;

            if(p->facet[n][5]==6)
            pt[4]=0;
			
        check1=1;
        }


        if(px > cx-0.5*p->DXN[IP]+eps && px < cx+0.5*p->DXN[IP]-eps
        && py > cy-0.5*p->DYN[JP]-eps && py < cy-0.5*p->DYN[JP]+eps
        && pz > cz-0.5*p->DZN[KP]-eps && pz < cz-0.5*p->DZN[KP]+eps
        && check1==0)
        {
            if(p->facet[n][3]==1)
            {
            pt[0]=0;
            pt[1]=1;
            cl[0]=-p->facet[n][7+q]-p->pointnum;
            }

            if(p->facet[n][3]==4)
            {
            pt[0]=1;
            pt[1]=0;
            cl[0]=-p->facet[n][7+q]-p->pointnum;
            }
        }

        // 1
        check2=0;
        if(px > cx+0.5*p->DXN[IP]-eps && px < cx+0.5*p->DXN[IP]+eps
        && py > cy-0.5*p->DYN[JP]-eps && py < cy-0.5*p->DYN[JP]+eps
        && pz > cz-0.5*p->DZN[KP]-eps && pz < cz-0.5*p->DZN[KP]+eps)
        {
			pt[1]=1;
			
			if(p->facet[n][3]==1)
            pt[0]=0;

            if(p->facet[n][3]==4)
            pt[0]=1;
			
            if(p->facet[n][4]==3)
            pt[2]=1;

            if(p->facet[n][4]==2)
            pt[2]=0;

			if(p->facet[n][5]==5)
            pt[5]=1;

            if(p->facet[n][5]==6)
            pt[5]=0;
			
        check2=1;
        }

        if(px > cx+0.5*p->DXN[IP]-eps && px < cx+0.5*p->DXN[IP]+eps
        && py > cy-0.5*p->DYN[JP]+eps && py < cy+0.5*p->DYN[JP]-eps
        && pz > cz-0.5*p->DZN[KP]-eps && pz < cz-0.5*p->DZN[KP]+eps
        && check2==0)
        {
            if(p->facet[n][4]==3)
            {
            pt[1]=0;
            pt[2]=1;
            cl[1]=-p->facet[n][7+q]-p->pointnum;
            }

            if(p->facet[n][4]==2)
            {
            pt[1]=1;
            pt[2]=0;
            cl[1]=-p->facet[n][7+q]-p->pointnum;
            }
        }

        // 2
        check3=0;
        if(px > cx+0.5*p->DXN[IP]-eps && px < cx+0.5*p->DXN[IP]+eps
        && py > cy+0.5*p->DYN[JP]-eps && py < cy+0.5*p->DYN[JP]+eps
        && pz > cz-0.5*p->DZN[KP]-eps && pz < cz-0.5*p->DZN[KP]+eps)
        {
            pt[2]=1;
			
			if(p->facet[n][3]==1)
            pt[3]=0;

            if(p->facet[n][3]==4)
            pt[3]=1;
			
            if(p->facet[n][4]==3)
            pt[1]=0;

            if(p->facet[n][4]==2)
            pt[1]=1;

			if(p->facet[n][5]==5)
            pt[6]=1;

            if(p->facet[n][5]==6)
            pt[6]=0;
			
        check3=1;
        }

        if(px > cx-0.5*p->DXN[IP]+eps && px < cx+0.5*p->DXN[IP]-eps
        && py > cy+0.5*p->DYN[JP]-eps && py < cy+0.5*p->DYN[JP]+eps
        && pz > cz-0.5*p->DZN[KP]-eps && pz < cz-0.5*p->DZN[KP]+eps
        && check3==0)
        {
            if(p->facet[n][3]==1)
            {
            pt[2]=1;
            pt[3]=0;
            cl[2]=-p->facet[n][7+q]-p->pointnum;
            }

            if(p->facet[n][3]==4)
            {
            pt[2]=0;
            pt[3]=1;
            cl[2]=-p->facet[n][7+q]-p->pointnum;
            }
        }


        // 3
        check4=0;
        if(px > cx-0.5*p->DXN[IP]-eps && px < cx-0.5*p->DXN[IP]+eps
        && py > cy+0.5*p->DYN[JP]-eps && py < cy+0.5*p->DYN[JP]+eps
        && pz > cz-0.5*p->DZN[KP]-eps && pz < cz-0.5*p->DZN[KP]+eps)
        {
            pt[3]=1;
			
			if(p->facet[n][3]==1)
            pt[2]=1;

            if(p->facet[n][3]==4)
            pt[2]=0;
			
            if(p->facet[n][4]==3)
            pt[0]=0;

            if(p->facet[n][4]==2)
            pt[0]=1;

			if(p->facet[n][5]==5)
            pt[7]=1;

            if(p->facet[n][5]==6)
            pt[7]=0;
			
        check4=1;
        }


        if(px > cx-0.5*p->DXN[IP]-eps && px < cx-0.5*p->DXN[IP]+eps
        && py > cy-0.5*p->DYN[JP]+eps && py < cy+0.5*p->DYN[JP]-eps
        && pz > cz-0.5*p->DZN[KP]-eps && pz < cz-0.5*p->DZN[KP]+eps
        && check4==0)
        {
            if(p->facet[n][4]==3)
            {
            pt[3]=1;
            pt[0]=0;
            cl[3]=-p->facet[n][7+q]-p->pointnum;
            }

            if(p->facet[n][4]==2)
            {
            pt[3]=0;
            pt[0]=1;
            cl[3]=-p->facet[n][7+q]-p->pointnum;
            }
        }

//-----------------------------------

        // 4
        if(px > cx-0.5*p->DXN[IP]-eps && px < cx-0.5*p->DXN[IP]+eps
        && py > cy-0.5*p->DYN[JP]-eps && py < cy-0.5*p->DYN[JP]+eps
        && pz > cz-0.5*p->DZN[KP]+eps && pz < cz+0.5*p->DZN[KP]-eps)
        {
            if(p->facet[n][5]==5)
            {
            pt[0]=0;
            pt[4]=1;
            cl[4]=-p->facet[n][7+q]-p->pointnum;
            }

            if(p->facet[n][5]==6)
            {
            pt[0]=1;
            pt[4]=0;
            cl[4]=-p->facet[n][7+q]-p->pointnum;
            }
        }

        // 5
        if(px > cx+0.5*p->DXN[IP]-eps && px < cx+0.5*p->DXN[IP]+eps
        && py > cy-0.5*p->DYN[JP]-eps && py < cy-0.5*p->DYN[JP]+eps
        && pz > cz-0.5*p->DZN[KP]+eps && pz < cz+0.5*p->DZN[KP]-eps)
        {
            if(p->facet[n][5]==5)
            {
            pt[1]=0;
            pt[5]=1;
            cl[5]=-p->facet[n][7+q]-p->pointnum;
            }

            if(p->facet[n][5]==6)
            {
            pt[1]=1;
            pt[5]=0;
            cl[5]=-p->facet[n][7+q]-p->pointnum;
            }
        }

        // 6
        if(px > cx+0.5*p->DXN[IP]-eps && px < cx+0.5*p->DXN[IP]+eps
        && py > cy+0.5*p->DYN[JP]-eps && py < cy+0.5*p->DYN[JP]+eps
        && pz > cz-0.5*p->DZN[KP]+eps && pz < cz+0.5*p->DZN[KP]-eps)
        {
            if(p->facet[n][5]==5)
            {
            pt[2]=0;
            pt[6]=1;
            cl[6]=-p->facet[n][7+q]-p->pointnum;
            }

            if(p->facet[n][5]==6)
            {
            pt[2]=1;
            pt[6]=0;
            cl[6]=-p->facet[n][7+q]-p->pointnum;
            }
        }


        // 7
        if(px > cx-0.5*p->DXN[IP]-eps && px < cx-0.5*p->DXN[IP]+eps
        && py > cy+0.5*p->DYN[JP]-eps && py < cy+0.5*p->DYN[JP]+eps
        && pz > cz-0.5*p->DZN[KP]+eps && pz < cz+0.5*p->DZN[KP]-eps)
        {
            if(p->facet[n][5]==5)
            {
            pt[3]=0;
            pt[7]=1;
            cl[7]=-p->facet[n][7+q]-p->pointnum;
            }

            if(p->facet[n][5]==6)
            {
            pt[3]=1;
            pt[7]=0;
            cl[7]=-p->facet[n][7+q]-p->pointnum;
            }
        }

//-----------------------------------

        // 8
        check9=0;
        if(px > cx-0.5*p->DXN[IP]-eps && px < cx-0.5*p->DXN[IP]+eps
        && py > cy-0.5*p->DYN[JP]-eps && py < cy-0.5*p->DYN[JP]+eps
        && pz > cz+0.5*p->DZN[KP]-eps && pz < cz+0.5*p->DZN[KP]+eps)
        {
            pt[4]=1;
			
			if(p->facet[n][3]==1)
            pt[5]=1;

            if(p->facet[n][3]==4)
            pt[5]=0;
			
            if(p->facet[n][4]==3)
            pt[7]=1;

            if(p->facet[n][4]==2)
            pt[7]=0;

			if(p->facet[n][5]==5)
            pt[0]=0;

            if(p->facet[n][5]==6)
            pt[0]=1;
			
        check9=1;
        }

        if(px > cx-0.5*p->DXN[IP]+eps && px < cx+0.5*p->DXN[IP]-eps
        && py > cy-0.5*p->DYN[JP]-eps && py < cy-0.5*p->DYN[JP]+eps
        && pz > cz+0.5*p->DZN[KP]-eps && pz < cz+0.5*p->DZN[KP]+eps
        && check9==0)
        {
            if(p->facet[n][3]==1)
            {
            pt[4]=0;
            pt[5]=1;
            cl[8]=-p->facet[n][7+q]-p->pointnum;
            }

            if(p->facet[n][3]==4)
            {
            pt[4]=1;
            pt[5]=0;
            cl[8]=-p->facet[n][7+q]-p->pointnum;
            }
        }


        // 9
        check10=0;
        if(px > cx+0.5*p->DXN[IP]-eps && px < cx+0.5*p->DXN[IP]+eps
        && py > cy-0.5*p->DYN[JP]-eps && py < cy-0.5*p->DYN[JP]+eps
        && pz > cz+0.5*p->DZN[KP]-eps && pz < cz+0.5*p->DZN[KP]+eps)
        {
			pt[5]=1;
			
			if(p->facet[n][3]==1)
            pt[4]=0;

            if(p->facet[n][3]==4)
            pt[4]=1;
			
            if(p->facet[n][4]==3)
            pt[6]=1;

            if(p->facet[n][4]==2)
            pt[6]=0;

			if(p->facet[n][5]==5)
            pt[1]=0;

            if(p->facet[n][5]==6)
            pt[1]=1;

        check10=1;
        }

        if(px > cx+0.5*p->DXN[IP]-eps && px < cx+0.5*p->DXN[IP]+eps
        && py > cy-0.5*p->DYN[JP]+eps && py < cy+0.5*p->DYN[JP]-eps
        && pz > cz+0.5*p->DZN[KP]-eps && pz < cz+0.5*p->DZN[KP]+eps
        && check10==0)
        {
            if(p->facet[n][4]==3)
            {
            pt[5]=0;
            pt[6]=1;
            cl[9]=-p->facet[n][7+q]-p->pointnum;
            }

            if(p->facet[n][4]==2)
            {
            pt[5]=1;
            pt[6]=0;
            cl[9]=-p->facet[n][7+q]-p->pointnum;
            }
        }

        // 10
        check11=0;
        if(px > cx+0.5*p->DXN[IP]-eps && px < cx+0.5*p->DXN[IP]+eps
        && py > cy+0.5*p->DYN[JP]-eps && py < cy+0.5*p->DYN[JP]+eps
        && pz > cz+0.5*p->DZN[KP]-eps && pz < cz+0.5*p->DZN[KP]+eps)
        {
            pt[6]=1;
			
			if(p->facet[n][3]==1)
            pt[7]=0;

            if(p->facet[n][3]==4)
            pt[7]=1;
			
            if(p->facet[n][4]==3)
            pt[5]=0;

            if(p->facet[n][4]==2)
            pt[5]=1;

			if(p->facet[n][5]==5)
            pt[2]=0;

            if(p->facet[n][5]==6)
            pt[2]=1;
			
        check11=1;
        }

        if(px > cx-0.5*p->DXN[IP]+eps && px < cx+0.5*p->DXN[IP]-eps
        && py > cy+0.5*p->DYN[JP]-eps && py < cy+0.5*p->DYN[JP]+eps
        && pz > cz+0.5*p->DZN[KP]-eps && pz < cz+0.5*p->DZN[KP]+eps
        && check11==0)
        {
            if(p->facet[n][3]==1)
            {
            pt[6]=1;
            pt[7]=0;
            cl[10]=-p->facet[n][7+q]-p->pointnum;
            }

            if(p->facet[n][3]==4)
            {
            pt[6]=0;
            pt[7]=1;
            cl[10]=-p->facet[n][7+q]-p->pointnum;
            }
        }


        // 11
        check12=0;
        if(px > cx-0.5*p->DXN[IP]-eps && px < cx-0.5*p->DXN[IP]+eps
        && py > cy+0.5*p->DYN[JP]-eps && py < cy+0.5*p->DYN[JP]+eps
        && pz > cz+0.5*p->DZN[KP]-eps && pz < cz+0.5*p->DZN[KP]+eps)
        {
            pt[7]=1;
			
			if(p->facet[n][3]==1)
            pt[6]=1;

            if(p->facet[n][3]==4)
            pt[6]=0;
			
            if(p->facet[n][4]==3)
            pt[4]=0;

            if(p->facet[n][4]==2)
            pt[4]=1;

			if(p->facet[n][5]==5)
            pt[3]=0;

            if(p->facet[n][5]==6)
            pt[3]=1;
			
        check12=1;
        }

        if(px > cx-0.5*p->DXN[IP]-eps && px < cx-0.5*p->DXN[IP]+eps
        && py > cy-0.5*p->DYN[JP]+eps && py < cy+0.5*p->DYN[JP]-eps
        && pz > cz+0.5*p->DZN[KP]-eps && pz < cz+0.5*p->DZN[KP]+eps
        && check12==0)
        {
            if(p->facet[n][4]==3)
            {
            pt[7]=1;
            pt[4]=0;
            cl[11]=-p->facet[n][7+q]-p->pointnum;
            }

            if(p->facet[n][4]==2)
            {
            pt[7]=0;
            pt[4]=1;
            cl[11]=-p->facet[n][7+q]-p->pointnum;
            }
        }
    }
}
