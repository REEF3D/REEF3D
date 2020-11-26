/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

void pvccparse::cell_divide4(lexer* p, fdm* a, ghostcell *pgc)
{
    //Paraview
    if(pcount>8 && pointcount==6 && clcount==4)
    {
// Z
        // 0 --------------------------------------------
        if(pt[0]==0 && pt[4]==0)
        {

        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[2];

        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        count++;

        // 2nd cell
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[4];

        a->pvccnode[count][4]=nd[5];
        a->pvccnode[count][5]=nd[7];
        a->pvccnode[count][6]=nd[8];
        a->pvccnode[count][7]=nd[9];

        a->ccedge[count]=8;

        count++;
        }

        // 1 --------------------------------------------
        if(pt[1]==0 && pt[5]==0)
        {

        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[4];

        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[9];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        count++;

        // 2nd cell
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[4];

        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[7];
        a->pvccnode[count][6]=nd[8];
        a->pvccnode[count][7]=nd[9];

        a->ccedge[count]=8;

        count++;
        }


        // 2 --------------------------------------------
        if(pt[2]==0 && pt[6]==0)
        {

        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[2];

        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        count++;

        // 2nd cell
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[4];

        a->pvccnode[count][4]=nd[5];
        a->pvccnode[count][5]=nd[7];
        a->pvccnode[count][6]=nd[8];
        a->pvccnode[count][7]=nd[9];

        a->ccedge[count]=8;

        count++;
        }


        // 3 --------------------------------------------
        if(pt[3]==0 && pt[7]==0)
        {

        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[4];

        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[9];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;


        count++;

        // 2nd cell
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[3];

        a->pvccnode[count][4]=nd[9];
        a->pvccnode[count][5]=nd[6];
        a->pvccnode[count][6]=nd[7];
        a->pvccnode[count][7]=nd[8];

        a->ccedge[count]=8;

        count++;
        }

 // Y
        // 0 --------------------------------------------
        if(pt[0]==0 && pt[3]==0)
        {

        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[7];

        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[3];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        count++;

        // 2nd cell
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[5];

        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[7];
        a->pvccnode[count][6]=nd[8];
        a->pvccnode[count][7]=nd[9];

        a->ccedge[count]=8;

        count++;
        }

        // 1 --------------------------------------------
        if(pt[1]==0 && pt[2]==0)
        {

        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[7];

        a->pvccnode[count][3]=nd[9];
        a->pvccnode[count][4]=nd[3];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        count++;

        // 2nd cell
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[5];

        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[2];
        a->pvccnode[count][6]=nd[3];
        a->pvccnode[count][7]=nd[9];

        a->ccedge[count]=8;

        count++;
        }


        // 2 --------------------------------------------
        if(pt[5]==0 && pt[6]==0)
        {

        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[4];

        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[8];
        a->pvccnode[count][5]=nd[9];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        count++;

        // 2nd cell
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[3];

        a->pvccnode[count][4]=nd[5];
        a->pvccnode[count][5]=nd[6];
        a->pvccnode[count][6]=nd[7];
        a->pvccnode[count][7]=nd[8];

        a->ccedge[count]=8;

        count++;
        }


        // 3 --------------------------------------------
        if(pt[4]==0 && pt[7]==0)
        {

        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[4];

        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[9];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;


        count++;

        // 2nd cell
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[9];

        a->pvccnode[count][4]=nd[5];
        a->pvccnode[count][5]=nd[6];
        a->pvccnode[count][6]=nd[7];
        a->pvccnode[count][7]=nd[8];

        a->ccedge[count]=8;

        count++;
        }

    // X
        // 0 --------------------------------------------
        if(pt[0]==0 && pt[1]==0)
        {

        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[9];

        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        count++;

        // 2nd cell
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[5];

        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[7];
        a->pvccnode[count][6]=nd[8];
        a->pvccnode[count][7]=nd[9];

        a->ccedge[count]=8;

        count++;
        }

        // 1 --------------------------------------------
        if(pt[4]==0 && pt[5]==0)
        {

        a->pvccnode[count][0]=nd[3];
        a->pvccnode[count][1]=nd[0];
        a->pvccnode[count][2]=nd[4];

        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=nd[1];
        a->pvccnode[count][5]=nd[5];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        count++;

        // 2nd cell
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[3];

        a->pvccnode[count][4]=nd[9];
        a->pvccnode[count][5]=nd[6];
        a->pvccnode[count][6]=nd[7];
        a->pvccnode[count][7]=nd[8];

        a->ccedge[count]=8;

        count++;
        }


        // 2 --------------------------------------------
        if(pt[7]==0 && pt[6]==0)
        {

        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[9];

        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=nd[5];
        a->pvccnode[count][5]=nd[6];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        count++;

        // 2nd cell
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[3];

        a->pvccnode[count][4]=nd[9];
        a->pvccnode[count][5]=nd[6];
        a->pvccnode[count][6]=nd[7];
        a->pvccnode[count][7]=nd[8];

        a->ccedge[count]=8;

        count++;
        }


        // 3 --------------------------------------------
        if(pt[3]==0 && pt[2]==0)
        {

        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[9];

        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;


        count++;

        // 2nd cell
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[5];

        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[7];
        a->pvccnode[count][6]=nd[3];
        a->pvccnode[count][7]=nd[4];

        a->ccedge[count]=8;

        count++;
        }

    }

}


