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

void pvccparse::cell_divide6(lexer* p, fdm* a, ghostcell *pgc)
{
        // cell 1
        if(pt[0]==1 && pt[2]==1 && pt[3]==1 && pt[7]==1
        && cl[0]<=0 && cl[1]<=0 && cl[4]<=0 && cl[6]<=0 && cl[10]<=0 && cl[11]<=0)
        {
        p->ccpoint[p->ccptnum][0]=p->XN[IP];
        p->ccpoint[p->ccptnum][1]=p->YN[JP1];
        p->ccpoint[p->ccptnum][2]=p->ZP[KP];
        ++p->ccptnum;

        // wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[6];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 2
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][5]=nd[6];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 3
        a->pvccnode[count][0]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[6];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=nd[9];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
        }


        // cell 2
        if(pt[0]==1 && pt[1]==1 && pt[3]==1 && pt[4]==1
        && cl[1]<=0 && cl[2]<=0 && cl[5]<=0 && cl[7]<=0 && cl[8]<=0 && cl[11]<=0)
        {
        p->ccpoint[p->ccptnum][0]=p->XN[IP];
        p->ccpoint[p->ccptnum][1]=p->YN[JP];
        p->ccpoint[p->ccptnum][2]=p->ZP[KP];
        ++p->ccptnum;

        // wedge 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[8];
        a->pvccnode[count][5]=nd[3];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 2
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[0];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 3
        a->pvccnode[count][0]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[8];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[9];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
        }


        // cell 3
        if(pt[0]==1 && pt[1]==1 && pt[2]==1 && pt[5]==1
        && cl[2]<=0 && cl[3]<=0 && cl[4]<=0 && cl[6]<=0 && cl[8]<=0 && cl[9]<=0)
        {
        p->ccpoint[p->ccptnum][0]=p->XN[IP1];
        p->ccpoint[p->ccptnum][1]=p->YN[JP];
        p->ccpoint[p->ccptnum][2]=p->ZP[KP];
        ++p->ccptnum;

        // wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=nd[9];
        a->pvccnode[count][5]=nd[3];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 2
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][5]=nd[9];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 3
        a->pvccnode[count][0]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][1]=nd[9];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[7];
        a->pvccnode[count][4]=nd[8];
        a->pvccnode[count][5]=nd[6];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
        }


        // cell 4
        if(pt[1]==1 && pt[2]==1 && pt[3]==1 && pt[6]==1
        && cl[0]<=0 && cl[3]<=0 && cl[5]<=0 && cl[7]<=0 && cl[9]<=0 && cl[10]<=0)
        {
        p->ccpoint[p->ccptnum][0]=p->XN[IP1];
        p->ccpoint[p->ccptnum][1]=p->YN[JP1];
        p->ccpoint[p->ccptnum][2]=p->ZP[KP];
        ++p->ccptnum;

        // wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[3];
        a->pvccnode[count][5]=nd[9];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 2
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][5]=nd[9];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 3
        a->pvccnode[count][0]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][1]=nd[9];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[7];
        a->pvccnode[count][4]=nd[8];
        a->pvccnode[count][5]=nd[6];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
        }

  // -----------------------------------------------

        // cell 5
        if(pt[1]==1 && pt[4]==1 && pt[5]==1 && pt[6]==1
        && cl[0]<=0 && cl[1]<=0 && cl[4]<=0 && cl[6]<=0 && cl[10]<=0 && cl[11]<=0)
        {
        p->ccpoint[p->ccptnum][0]=p->XN[IP1];
        p->ccpoint[p->ccptnum][1]=p->YN[JP];
        p->ccpoint[p->ccptnum][2]=p->ZP[KP];
        ++p->ccptnum;

        // wedge 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[0];
        a->pvccnode[count][4]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][5]=nd[4];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 2
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 3
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[9];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
        }


        // cell 6
        if(pt[2]==1 && pt[5]==1 && pt[6]==1 && pt[7]==1
        && cl[1]<=0 && cl[2]<=0 && cl[5]<=0 && cl[7]<=0 && cl[8]<=0 && cl[11]<=0)
        {
        p->ccpoint[p->ccptnum][0]=p->XN[IP1];
        p->ccpoint[p->ccptnum][1]=p->YN[JP1];
        p->ccpoint[p->ccptnum][2]=p->ZP[KP];
        ++p->ccptnum;

        // wedge 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[0];
        a->pvccnode[count][4]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][5]=nd[4];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 2
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[6];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 3
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[8];
        a->pvccnode[count][5]=nd[9];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
        }


        // cell 7
        if(pt[3]==1 && pt[4]==1 && pt[6]==1 && pt[7]==1
        && cl[2]<=0 && cl[3]<=0 && cl[4]<=0 && cl[6]<=0 && cl[8]<=0 && cl[9]<=0)
        {
        p->ccpoint[p->ccptnum][0]=p->XN[IP];
        p->ccpoint[p->ccptnum][1]=p->YN[JP1];
        p->ccpoint[p->ccptnum][2]=p->ZP[KP];
        ++p->ccptnum;

        // wedge 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][5]=nd[0];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 2
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][2]=nd[0];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=nd[9];
        a->pvccnode[count][5]=nd[5];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 3
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
        }


        // cell 8
        if(pt[0]==1 && pt[4]==1 && pt[5]==1 && pt[7]==1
        && cl[0]<=0 && cl[3]<=0 && cl[5]<=0 && cl[7]<=0 && cl[9]<=0 && cl[10]<=0)
        {
        p->ccpoint[p->ccptnum][0]=p->XN[IP];
        p->ccpoint[p->ccptnum][1]=p->YN[JP];
        p->ccpoint[p->ccptnum][2]=p->ZP[KP];
        ++p->ccptnum;

        // wedge 1
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[0];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][5]=nd[2];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 2
        a->pvccnode[count][0]=nd[3];
        a->pvccnode[count][1]=-p->ccptnum+1-p->pointnum;
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[9];
        a->pvccnode[count][4]=nd[5];
        a->pvccnode[count][5]=nd[6];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 3
        a->pvccnode[count][0]=nd[3];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[9];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[6];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
        }

}
