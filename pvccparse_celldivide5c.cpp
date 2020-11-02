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

void pvccparse::cell_divide5c(lexer* p, fdm* a, ghostcell *pgc)
{
		// cell 1
        if(pt[1]==1 && pt[2]==1 && pt[3]==1
        && cl[0]<=0 && cl[3]<=0 && cl[5]<=0 && cl[7]<=0 && (cl[6]<=0 || pt[6]==1))
        {
        // wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[3];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;

        // wedge 2
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
        }
	
}
