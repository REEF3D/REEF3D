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

void pvccparse::cell_wedge(lexer* p, fdm* a, ghostcell *pgc)
{	

		if(fcount[0]==3 && fcount[1]==3)
        {
        a->pvccnode[count][0]=fd[0][0];
        a->pvccnode[count][1]=fd[0][1];
        a->pvccnode[count][2]=fd[0][2];

		a->pvccnode[count][3]=fd[1][0];
		a->pvccnode[count][4]=fd[1][1];
		a->pvccnode[count][5]=fd[1][2];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        count++;
        }
		
		if(fcount[2]==3 && fcount[3]==3)
        {
        a->pvccnode[count][0]=fd[2][0];
        a->pvccnode[count][1]=fd[2][1];
        a->pvccnode[count][2]=fd[2][2];

		a->pvccnode[count][3]=fd[3][0];
		a->pvccnode[count][4]=fd[3][1];
		a->pvccnode[count][5]=fd[3][2];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        count++;
        }
		
		if(fcount[4]==3 && fcount[5]==3)
        {
        a->pvccnode[count][0]=fd[4][0];
        a->pvccnode[count][1]=fd[4][1];
        a->pvccnode[count][2]=fd[4][2];

		a->pvccnode[count][3]=fd[5][0];
		a->pvccnode[count][4]=fd[5][1];
		a->pvccnode[count][5]=fd[5][2];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        count++;
        }

}

