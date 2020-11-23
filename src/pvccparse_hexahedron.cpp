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

// see p. 372

void pvccparse::cell_hexahedron(lexer* p, fdm* a, ghostcell *pgc)
{
		int check = 0;
	
	
		if(fcount[0]==4 && fcount[1]==4)
        {
        a->pvccnode[count][0]=fd[0][0];
        a->pvccnode[count][1]=fd[0][1];
        a->pvccnode[count][2]=fd[0][2];
		a->pvccnode[count][3]=fd[0][3];

		a->pvccnode[count][4]=fd[1][0];
		a->pvccnode[count][5]=fd[1][1];
		a->pvccnode[count][6]=fd[1][2];
		a->pvccnode[count][7]=fd[1][3];

        a->ccedge[count]=8;

        count++;
		check=1;
        }
		
		if(fcount[2]==4 && fcount[3]==4 && check==0)
        {
        a->pvccnode[count][0]=fd[2][0];
        a->pvccnode[count][1]=fd[2][1];
        a->pvccnode[count][2]=fd[2][2];
		a->pvccnode[count][3]=fd[2][3];

		a->pvccnode[count][4]=fd[3][0];
		a->pvccnode[count][5]=fd[3][1];
		a->pvccnode[count][6]=fd[3][2];
		a->pvccnode[count][7]=fd[3][3];

        a->ccedge[count]=8;

        count++;
		check=1;
        }
		
		if(fcount[4]==4 && fcount[5]==4 && check==0)
        {
        a->pvccnode[count][0]=fd[4][0];
        a->pvccnode[count][1]=fd[4][1];
        a->pvccnode[count][2]=fd[4][2];
		a->pvccnode[count][3]=fd[4][3];

		a->pvccnode[count][4]=fd[5][0];
		a->pvccnode[count][5]=fd[5][1];
		a->pvccnode[count][6]=fd[5][2];
		a->pvccnode[count][7]=fd[5][3];

        a->ccedge[count]=8;

        count++;
		check=1;
        }
	

}
