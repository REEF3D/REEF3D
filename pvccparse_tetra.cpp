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

#include"pvccparse.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

// see p. 372

void pvccparse::cell_tetra(lexer* p, fdm* a, ghostcell *pgc)
{
		for(int qn=0;qn<6;++qn)
		if(fcount[qn]==3)
        {
        a->pvccnode[count][0]=fd[qn][0];
        a->pvccnode[count][1]=fd[qn][1];
        a->pvccnode[count][2]=fd[qn][2];
			
			for(int qnn=qn+1;qnn<6;++qnn)
			if(fcount[qnn]==3)
			for(int qqnn=0;qqnn<3;++qqnn)
			if(fd[qnn][qqnn]!=fd[qn][0] && fd[qnn][qqnn]!=fd[qn][1] && fd[qnn][qqnn]!=fd[qn][2])
			{
			a->pvccnode[count][3]=fd[qnn][qqnn];
			break;
			}
			
			
		a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        count++;
		
		break;
        }
}
