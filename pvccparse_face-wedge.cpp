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

void pvccparse::face_wedge(lexer* p, fdm* a, ghostcell *pgc)
{		
	for(int qn=0;qn<6;++qn)
	{
    fcount[qn]=0;
	for(int qnn=0;qnn<6;++qnn)
	fd[qn][qnn]=0;
	}

// Face 0

    // 0
    if(pt[0]==1)
    {
    fd[0][fcount[0]]=a->nodeval(i-1,j-1,k-1)-1;
    ++fcount[0];
    }

    if(cl[0]<=0)
    {
    fd[0][fcount[0]]=cl[0];
    ++fcount[0];
    }
	
	if(pt[0]==0 && cl[0]>0 && cl[4]<=0)
    {
    fd[0][fcount[0]]=cl[4];
    ++fcount[0];
    }

    // 1
    if(pt[1]==1)
    {
    fd[0][fcount[0]]=a->nodeval(i,j-1,k-1)-1;
    ++fcount[0];
    }

    if(cl[5]<=0)
    {
    fd[0][fcount[0]]=cl[5];
    ++fcount[0];
    }
	
	// 5
    if(pt[5]==1)
    {
    fd[0][fcount[0]]=a->nodeval(i,j-1,k)-1;
    ++fcount[0];
    }

    if(cl[8]<=0)
    {
    fd[0][fcount[0]]=cl[8];
    ++fcount[0];
    }

	// 4
    if(pt[4]==1)
    {
    fd[0][fcount[0]]=a->nodeval(i-1,j-1,k)-1;
    ++fcount[0];
    }

    if(cl[4]<=0 && (pt[0]==1 || cl[0]<=0))
    {
    fd[0][fcount[0]]=cl[4];
    ++fcount[0];
    }
	
// Face 1

	// 3
    if(pt[3]==1)
    {
    fd[1][fcount[1]]=a->nodeval(i-1,j,k-1)-1;
    ++fcount[1];
    }

    if(cl[2]<=0)
    {
    fd[1][fcount[1]]=cl[2];
    ++fcount[1];
    }
	
	if(pt[3]==0 && cl[2]>0 && cl[7]<=0)
    {
    fd[1][fcount[1]]=cl[7];
    ++fcount[1];
    }
	
	// 2
    if(pt[2]==1)
    {
    fd[1][fcount[1]]=a->nodeval(i,j,k-1)-1;
    ++fcount[1];
    }
	
    if(cl[6]<=0)
    {
    fd[1][fcount[1]]=cl[6];
    ++fcount[1];
    }

	// 6
    if(pt[6]==1)
    {
    fd[1][fcount[1]]=a->nodeval(i,j,k)-1;
    ++fcount[1];
    }

    if(cl[10]<=0)
    {
    fd[1][fcount[1]]=cl[10];
    ++fcount[1];
    }

	// 7
    if(pt[7]==1)
    {
    fd[1][fcount[1]]=a->nodeval(i-1,j,k)-1;
    ++fcount[1];
    }

    if(cl[7]<=0 && (pt[3]==1 || cl[2]<=0))
    {
    fd[1][fcount[1]]=cl[7];
    ++fcount[1];
    }
	
// Face 2
	
	// 1
    if(pt[1]==1)
    {
    fd[2][fcount[2]]=a->nodeval(i,j-1,k-1)-1;
    ++fcount[2];
    }

    if(cl[1]<=0)
    {
    fd[2][fcount[2]]=cl[1];
    ++fcount[2];
    }
	
	if(pt[1]==0 && cl[1]>0 && cl[5]<=0)
    {
    fd[2][fcount[2]]=cl[5];
    ++fcount[2];
    }
	
	// 2
    if(pt[2]==1)
    {
    fd[2][fcount[2]]=a->nodeval(i,j,k-1)-1;
    ++fcount[2];
    }

    if(cl[6]<=0)
    {
    fd[2][fcount[2]]=cl[6];
    ++fcount[2];
    }
	
	// 6
    if(pt[6]==1)
    {
    fd[2][fcount[2]]=a->nodeval(i,j,k)-1;
    ++fcount[2];
    }

    if(cl[9]<=0)
    {
    fd[2][fcount[2]]=cl[9];
    ++fcount[2];
    }
	
	// 5
    if(pt[5]==1)
    {
    fd[2][fcount[2]]=a->nodeval(i,j-1,k)-1;
    ++fcount[2];
    }

    if(cl[5]<=0 && (pt[1]==1 || cl[1]<=0))
    {
    fd[2][fcount[2]]=cl[5];
    ++fcount[2];
    }
	

// Face 3

	// 0
    if(pt[0]==1)
    {
    fd[3][fcount[3]]=a->nodeval(i-1,j-1,k-1)-1;
    ++fcount[3];
    }

    if(cl[3]<=0)
    {
    fd[3][fcount[3]]=cl[3];
    ++fcount[3];
    }
	
	if(pt[0]==0 && cl[3]>0 && cl[4]<=0)
    {
    fd[3][fcount[3]]=cl[4];
    ++fcount[3];
    }
	
	// 3
    if(pt[3]==1)
    {
    fd[3][fcount[3]]=a->nodeval(i-1,j,k-1)-1;
    ++fcount[3];
    }

    if(cl[7]<=0)
    {
    fd[3][fcount[3]]=cl[7];
    ++fcount[3];
    }

	// 7
    if(pt[7]==1)
    {
    fd[3][fcount[3]]=a->nodeval(i-1,j,k)-1;
    ++fcount[3];
    }

    if(cl[11]<=0)
    {
    fd[3][fcount[3]]=cl[11];
    ++fcount[3];
    }
	
	// 4
    if(pt[4]==1)
    {
    fd[3][fcount[3]]=a->nodeval(i-1,j-1,k)-1;
    ++fcount[3];
    }

    if(cl[4]<=0 && (pt[0]==1 || cl[3]<=0)) 
    {
    fd[3][fcount[3]]=cl[4];
    ++fcount[3];
    }
	
// Face 4
	
	// 0
    if(pt[0]==1)
    {
    fd[4][fcount[4]]=a->nodeval(i-1,j-1,k-1)-1;
    ++fcount[4];
    }

    if(cl[0]<=0)
    {
    fd[4][fcount[4]]=cl[0];
    ++fcount[4];
    }
	
	if(pt[0]==0 && cl[0]>0 && cl[3]<=0)
    {
    fd[4][fcount[4]]=cl[3];
    ++fcount[4];
    }

    // 1
    if(pt[1]==1)
    {
    fd[4][fcount[4]]=a->nodeval(i,j-1,k-1)-1;
    ++fcount[4];
    }

    if(cl[1]<=0)
    {
    fd[4][fcount[4]]=cl[1];
    ++fcount[4];
    }
	
	// 2
    if(pt[2]==1)
    {
    fd[4][fcount[4]]=a->nodeval(i,j,k-1)-1;
    ++fcount[4];
    }
	
    if(cl[2]<=0)
    {
    fd[4][fcount[4]]=cl[2];
    ++fcount[4];
    }
	
	// 3
    if(pt[3]==1)
    {
    fd[4][fcount[4]]=a->nodeval(i-1,j,k-1)-1;
    ++fcount[4];
    }

    if(cl[3]<=0 && (pt[0]==1 || cl[0]<=0))
    {
    fd[4][fcount[4]]=cl[3];
    ++fcount[4];
    }


// Face 5

	// 4
    if(pt[4]==1)
    {
    fd[5][fcount[5]]=a->nodeval(i-1,j-1,k)-1;
    ++fcount[5];
    }

    if(cl[8]<=0)
    {
    fd[5][fcount[5]]=cl[8];
    ++fcount[5];
    }
	
	if(pt[4]==0 && cl[6]>0 && cl[11]<=0)
    {
    fd[5][fcount[5]]=cl[11];
    ++fcount[5];
    }
	
	// 5
    if(pt[5]==1)
    {
    fd[5][fcount[5]]=a->nodeval(i,j-1,k)-1;
    ++fcount[5];
    }

    if(cl[9]<=0)
    {
    fd[5][fcount[5]]=cl[9];
    ++fcount[5];
    }

	// 6
    if(pt[6]==1)
    {
    fd[5][fcount[5]]=a->nodeval(i,j,k)-1;
    ++fcount[5];
    }

    if(cl[10]<=0)
    {
    fd[5][fcount[5]]=cl[10];
    ++fcount[5];
    }

	// 7
    if(pt[7]==1)
    {
    fd[5][fcount[5]]=a->nodeval(i-1,j,k)-1;
    ++fcount[5];
    }

    if(cl[11]<=0 && (pt[4]==1 || cl[6]<=0))
    {
    fd[5][fcount[5]]=cl[11];
    ++fcount[5];
    }
}

