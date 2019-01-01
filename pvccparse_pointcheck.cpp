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

void pvccparse::pointcheck(lexer* p, fdm* a, ghostcell *pgc)
{

// Point check (non assigned pt's are evaluated)
	
    for(int qn=0; qn<8; ++qn)
    {
    if(pt[0]==-1 && ((pt[1]==1 && pt[3]==1) || (pt[1]==1 && pt[4]==1) || (pt[3]==1 && pt[4]==1)))
    pt[0]=1;

    if(pt[1]==-1 && ((pt[0]==1 && pt[2]==1) || (pt[0]==1 && pt[5]==1) || (pt[2]==1 && pt[5]==1)))
    pt[1]=1;

    if(pt[2]==-1 && ((pt[1]==1 && pt[3]==1) || (pt[1]==1 && pt[6]==1) || (pt[3]==1 && pt[6]==1)))
    pt[2]=1;

    if(pt[3]==-1 && ((pt[0]==1 && pt[2]==1) ||(pt[0]==1 && pt[7]==1) || (pt[2]==1 && pt[7]==1)))
    pt[3]=1;

    if(pt[4]==-1 && ((pt[0]==1 && pt[5]==1) || (pt[0]==1 && pt[7]==1) || (pt[5]==1 && pt[7]==1)))
    pt[4]=1;

    if(pt[5]==-1 && ((pt[1]==1 && pt[4]==1) || (pt[1]==1 && pt[6]==1) || (pt[4]==1 && pt[6]==1)))
    pt[5]=1;

    if(pt[6]==-1 && ((pt[2]==1 && pt[5]==1) || (pt[2]==1 && pt[7]==1) || (pt[5]==1 && pt[7]==1)))
    pt[6]=1;

    if(pt[7]==-1 && ((pt[3]==1 && pt[4]==1) || (pt[3]==1 && pt[6]==1) || (pt[4]==1 && pt[6]==1)))
    pt[7]=1;
    }
	
    for(int qn=0; qn<8; ++qn)
    if(pt[qn]==-1)
    pt[qn]=0;
}
