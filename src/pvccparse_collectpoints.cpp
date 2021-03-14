/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

void pvccparse::collectpoints(lexer* p, fdm* a, ghostcell *pgc)
{

    pcount=0;
    clcount=0;
	
	pcount0=0;
	pcount1=0;

    // 0
    if(pt[0]==1)
    {
    nd[pcount]=a->nodeval(i-1,j-1,k-1)-1;
    ++pcount;
	++pcount0;
    }

    if(cl[4]<=0 && pt[0]==0)
    {
    nd[pcount]=cl[4];
    ++pcount;
	++pcount0;
    ++clcount;
    }

    if(cl[0]<=0)
    {
    nd[pcount]=cl[0];
    ++pcount;
	++pcount0;
    ++clcount;
    }

    // 1
    if(pt[1]==1)
    {
    nd[pcount]=a->nodeval(i,j-1,k-1)-1;
    ++pcount;
	++pcount0;
    }

    if(cl[5]<=0 && pt[1]==0)
    {
    nd[pcount]=cl[5];
    ++pcount;
	++pcount0;
    ++clcount;
    }

    if(cl[1]<=0)
    {
    nd[pcount]=cl[1];
    ++pcount;
	++pcount0;
    ++clcount;
    }

    // 2
    if(pt[2]==1)
    {
    nd[pcount]=a->nodeval(i,j,k-1)-1;
    ++pcount;
	++pcount0;
    }

    if(cl[6]<=0 && pt[2]==0)
    {
    nd[pcount]=cl[6];
    ++pcount;
	++pcount0;
    ++clcount;
    }

    if(cl[2]<=0)
    {
    nd[pcount]=cl[2];
    ++pcount;
	++pcount0;
    ++clcount;
    }

    // 3
    if(pt[3]==1)
    {
    nd[pcount]=a->nodeval(i-1,j,k-1)-1;
    ++pcount;
	++pcount0;
    }

    if(cl[7]<=0 && pt[3]==0)
    {
    nd[pcount]=cl[7];
    ++pcount;
	++pcount0;
    ++clcount;
    }

    if(cl[3]<=0)
    {
    nd[pcount]=cl[3];
    ++pcount;
	++pcount0;
    ++clcount;
    }

    //-------
    //-------

    // 4
    if(pt[4]==1)
    {
    nd[pcount]=a->nodeval(i-1,j-1,k)-1;
    ++pcount;
	++pcount1;
    }

    if(cl[4]<=0 && pt[4]==0)
    {
    nd[pcount]=cl[4];
    ++pcount;
	++pcount1;
    ++clcount;
    }


    if(cl[8]<=0)
    {
    nd[pcount]=cl[8];
    ++pcount;
	++pcount1;
    ++clcount;
    }


    // 5
    if(pt[5]==1)
    {
    nd[pcount]=a->nodeval(i,j-1,k)-1;
    ++pcount;
	++pcount1;
    }

    if(cl[5]<=0 && pt[5]==0)
    {
    nd[pcount]=cl[5];
    ++pcount;
	++pcount1;
    ++clcount;
    }


    if(cl[9]<=0)
    {
    nd[pcount]=cl[9];
    ++pcount;
	++pcount1;
    ++clcount;
    }

    // 6
    if(pt[6]==1)
    {
    nd[pcount]=a->nodeval(i,j,k)-1;
    ++pcount;
	++pcount1;
    }

    if(cl[6]<=0 && pt[6]==0)
    {
    nd[pcount]=cl[6];
    ++pcount;
	++pcount1;
    ++clcount;
    }


    if(cl[10]<=0)
    {
    nd[pcount]=cl[10];
    ++pcount;
	++pcount1;
    ++clcount;
    }

    // 7
    if(pt[7]==1)
    {
    nd[pcount]=a->nodeval(i-1,j,k)-1;
    ++pcount;
	++pcount1;
    }

    if(cl[7]<=0 && pt[7]==0)
    {
    nd[pcount]=cl[7];
    ++pcount;
	++pcount1;
    ++clcount;
    }


    if(cl[11]<=0)
    {
    nd[pcount]=cl[11];
    ++pcount;
	++pcount1;
    ++clcount;
    }
}
