/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"mgcslice1.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void mgcslice1::gcb_seed(lexer *p)
{	
    // count gcbsl
	count=0;
	SLICELOOP1
    {   
        if(p->flagslice1[Im1J]<0)
        ++count;
	
        if(p->flagslice1[IJp1]<0)
        ++count;

        if(p->flagslice1[IJm1]<0)
        ++count;

        if(p->flagslice1[Ip1J]<0)
        ++count;
    }
    
	p->Iresize(p->gcbsl1,p->gcbsl1_count, count, 6, 6); 
    
    // solid depth behind each boundary cell, column [2]: the number of cells (at most 3) in
    // the direction of the side cs until the next fluid cell. The ghost cell routines write at
    // most that deep, so a solid one or two cells thick does not overwrite the fluid behind it.
    // The three ghost layers outside the domain have flagslice1 -10: depth 3 there, as before.
    auto depth = [p](int ii, int jj, int di, int dj)
    {
        int d=0;
        for(int q=1;q<=3;++q)
        {
        if(p->flagslice1[(ii+q*di-p->imin)*p->jmax + (jj+q*dj-p->jmin)]>0)
        break;
        ++d;
        }
        return d;
    };
    
    // find gcbsl
	count=0;
	SLICELOOP1
    {
        if(p->flagslice1[Im1J]<0)
        {
        p->gcbsl1[count][0]=i;
        p->gcbsl1[count][1]=j;
        p->gcbsl1[count][3]=1;
        p->gcbsl1[count][2]=depth(i,j,-1,0);
        p->gcbsl1[count][4]=21;
        ++count;
        }

        if(p->flagslice1[IJp1]<0)
        {
        p->gcbsl1[count][0]=i;
        p->gcbsl1[count][1]=j;
        p->gcbsl1[count][3]=2;
        p->gcbsl1[count][2]=depth(i,j,0,1);
        p->gcbsl1[count][4]=21;
        ++count;
        }

        if(p->flagslice1[IJm1]<0)
        {
        p->gcbsl1[count][0]=i;
        p->gcbsl1[count][1]=j;
        p->gcbsl1[count][3]=3;
        p->gcbsl1[count][2]=depth(i,j,0,-1);
        p->gcbsl1[count][4]=21;
        ++count;
        }

        if(p->flagslice1[Ip1J]<0)
        {
        p->gcbsl1[count][0]=i;
        p->gcbsl1[count][1]=j;
        p->gcbsl1[count][3]=4;
        p->gcbsl1[count][2]=depth(i,j,1,0);
        p->gcbsl1[count][4]=21;
        ++count;
        }

    } 
    p->gcbsl1_count=count;
}
