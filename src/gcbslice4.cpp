/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"mgcslice4.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void mgcslice4::gcb_seed(lexer *p)
{	
    // count gcbsl
	count=0;
	SLICELOOP4
    {   
        if(p->flagslice4[Im1J]<0)
        ++count;
	
        if(p->flagslice4[IJp1]<0)
        ++count;

        if(p->flagslice4[IJm1]<0)
        ++count;

        if(p->flagslice4[Ip1J]<0)
        ++count;
    }
	
	p->Iresize(p->gcbsl4,p->gcbsl4_count, count, 6, 6); 
	p->Dresize(p->gcdsl4,p->gcbsl4_count, count); 
	
    // find gcbsl
	count=0;
	SLICELOOP4
    {
        if(p->flagslice4[Im1J]<0)
        {
        p->gcbsl4[count][0]=i;
        p->gcbsl4[count][1]=j;
        p->gcbsl4[count][3]=1;
        p->gcbsl4[count][4]=21;
        ++count;
        }

        if(p->flagslice4[IJp1]<0)
        {
        p->gcbsl4[count][0]=i;
        p->gcbsl4[count][1]=j;
        p->gcbsl4[count][3]=2;
        p->gcbsl4[count][4]=21;
        ++count;
        }

        if(p->flagslice4[IJm1]<0)
        {
        p->gcbsl4[count][0]=i;
        p->gcbsl4[count][1]=j;
        p->gcbsl4[count][3]=3;
        p->gcbsl4[count][4]=21;
        ++count;
        }

        if(p->flagslice4[Ip1J]<0)
        {
        p->gcbsl4[count][0]=i;
        p->gcbsl4[count][1]=j;
        p->gcbsl4[count][3]=4;
        p->gcbsl4[count][4]=21;
        ++count;
        }

    }
    p->gcbsl4_count=count;
	
}
