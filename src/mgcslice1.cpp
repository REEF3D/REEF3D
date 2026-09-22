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
#include"lexer.h"
#include"ghostcell.h"


mgcslice1::mgcslice1(lexer *p)
{
	imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
}

mgcslice1::~mgcslice1()
{
}

void mgcslice1::makemgc(lexer* p)
{
    
//flagslice1
    for(i=0;i<p->imax*p->jmax; ++i)
	{
	p->flagslice1[i]=p->flagslice4[i];
	}

    SLICELOOP4
    {

        if(p->flagslice4[Ip1J]<0)
        p->flagslice1[IJ]=p->flagslice4[Ip1J];
    }
}
