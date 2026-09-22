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

#include"mgcslice2.h"
#include"lexer.h"
#include"ghostcell.h"


mgcslice2::mgcslice2(lexer *p)
{
	imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
}

mgcslice2::~mgcslice2()
{
}

void mgcslice2::makemgc(lexer* p)
{
    
//flag2
    for(i=0;i<p->imax*p->jmax; ++i)
	{
	p->flagslice2[i]=p->flagslice4[i];
	}

    SLICELOOP4
    {

        if(p->flagslice4[IJp1]<0)
        p->flagslice2[IJ]=p->flagslice4[IJp1];
    }
}

