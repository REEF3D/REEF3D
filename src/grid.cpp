/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"grid.h"
#include"lexer.h"
#include"ghostcell.h"

grid::grid(lexer *p)
{
    imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
    kmin=p->kmin;
    kmax=p->kmax;
}

grid::~grid()
{
}

void grid::make_dgc(lexer* p)
{
    p->dgc1_count=1;
	p->dgc2_count=1;
	p->dgc3_count=1;
	p->dgc4_count=1;
	
	p->Iarray(p->dgc1,p->dgc1_count,8);
	p->Iarray(p->dgc2,p->dgc2_count,8);
	p->Iarray(p->dgc3,p->dgc3_count,8);
	p->Iarray(p->dgc4,p->dgc4_count,8);
    
    
    p->Iarray(hgc,imax*jmax*kmax);
    
    for(i=0;i<imax*jmax*kmax;++i)
    hgc[i]=0;
}

void grid::unmake_dgc(lexer* p)
{
    p->del_Iarray(hgc,imax*jmax*kmax);
}
