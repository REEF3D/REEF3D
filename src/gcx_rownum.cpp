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

#include"ghostcell.h"
#include"lexer.h"
#include"fieldint.h"

void ghostcell::rownum4_update(lexer* p,fieldint &rownum4)
{
	p->N4_row=0;
	p->N4_col=0;

    FLUIDLOOP
	{
    rownum4(i,j,k)=p->N4_row;
    ++p->N4_row;
	++p->N4_col;
	}

    rangex(p,p->range_row4,p->N4_row);

	FLUIDLOOP
    rownum4(i,j,k)+=p->range_row4[p->mpirank];
}

void ghostcell::rownum7_update(lexer* p, int *rownum7)
{
	p->N7_row=0;
	p->N7_col=0;

    FLOOP
	{
    rownum7[FIJK]=p->N7_row;
    ++p->N7_row;
	++p->N7_col;
	}

    rangex(p,p->range_row7,p->N7_row);

	FLOOP
    rownum7[FIJK]+=p->range_row7[p->mpirank];
}
