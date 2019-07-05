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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm2D.h"

void ghostcell::cval2D_gcb1(lexer* p, sliceint &cval)
{
	GCSL1LOOP
    {
    i=p->gcbsl1[n][0];
    j=p->gcbsl1[n][1];
	
	p->gcbsl1[n][5]=cval(i,j);
	}
}

void ghostcell::cval2D_gcb2(lexer* p, sliceint &cval)
{
	GCSL2LOOP
    {
    i=p->gcbsl2[n][0];
    j=p->gcbsl2[n][1];
	
	p->gcbsl2[n][5]=cval(i,j);
	}
}

void ghostcell::cval2D_gcb4(lexer* p, sliceint &cval)
{
	GCSL4LOOP
    {
    i=p->gcbsl4[n][0];
    j=p->gcbsl4[n][1];
	
	p->gcbsl4[n][5]=cval(i,j);
	}
}

