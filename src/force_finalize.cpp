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

#include"force.h"
#include"lexer.h"
#include"fdm.h"

void force::finalize(lexer *p, fdm* a)
{
	p->del_Iarray(tri,numtri_mem,4);
    p->del_Darray(pt,numvert_mem,3);
    p->del_Darray(ls,numvert_mem);
    p->del_Iarray(facet,numtri_mem,4);
    p->del_Iarray(confac,numtri_mem);
    p->del_Iarray(numfac,numtri_mem);
	p->del_Iarray(numpt,numtri_mem);
    p->del_Darray(ccpt,numtri_mem*4,3);
}
