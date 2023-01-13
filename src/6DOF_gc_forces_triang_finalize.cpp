/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void sixdof_gc::forces_triang_finalize(lexer* p, fdm *a, ghostcell *pgc)
{
    
    p->del_Iarray(tri,numtri,4);
    p->del_Darray(pt,numvert,3);
    p->del_Darray(ls,numvert);
    p->del_Iarray(facet,numtri,4);
    p->del_Iarray(confac,numtri);
    p->del_Iarray(numfac,numtri);
	p->del_Iarray(numpt,numtri);
    p->del_Darray(ccpt,numtri*4,3);
    
}
