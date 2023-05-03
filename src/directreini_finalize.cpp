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

#include"directreini.h"
#include"lexer.h"
#include"fdm.h"

void directreini::finalize(lexer *p,fdm* a)
{
    del_Iarray(tri,numtri_mem,4);

    del_Darray(pt,numvert_mem,3);
    del_Iarray(ijk,numvert_mem,3);
    del_Darray(ls,numvert_mem);
	del_Darray(ls0,numvert_mem);
	del_Darray(ls1,numvert_mem);
	del_Darray(lsvert,numvert_mem);
	del_Darray(lsfac,numvert_mem);
	del_Iarray(reiniflag,numvert_mem);

    del_Iarray(facet,numtri_mem,4);
    del_Iarray(confac,numtri_mem);
    del_Iarray(numfac,numtri_mem);
    del_Darray(ccpt,numtri_mem*4,3);
}
