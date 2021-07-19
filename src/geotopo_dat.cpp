/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"lexer.h"
#include"geotopo.h"
#include"fdm.h"
#include"ghostcell.h"

void geotopo::dat(lexer* p, fdm* a, ghostcell* pgc)
{
    if(p->G51>0)
    ALOOP
    a->topo(i,j,k)=-p->geodat[(i-p->imin)*p->jmax + (j-p->jmin)]+p->pos_z();
    
    pgc->start4a(p,a->topo,154);
    
}






