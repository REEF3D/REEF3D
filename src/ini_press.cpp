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

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

void initialize::pressini(lexer *p, fdm *a, ghostcell *pgc)
{
    double depth=0.0;

    LOOP
    depth=MAX(depth, p->ZN[KP]);

    LOOP
    a->press(i,j,k)=a->ro(i,j,k)*fabs(p->W22)*(depth-p->ZN[KP]);


}
