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

#include"initialise_ptf.h"
#include"fdm_ptf.h"
#include"lexer.h"
#include"ghostcell.h"

void initialise_ptf::pressini(lexer *p, fdm_ptf *e, ghostcell *pgc)
{
    double depth=0.0;

    LOOP
    depth=MAX(depth, p->ZN[KP]);

    LOOP
    e->press(i,j,k)=e->ro(i,j,k)*fabs(p->W22)*(depth-p->ZN[KP]);


}
