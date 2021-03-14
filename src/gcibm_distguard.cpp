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

#include"ghostcell.h"
#include"lexer.h"

void ghostcell::distguard(lexer* p)
{
    /*
    for(n=0;n<p->gcb4_count;n++)
    p->gcd4[n]=MAX(p->gcd4[n],0.05*p->DXM);

    for(n=0;n<p->gcb4_count;n++)
    p->gcd4[n]=MIN(p->gcd4[n],0.95*p->DXM);
    */
}
