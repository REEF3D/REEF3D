/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"part.h"
#include"lexer.h"
#include"ghostcell.h"

void part::xchange_sendid(lexer *p, ghostcell *pgc, int mode)
{
    for(q=0;q<6;++q)
        sendnum[q]=0;

    // target sides from xchange_count (the positions have not changed since)
    for(n=0;n<index;++n)
    if(Flag[n]==ACTIVE && xdest[n]>0)
    {
        q = xdest[n]-1;
        sendid[q][sendnum[q]]=n;
        ++sendnum[q];
    }
}
