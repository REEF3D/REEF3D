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

#include"particle_pls.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void particle_pls::picardmove(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->F46>0)
    {
        for(n=0;n<posactive;++n)
        if(posflag[n]>0)
        pos[n][3]=phipol(p,a,pos[n][0],pos[n][1],pos[n][2]);

        for(n=0;n<negactive;++n)
        if(negflag[n]>0)
        neg[n][3]=phipol(p,a,neg[n][0],neg[n][1],neg[n][2]);
    }
}
