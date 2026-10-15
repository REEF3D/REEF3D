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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"sflow_f.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"river_ini.h"

// River free-surface initialisation (F 65, river_ini.h): eta from ini_fsf (vel=0), the velocities after
// the potential-flow initialisation (vel=1, F 65 2 only)
void sflow_f::ini_river(lexer *p, fdm2D *b, ghostcell *pgc, int vel)
{
    if(p->F65<1 || p->F65>2)
    return;

    if(vel==1 && p->F65!=2)
    return;

    river_ini river;
    if(river.start(p,pgc,b->bed,vel==0?1:0,0)==0)
    return;

    if(vel==0)
    {
        SLICELOOP4
        b->eta(i,j) = river.eta(p,i,j,b->bed(i,j));
        return;
    }

    SLICELOOP4
    river.velocity(p,i,j,b->bed(i,j),b->U(i,j),b->V(i,j));
}
