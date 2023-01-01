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

#include"weno3_nug_func.h"
#include<math.h>

// IS ----
void weno3_nug_func::is_min_x()
{
    is1x = isfx[IP][uf][0]*pow(q3-q2,2.0);
    is2x = isfx[IP][uf][1]*pow(q2-q1,2.0);
}

void weno3_nug_func::is_min_y()
{
    is1y = isfy[JP][vf][0]*pow(q3-q2,2.0);
    is2y = isfy[JP][vf][1]*pow(q2-q1,2.0);
}

void weno3_nug_func::is_min_z()
{
    is1z = isfz[KP][wf][0]*pow(q3-q2,2.0);
    is2z = isfz[KP][wf][1]*pow(q2-q1,2.0);
}

void weno3_nug_func::is_max_x()
{
    is1x = isfx[IP][uf][2]*pow(q3-q2,2.0);
    is2x = isfx[IP][uf][3]*pow(q2-q1,2.0);
}

void weno3_nug_func::is_max_y()
{
    is1y = isfy[JP][vf][2]*pow(q3-q2,2.0);
    is2y = isfy[JP][vf][3]*pow(q2-q1,2.0);
}

void weno3_nug_func::is_max_z()
{
    is1z = isfz[KP][wf][2]*pow(q3-q2,2.0);
    is2z = isfz[KP][wf][3]*pow(q2-q1,2.0);
}
