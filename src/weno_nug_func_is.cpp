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

#include"weno_nug_func.h"
#include"lexer.h"
#include<math.h>

// IS ----
// x
void weno_nug_func::is_min_x()
{
    is1x = isfx[IP][uf][0][0]*pow(q5-q4,2.0) + isfx[IP][uf][0][1]*(q5-q4)*(q3-q4) + isfx[IP][uf][0][2]*pow(q3-q4,2.0);
    is2x = isfx[IP][uf][1][0]*pow(q2-q3,2.0) + isfx[IP][uf][1][1]*(q4-q3)*(q2-q3) + isfx[IP][uf][1][2]*pow(q4-q3,2.0);
    is3x = isfx[IP][uf][2][0]*pow(q1-q2,2.0) + isfx[IP][uf][2][1]*(q3-q2)*(q1-q2) + isfx[IP][uf][2][2]*pow(q3-q2,2.0);
}

void weno_nug_func::is_max_x()
{
    is1x = isfx[IP][uf][3][0]*pow(q5-q4,2.0) + isfx[IP][uf][3][1]*(q5-q4)*(q3-q4) + isfx[IP][uf][3][2]*pow(q3-q4,2.0);
    is2x = isfx[IP][uf][4][0]*pow(q2-q3,2.0) + isfx[IP][uf][4][1]*(q4-q3)*(q2-q3) + isfx[IP][uf][4][2]*pow(q4-q3,2.0);
    is3x = isfx[IP][uf][5][0]*pow(q1-q2,2.0) + isfx[IP][uf][5][1]*(q3-q2)*(q1-q2) + isfx[IP][uf][5][2]*pow(q3-q2,2.0);
}

// y
void weno_nug_func::is_min_y()
{
    is1y = isfy[JP][vf][0][0]*pow(q5-q4,2.0) + isfy[JP][vf][0][1]*(q5-q4)*(q3-q4) + isfy[JP][vf][0][2]*pow(q3-q4,2.0);
    is2y = isfy[JP][vf][1][0]*pow(q2-q3,2.0) + isfy[JP][vf][1][1]*(q4-q3)*(q2-q3) + isfy[JP][vf][1][2]*pow(q4-q3,2.0);
    is3y = isfy[JP][vf][2][0]*pow(q1-q2,2.0) + isfy[JP][vf][2][1]*(q3-q2)*(q1-q2) + isfy[JP][vf][2][2]*pow(q3-q2,2.0);
}

void weno_nug_func::is_max_y()
{
    is1y = isfy[JP][vf][3][0]*pow(q5-q4,2.0) + isfy[JP][vf][3][1]*(q5-q4)*(q3-q4) + isfy[JP][vf][3][2]*pow(q3-q4,2.0);
    is2y = isfy[JP][vf][4][0]*pow(q2-q3,2.0) + isfy[JP][vf][4][1]*(q4-q3)*(q2-q3) + isfy[JP][vf][4][2]*pow(q4-q3,2.0);
    is3y = isfy[JP][vf][5][0]*pow(q1-q2,2.0) + isfy[JP][vf][5][1]*(q3-q2)*(q1-q2) + isfy[JP][vf][5][2]*pow(q3-q2,2.0);
}

// z
void weno_nug_func::is_min_z()
{
    is1z = isfz[KP][wf][0][0]*pow(q5-q4,2.0) + isfz[KP][wf][0][1]*(q5-q4)*(q3-q4) + isfz[KP][wf][0][2]*pow(q3-q4,2.0);
    is2z = isfz[KP][wf][1][0]*pow(q2-q3,2.0) + isfz[KP][wf][1][1]*(q4-q3)*(q2-q3) + isfz[KP][wf][1][2]*pow(q4-q3,2.0);
    is3z = isfz[KP][wf][2][0]*pow(q1-q2,2.0) + isfz[KP][wf][2][1]*(q3-q2)*(q1-q2) + isfz[KP][wf][2][2]*pow(q3-q2,2.0);
}

void weno_nug_func::is_max_z()
{
    is1z = isfz[KP][wf][3][0]*pow(q5-q4,2.0) + isfz[KP][wf][3][1]*(q5-q4)*(q3-q4) + isfz[KP][wf][3][2]*pow(q3-q4,2.0);
    is2z = isfz[KP][wf][4][0]*pow(q2-q3,2.0) + isfz[KP][wf][4][1]*(q4-q3)*(q2-q3) + isfz[KP][wf][4][2]*pow(q4-q3,2.0);
    is3z = isfz[KP][wf][5][0]*pow(q1-q2,2.0) + isfz[KP][wf][5][1]*(q3-q2)*(q1-q2) + isfz[KP][wf][5][2]*pow(q3-q2,2.0);
}
