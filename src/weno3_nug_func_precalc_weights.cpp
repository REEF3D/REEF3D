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
#include"lexer.h"
#include"fdm.h"

// Weights ----
void weno3_nug_func::weight_min_x()
{
    w1x = cfx[IP][uf][0]/(epsilon + pow(is1x+psi,2.0)*(cfx[IP][uf][0]/pow(is1x+psi,2.0) + cfx[IP][uf][1]/pow(is2x+psi,2.0)));
    w2x = cfx[IP][uf][1]/(epsilon + pow(is2x+psi,2.0)*(cfx[IP][uf][0]/pow(is1x+psi,2.0) + cfx[IP][uf][1]/pow(is2x+psi,2.0)));
}

void weno3_nug_func::weight_max_x()
{
    w1x = cfx[IP][uf][2]/(epsilon + pow(is1x+psi,2.0)*(cfx[IP][uf][2]/pow(is1x+psi,2.0) + cfx[IP][uf][3]/pow(is2x+psi,2.0)));
    w2x = cfx[IP][uf][3]/(epsilon + pow(is2x+psi,2.0)*(cfx[IP][uf][2]/pow(is1x+psi,2.0) + cfx[IP][uf][3]/pow(is2x+psi,2.0)));
}

void weno3_nug_func::weight_min_y()
{
    w1y = cfy[JP][vf][0]/(epsilon + pow(is1y+psi,2.0)*(cfy[JP][vf][0]/pow(is1y+psi,2.0) + cfy[JP][vf][1]/pow(is2y+psi,2.0)));
    w2y = cfy[JP][vf][1]/(epsilon + pow(is2y+psi,2.0)*(cfy[JP][vf][0]/pow(is1y+psi,2.0) + cfy[JP][vf][1]/pow(is2y+psi,2.0)));
}

void weno3_nug_func::weight_max_y()
{
    w1y = cfy[JP][vf][2]/(epsilon + pow(is1y+psi,2.0)*(cfy[JP][vf][2]/pow(is1y+psi,2.0) + cfy[JP][vf][3]/pow(is2y+psi,2.0)));
    w2y = cfy[JP][vf][3]/(epsilon + pow(is2y+psi,2.0)*(cfy[JP][vf][2]/pow(is1y+psi,2.0) + cfy[JP][vf][3]/pow(is2y+psi,2.0)));
}

void weno3_nug_func::weight_min_z()
{
    w1z = cfz[KP][wf][0]/(epsilon + pow(is1z+psi,2.0)*(cfz[KP][wf][0]/pow(is1z+psi,2.0) + cfz[KP][wf][1]/pow(is2z+psi,2.0)));
    w2z = cfz[KP][wf][1]/(epsilon + pow(is2z+psi,2.0)*(cfz[KP][wf][0]/pow(is1z+psi,2.0) + cfz[KP][wf][1]/pow(is2z+psi,2.0)));
}

void weno3_nug_func::weight_max_z()
{
    w1z = cfz[KP][wf][2]/(epsilon + pow(is1z+psi,2.0)*(cfz[KP][wf][2]/pow(is1z+psi,2.0) + cfz[KP][wf][3]/pow(is2z+psi,2.0)));
    w2z = cfz[KP][wf][3]/(epsilon + pow(is2z+psi,2.0)*(cfz[KP][wf][2]/pow(is1z+psi,2.0) + cfz[KP][wf][3]/pow(is2z+psi,2.0)));
}

