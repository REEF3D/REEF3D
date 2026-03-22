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

#include"position.h"
#include"lexer.h"

double position::pos_x()
{
    return p->XP[IP];
}

double position::pos_y()
{
    return p->YP[JP];
}

double position::pos_z()
{
    if(p->G2==0)
    return p->ZP[KP];

    if(p->G2==1)
    return p->ZSP[IJK];
}

double position::pos1_x()
{
    return p->XN[IP1];
}

double position::pos1_y()
{
    return pos_y();
}

double position::pos1_z()
{
    return pos_z();
}

double position::pos2_x()
{
    return pos_x();
}

double position::pos2_y()
{
    return p->YN[JP1];
}

double position::pos2_z()
{
    return pos_z();
}

double position::pos3_x()
{
    return pos_x();
}

double position::pos3_y()
{
    return pos_y();
}

double position::pos3_z()
{
    if(p->G2==0)
    return p->ZN[KP1];

    if(p->G2==1)
    return p->ZSN[FIJKp1];
}

double position::posnode_x()
{
    return p->XN[IP1];
}

double position::posnode_y()
{
    return p->YN[JP1];
}

double position::posnode_z()
{
    if(p->G2==0)
    return p->ZN[KP1];

    if(p->G2==1)
    return p->ZSN[FIJKp1];
}
