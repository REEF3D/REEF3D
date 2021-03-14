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

#include"position.h"
#include"lexer.h"

double position::pos_x()
{
    double pos = p->XP[IP];

    return pos;
}

double position::pos_y()
{
    double pos = p->YP[JP];

    return pos;
}

double position::pos_z()
{
    double pos = p->ZP[KP];

    return pos;
}

double position::pos1_x()
{
    double pos = p->XN[IP1];

    return pos;
}

double position::pos1_y()
{
    double pos = p->YP[JP];

    return pos;
}

double position::pos1_z()
{
    double pos = p->ZP[KP];

    return pos;
}

double position::pos2_x()
{
    double pos = p->XP[IP];

    return pos;
}

double position::pos2_y()
{
    double pos = p->YN[JP1];

    return pos;
}

double position::pos2_z()
{
    double pos = p->ZP[KP];

    return pos;
}

double position::pos3_x()
{
    double pos = p->XP[IP];

    return pos;
}

double position::pos3_y()
{
    double pos = p->YP[JP];

    return pos;
}

double position::pos3_z()
{
    double pos = p->ZN[KP1];

    return pos;
}


double position::posnode_x()
{
    double pos = p->XN[IP1];

    return pos;
}

double position::posnode_y()
{
    double pos = p->YN[JP1];

    return pos;
}

double position::posnode_z()
{
    double pos = p->ZN[KP1];

    return pos;
}
