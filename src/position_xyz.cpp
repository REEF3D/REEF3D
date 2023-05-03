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

#include"position.h"
#include"lexer.h"

double position::pos_x()
{
    pos = p->XP[IP];

    return pos;
}

double position::pos_y()
{
    pos = p->YP[JP];

    return pos;
}

double position::pos_z()
{
    if(p->G2==0)
    pos = p->ZP[KP];
    
    if(p->G2==1)
    pos = p->ZSP[IJK];

    return pos;
}

double position::pos1_x()
{
    pos = p->XN[IP1];

    return pos;
}

double position::pos1_y()
{
    pos = p->YP[JP];

    return pos;
}

double position::pos1_z()
{
    if(p->G2==0)
    pos = p->ZP[KP];
    
    if(p->G2==1)
    pos = p->ZSP[IJK];

    return pos;
}

double position::pos2_x()
{
    pos = p->XP[IP];

    return pos;
}

double position::pos2_y()
{
    pos = p->YN[JP1];

    return pos;
}

double position::pos2_z()
{
    if(p->G2==0)
    pos = p->ZP[KP];
    
    if(p->G2==1)
    pos = p->ZSP[IJK];
    
    return pos;
}

double position::pos3_x()
{
    pos = p->XP[IP];

    return pos;
}

double position::pos3_y()
{
    pos = p->YP[JP];

    return pos;
}

double position::pos3_z()
{
    if(p->G2==0)
    pos = p->ZN[KP1];
    
    if(p->G2==1)
    pos = p->ZSN[FIJKp1];
    
    return pos;
}


double position::posnode_x()
{
    pos = p->XN[IP1];

    return pos;
}

double position::posnode_y()
{
    pos = p->YN[JP1];

    return pos;
}

double position::posnode_z()
{
    if(p->G2==0)
    pos = p->ZN[KP1];
    
    if(p->G2==1)
    pos = p->ZSN[FIJKp1];

    return pos;
}
