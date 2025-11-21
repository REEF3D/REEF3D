/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"geo_container.h"
#include"lexer.h"

geo_container::geo_container(lexer *p)
{
}

geo_container::~geo_container()
{
}

void geo_container::create_obj(lexer *p, int in_ID, int in_type, int in_trinum)
{
    ID = in_ID;
    type = in_type;
    trinum = in_trinum;
    
    p->Darray(tri,trinum,3);
}

void geo_container::delete_obj(lexer *p, int ID)
{
    p->del_Darray(tri,trinum,3);
}