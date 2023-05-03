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

#include"vec2D.h"
#include"lexer.h"

vec2D::vec2D(lexer* p)
{
    p->Darray(V,p->vec2Dlength);
}

vec2D::~vec2D()
{
    delete [] V;
}

void vec2D::resize(lexer *pp, int size_old, int size_new)
{
    pp->Dresize(V,size_old,size_new);
}
