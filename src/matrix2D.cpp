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
#include"matrix2D.h"
#include"lexer.h"

matrix2D::matrix2D(lexer *pp)
{
	pp->Darray(n,pp->vec2Dlength);
	pp->Darray(s,pp->vec2Dlength);
	pp->Darray(e,pp->vec2Dlength);
	pp->Darray(w,pp->vec2Dlength);
	pp->Darray(p,pp->vec2Dlength);
}

matrix2D::~matrix2D()
{
    delete [] n;
	delete [] s;
	delete [] w;
	delete [] e;
	delete [] p;
}

void matrix2D::resize(lexer *pp, int size_old, int size_new)
{
    pp->Dresize(n,size_old,size_new);
    pp->Dresize(s,size_old,size_new);
    pp->Dresize(e,size_old,size_new);
    pp->Dresize(w,size_old,size_new);
    pp->Dresize(p,size_old,size_new);
}
