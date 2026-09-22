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

#include"matrix_diag.h"
#include"lexer.h"

matrix_diag::matrix_diag(lexer *pp)
{
	pp->Darray(n,pp->veclength);
	pp->Darray(s,pp->veclength);
	pp->Darray(e,pp->veclength);
	pp->Darray(w,pp->veclength);
	pp->Darray(t,pp->veclength);
	pp->Darray(b,pp->veclength);
	pp->Darray(p,pp->veclength);
    
}

matrix_diag::~matrix_diag()
{
    delete [] n;
	delete [] s;
	delete [] w;
	delete [] e;
	delete [] t;
	delete [] b;
	delete [] p;
}

void matrix_diag::resize(lexer *pp, int size_old, int size_new)
{
    pp->Dresize(n,size_old,size_new);
    pp->Dresize(s,size_old,size_new);
    pp->Dresize(e,size_old,size_new);
    pp->Dresize(w,size_old,size_new);
    pp->Dresize(t,size_old,size_new);
    pp->Dresize(b,size_old,size_new);
    pp->Dresize(p,size_old,size_new);
}
