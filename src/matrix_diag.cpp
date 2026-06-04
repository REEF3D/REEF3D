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

#include "matrix_diag.h"
#include "lexer.h"

#include <algorithm>

matrix_diag::matrix_diag(lexer *pp)
{
	n.resize(pp->veclength);
	s.resize(pp->veclength);
	e.resize(pp->veclength);
	w.resize(pp->veclength);
	t.resize(pp->veclength);
	b.resize(pp->veclength);
	p.resize(pp->veclength);
}

void matrix_diag::resize(lexer *pp, int size_new)
{
    n.resize(size_new);
    s.resize(size_new);
    e.resize(size_new);
    w.resize(size_new);
    t.resize(size_new);
    b.resize(size_new);
    p.resize(size_new);
}
