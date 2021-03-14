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

#include"diff_void.h"

diff_void::diff_void()
{
}

diff_void::~diff_void()
{
}


void diff_void::diff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &b, field &visc, double sig, double alpha)
{
}

void diff_void::idiff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &b, field &visc, double sig, double alpha)
{
}

void diff_void::diff_u(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
}

void diff_void::diff_v(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
}

void diff_void::diff_w(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
}

