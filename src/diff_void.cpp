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

#include"diff_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"

diff_void::diff_void()
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
}

diff_void::~diff_void()
{
}


void diff_void::diff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &b, field &visc, field &eddyv, double sig, double alpha)
{
}

void diff_void::diff_scalar(lexer* p, fdm* a, ghostcell* pgc, solver* psolv, field &diff, field &b, field &visc, field &eddyv, double sig, double alpha)
{
    LOOP
	diff(i,j,k) = b(i,j,k);

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

void diff_void::diff_u(lexer* p, fdm* a, ghostcell* pgc, solver* psolv, field &udiff, field &u_in, field &u, field &v, field &w, double alpha)
{
	ULOOP
	udiff(i,j,k) = u_in(i,j,k);

    pgc->start1(p,udiff,gcval_u);
}

void diff_void::diff_v(lexer* p, fdm* a, ghostcell* pgc, solver* psolv, field &vdiff, field &v_in, field &u, field &v, field &w, double alpha)
{
	VLOOP
	vdiff(i,j,k) = v_in(i,j,k);

    pgc->start2(p,vdiff,gcval_v);
}

void diff_void::diff_w(lexer* p, fdm* a, ghostcell* pgc, solver* psolv, field &wdiff, field &w_in, field &u, field &v, field &w, double alpha)
{
	WLOOP
	wdiff(i,j,k) = w_in(i,j,k);

    pgc->start3(p,wdiff,gcval_w);
}
