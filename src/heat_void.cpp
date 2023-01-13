/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"heat_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"

heat_void::heat_void(lexer* p, fdm* a, ghostcell *pgc)
{
}

heat_void::~heat_void()
{
}

void heat_void::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow)
{
}

void heat_void::ttimesave(lexer *p, fdm* a)
{
}

void heat_void::diff_update(lexer *p, fdm *a, ghostcell *pgc)
{
}

void heat_void::print_3D(lexer *p, fdm *a, ghostcell *pgc, ofstream& r)
{
}

double heat_void::val(int ii, int jj, int kk)
{
    double val=0.0;

    return val;
}

void heat_void::heat_ini(lexer* p, fdm *a, ghostcell* pgc, heat *pheat)
{
}

void heat_void::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
}

void heat_void::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void heat_void::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}



