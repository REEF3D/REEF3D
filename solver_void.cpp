/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"solver_void.h"

solver_void::solver_void(lexer* p,fdm* a,ghostcell *pgc)
{
}

void solver_void::start(lexer* p,fdm* a, ghostcell* pgc, field& xfield, vec& xvec, vec& rhsvec, int var, int gcv, double stop_crit)
{
}

void solver_void::startF(lexer* p, ghostcell* pgc, double *f, vec& rhsvec, matrix_diag &M, int var, int gcv, double stop_crit)
{
}

solver_void::~solver_void()
{
}

void solver_void::setup(lexer* p,fdm* a, ghostcell* pgc, int var, cpt &C)
{
}

void solver_void::solve(lexer* p,fdm* a, ghostcell* pgc, vec& xvec, vec& rhsvec, int var, int gcv, int &solveriter, int maxiter, double stop_crit, cpt &C)
{	
}

void solver_void::fillxvec1(lexer* p, fdm* a, field& f)
{
}

void solver_void::fillxvec2( lexer* p, fdm* a, field& f)
{
}

void solver_void::fillxvec3( lexer* p, fdm* a, field& f)
{
}

void solver_void::fillxvec4( lexer* p, fdm* a, field& f)
{
}




