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

#include"solver_void.h"

solver_void::solver_void(lexer* p,fdm* a,ghostcell *pgc)
{
}

solver_void::~solver_void()
{
}

void solver_void::start(lexer* p,fdm* a, ghostcell* pgc, field& xfield, vec& rhsvec, int var)
{
}

void solver_void::startf(lexer* p, ghostcell* pgc, field &f, vec& rhs, matrix_diag &M, int var)
{
    
}

void solver_void::startF(lexer* p, ghostcell* pgc, double *f, vec& rhsvec, matrix_diag &M, int var)
{
}

void solver_void::startV(lexer* p, ghostcell* pgc, double *f, vec& rhs, matrix_diag &M, int var)
{
    
}

void solver_void::startM(lexer* p, ghostcell* pgc, double *x, double *rhs, double *M, int var)
{
}




