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

#include"solver_void_ptf.h"

solver_void_ptf::solver_void_ptf(lexer* p,fdm_ptf *e,ghostcell *pgc)
{
}

solver_void_ptf::~solver_void_ptf()
{
}

void solver_void_ptf::start(lexer* p,fdm_ptf *e, ghostcell* pgc, field& xfield, vec& rhsvec, int var)
{
}

void solver_void_ptf::startf(lexer* p, ghostcell* pgc, field &f, vec& rhs, matrix_diag &M, int var)
{
    
}

void solver_void_ptf::startF(lexer* p, ghostcell* pgc, double *f, vec& rhsvec, matrix_diag &M, int var)
{
}

void solver_void_ptf::startM(lexer* p, ghostcell* pgc, double *x, double *rhs, double *M, int var)
{
}




