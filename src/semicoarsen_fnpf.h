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

#ifndef SEMICOARSEN_FNPF_H_
#define SEMICOARSEN_FNPF_H_

#include "solver.h"
#include "increment.h"
#include "semicoarsen_core.h"

class lexer;
class fdm;
class ghostcell;
class field;
class vec;
class matrix_diag;

using namespace std;

//  x-y semicoarsening multigrid with vertical line relaxation for the
//  REEF3D::FNPF Laplace equation.
//
//  The sigma direction is never coarsened and is solved exactly by a
//  tridiagonal sweep, so vertical stretching does not enter the convergence
//  rate and the hierarchy does not run out of levels when there are only a
//  few layers.  Coarse operators are built from the fine coefficients, so
//  variable depth and grid stretching are inherited rather than rederived.

class semicoarsen_fnpf final : public solver, public increment
{
public:

    semicoarsen_fnpf(lexer*, ghostcell*, int, int);
    virtual ~semicoarsen_fnpf();

    void start(lexer*, fdm*, ghostcell*, field&, vec&, int) override final;
    void startf(lexer*, ghostcell*, field&, vec&, matrix_diag&, int) override final;
    void startF(lexer*, ghostcell*, double*, vec&, matrix_diag&, int) override final;
    void startV(lexer*, ghostcell*, double*, vec&, matrix_diag&, int) override final;
    void startM(lexer*, ghostcell*, double*, double*, double*, int) override final;

private:

    void topology(lexer*, ghostcell*);
    void start_solver8(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    void fill_matrix8(lexer*, double*, vec&, matrix_diag&);
    void fillbackvec8(lexer*, double*);

    semicoarsen_core mg;

    int *CVAL4;
    int count;

    int solve_mode, presweep, postsweep;
    int npx,npy,cx,cy;

    double starttime;
};

#endif
