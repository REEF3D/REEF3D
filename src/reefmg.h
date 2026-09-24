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

#ifndef REEFMG_H_
#define REEFMG_H_

#include "solver.h"
#include "increment.h"
#include "reefmg_core.h"

class lexer;
class fdm;
class ghostcell;
class field;
class vec;
class matrix_diag;

using namespace std;

//  x-y semicoarsening multigrid with vertical line relaxation for the
//  REEF3D::FNPF Laplace equation, the NHFLOW pressure and potential
//  problems, and the REEF3D::CFD pressure Poisson equation.
//
//  The sigma direction is never coarsened and is solved exactly by a
//  tridiagonal sweep, so vertical stretching does not enter the convergence
//  rate and the hierarchy does not run out of levels when there are only a
//  few layers.  Coarse operators are built from the fine coefficients, so
//  variable depth and grid stretching are inherited rather than rederived.

class reefmg final : public solver, public increment, public sc_operator
{
public:

    reefmg(lexer*, ghostcell*, int, int);
    virtual ~reefmg();

    void start(lexer*, fdm*, ghostcell*, field&, vec&, int) override final;
    void startf(lexer*, ghostcell*, field&, vec&, matrix_diag&, int) override final;
    void startF(lexer*, ghostcell*, double*, vec&, matrix_diag&, int) override final;
    void startV(lexer*, ghostcell*, double*, vec&, matrix_diag&, int) override final;
    void startM(lexer*, ghostcell*, double*, double*, double*, int) override final;

    //  exact fine operator for the Krylov iteration in fp32 mode (sc_operator)
    void fine_apply(const sc_level &L,const double *x,double *y) override;

private:

    void topology(lexer*, ghostcell*);
    void start_solver8(lexer*, ghostcell*, double*, vec&, matrix_diag&, int);
    void fill_matrix8(lexer*, double*, vec&, matrix_diag&);
    void fillbackvec8(lexer*, double*);

    //  CFD pressure Poisson equation, start(...,5): pjm, pjm_corr
    void start_solver5(lexer*, ghostcell*, field&, vec&, matrix_diag&);
    void fill_matrix5(lexer*, field&, vec&, matrix_diag&);
    void fillbackvec5(lexer*, field&);

    //  first row of every fully active column with consecutive rows (fp32)
    void build_colrow0(const sc_level&);

    //  Potential-flow initialisation, var 44: NHFLOW through startV (double*,
    //  ff==0) and CFD through start (field, ff!=0).
    void start_solver44(lexer*, fdm*, ghostcell*, field*, double*, vec&, matrix_diag&);
    void fill_matrix44(lexer*, reefmg_core&, field*, double*, vec&, matrix_diag&);
    void fillbackvec44(lexer*, reefmg_core&, field*, double*);

    reefmg_core mg;

    int *CVAL4;

    matrix_diag *Mcur;          // REEF3D's matrix for the current solve
    std::vector<int> rowmap;    // reefmg fine cell -> matrix_diag row, -1 inactive
    std::vector<int> colrow0;   // first row of a fully active column with
                                // consecutive rows, -1 otherwise
    bool fill_ini=false;        // CVAL4 and the zero halo set up (first fill)
    std::vector<int> wetsig;    // p->wet at the last rowmap/colrow0 build
    int count;

    int solve_mode, presweep, postsweep;
    int npx,npy,cx,cy;

    double starttime;
};

#endif
