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

#ifndef REEFMG2D_H_
#define REEFMG2D_H_

#include "solver2D.h"
#include "increment.h"
#include "reefmg_core.h"
#include <vector>

class lexer;
class ghostcell;
class slice;
class matrix2D;
class vec2D;

using namespace std;

//  REEFMG for the depth-averaged solvers (REEF3D::SFLOW): the non-hydrostatic
//  pressure (sflow_pjm_lin, sflow_pjm_quad) and the potential-flow
//  initialisation (sflow_potential_f), selected with N 10 1 like the 3D
//  solver.
//
//  A connector, not a second multigrid: the 5-point slice matrix is handed to
//  the same reefmg_core with a single layer (nz = 1).  The column solve then
//  reduces to point Gauss-Seidel, while the x-y coarsening, the geometry-
//  weighted coarse operators, agglomeration, fp32 storage and the fused
//  BiCGStab are shared with FNPF, NHFLOW and CFD.
//
//  Two things are specific to the depth-averaged problems:
//   - the pressure operator -div(h grad q) + c/h q is screened: error modes
//     longer than l ~ h/2 are damped by the reaction term alone, so levels
//     much coarser than l only cost time.  With N 13 0 the hierarchy is
//     truncated per solve from the matrix itself (set_active_levels).
//   - hydrostatic, breaking and dry cells arrive as decoupled rows (no
//     off-diagonals).  They are solved directly, kept out of the multigrid,
//     and their couplings from active neighbours become Dirichlet data, so the
//     coarse operators see them as boundaries rather than as fluid.

class reefmg2D final : public solver2D, public increment, public sc_operator
{
public:

    reefmg2D(lexer*, ghostcell*);
    virtual ~reefmg2D();

    void start(lexer*, ghostcell*, slice&, matrix2D&, vec2D&, vec2D&, int) override final;

    //  exact fine operator for the Krylov iteration in fp32 mode (sc_operator)
    void fine_apply(const sc_level &L,const double *x,double *y) override;

private:

    void topology(lexer*, ghostcell*);

    //  Copies matrix2D into the fine level of core.  Returns the global
    //  maximum over active rows of the screening length in cells, l/dx
    //  (a large value for unscreened rows), and in singular whether no row
    //  anywhere carries a reaction or Dirichlet term (pure Neumann).
    double fill(lexer*, ghostcell*, reefmg_core&, slice&, matrix2D&, vec2D&, bool &singular);
    void fillback(lexer*, reefmg_core&, slice&);

    void solve_neumann(lexer*, ghostcell*, slice&, matrix2D&, vec2D&, vec2D&, int);
    void solve_hypre(lexer*, ghostcell*, slice&, matrix2D&, vec2D&, vec2D&, int);

    int  screened_levels(double ldx, int nlev) const;

    reefmg_core mg;

    matrix2D *Mcur;               // matrix of the current solve (fine_apply)
    std::vector<int> cval;        // slice index IJ -> matrix2D row, -1 if none
    std::vector<int> rowmap;      // fine cell -> matrix2D row, -1 inactive
    std::vector<unsigned char> nbm; // fine cell: couplings kept (1 n, 2 s, 4 w, 8 e)
    std::vector<double> dval;     // fine cell: value of a decoupled row, 0 otherwise
    std::vector<char> dflag;      // fine cell: 1 decoupled row, solved directly
    std::vector<double> amask;    // fine cell incl. halo: 1 active, 0 not
    std::vector<double> dhalo;    // fine cell incl. halo: decoupled values

    int solve_mode, presweep, postsweep;
    int npx,npy,ndir;
    double starttime;
};

#endif
