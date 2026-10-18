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

#ifndef SI_ETA2D_H_
#define SI_ETA2D_H_

#include"increment.h"
#include"matrix2D.h"
#include"vec2D.h"
#include"slice4.h"

class lexer;
class ghostcell;
class solver2D;
class slice;

using namespace std;

// Semi-implicit free surface (theta method): the 2D eta system, shared by the
// semi-implicit SFLOW (A 226 1) and the hydrostatic semi-implicit NHFLOW (A 534 1).
//
//   eta^{n+1} - div_h( a_f grad_f eta^{n+1} ) = rhs
//
//   a_f   = g theta^2 dt^2 Ht_f / dxp_f     (face coefficient, per face, 0: closed face)
//   Ht_f  : conveyance of the face:  SFLOW  H_f/(1 + dt gamma_f)  (implicit bed friction)
//                                   NHFLOW dz^T A_f^-1 dz          (column operator)
//   grad_f: compact face gradient (eta_{i+1} - eta_i)/dxp_f, div_h over the cell widths
//           (curvilinear grid G 50: dxp_f the normal distance of the centres, div_h with the face
//           lengths over the cell area)
//
// The caller supplies the face coefficients (ax on the x-faces, slice1; ay on the
// y-faces, slice2) and the right-hand side (slice4).  Cells without an open face
// become decoupled rows eta = rhs.  Rows are numbered in SLICELOOP4 order, stencil
// n i+1, s i-1, w j+1, e j-1, as in every solver2D assembly.  The matrix is symmetric
// positive definite after scaling with the cell area; the 5-point operator with the
// reaction term 1 is the screened form reefmg2D is built for.

class si_eta2D : public increment
{
public:
    si_eta2D(lexer*);
    virtual ~si_eta2D();

    // eta: initial guess on entry, solution on exit (interior cells, no ghost update)
    void solve(lexer*, ghostcell*, solver2D*, slice &ax, slice &ay, slice &rhs, slice &eta);

    // general form:  r eta - div( a_f grad_f eta ) = rhs
    //  r      : reaction term per cell (1: wet; 0: dry in the wetting-drying iteration; plus the
    //           implicit part of a radiating open boundary)
    //  a face coefficient towards a cell outside the domain (flagslice4 < 0) is a Dirichlet
    //  condition with the value of eta in that ghost cell (open water level boundary)
    void solve(lexer*, ghostcell*, solver2D*, slice &ax, slice &ay, slice &rhs, slice &eta, slice &r);

    // wetting and drying (Casulli 2009): the water depth is WL(eta) = max(eta + d, 0), the volume
    // equation  WL(eta) + bdiag eta - div( a_f grad_f eta ) = rhsW  is solved by a Newton
    // iteration on the wet pattern pw (1 wet, 0 dry), starting from all wet (the linear system):
    //     (pw + bdiag) eta - div( a_f grad_f eta ) = rhsW - pw d
    // a cell without open faces and boundaries stays wet (its row is WL = rhsW >= 0)
    //  eta    : initial guess on entry, solution on exit (interior cells and MPI halo)
    //  itmax  : 1 gives the linear system of the fully wet domain
    //  returns the number of solves; dry: number of dry cells (all ranks)
    int solve_wetdry(lexer*, ghostcell*, solver2D*, slice &ax, slice &ay, slice &rhsW, slice &depth, slice &bdiag, slice &eta, int itmax, int &dry);

    // the parts of solve / solve_wetdry for a composite solve over several grids (multi-block NHFLOW,
    // nhflow_amr_si.cpp): the rows M, rhsvec of the system (build), the wet pattern pw (wd_start, the
    // reaction and right-hand side of the pattern wd_rows, the new pattern wd_update: local changes)
    void build(lexer*, slice &ax, slice &ay, slice &rhs, slice &eta, slice &r);
    void wd_start(lexer*);
    void wd_rows(lexer*, slice &rhsW, slice &depth, slice &bdiag);
    int wd_update(lexer*, slice &ax, slice &ay, slice &bdiag, slice &eta, slice &depth);
    const matrix2D& mat() const { return M; }
    const vec2D& rhsv() const { return rhsvec; }
    slice& wd_react() { return react; }
    slice& wd_rhs() { return rhs; }
    double wd_pattern(int ii, int jj) { return pw(ii,jj); }

    int iter;
    double time;

private:
    void assemble(lexer*, slice&, slice&, slice&, slice&, double);
    void run(lexer*, ghostcell*, solver2D*, slice&);

    matrix2D M;
    vec2D xvec, rhsvec;
    
    // wet pattern, reaction term and right-hand side of the wetting and drying iteration
    slice4 pw, react, rhs;
};

#endif
