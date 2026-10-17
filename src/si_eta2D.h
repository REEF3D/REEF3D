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

class lexer;
class ghostcell;
class solver2D;
class slice;

using namespace std;

// Semi-implicit free surface (theta method): the 2D eta system, shared by the
// semi-implicit SFLOW (A 226 1) and, later, the hydrostatic semi-implicit NHFLOW.
//
//   eta^{n+1} - div_h( a_f grad_f eta^{n+1} ) = rhs
//
//   a_f   = g theta^2 dt^2 Ht_f / dxp_f     (face coefficient, per face, 0: closed face)
//   Ht_f  : conveyance of the face:  SFLOW  H_f/(1 + dt gamma_f)  (implicit bed friction)
//                                   NHFLOW dz^T A_f^-1 dz          (column operator)
//   grad_f: compact face gradient (eta_{i+1} - eta_i)/dxp_f, div_h over the cell widths
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

    int iter;
    double time;

private:
    matrix2D M;
    vec2D xvec, rhsvec;
};

#endif
