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

#ifndef SFLOW_CURV_H_
#define SFLOW_CURV_H_

#include"increment.h"

class lexer;
class fdm2D;
class slice;
class ghostcell;

using namespace std;

// SFLOW on a curvilinear horizontal grid (G 50, p->geo_curv==1).
//
// The state stays cell-centred with Cartesian velocity components (U,V) and conserved variables
// (WL, UH, VH, WH). The faces of the logical grid get the normal vectors of the geometry layer:
// face i+1/2 of cell (i,j) the unit normal n1 = (SX1,SY1)/LEN1 (towards +i), face j+1/2 the unit
// normal n2 = (SX2,SY2)/LEN2 (towards +j). The reconstruction stays in index space, the Riemann
// problem of a face is solved in its normal direction:
//   - the face velocities Us/Un (i-faces) and Ve/Vw (j-faces) hold the normal velocity u.n
//   - momentum flux of UH: UH (u.n) + P n_x, of VH: VH (u.n) + P n_y, P = 0.5 g eta^2 + g eta d_face
//   - mass flux: UH n_x + VH n_y
//   - the HLL formula is unchanged (fluxes per unit length, Cartesian conserved variables)
//   - divergence: sum of the face fluxes times the face lengths over the cell area
//   - bed-slope source: g eta sum(d_face S_face)/A, so a lake at rest stays at rest for any grid
//     (the face vectors of every cell close)
// On a Cartesian grid (n1=(1,0), n2=(0,1)) this is the Cartesian scheme; that path is not used
// there, the Cartesian kernels stay untouched.

class sflow_curv : public increment
{
public:
    // unsupported options are refused (P1a: hydrostatic, single level)
    static void check(lexer*);

    // Us/Un, Ve/Vw: normal velocity from (U,V) at both sides of the face (after the reconstruction
    // of U and V on the i-faces (Us/Un, Vs/Vn) and j-faces (Ue/Uw, Ve/Vw))
    static void normal_velocities(lexer*, fdm2D*);

    // face fluxes of the states (sflow_flux_build_f): ipol 1 UH, 2 VH, 3 WH, 4 mass
    static void flux(lexer*, fdm2D*, int ipol);

    // boundary faces: walls (mirrored state), inflow and outflow (physical flux of the ghost state)
    static void flux_bc(lexer*, fdm2D*, int ipol, slice &Fx, slice &Fy, int inflow, int outflow);

    // f -= div(F) over the cell area
    static void divergence(lexer*, fdm2D*, slice &f, slice &Fx, slice &Fy);

    // bed-slope source of the pre-balanced fluxes: comp 1 x (to b->F), comp 2 y (to b->G)
    static void pgrad(lexer*, fdm2D*, slice &eta, slice &f, int comp);

    // 1/sigma of the CFL condition of a cell: sum over the two directions of (|u.n|+c) L / A
    static double sigma(lexer*, fdm2D*, int i, int j, double c);

    // free-slip wall: ghost velocity = interior velocity with the normal component reversed
    static void wall_mirror(lexer*, slice &u, slice &v, int i, int j, int cs);

    // unit normal and length of face i+1/2 (dir 1) or j+1/2 (dir 2) of cell (i,j)
    static void normal(lexer*, int dir, int i, int j, double &nx, double &ny, double &len);
};

#endif
