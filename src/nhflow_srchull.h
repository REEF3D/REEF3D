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
Author: Ahmet Soydan
--------------------------------------------------------------------*/

#ifndef NHFLOW_SRCHULL_H_
#define NHFLOW_SRCHULL_H_

#include<vector>

// EXPERIMENT (branch ahmet_dev_ship; switched on with environment variables, no input parameters yet).
//
// NHFLOW floating body (X 10 1) as a distribution of mass sources instead of forced body cells.
//
// The direct forcing holds the water in the cells of the hull at the body velocity. Next to them the water loses
// its tangential momentum (every variant X 15 0-11): the inviscid solver then has a boundary layer of numerical
// viscosity that separates at the stern, a wake with 10 times the momentum deficit of the real ship, and no stern
// and transverse waves (KCS, 2D ellipse; Cases/forcing_tests). Without body cells NHFLOW makes the waves of a
// surface pressure to 0.5 % of the linear theory.
//
// Here the hull is replaced by the sources that make the water flow around it (Hogner; Neumann-Michell, Noblesse
// et al. 2013): with the hull indicator chi (1 inside, the smoothed Heaviside of the body level set) and the
// stream Us relative to the body in x,
//
//      div u = q,     q = (1 + k1) Us d(chi)/dx      [1/s]
//
// (surface sources of density -(1 + k1) Us n_x, n the outward normal of the hull; grad(chi) = -n delta)
//
// sources where the hull opens (bow), sinks where it closes (stern), sum 0. The water moves through the hull
// region as through the stream; there is no wall, no forcing and no wake. For an ellipsoid the flow outside is
// exact with k1 its added mass coefficient in x (ellipse: b/a, circle: 1); for a slender ship k1 ~ 0.03.
// The sources enter
//      the pressure projection (right-hand side, nhflow_pjm::rhs, nhflow_pjm_corr::rhs),
//      the water level (nhflow_fsf_f::rk2_step1/2), the sigma velocity (nhflow_sigma::omega_update),
//      the momentum equations (the water enters with the local velocity; nhflow_HLL::aij_U/V/W).
// Resistance: the wave pattern behind the ship (wave cut), and Lagally's theorem, Fx = -rho sum q u dV
// (REEF3D_NHFLOW_6DOF/REEF3D_srchull.dat).
//
//      NHF_SRC    = Us [m/s]   stream in +x relative to the (fixed) body; set: sources instead of forcing
//      NHF_SRC_K  = k1         default 0
//      NHF_SRC_FB = 1          the strength from the flow: q <- (1 - chi) q + (u - u_body).grad(chi) every time step,
//                              no flow through the hull from outside (6DOF_obj_forcing_nhflow.cpp); NHF_SRC then
//                              only switches the sources on
//      NHF_SRC_WL = 1          water line term (line sources -(Us^2/g) sigma n_x t_y in the top layer), see
//                              6DOF_obj_forcing_nhflow.cpp; factor
//      NHF_SRC_STL = file      the sources from the hull facets and the water line (Cases/WIGLEY/srcpts.py) instead of
//                              the level set; the water line term then needs NHF_SRC_WL
//      NHF_SRC_TOP = 1         with NHF_SRC_STL: facets above the free surface of their column keep their source (top layer)
//      NHF_SRC_FORCE = 1       transparent interior: the cells of the hull forced to the stream (NHF_SRC, 0, 0); with the
//                              sources no flow through the hull for any shape (k1 = 0), see 6DOF_obj_forcing_nhflow.cpp
//      NHF_SRC_FORCE = 3       mirror hull, NO sources: the cells inside the hull are ghost cells of a free-slip wall (velocity of a
//                              probe point outside, tangential part unchanged, normal part reflected), see
//                              6DOF_obj_forcing_nhflow.cpp; NHF_SRC_PROBE (cells, 1.5), NHF_SRC_MDEPTH (cells, 2.6)
//      NHF_SRC_ANISO = 1       band half width from the grid spacing normal to the hull (grids stretched to the hull)
//      NHF_SRC_WATER = 1       the cells of the hull stay water for the solver (DF = 1 instead of -2): the usual wave speeds
//                              of the HLL flux at and inside the hull (with DF < 0: the left state inside, dry-bed speeds at
//                              the faces to the water) and no mirrored no-slip wall in the diffusion
//      NHF_SRC_VMATCH = 1      the stream inside follows the outer flow: V_x(x, layer) = mean x-velocity of the water in the shell
//                              1.5-3 psi outside the hull, sources q = d(chi V_x)/dx (no slip between the stream inside and the
//                              water beside the hull, see 6DOF_obj_forcing_nhflow.cpp); NHF_SRC_VREL = relaxation per step (0.02)
//      NHF_SRC_VMATCH = 2      the field inside is the tangential extension of the outer flow (probe NHF_SRC_PROBE cells outside
//                              the hull, normal part linear through zero at the hull), sources q = div(chi V): no band of
//                              sources with a jump of the normal velocity (see 6DOF_obj_forcing_nhflow.cpp)
//      NHF_SRC_RAMP = T        the sources grow from 0 over the first T seconds (soft start, see 6DOF_obj_forcing_nhflow.cpp)
//      NHF_SRC_QV = theta      momentum of the source water: theta V + (1 - theta) u per mass, V = ((1 + k1) Us, 0, 0) the
//                              stream inside the hull, u the velocity of the cell (default theta = 0). With the transparent
//                              interior the sources make the water of the stream inside and the sinks take it away again:
//                              two fluids in the band, the one inside with V. With u (theta = 0) the sinks at the stern take
//                              less momentum than the stream brings, the rest leaves as a jet (Wigley double body: u = 1.04 U
//                              behind the stern, thrust 20 N = -2.9e-3 of 1/2 rho U^2 S, Cases/WIGLEY/db_n075)
// Explicit RK2 (A 510 2), Cartesian grid.

namespace nhflow_srchull
{
    bool active();
    double stream();
    double k1();
    bool water();                       // NHF_SRC_WATER=1: the body cells are water for the solver (DF = 1)
    double qv();                        // NHF_SRC_QV = theta: the source water carries theta V + (1 - theta) u, V the stream inside
    int feedback();                     // NHF_SRC_FB: 1 strength from the local flow (no k1), 2 without the w term, +10 every stage
    int opt();                          // NHF_SRC_OPT, tests: 1 no momentum term, 2 no sigma velocity term, 4 no water level term
    
    extern std::vector<double> q;       // source of the cells [1/s], index IJK; empty until the first stage
    
    // momentum budget of the hull representation (REEF3D_srchull.dat): x and z momentum that the level reset of the pierced
    // columns (X 13 1, nhflow_fsf_body) gave to the water in its last stage [kg m/s] of this rank, and alpha dt of that stage
    extern double x13_mx, x13_mz, x13_adt, x13_dv;
    
    inline double at(int n) {return q.empty() ? 0.0 : q[n];}
}

#endif
