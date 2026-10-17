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

#ifndef SFLOW_MOMENTUM_SI_H_
#define SFLOW_MOMENTUM_SI_H_

#include"sflow_momentum_func.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"

class si_eta2D;

using namespace std;

// ---------------------------------------------------------------------------
// SFLOW semi-implicit free surface (A 226 1), theta method (A 227), hydrostatic.
//
// For long simulations (river, tide): the gravity wave speed sqrt(gh) no longer limits
// the time step, only advection does (sflow_timestep_SI, cap A 228).
// Separate from the explicit HLL integrators (sflow_momentum_RK2/RK3), which are not
// touched; the shared stage functions of sflow_momentum_func (ghost cells, velocities,
// relaxation zones) are used as they are.
//
// Prognostic: eta, UH, VH (cell centres) and the face discharges qx, qy (slice1/slice2),
// which carry the mass flux of the continuity equation.
//
// One time step:
//  1. face depth H_f (upwind of q^n, max(eta) for q = 0), open faces, implicit bed friction
//     alpha = 1/(1 + dt gamma),  gamma = c_f |u| / H  (Manning, A 218 1)
//  2. predictor of the explicit terms, SSP-RK2/3 (A 210) with H^n and q^n frozen:
//     UH* = UH^n + dt L(UH),  L = -div(q^n u_upwind) + diffusion + sources
//     (no gravity term, no c in the fluxes: only the advective Courant number counts)
//  3. face predictor, incremental momentum interpolation:  q^ = q^n + avg(UH* - UH^n)
//  4. q^E = alpha_f ( q^ - dt g H_f (1-theta) grad_f eta^n )
//  5. eta^{n+1} - div( g theta^2 dt^2 alpha_f H_f grad_f eta^{n+1} )
//         = eta^n - dt div( theta q^E + (1-theta) q^n )                       (si_eta2D)
//  6. q^{n+1} = q^E - alpha_f dt g H_f theta grad_f eta^{n+1}
//  7. eta^{n+1} from the conservative divergence of the face fluxes (exact volume)
//  8. UH^{n+1} = alpha_c ( UH* - dt g H_c (theta G_c^{n+1} + (1-theta) G_c^n) ),
//     G_c: mean of the open face gradients of the cell
//
// The compact face gradient in the eta system and in q avoids the checkerboard mode of
// the collocated grid; the cells see the mean of the two face gradients.
// Lake at rest: grad eta = 0, exact for any bathymetry.
//
// Not yet (stops in check()): non-hydrostatic pressure and Boussinesq (A 220 > 0), mesh
// refinement and multi-block grids, curvilinear grids, open boundaries (wave generation,
// discharge, water level), 6DOF / direct forcing, SEASTATE.
// ---------------------------------------------------------------------------

class sflow_momentum_SI final : public sflow_momentum_func
{
public:
	sflow_momentum_SI(lexer*, fdm2D*, ghostcell*, sflow_HLL*, sflow_signal_speed*, sflow_reconstruct*, sflow_diffusion*, sflow_pressure*,
                        solver2D*, solver2D*, ioflow*, sflow_fsf*, sflow_forcing*, sixdof*);
	virtual ~sflow_momentum_SI();

	void start(lexer*, fdm2D*, ghostcell*) override final;

    static void check(lexer*);

private:
    void face_depth(lexer*, fdm2D*, ghostcell*);
    void friction(lexer*, fdm2D*);
    void predictor(lexer*, fdm2D*, ghostcell*);
    void rk_stage(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&, double);
    void advection(lexer*, fdm2D*, ghostcell*);
    void eta_solve(lexer*, fdm2D*, ghostcell*);
    void face_update(lexer*, fdm2D*, ghostcell*);
    void eta_update(lexer*, fdm2D*, ghostcell*);
    void cell_update(lexer*, fdm2D*, ghostcell*);

    void exchange1(lexer*, ghostcell*, slice&);
    void exchange2(lexer*, ghostcell*, slice&);
    void exchange4(lexer*, ghostcell*, slice&);
    double div(lexer*, slice&, slice&);

    si_eta2D *psi;

    // face discharges (prognostic) and the face quantities of the step:
    // q^n, q^E, depth H_f (0: closed face), friction factor alpha_f, eta-system
    // coefficient a_f, face gradient of eta^n
    slice1 qx, qxn, qxE, Hx, alx, ax, gx0;
    slice2 qy, qyn, qyE, Hy, aly, ay, gy0;

    // advective fluxes of UH and VH
    slice1 fxU, fxV;
    slice2 fyU, fyV;

    // cell state of the step: t^n, RK stages, predictor UH*, increments, implicit eta,
    // friction factor alpha_c, right-hand side of the eta system
    slice4 eta0, WL0, UH0, VH0, UHa, VHa, UHb, VHb, UHp, VHp, WHz, dUH, dVH, etaI, alc, rhs;

    // interior copies kept through the ghost cell update (thin solids, see rk_stage)
    slice4 tU, tV, tUH, tVH;

    // water depth at the end of the step, the start of the next step (see start)
    slice4 WLk;

    double theta, g;
    int rk;
    bool qini;
};

#endif
