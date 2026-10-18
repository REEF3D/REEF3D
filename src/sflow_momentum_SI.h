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
#include<vector>
#include<fstream>

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
// Prognostic: WL (eta), UH, VH (cell centres) and the face discharges qx, qy
// (slice1/slice2), which carry the mass flux of the continuity equation.
//
// One time step:
//  1. face depth H_f = eta_up + d_f (eta upwind of q^n, max(eta) for q = 0; d_f the mean
//     still water depth, the smaller one next to a dry cell), open faces H_f > A 244,
//     implicit bed friction alpha = 1/(1 + dt gamma), gamma = c_f |u| / H (Manning, A 218 1)
//  2. open boundaries from the ghost cells of ioflow / iowave (bc_setup):
//     discharge (flux given), water level (Dirichlet), radiating (Riemann, Flather, free
//     outflow: q_in = A - B eta^{n+1}, implicit)
//  3. predictor of the explicit terms, SSP-RK2/3 (A 210) with H^n and q^n frozen:
//     UH* = UH^n + dt L(UH),  L = -div(q^n u_upwind) + diffusion + sources
//     (no gravity term, no c in the fluxes: only the advective Courant number counts)
//  4. face predictor, momentum interpolation:  q^ = avg(UH*)  (the cell predictors carry no
//     gravity term of the step; an incremental form q^n + avg(UH* - UH^n) lets cells and faces
//     drift apart and oscillates in steady river flow for theta < 0.6)
//  5. q^E = alpha_f ( q^ - dt g H_f (1-theta) grad_f eta^n )
//  6. r(WL(eta^{n+1})) - div( a_f grad_f eta^{n+1} ) = WL^n - dt div( theta q^E + (1-theta) q^n )
//     a_f = theta dt (g theta dt alpha_f H_f + K_f)/dxp_f,  K_f: optional mass flux
//     dissipation (A 231) at faces with Froude number > 0.5
//     wetting and drying (A 229 1): WL(eta) = max(eta + d, 0), Newton iteration on the
//     wet/dry pattern (Casulli 2009): non-negative water depth, exact volume
//  7. q^{n+1} = q^E - (alpha_f dt g H_f theta + K_f) grad_f eta^{n+1}
//  8. WL^{n+1} from the conservative divergence of the face fluxes (exact volume)
//  9. UH^{n+1} = alpha_c ( UH* - dt g H_c (theta G_c^{n+1} + (1-theta) G_c^n) ),
//     G_c: mean of the open face gradients of the cell
// 10. relaxation zones once per step (eta, U, V), the face discharges follow the change
//
// The compact face gradient in the eta system and in q avoids the checkerboard mode of
// the collocated grid; the cells see the mean of the two face gradients.
// Lake at rest: grad eta = 0, exact for any bathymetry.
//
// Not yet (stops in check()): non-hydrostatic pressure and Boussinesq (A 220 > 0), mesh
// refinement and multi-block grids, curvilinear grids, patch boundary conditions (B 411-422),
// 6DOF / direct forcing, SEASTATE.
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
    // open boundary face of a domain edge
    enum bc_type {bc_flux=1, bc_level=2, bc_riemann=3, bc_flather=4};

    struct bc_face
    {
        int i,j;        // interior cell
        int cs;         // side: 1 x-, 2 y+, 3 y-, 4 x+
        int type;       // bc_type
        int fi,fj;      // face index (x-face for cs 1/4, y-face for cs 2/3)
        int gi,gj;      // first ghost cell
        double s;       // inward normal in the face direction (+1 / -1)
        double dxp;     // distance cell centre - ghost centre
        double dxn;     // cell width normal to the face
        double ds;      // face length
        double qin;     // flux type: inward discharge per width at t^{n+1}
        double A,B;     // radiating: q_in = A - B eta^{n+1}
    };

    void face_depth(lexer*, fdm2D*, ghostcell*);
    void friction(lexer*, fdm2D*);
    void bc_setup(lexer*, fdm2D*, ghostcell*);
    void bc_method(lexer*);
    void predictor(lexer*, fdm2D*, ghostcell*);
    void rk_stage(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&, double);
    void advection(lexer*, fdm2D*, ghostcell*);
    void eta_solve(lexer*, fdm2D*, ghostcell*);
    void face_update(lexer*, fdm2D*, ghostcell*);
    void eta_update(lexer*, fdm2D*, ghostcell*);
    void cell_update(lexer*, fdm2D*, ghostcell*);
    void relax(lexer*, fdm2D*, ghostcell*);
    void log(lexer*, fdm2D*, ghostcell*);

    void exchange1(lexer*, ghostcell*, slice&);
    void exchange2(lexer*, ghostcell*, slice&);
    void exchange4(lexer*, ghostcell*, slice&);
    double div(lexer*, slice&, slice&);
    double dissipation(lexer*, double, double, double);

    si_eta2D *psi;

    // face discharges (prognostic) and the face quantities of the step:
    // q^n, q^E, depth H_f (0: closed face), friction factor alpha_f, eta-system
    // coefficient a_f, dissipation coefficient K_f, face gradient of eta^n
    slice1 qx, qxn, qxE, Hx, alx, ax, kx, gx0;
    slice2 qy, qyn, qyE, Hy, aly, ay, ky, gy0;

    // advective fluxes of UH and VH
    slice1 fxU, fxV;
    slice2 fyU, fyV;

    // cell state of the step: t^n, RK stages, predictor UH*, increments, implicit eta,
    // friction factor alpha_c, right-hand side and reaction term of the eta system
    slice4 eta0, WL0, UH0, VH0, UHa, VHa, UHb, VHb, UHp, VHp, WHz, dUH, dVH, etaI, alc, rhs, react, rhsW;

    // wetting and drying: wet pattern of the Newton iteration, implicit diagonal of the open boundaries
    slice4 pw, bdiag;

    // interior copies kept through the ghost cell update (thin solids, see rk_stage)
    slice4 tU, tV, tUH, tVH;

    // water depth at the end of the step, the start of the next step (see start);
    // eta with the boundary values of ioflow / iowave (bc_setup)
    slice4 WLk, etaB;

    vector<bc_face> bcf;
    int method[5];          // iowave open edge method per side cs (B 520 3-6), 0 none

    double theta, g;
    int rk, wd_iter, wd_cells;
    bool qini;

    double Qin, Qout;       // discharge over the open boundaries of the step (theta weighted)
    double vrelax;          // volume change of the relaxation zones of the step
    ofstream silog;
};

#endif
