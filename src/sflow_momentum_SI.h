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
#include<functional>
#include<memory>

class si_eta2D;
class sflow_amr;
struct sflow_si_nh_lu;

using namespace std;

// ---------------------------------------------------------------------------
// SFLOW semi-implicit free surface (A 210 12/13), theta method (A 227), hydrostatic.
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
//     wetting and drying (A 243 > 0): WL(eta) = max(eta + d, 0), Newton iteration on the
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
// Grids: Cartesian and curvilinear (G 50 2, sflow_curv: face lengths and normals, the normal distance of the
// centres; the face discharges are normal to the faces, the gravity term of a cell from the mean normals of its open
// faces per index direction, as nhflow_momentum_SI); multi-block grids (G 50 1) and mesh refinement (G 1): every
// block and every patch has its own object, the stage runner sflow_amr::si_step runs the parts of the step on all
// grids, fills the cells around them in between (the coarse faces next to a patch take its fine fluxes) and solves
// one eta system over the leaf cells of all grids (reefamr_si.cpp).
//
// Non-hydrostatic quadratic pressure (A 220 2/3) (sflow_momentum_SI_nh.cpp): the pressure P^{n+1} fully implicit
// and eta^{n+1} (theta) from one coupled system of the continuity equation and the incompressibility constraint, so
// that neither the gravity nor the dispersive wave speed limits the time step; 1D and 2D, MPI, one Cartesian grid,
// without wetting and drying iteration (A 243 0).
//
// Not yet (stops in check()): Boussinesq (A 220 4); the non-hydrostatic pressure on curvilinear, multi-block and
// refined grids and with the wetting and drying iteration,
// patch boundary conditions (B 411-422), 6DOF / direct forcing, subcycling of the refinement (G 7 1), SEASTATE
// surfbeat.
// ---------------------------------------------------------------------------

class sflow_momentum_SI final : public sflow_momentum_func
{
public:

    void ckpt_restored() override { qini = true; }   // the face discharges of the checkpoint
    void water_column(lexer*, fdm2D*) override;      // REEF3D::HYDRO losses: WLk = WL
    bool boundary_flow(double &qin, double &qout) override { qin = Qin;  qout = -Qout;  return true; }
    void ckpt_fields(const std::function<void(const char*, char, double*)> &fd, const std::function<void(const char*, char, int*)> &fi) override
    {
        fd("momentum_SI_qx",'S',qx.V); fd("momentum_SI_qxn",'S',qxn.V); fd("momentum_SI_qxE",'S',qxE.V);
        fd("momentum_SI_Hx",'S',Hx.V); fd("momentum_SI_alx",'S',alx.V); fd("momentum_SI_ax",'S',ax.V);
        fd("momentum_SI_kx",'S',kx.V); fd("momentum_SI_gx0",'S',gx0.V); fd("momentum_SI_qy",'S',qy.V);
        fd("momentum_SI_qyn",'S',qyn.V); fd("momentum_SI_qyE",'S',qyE.V); fd("momentum_SI_Hy",'S',Hy.V);
        fd("momentum_SI_aly",'S',aly.V); fd("momentum_SI_ay",'S',ay.V); fd("momentum_SI_ky",'S',ky.V);
        fd("momentum_SI_gy0",'S',gy0.V); fd("momentum_SI_fxU",'S',fxU.V); fd("momentum_SI_fxV",'S',fxV.V);
        fd("momentum_SI_fyU",'S',fyU.V); fd("momentum_SI_fyV",'S',fyV.V); fd("momentum_SI_eta0",'S',eta0.V);
        fd("momentum_SI_WL0",'S',WL0.V); fd("momentum_SI_UH0",'S',UH0.V); fd("momentum_SI_VH0",'S',VH0.V);
        fd("momentum_SI_UHa",'S',UHa.V); fd("momentum_SI_VHa",'S',VHa.V); fd("momentum_SI_UHb",'S',UHb.V);
        fd("momentum_SI_VHb",'S',VHb.V); fd("momentum_SI_UHp",'S',UHp.V); fd("momentum_SI_VHp",'S',VHp.V);
        fd("momentum_SI_WHz",'S',WHz.V); fd("momentum_SI_dUH",'S',dUH.V); fd("momentum_SI_dVH",'S',dVH.V);
        fd("momentum_SI_etaI",'S',etaI.V); fd("momentum_SI_alc",'S',alc.V); fd("momentum_SI_rhs",'S',rhs.V);
        fd("momentum_SI_rhsW",'S',rhsW.V);
        fd("momentum_SI_bdiag",'S',bdiag.V); fd("momentum_SI_tU",'S',tU.V);
        fd("momentum_SI_tV",'S',tV.V); fd("momentum_SI_tUH",'S',tUH.V); fd("momentum_SI_tVH",'S',tVH.V);
        fd("momentum_SI_WLk",'S',WLk.V); fd("momentum_SI_etaB",'S',etaB.V);
    }
	sflow_momentum_SI(lexer*, fdm2D*, ghostcell*, sflow_HLL*, sflow_signal_speed*, sflow_reconstruct*, sflow_diffusion*, sflow_pressure*,
                        solver2D*, solver2D*, ioflow*, sflow_fsf*, sflow_forcing*, sixdof*, bool block=false);
	virtual ~sflow_momentum_SI();

	void start(lexer*, fdm2D*, ghostcell*) override final;

    static void check(lexer*);

    // the parts of one step (step_single runs them in this order on one grid; on a multi-block grid and with mesh
    // refinement the stage runner, sflow_amr::si_step, fills the cells around the grids between them and solves
    // the eta system of all grids together)
    void step_single(lexer*, fdm2D*, ghostcell*);
    void si_restore(lexer*, fdm2D*, ghostcell*);            // the water column of the step before
    void si_begin(lexer*, fdm2D*, ghostcell*);              // t^n, faces, friction, open boundaries
    void si_stage_vel(lexer*, fdm2D*, ghostcell*, int);     // velocities of predictor stage s
    void si_stage_rhs(lexer*, fdm2D*, ghostcell*, int);     // explicit terms of stage s
    void si_predictor_end(lexer*, fdm2D*, ghostcell*);      // the cell predictor to the faces
    void si_assemble(lexer*, fdm2D*, ghostcell*);           // face predictor and eta-system coefficients
    void si_rhs(lexer*, fdm2D*, ghostcell*);                // right-hand side of the eta system
    void si_end_faces(lexer*, fdm2D*, ghostcell*);          // face discharges
    void si_end_cells(lexer*, fdm2D*, ghostcell*);          // WL, cells, relaxation zones
    void ship_begin(lexer*, fdm2D*, ghostcell*, bool);      // X 10 3: the body to the end of the step (true), the pressure of the step
    int si_stages() const { return rk; }
    void stage_io_si(int, slice*&, slice*&, slice*&, slice*&);     // input and output UH, VH of stage s

    sflow_amr *prun = nullptr;  // level 0 of a multi-block grid or of the mesh refinement: the stage runner
    bool patch = false;         // a block 1.. or a refined patch (no log file)
    bool amrgrid = false;       // one grid of several: the state of t^n of all cells of the arrays (the filled ones)

    // mesh refinement: called by the advection with the face fluxes of UH, VH before the divergence (the coarse
    // faces next to a finer patch take the fine fluxes, a patch records its box faces)
    std::function<void(lexer*, slice&, slice&, slice&, slice&)> fhook;

    friend class sflow_amr;

    // geometry of the horizontal grid (Cartesian or curvilinear): face i+1/2 (dir 1) or j+1/2 (dir 2) of cell (i,j):
    // length, normal distance of the centres, unit normal towards +i / +j; the cell area
    double flen(lexer*, int, int, int);
    double fdn(lexer*, int, int, int);
    void fnorm(lexer*, int, int, int, double&, double&);
    double carea(lexer*, int, int);

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
        double nx,ny;   // unit normal of the face towards +i / +j
        double qin;     // flux type: inward discharge per width at t^{n+1}
        double A,B;     // radiating: q_in = A - B eta^{n+1}
    };

    void face_depth(lexer*, fdm2D*, ghostcell*);
    void face_weights(lexer*, fdm2D*, int, int, int, double&, double&);    // REEF3D::HYDRO: the cell weights of the face predictor
    void friction(lexer*, fdm2D*);
    void bc_setup(lexer*, fdm2D*, ghostcell*);
    void bc_method(lexer*);
    void stage_vel(lexer*, fdm2D*, ghostcell*, slice&, slice&);
    void stage_rhs(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&, double);
    void advection(lexer*, fdm2D*, ghostcell*);
    void advection_bcg(lexer*, fdm2D*, ghostcell*, slice&, slice&);
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
    // friction factor alpha_c, scratch, right-hand side of the eta system
    slice4 eta0, WL0, UH0, VH0, UHa, VHa, UHb, VHb, UHp, VHp, WHz, dUH, dVH, etaI, alc, rhs, rhsW;

    // implicit diagonal of the open boundaries
    slice4 bdiag;
    slice4 wsrc;        // SEASTATE vortex force (A 751 2): the Stokes mass source -div(M) of the step [m/s]
    slice4 hsrc;        // REEF3D::HYDRO: rain and point sources of the step [m/s]
    bool mass_source = false;

    // interior copies kept through the ghost cell update (thin solids, see rk_stage)
    slice4 tU, tV, tUH, tVH;

    // water depth at the end of the step, the start of the next step (see start);
    // eta with the boundary values of ioflow / iowave (bc_setup)
    slice4 WLk, etaB;

    vector<bc_face> bcf;
    vector<bc_face> bco;    // open boundary faces closed for the step (dry cell or face): no discharge of t^n either
    int method[5];          // iowave open edge method per side cs (B 520 3-6), 0 none

    double theta, g;
    int rk, wd_iter, wd_cells;
    bool qini;

    // A 210 11 (inside A210 = 1): one predictor stage, the advection by the BCG scheme with its source term (as NHFLOW
    // A 510 11), the time levels of nh_levels also without the non-hydrostatic pressure (tlev); bcA: the cells the BCG
    // slopes may use
    bool onestep = false, tlev = false;
    slice4 bcA;
    // A 210 11: a wall face (solid, outside the domain, no open boundary) counts in the cell momentum with its zero
    // gradient and zero non-hydrostatic face term (cell_update, nh_cell, the BCG source term); bcS: the open sides
    slice4 bcS;
    void wall_sides(lexer*);
    double wall(lexer*, int, int, int);

    double Qin, Qout;       // discharge over the open boundaries of the step (theta weighted)
    double vrelax;          // volume change of the relaxation zones of the step
    ofstream silog;

    // non-hydrostatic quadratic pressure (A 220 2/3) in the coupled eta-P system (sflow_momentum_SI_nh.cpp)
    //  qxN, qyN: the non-hydrostatic part of the face discharges of the step (pressure gradient, bed pressure,
    //            A 220 3 dispersion term), added to q in face_update and, as the mean of the faces, to UH, VH
    //  phiN: bed acceleration Dw_b/Dt of the step (the acceleration of the previous step), WHp: WH predictor
    //  Lc0 : h^3 lap(eta^n) of the A 220 3 dispersion term, Unm1, Vnm1: U^{n-1}, V^{n-1} of the cells
    //  time levels of the nonlinear terms (nh_levels): face and cell depths of the gravity and A 220 3 terms at
    //  t^{n+1/2} (Hx, Hy, WLh), the depth of the pressure terms and of the constraint at t^{n+1} (WLq), the discharge
    //  of the advection at t^{n+1/2} (qxa, qya), all extrapolated from t^n and t^{n-1} (Hxo, Hyo, WLo, qxo, qyo);
    //  W advection Adams-Bashforth 2 (AWo: the tendency of the previous step)
    //  face coefficients of the step (nh_solve): q_NH = -R P_R + L P_L + B (L_R - L_L) + c0
    bool nh = false;
    double cb, cwn, adisp, Bdisp;
    slice1 qxN, Hxo, qxo, qxa, nxR, nxL, nxB, nx0, kxEL, kxER, kxPL, kxPR, kx1, kx2;
    slice2 qyN, Hyo, qyo, qya, nyR, nyL, nyB, ny0, kyEL, kyER, kyPL, kyPR, ky1, ky2;
    slice4 phiN, WHp, Lc0, Unm1, Vnm1, WLh, WLo, AWo, WLq;
    slice4 nPW, nEb, nEo, nPo, nU, nV, nUs, nVs, nA, nsx, nsy, nE, nP, nL, nZ, nAE, nAP, nBE, nBP;   // scratch of the solve
    slice4 nGx, nGy;    // curvilinear grid: the depth gradient of the cells (bed acceleration phi)
    bool nh_first = true, nh_krylov = false, nh_x0 = false;
    int nh_cells_active = 0, nh_iter = 0, nh_wd_iter = 0, nh_wd_dry = 0;
    double nh_time = 0.0;
    std::shared_ptr<sflow_si_nh_lu> nhlu;
    void nh_ini(lexer*);
    void nh_levels(lexer*, fdm2D*, ghostcell*);
    int nh_active(lexer*, fdm2D*, int, int);
    void nh_predictor(lexer*, fdm2D*, ghostcell*);
    void nh_solve(lexer*, fdm2D*, ghostcell*);
    void nh_assemble(lexer*, fdm2D*, ghostcell*);    // nA, face and cell coefficients of the coupled system
    void nh_cell(lexer*, fdm2D*, ghostcell*);

    // ship waves (X 10 3): the pressure of the moving body as the surface level it holds, in the face predictor and
    // the cell momentum with the gravity term
    bool ship = false;
    slice4 sps;
    void ship_level(lexer*, fdm2D*, ghostcell*);

    // the coupled system on the grids of the mesh refinement (sflow_momentum_SI_nh_amr.cpp; sflow_amr::si_nh_solve
    // runs the composite solve over the leaf cells of all grids): per grid the cell terms, the h^3 lap(eta) term and
    // the face discharges of (nE, nP), the rows, the preconditioner of the grid and the face terms of the solution
    //  nF1x/y: face discharges with the constants of the mass row, nF2x/y: with all constants (the constraint)
    slice1 nF1x, nF2x;
    slice2 nF1y, nF2y;
    slice4 nH3, nW2, nCs, nRm, nRc, nZE, nZP;
    std::vector<char> nh_own;              // the own cells of the solve (nh_amr_prec_setup)
    std::vector<int> nh_cq;                // the cells in SLICELOOP4 order (nh_amr_prec_setup)
    std::vector<double> nh_dxn, nh_dyn;    // their widths
    std::vector<double> nh_lc;             // the coefficients of lap(nE) of the cells (nh_amr_cells)
    std::vector<int> nh_lq;                // the leaf cells of the solve (nh_amr_leaf)
    std::vector<double> nh_ldx, nh_ldy;    // their widths
    bool nh_rb1 = true, nh_rb2 = true;     // the V-cycle hierarchies of the preconditioner: rebuilt after nh_amr_prec_setup
    void nh_amr_cells(lexer*, fdm2D*);
    void nh_amr_guess(lexer*, fdm2D*);
    void nh_amr_lap(lexer*, ghostcell*);
    void nh_amr_faces(lexer*, bool affine);
    void nh_amr_row(lexer*, fdm2D*, int, int, bool affine, double &ym, double &yc);
    void nh_amr_leaf(lexer*, const std::vector<int> &leaf);   // the leaf cells (i, j pairs) for nh_amr_rows
    void nh_amr_rows(lexer*, bool affine, double *y);           // nh_amr_row of all leaf cells
    void nh_amr_prec_setup(lexer*, ghostcell*, const std::function<bool(int,int)> &own);
    void nh_amr_prec(lexer*, ghostcell*, const std::function<bool(int,int)> &own);
    void nh_amr_finish(lexer*, fdm2D*, ghostcell*);
    int nh_amr_wetdry(lexer*, fdm2D*, const std::function<bool(int,int)> &leaf);
};

#endif
