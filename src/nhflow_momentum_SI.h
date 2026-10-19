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

#ifndef NHFLOW_MOMENTUM_SI_H_
#define NHFLOW_MOMENTUM_SI_H_

#include"nhflow_momentum_func.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include<vector>
#include<fstream>

class si_eta2D;
class solver2D;
class wind;

using namespace std;

// ---------------------------------------------------------------------------
// NHFLOW hydrostatic semi-implicit free surface (A 534 1), theta method (A 535).
//
// For long simulations (tide, river, coastal ocean): the gravity wave speed sqrt(gh) no longer
// limits the time step, only advection does (nhflow_timestep_SI, cap A 536).
// Hydrostatic only (A 520 0). Separate from the explicit integrators (nhflow_momentum_RK2/RK3),
// which are not touched; the shared parts of nhflow_momentum_func (sigma metrics, omega,
// velocities, bed friction) and the 2D eta system of the semi-implicit SFLOW (si_eta2D) are used.
//
// Prognostic: WL, UH, VH, WH per sigma layer (cell centres) and the layer face discharges
// Qx, Qy (per unit sigma and unit face width, as d->FEx/FEy), which carry the mass flux.
// The vertical structure enters the eta system through the column operator of a face (and a cell)
//
//     A = I - dt L_v + dt V_omega + dt gamma e_0 e_0^T,   r = A^-1 1,   alpha_f = sum_k dsig_k r_k
//
// L_v the vertical diffusion (nu = VISC + EV, sigma metric 1/h^2, the vertical rows of nhflow_idiff),
// V_omega the vertical advection (upwind, omega of the step before; conservative, so alpha_f does not
// change and omega does not limit the time step), gamma the implicit log-law bed friction of the bottom
// layer (A 519 1). The hydrostatic pressure
// gradient -g WL grad(eta) is the same in every layer, so the column reduces to one 2D Helmholtz
// problem for eta^{n+1} (Casulli and Cattani 1994). The horizontal diffusion is explicit in the
// predictor (along the sigma layers, compact Laplacian); the diffusion is active as nhflow_idiff
// would be (A 512 > 0 or A 560 > 0). Vertical diffusion, bed friction and the pressure gradient act
// together in one implicit column step, so that the steady state does not depend on dt.
//
// Grids: Cartesian and curvilinear (G 50 2, nhflow_curv: face lengths and normals, the normal
// distance of the centres; with A 539 1 the non-orthogonal part of the face gradient of eta, explicit;
// the gravity term of a cell from the mean normals of its open faces per index direction). Multi-block grids (G 50 1):
// every block has its own object, the stage runner nhflow_amr::si_step (nhflow_amr_si.cpp) runs the
// parts of the step on all blocks, fills the cells around them in between and solves one eta system
// over the cells of all blocks.
//
// One time step:
//  1. face depth H_f = eta_up + d_f (eta upwind of the column discharge, d_f the mean still
//     water depth; next to a dry cell the eta of the wet cell and the smaller still water depth,
//     closed between two dry cells), open faces H_f > A 544; bed friction rates gamma of the
//     bottom layer on the faces and in the cells
//  2. open boundaries from the ghost cells of ioflow / iowave (bc_setup), as the semi-implicit
//     SFLOW: discharge (layer fluxes given), water level (Dirichlet), radiating (Riemann,
//     Flather, free outflow: q_in = A - B eta^{n+1}, implicit; layer profile of the upwind side)
//  3. predictor of the explicit terms, SSP-RK2/3 (A 510) on the frozen sigma grid of t^n:
//     UH* = UH^n + dt L(UH),  L = -div_h(T u_up) + div_h(nu_h grad u) WL + sources,
//     T the layer transport of the previous step (the vertical advection is in the column operator)
//     (no gravity term, no c in the fluxes: only the advective Courant number counts)
//  4. layer face predictor, momentum interpolation:  Q^_k = avg(UH*_k).n,  Q~ = A_f^-1 Q^
//  5. column: q^E = sum_k dsig_k Q~_k - dt g H_f alpha_f (1-theta) grad_f eta^n
//  6. r(WL(eta^{n+1})) - div( a_f grad_f eta^{n+1} ) = WL^n - dt div( theta q^E + (1-theta) q^n )
//     a_f = theta dt (g theta dt alpha_f H_f + K_f)/dxp_f (K_f: mass flux dissipation, A 538)
//     wetting and drying (A 537 1): Newton iteration on WL(eta) = max(eta + d, 0) (si_eta2D)
//  7. layer face discharges  Q_k^{n+1} = Q~_k - dt g H_f (theta G^{n+1} + (1-theta) G^n) r_k - K_f G^{n+1}
//  8. WL^{n+1} from the conservative divergence of the layer transport T_k = theta Q_k^{n+1} +
//     (1-theta) Q_k^n (exact volume); omega from the layer continuity with the same T_k
//  9. UH^{n+1} = A_c^-1 UH* - dt g [H G]_c r_c, [H G]_c the mean of H_f (theta G^{n+1} + (1-theta) G^n)
//     over the open faces; WH: the predictor with the implicit vertical advection (no pressure)
// 10. relaxation zones once per step, the layer face discharges follow the change
//
// Not yet (stops in check()): non-hydrostatic pressure (A 520 > 0), mesh refinement (G 1), VRANS,
// 6DOF / floating bodies, SEASTATE, side-wall friction A 519 2, membranes, FEM coupling, patch
// boundary conditions.
// ---------------------------------------------------------------------------

class nhflow_momentum_SI final : public nhflow_momentum_func
{
public:
	nhflow_momentum_SI(lexer*, fdm_nhf*, ghostcell*, sixdof*, vrans_nhflow*, nhflow_forcing*, bool block=false);
	virtual ~nhflow_momentum_SI();

	void start(lexer*, fdm_nhf*, ghostcell*, ioflow*, nhflow_signal_speed*, nhflow_reconstruct*, nhflow_convection*,
                nhflow_diffusion*, nhflow_pressure*, solver*, solver*, nhflow*, nhflow_fsf*, nhflow_turbulence*, vrans_nhflow*) override final;

    static void check(lexer*);

    // the parts of one step (start runs them in this order on one grid; on a multi-block grid the
    // stage runner, nhflow_amr::si_step, fills the cells around the blocks between them and solves
    // the eta system of all blocks together)
    void si_restore(lexer*, fdm_nhf*, ghostcell*);                      // the water column of the step before
    void si_begin(lexer*, fdm_nhf*, ghostcell*, nhflow_stage_obj&);     // t^n, faces, friction, open boundaries
    void si_stage_vel(lexer*, fdm_nhf*, ghostcell*, nhflow_stage_obj&, int);   // velocities of predictor stage s
    void si_stage_rhs(lexer*, fdm_nhf*, ghostcell*, nhflow_stage_obj&, int);   // explicit terms, stage s
    void si_predictor_end(lexer*, fdm_nhf*, ghostcell*);                // the predictor to the faces (UHD, VHD)
    void si_assemble(lexer*, fdm_nhf*, ghostcell*);                     // column operator, eta system
    void si_end(lexer*, fdm_nhf*, ghostcell*, nhflow_stage_obj&);       // faces, WL, cells, omega
    int si_stages() const { return rk; }

    bool patch;             // a block 1.. of a multi-block grid (no log file, no global output)

    friend class nhflow_amr;

    // the stage interface of the mesh refinement (nhflow_amr) is not used: one theta step
    int stages() const override final { return 1; }
    void step_begin(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&) override final;
    void phase_F(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) override final;
    void phase_M(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) override final;
    void phase_D(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int,int) override final;
    void phase_P1(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) override final;
    void phase_P2(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) override final;
    void phase_E(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) override final;
    double stage_alpha(int) const override final { return 1.0; }
    slice& stage_WL(fdm_nhf*,int) override final;

    // checkpoints (P 48, I 45): the face discharges and the state of the step
    void ckpt_fields(const std::function<void(const char*, char, double*)> &fd, const std::function<void(const char*, char, int*)> &fi) override
    {
        fsfbody_ckpt_fields(fd);
        fd("si_qx",'S',qx.V); fd("si_qxn",'S',qxn.V); fd("si_qxE",'S',qxE.V); fd("si_Hx",'S',Hx.V);
        fd("si_alx",'S',alx.V); fd("si_gmx",'S',gmx.V); fd("si_ax",'S',ax.V); fd("si_kx",'S',kx.V);
        fd("si_gx0",'S',gx0.V); fd("si_cx0",'S',cx0.V); fd("si_qy",'S',qy.V); fd("si_qyn",'S',qyn.V);
        fd("si_qyE",'S',qyE.V); fd("si_Hy",'S',Hy.V); fd("si_aly",'S',aly.V); fd("si_gmy",'S',gmy.V);
        fd("si_ay",'S',ay.V); fd("si_ky",'S',ky.V); fd("si_gy0",'S',gy0.V); fd("si_cy0",'S',cy0.V);
        fd("si_eta0",'S',eta0.V); fd("si_WL0",'S',WL0.V); fd("si_etaI",'S',etaI.V); fd("si_gmc",'S',gmc.V);
        fd("si_rhsW",'S',rhsW.V); fd("si_bdiag",'S',bdiag.V); fd("si_tS",'S',tS.V); fd("si_tE",'S',tE.V);
        fd("si_WLk",'S',WLk.V); fd("si_etaB",'S',etaB.V); fd("si_WLb",'S',WLb.V); fd("si_Qx",'C',Qx);
        fd("si_Qy",'C',Qy); fd("si_Qxn",'C',Qxn); fd("si_Qyn",'C',Qyn); fd("si_Tx",'C',Tx);
        fd("si_Ty",'C',Ty); fd("si_QAx",'C',QAx); fd("si_QAy",'C',QAy); fd("si_Rx",'C',Rx);
        fd("si_Ry",'C',Ry); fd("si_UH0",'C',UH0); fd("si_VH0",'C',VH0); fd("si_WH0",'C',WH0);
        fd("si_UHa",'C',UHa); fd("si_VHa",'C',VHa); fd("si_WHa",'C',WHa); fd("si_UHb",'C',UHb);
        fd("si_VHb",'C',VHb); fd("si_WHb",'C',WHb); fd("si_UHp",'C',UHp); fd("si_VHp",'C',VHp);
        fd("si_WHp",'C',WHp); fd("si_UHD",'C',UHD); fd("si_VHD",'C',VHD); fd("si_WHD",'C',WHD);
        fd("si_dUH",'C',dUH); fd("si_dVH",'C',dVH);
    }
    void ckpt_restored() override { qini = true; }   // the face discharges of the checkpoint
    double* stage_UH(fdm_nhf*,int,int) override final;

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
        double qin;     // flux type: inward column discharge per width
        double A,B;     // radiating: q_in = A - B eta^{n+1}
        vector<double> Qg;          // flux type: layer discharges (face direction); radiating: layer profile of the inflow
        vector<double> wo;          // radiating: layer profile of the outflow (interior)
        vector<double> ug,vg,wg;    // velocities of the ghost cell per layer (advection of the inflow)
    };

    void face_depth(lexer*, fdm_nhf*, ghostcell*);
    void friction(lexer*, fdm_nhf*);
    void bc_setup(lexer*, fdm_nhf*, ghostcell*);
    void bc_method(lexer*);
    void stage_arrays(int, double*&, double*&, double*&, double*&, double*&, double*&, double&);
    solver2D* new_solver2D(lexer*, ghostcell*);
    void advection(lexer*, fdm_nhf*, ghostcell*, nhflow_reconstruct*);
    void hdiffusion(lexer*, fdm_nhf*);
    void face_update(lexer*, fdm_nhf*, ghostcell*);
    void eta_update(lexer*, fdm_nhf*, ghostcell*, nhflow_fsf*);
    void cell_update(lexer*, fdm_nhf*, ghostcell*);
    void relax(lexer*, fdm_nhf*, ghostcell*);
    void finish(lexer*, fdm_nhf*, ghostcell*);
    void level_ghosts(lexer*, fdm_nhf*, ghostcell*);
    void log(lexer*, fdm_nhf*, ghostcell*);
    void write_log(lexer*, double, double, double, double, int, int, int);

    double bed_gamma(lexer*, double, double, double);

    // geometry of the horizontal grid (Cartesian or curvilinear, G 50 with nhflow_curv): face i+1/2
    // (dir 1) or j+1/2 (dir 2) of cell (i,j): length, normal distance of the centres, unit normal
    // towards +i / +j, the non-orthogonal part of a face normal gradient; the cell area
    double flen(lexer*, int, int, int);
    double fdn(lexer*, int, int, int);
    void fnorm(lexer*, int, int, int, double&, double&);
    double fcorr(lexer*, slice&, int, int, int);
    void cgrad(lexer*, slice&, int, int, double&, double&);
    double carea(lexer*, int, int);
    
    // column operator A = I - dt L_v + dt gamma e_0 e_0^T of a column of depth h with the layer
    // viscosities nu: the off-diagonal parts cu (k+1), cl (k-1); solve A x = b
    void col_coef(lexer*, fdm_nhf*, double h, int, int, int, int, const double *nu, double *cu, double *cl, double *cd);
    void col_solve(lexer*, const double *cu, const double *cl, const double *cd, double dtgam, const double *b, double *x);
    void col_nu(lexer*, fdm_nhf*, int, int, int, int, double*);
    bool visc;
    vector<double> wcu, wcl, wcd, wnu, wa, wb, wr, wx;   // column work arrays (knoz)
    double dissipation(lexer*, double, double, double);
    double div(lexer*, slice&, slice&);
    double layer_div(lexer*, double*, double*);

    void exchange1(lexer*, ghostcell*, slice&);
    void exchange2(lexer*, ghostcell*, slice&);
    void exchange4(lexer*, ghostcell*, slice&);
    void exchangeF(lexer*, ghostcell*, double*);
    void exchangeC(lexer*, ghostcell*, double*);

    si_eta2D *psi;
    solver2D *psolv2D;
    int psi_iter;
    double psi_time;
    wind *pwind;
    sixdof *p6dof;
    nhflow_forcing *pnhfdf;
    vrans_nhflow *pvrans;
    ioflow *pflow;
    nhflow_turbulence *pturb;
    solver *psolv;

    // column face quantities of the step: discharge q^{n+1}, q^n, q^E, depth H_f (0: closed face),
    // conveyance alpha_f, bed friction rate gamma, eta-system coefficient a_f, dissipation K_f,
    // face gradient of eta^n, its non-orthogonal part (curvilinear grid)
    slice1 qx, qxn, qxE, Hx, alx, gmx, ax, kx, gx0, cx0;
    slice2 qy, qyn, qyE, Hy, aly, gmy, ay, ky, gy0, cy0;

    // cell state of the step: t^n, implicit eta, bed friction rate, right-hand side, scratch
    slice4 eta0, WL0, etaI, gmc, rhsW, bdiag, tS, tE;

    // water depth at the end of the step, the start of the next (see start); boundary copies
    slice4 WLk, etaB, WLb;

    // layer face discharges: prognostic Q^{n+1}, Q^n, layer transport of the step (theta weighted),
    // column parts Q~ = A^-1 (I - dt L_v) Q^ and r = A^-1 1
    double *Qx, *Qy, *Qxn, *Qyn, *Tx, *Ty, *QAx, *QAy, *Rx, *Ry;

    // layer state: t^n, RK stages, predictor, implicit diffusion
    double *UH0, *VH0, *WH0, *UHa, *VHa, *WHa, *UHb, *VHb, *WHb, *UHp, *VHp, *WHp;
    double *UHD, *VHD, *WHD;
    double *dUH, *dVH;

    vector<bc_face> bcf;
    int method[5];          // iowave open edge method per side cs (B 520 3-6), 0 none

    double theta, g;
    int rk, wd_iter, wd_cells, gcval_eta;
    bool qini;

    double Qin, Qout;       // discharge over the open boundaries of the step (theta weighted)
    double vrelax;          // volume change of the relaxation zones of the step
    ofstream silog;
};

#endif
