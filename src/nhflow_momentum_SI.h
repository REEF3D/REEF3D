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
#include<functional>
#include<memory>

class si_eta2D;
struct sflow_si_nh_lu;
class nhflow_si_nhmg;
class nhflow_tracer_SI;
class solver2D;
class wind;
class met_forcing;
class patchBC_interface;
class patch_obj;

using namespace std;

// ---------------------------------------------------------------------------
// Layout of the SI-private 3D arrays (the work arrays of the step: Qx .. dVH, the non-hydrostatic NQx .. NSY, PHI,
// NQcx/NQcy, the time levels accU .. FSy, nUZx/nUZy, hoPS): cell (i, j, k) at lso + i lsi + j lsj + k (layout_ini).
// The full layout (lsi = jmax kmax, lsj = kmax, lso = -(imin lsi + jmin lsj + kmin)) is that of IJK: the same index.
// The compact layout (one Cartesian grid, not a block, no curvilinear or refined grid): the cells i in [-3, knox+3),
// j in [-3, knoy+3) (one row j = 0 with j_dir 0) and k in [0, knoz] (the cells and the top face), about 2/5 of the
// entries in 3D and 1/15 in 1D. Environment REEF3D_SI_LAYOUT=full keeps the full layout.
// The shared fields of fdm_nhf (U, UH, F, ...) and the other arrays keep IJK. P-macros index the private arrays;
// code that indexes a private array with IJK (ghost cells, the open boundaries, MPI exchange, checkpoints, tracers)
// gets it through lfull / lback (a copy in the full layout with the compact layout), lset / lget (from / to a full
// array). Environment REEF3D_SI_LSHIFT (test, full layout): the private arrays are shifted by lshift entries, so that
// an access with IJK instead of a P-macro reads another cell. Compiled with -DSI_LCHECK: every P-macro index is
// checked against the layout (lchk).
// ---------------------------------------------------------------------------
#ifndef SI_LCHECK
#define PIJK (lso + i*lsi + j*lsj + k)
#define PIp1JK (PIJK + lsi)
#define PIm1JK (PIJK - lsi)
#define PIp2JK (PIJK + 2*lsi)
#define PIm2JK (PIJK - 2*lsi)
#define PIJp1K (PIJK + lsj)
#define PIJm1K (PIJK - lsj)
#define PIJp2K (PIJK + 2*lsj)
#define PIJm2K (PIJK - 2*lsj)
#define PIJKp1 (PIJK + 1)
#define PIJKm1 (PIJK - 1)
#define PIX(a,b,c) (lso + (a)*lsi + (b)*lsj + (c))
#else
#define PIJK lchk(i,j,k,__FILE__,__LINE__)
#define PIp1JK lchk(i+1,j,k,__FILE__,__LINE__)
#define PIm1JK lchk(i-1,j,k,__FILE__,__LINE__)
#define PIp2JK lchk(i+2,j,k,__FILE__,__LINE__)
#define PIm2JK lchk(i-2,j,k,__FILE__,__LINE__)
#define PIJp1K lchk(i,j+1,k,__FILE__,__LINE__)
#define PIJm1K lchk(i,j-1,k,__FILE__,__LINE__)
#define PIJp2K lchk(i,j+2,k,__FILE__,__LINE__)
#define PIJm2K lchk(i,j-2,k,__FILE__,__LINE__)
#define PIJKp1 lchk(i,j,k+1,__FILE__,__LINE__)
#define PIJKm1 lchk(i,j,k-1,__FILE__,__LINE__)
#define PIX(a,b,c) lchk((a),(b),(c),__FILE__,__LINE__)
#endif
#define PEXT(x) ((x) + lshift)

// ---------------------------------------------------------------------------
// NHFLOW hydrostatic semi-implicit free surface (A 510 11/12/13), theta method (A 535).
//
// For long simulations (tide, river, coastal ocean): the gravity wave speed sqrt(gh) does not
// limit the time step, only advection does (nhflow_timestep_SI, cap A 536).
// Separate from the explicit integrators (nhflow_momentum_RK2/RK3); the shared parts of
// nhflow_momentum_func (sigma metrics, omega, velocities, bed friction) and the 2D eta system of
// the semi-implicit SFLOW (si_eta2D) are used. The description below is the hydrostatic step;
// with A 520 > 0 the non-hydrostatic pressure joins eta in one coupled eta-phi system
// (nhflow_momentum_SI_nh.cpp, Cartesian, one grid).
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
// problem for eta^{n+1} (Casulli and Cattani 1994). The horizontal diffusion (along the sigma layers, compact
// Laplacian) is explicit in the predictor, with A 512 2 on one grid implicit after the predictor
// (nhflow_momentum_SI_hdiff.cpp); the diffusion is active as nhflow_idiff would be (A 512 > 0 or A 560 > 0). Vertical diffusion, bed friction and the pressure gradient act
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
// Mesh refinement (G 1): every patch has its own object; the stage runner runs the parts on all grids, finest
// first where the coarse faces next to a patch take the fine fluxes (advection, discharges), and one eta system over
// the leaf cells of all grids (reefamr_si.cpp); not with subcycling (G 7 1), vertical refinement (G 6 1) or tracers.
//
// Fixed solids with direct forcing (A 580 1, A 581-590; one grid): the explicit forcing of nhflow_forcing on the
// cell predictor and on the cells at the end of the step (df_*); the layer face discharges are not forced.
//
// Not available (stops in check()): VRANS, 6DOF / floating bodies, side-wall
// friction A 519 2, membranes, FEM coupling, patch boundary conditions, SEASTATE surfbeat (A 770 1);
// the non-hydrostatic pressure and wave breaking on one Cartesian grid only.
// ---------------------------------------------------------------------------

class nhflow_breaking;

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
    void si_assemble(lexer*, fdm_nhf*, ghostcell*);                     // column operator, face coefficients
    void si_rhs(lexer*, fdm_nhf*, ghostcell*);                          // right-hand side of the eta system
    void si_old_closed(lexer*, ghostcell*);                             // no old discharge on the closed faces
    void si_density(lexer*, fdm_nhf*);                                  // A 620: the density, its pressure gradient
    void si_end(lexer*, fdm_nhf*, ghostcell*, nhflow_stage_obj&);       // faces, WL, cells, omega
    void si_end_faces(lexer*, fdm_nhf*, ghostcell*);                    // the face discharges (si_end in two parts)
    void si_end_cells(lexer*, fdm_nhf*, ghostcell*, nhflow_stage_obj&); // WL, cells, omega
    void step_single(lexer*, fdm_nhf*, ghostcell*, nhflow_stage_obj&);  // the step of one grid

    // mesh refinement: called by the advection of a component (0 U, 1 V, 2 W) with its layer face fluxes before the
    // divergence (the coarse faces next to a finer patch take the fine fluxes, a patch records its box faces)
    std::function<void(lexer*, int, double*, double*)> fhook;
    int si_stages() const { return rk; }

    bool patch;             // a block 1.. of a multi-block grid (no log file, no global output)

    // patch boundary conditions (B 411-422 on B 440-442; patchBC_2D on the slice faces gcbsl4, set by the driver): the
    // inlets as discharge faces, the outlets with a water level as level faces, the others as free outflow (bc_setup);
    // the ghost cells of the patch faces (patch_ghosts): the inlet velocity on all layers, zero gradient at the outlets
    patchBC_interface *pBC = nullptr;
    void patch_ghosts(lexer*, fdm_nhf*, double*, double*, double*, double*, double*, double*);
    const patch_obj* patch_of(int, int, int) const;

    friend class nhflow_amr;
    friend class nhflow_tracer_SI;

    // tracers (A 600): salinity, temperature (nhflow_tracer_SI), after the momentum step
    nhflow_tracer_SI *ptrc = nullptr;

    // density from the tracers and its pressure gradient (A 620, nhflow_momentum_SI_density.cpp)
    static double eos(lexer*, double S, double T);
    static double internal_speed(lexer*, fdm_nhf*, int, int);
    // A 512 2 on one grid without wave breaking: the horizontal diffusion is implicit (no time step limit of it)
    static bool hdiff_implicit(lexer*);
    // the advective time step limit of a cell (cmin of dt = 2 cfl cmin), A 549 0 and 1
    static double adv_cmin(lexer*, double u, double v, double ci, int i, int j, double cfl);
    static int elm_mode(lexer *p) { return p->A549; }      // A 549 m C: the mode and the Courant number
    static double elm_c(lexer *p) { return p->A549_c; }
    bool trc_mask = false;

    // A 623 with A 547: the vertical eddy diffusivity of the tracers K = nu_n/sigma_t f_h(Ri) + nu0/sigma_t
    // (nullptr: EV/sigma_t), the buoyancy frequency
    double *KV = nullptr, *N2 = nullptr;   // N2 also with A 627 + A 623 (the floor damped by f_m(Ri))
    void tracer_step(lexer*, fdm_nhf*, ghostcell*, int phase=0);

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
        // the private arrays in the layout of IJK (px, nhflow_momentum_SI.h)
        // compact layout: each private array through a full copy (gather, write or read, scatter)
        auto px = [&](const char *n, double *a)
        {
            if(a==nullptr) return;
            double *F = lfull(nullptr,a,0);
            fd(n,'C',F);
            lback(nullptr,a,0);
        };
        fsfbody_ckpt_fields(fd);
        fd("si_qx",'S',qx.V); fd("si_qxn",'S',qxn.V); fd("si_qxE",'S',qxE.V); fd("si_Hx",'S',Hx.V);
        fd("si_alx",'S',alx.V); fd("si_gmx",'S',gmx.V); fd("si_ax",'S',ax.V); fd("si_kx",'S',kx.V);
        fd("si_gx0",'S',gx0.V); fd("si_cx0",'S',cx0.V); fd("si_qy",'S',qy.V); fd("si_qyn",'S',qyn.V);
        fd("si_qyE",'S',qyE.V); fd("si_Hy",'S',Hy.V); fd("si_aly",'S',aly.V); fd("si_gmy",'S',gmy.V);
        fd("si_ay",'S',ay.V); fd("si_ky",'S',ky.V); fd("si_gy0",'S',gy0.V); fd("si_cy0",'S',cy0.V);
        fd("si_eta0",'S',eta0.V); fd("si_WL0",'S',WL0.V); fd("si_etaI",'S',etaI.V); fd("si_gmc",'S',gmc.V);
        fd("si_rhsW",'S',rhsW.V); fd("si_bdiag",'S',bdiag.V); fd("si_tS",'S',tS.V); fd("si_tE",'S',tE.V);
        fd("si_WLk",'S',WLk.V); fd("si_etaB",'S',etaB.V); fd("si_WLb",'S',WLb.V); px("si_Qx",Qx);
        px("si_Qy",Qy); px("si_Qxn",Qxn); px("si_Qyn",Qyn); px("si_Tx",Tx);
        px("si_Ty",Ty); px("si_QAx",QAx); px("si_QAy",QAy); px("si_Rx",Rx);
        px("si_Ry",Ry); px("si_UH0",UH0); px("si_VH0",VH0); px("si_WH0",WH0);
        px("si_UHa",UHa); px("si_VHa",VHa); px("si_WHa",WHa); px("si_UHb",UHb);
        px("si_VHb",VHb); px("si_WHb",WHb); px("si_UHp",UHp); px("si_VHp",VHp);
        px("si_WHp",WHp); px("si_UHD",UHD); px("si_VHD",VHD); px("si_WHD",WHD);
        px("si_dUH",dUH); px("si_dVH",dVH);
        // the time levels of the step before (second order in time, tl2; nh_begin, nh_levels)
        px("si_accU",accU); px("si_accV",accV); px("si_accW",accW);
        px("si_vaU",vaU); px("si_vaV",vaV); px("si_vaW",vaW);
        fd("si_Bx",'C',Bx); fd("si_By",'C',By); fd("si_BFx",'C',BFx); fd("si_BFy",'C',BFy);
        fd("si_nOmP",'F',nOmP.data()); fd("si_nHxo",'S',nHxo.data()); fd("si_nHyo",'S',nHyo.data());
        if(!nWLo.empty()) fd("si_nWLo",'S',nWLo.data());
        tracer_ckpt_fields(fd);
    }
    void ckpt_restored() override { qini = true; }   // the face discharges of the checkpoint
    double* stage_UH(fdm_nhf*,int,int) override final;
    void tracer_ckpt_fields(const std::function<void(const char*, char, double*)> &fd);

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
    void mlflux_faces(lexer*, fdm_nhf*, ghostcell*);
    void friction(lexer*, fdm_nhf*);
    void density(lexer*, fdm_nhf*);
    void baroclinic(lexer*, fdm_nhf*);
    void baroclinic_face(lexer*, fdm_nhf*, int dir, int fi, int fj, double *bn);
    void baroclinic_column(lexer*, fdm_nhf*, int, int, double *P);
    double *Bx = nullptr, *By = nullptr;     // (1/rho0) grad_z p' of the cells (A 620)
    double *BFx = nullptr, *BFy = nullptr;   // its normal component on the faces (the layer face discharges)
    vector<double> bpl, bpr, bcr, bcz, bcd, bce;
    void eddy_algebraic(lexer*, fdm_nhf*);     // A 547: parabolic eddy viscosity from the bed friction velocity
    void eddy_floor(lexer*, fdm_nhf*);         // A 627: turbulence model, eddy viscosity >= c times that profile (thin water)
    double col_ustar(lexer*, fdm_nhf*, double h);   // bed (+ wind) friction velocity of the column (i,j), 0 dry
    // A 549 1: semi-Lagrangian momentum advection (nhflow_momentum_SI_elm.cpp), the velocities at the feet of the
    // characteristics of the cells
    void elm_departure(lexer*, fdm_nhf*);
    bool elm_fluid(lexer*, int, int);
    double *ELU = nullptr, *ELV = nullptr, *ELW = nullptr;
    void bc_setup(lexer*, fdm_nhf*, ghostcell*);
    void bc_method(lexer*);
    void stage_arrays(int, double*&, double*&, double*&, double*&, double*&, double*&, double&);
    solver2D* new_solver2D(lexer*, ghostcell*);
    void advection(lexer*, fdm_nhf*, ghostcell*, nhflow_reconstruct*);
    void si_corrector(lexer*, fdm_nhf*, ghostcell*);
    void advection_bcg(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*);
    std::vector<double> bcgF;     // advection_bcg: the half-step field f + dt/2 S (full layout)
    void hdiffusion(lexer*, fdm_nhf*, ghostcell*);
    void hdiffusion_brk(lexer*, fdm_nhf*, ghostcell*);
    // implicit horizontal diffusion (nhflow_momentum_SI_hdiff.cpp): (I - dt L_h) u = u* per layer and velocity
    // component after the predictor, the 2D systems with si_eta2D
    void hdiffusion_implicit(lexer*, fdm_nhf*, ghostcell*);
    bool hdimp = false;
    // the work arrays of the implicit horizontal diffusion (full layout, IJK): face coefficients x, y, diagonal; per
    // component the solution, residual, search direction and A p
    std::vector<double> hiF[3], hdX[3], hdR[3], hdP[3], hdA[3];
    int hd_iter = 0, hd_sys = 0;
    double hd_time = 0.0;
    std::vector<double> hdU, hdV, hdW, hdN;
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
    void col_vadv_central(lexer*, int, double, double, double);   // central vertical advection of one layer of col_coef
    void col_solve(lexer*, const double *cu, const double *cl, const double *cd, double dtgam, const double *b, double *x);
    void col_nu(lexer*, fdm_nhf*, int, int, int, int, double*);
    bool visc;
    vector<double> wcu, wcl, wcd, wnu, wa, wb, wr, wx;   // column work arrays (knoz)
    vector<double> wnx, wny;                             // the phi face terms of a cell (nh_cell_faces)
    double dissipation(lexer*, double, double, double);
    double div(lexer*, slice&, slice&);
    double layer_div(lexer*, double*, double*);
    double layer_div_p(lexer*, double*, double*);   // of the private arrays

    void exchange1(lexer*, ghostcell*, slice&);
    void exchange2(lexer*, ghostcell*, slice&);
    void exchange4(lexer*, ghostcell*, slice&);
    void exchangeF(lexer*, ghostcell*, double*);
    void exchangeC(lexer*, ghostcell*, double*);
    void pexchF(lexer*, ghostcell*, double*);    // a private array (lfull)
    void pexchC(lexer*, ghostcell*, double*);

    // non-hydrostatic pressure (A 520 > 0) coupled with eta^{n+1}, after the multilayer solver of Popinet 2020
    // (nhflow_momentum_SI_nh.cpp): phi = p/rho on the layer interfaces (PHI, k: the lower interface of layer k), the
    // layer face discharges of phi (NQx, NQy, per unit sigma), w* of the cells (WHS), the interface slopes of the
    // cells (NSX, NSY, nodes 0 .. K), the layer thicknesses at t^{n+1/2} (nWLh), the non-hydrostatic cells (nAct)
    bool nh = false;
    double *NQx = nullptr, *NQy = nullptr, *WHS = nullptr, *PHI = nullptr, *PHIo = nullptr, *NSX = nullptr, *NSY = nullptr;
    // hMx, hMy: the face depths of the mass flux of the predictor (mlflux_faces, nh_begin)
    std::vector<double> nWLh, nWLo, nHxo, nHyo, hMx, hMy;
    // [u . grad z] with the geometry of t^{n+1} from the solve (nh_prepare, nh_residual, nh_assemble): the cells of the
    // centred interface slopes (offsets east, west / north, south; -9 none) and their distance; NSFx, NSFy: 1 where
    // the slope takes the change of eta (not at its limit)
    std::vector<signed char> nSx, nSy;
    std::vector<double> nSdx, nSdy;
    double *NSFx = nullptr, *NSFy = nullptr;
    std::vector<double> nUZx, nUZy;

    // fourth-order horizontal operators (nhflow_momentum_SI_ho.cpp): hoC the fluid cells of the stencils (with
    // those of the neighbouring ranks), hoE eta of the estimate pass, NQcx, NQcy the explicit part of the phi face terms
    bool ho = false, ho_ux = false, ho_uy = false, ho_phi = false;
    double *NQcx = nullptr, *NQcy = nullptr;
    std::vector<double> hoPS;
    void ho_begin(lexer*,fdm_nhf*,ghostcell*);
    void ho_estimate(lexer*,ghostcell*);
    bool ho_face(lexer*,int,int,int);
    bool ho_cell(lexer*,fdm_nhf*,int,int,int);
    void ho_gradient(lexer*,fdm_nhf*,int,double&,double&);
    void ho_predictor(lexer*,int,const double*,double*);
    void ho_cell_gravity(lexer*,fdm_nhf*,double&,double&);
    void ho_phi_corr(lexer*,fdm_nhf*,ghostcell*);
    void ho_phi_add(lexer*);            // tl2: UH^{n+1}, VH^{n+1} estimated for [u . grad z] (nh_prepare)
    std::vector<double> nOmP;                  // omega of the step before (extrapolation to t^{n+1/2}), NaN: none yet
    std::vector<char> nAct;
    std::shared_ptr<sflow_si_nh_lu> nhlu;
    std::shared_ptr<nhflow_si_nhmg> nhmg;
    static std::vector<char> nh_mask(int B, bool yd);    // the stored entries of the blocks (nh_assemble)
    std::shared_ptr<nhflow_si_nhmg> nhme, nhmp;     // A 507 1: the eta system (one unknown per column) and the phi system (K)
    int nh_lagc = 0, nh_lagit = 0, nh_lagnact = -1; // the preconditioner kept for some steps: steps since its refresh, the
    double nh_lagdt = 0.0;                          // iterations of that step, the active columns and dt then
    std::vector<char> nh_lagpat;                    // the wet pattern (pw) and active columns of the preconditioner
    int nh_reb = 0;                                 // the preconditioner rebuilt inside FGMRES (stall) in this solve
    int nh_iter = 0, nh_iter1 = 0;   // the iterations of the last solve and of the first pass of the step
    double nh_time = 0.0, nh_res = 0.0, nh_brk = 0.0;
    // wetting and drying of the coupled solve (A 537 1): pw of the mass rows, the iterations of the pattern and the dry
    // columns
    std::vector<double> nPw;
    int nh_wdit = 1, nh_dry = 0;
    // a further pass of the pattern (A 537 1, nh_wdpass); with tl2 (nh_warm) FGMRES starts from the solution of the pass
    // before and stops at the smaller of the residual of the first pass (nh_tolabs, metres) and A 506 of its own start
    bool nh_wdpass = false, nh_warm = false;
    double nh_tolabs = 0.0;
    void nh_ini(lexer*);
    void nh_free();
    void nh_levels(lexer*, fdm_nhf*);
    double nh_wl1(lexer*, int, int);      // the depth of t^{n+1}, extrapolated (the constraint)
    void nh_begin(lexer*, fdm_nhf*, ghostcell*);
    void nh_prepare(lexer*, fdm_nhf*, ghostcell*);
    void nh_faces(lexer*, fdm_nhf*);
    // the coupled system assembled directly (instead of probing nh_residual): the derivatives of the layer discharges
    // of the faces, d Q_k / d eta (the cells left and right) and d Q_k / d phi (phi_k and phi_{k+1} left and right)
    void nh_dfaces(lexer*, fdm_nhf*, ghostcell*);
    void nh_assemble(lexer*, const std::function<void(int,int,int,int,int,double)>&);
    // the layers of the columns only: index dqi
    size_t dqi(lexer *p, int a, int c, int m) const { return ((size_t)(a-p->imin)*p->jmax + (c-p->jmin))*p->knoz + m; }
    std::vector<double> dQxL, dQxR, dQyL, dQyR;                              // d Q / d eta (per layer)
    std::vector<double> pQxL0, pQxR0, pQxL1, pQxR1, pQyL0, pQyR0, pQyL1, pQyR1;  // d Q / d phi_k, phi_{k+1}
    std::vector<double> nActG;                                               // nAct with the columns of the neighbouring ranks
    bool nh_fixed = false;
    std::vector<std::vector<double>> fgV, fgZ;   // FGMRES of the coupled system: the Krylov vectors (kept, allocated as needed)
    std::vector<std::vector<double>> fgVe, fgZe, fgVp, fgZp;   // the same for the eta and phi systems (A 507 1)
    std::vector<double> nEtaP, nPhiP;            // eta and phi of the step before (the initial guess)
    std::vector<double> nsx, nsF0, nsy, nmb, nmx, nmw, nmv;     // work vectors of nh_solve and nh_solve_mg (kept)
    // zero: x = 0 (F(0); no phi faces)
    void nh_residual(lexer*, fdm_nhf*, ghostcell*, const std::vector<double>&, std::vector<double>&, bool zero=false, bool rows=true);
    void nh_solve(lexer*, fdm_nhf*, ghostcell*);
    void nh_solve_mg(lexer*, fdm_nhf*, ghostcell*, const std::vector<double>&, std::vector<double>&);
    void nh_cell(lexer*, fdm_nhf*, ghostcell*);
    void nh_cell_faces(lexer*, fdm_nhf*, double*, double*);

    // a cell side towards a wall (a solid or the outside of the domain, no open boundary face) counts in the
    // cell momentum with its zero face term (cell_update, nh_cell); wSide: the sides of the cells with an open boundary
    // face (bits 1 x-, 2 x+, 4 y-, 8 y+), from bcf (wall_sides)
    std::vector<unsigned char> wSide;
    void wall_sides(lexer*);
    double wall(lexer*, int, int, int);

    si_eta2D *psi;
    solver2D *psolv2D;
    int psi_iter;
    double psi_time;
    wind *pwind;
    met_forcing *pmet = nullptr;     // atmospheric forcing (W 200)
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
    slice4 eta0, WL0, etaI, gmc, rhsW, bdiag, tS, tE, tA, hoC, hoE;

    // water depth at the end of the step, the start of the next (see start); boundary copies
    slice4 WLk, etaB, WLb;

    // SEASTATE vortex force (A 751 2): the Stokes mass source -div(M) of the step [m/s], in rhsW and eta_update
    slice4 wsrc;

    // implicit horizontal diffusion: face coefficients dt nu_f/dn_f, the velocity of a layer and the right-hand side
    bool mass_source = false;

    // direct forcing of fixed solids (A 580 1 grid solids, A 581-590), as the explicit NHFLOW: the forcing of
    // nhflow_forcing on the cell predictor (phase_P1 of RK) and on the cells at the end of the step (phase_P2);
    // the layer face discharges are not forced (df_*, nhflow_momentum_SI.cpp)
    bool df_on() const;
    void df_predictor(lexer*, fdm_nhf*, ghostcell*);
    void df_reforce(lexer*, fdm_nhf*, ghostcell*);
    void df_faces(lexer*, fdm_nhf*, ghostcell*, bool);
    // the incremental force of the layer faces in a solid (private layout as Qx, Qy): added to the face predictor of
    // the next step, corrected by the discharge the solve leaves in the face (df_faces)
    double *DFx = nullptr, *DFy = nullptr;
    // columns inside a solid over the whole water depth (DF < 0 in every layer): walls of the implicit solve (their
    // faces closed in face_depth, wall sides of the gravity and phi face terms, hydrostatic, no fourth-order stencil)
    void df_columns(lexer*, fdm_nhf*, ghostcell*);
    bool df_col(lexer*, int, int);
    bool dfcols = false;
    void df_apply(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*, slice&, bool);
    std::vector<double> dfU, dfV, dfW, dfSU, dfSV, dfSW;
    slice4 dfWL, dfEta, dfC;

    // layer face discharges: prognostic Q^{n+1}, Q^n, layer transport of the step (theta weighted),
    // column parts Q~ = A^-1 (I - dt L_v) Q^ and r = A^-1 1
    // layout of the private arrays (see the top of the file): strides, offset, shift (test), number of entries
    int lsi = 0, lsj = 0, lso = 0, lshift = 0;
    size_t lN = 0;
    void layout_ini(lexer*, bool);
    bool lcomp = false;                 // compact layout
    int lmx = 0, lmy = 0, lnk = 0;      // compact: ghost cells in i, j, the entries of a column (knoz+1)
    int lki = 0, lkj = 0;               // compact: the cells in i, j incl. ghost cells
    int lrng[4] = {0,0,0,0};            // the columns of the layout: i in [lrng[0], lrng[1]), j in [lrng[2], lrng[3])
    size_t fN = 0;                      // entries of a full array (IJK)
    int fim = 0, fjm = 0, fkm = 0, fjx = 0, fkx = 0;    // the full layout: imin, jmin, kmin, jmax, kmax
    double *lscr[3] = {nullptr,nullptr,nullptr};   // compact: full copies of private arrays (lfull)
    void lgather(lexer*, const double *P, double *F) const;     // compact -> full (the cells of the layout)
    void lscatter(lexer*, const double *F, double *P) const;    // full -> compact
public:
    double* lfull(lexer*, double *P, int s);    // the array in the layout of IJK: P itself (full) or scratch s (compact)
    void lback(lexer*, double *P, int s);       // compact: scratch s (changed by IJK code) back to P
    void lset(lexer*, double *P, const double *F);          // P = F (F a full array)
    void lget(lexer*, const double *P, double *F);          // F = P (F a full array; compact: the cells of the layout)
private:
#ifdef SI_LCHECK
    int lchk(int a, int b, int c, const char *file, int line) const;
#endif
    double *Qx, *Qy, *Qxn, *Qyn, *Tx, *Ty, *QAx, *QAy, *Rx, *Ry;

    // layer state: t^n, RK stages, predictor, implicit diffusion
    double *UH0, *VH0, *WH0, *UHa, *VHa, *WHa, *UHb, *VHb, *WHb, *UHp, *VHp, *WHp;
    double *UHD, *VHD, *WHD;
    double *dUH, *dVH;

    vector<bc_face> bcf;
    int method[5];          // iowave open edge method per side cs (B 520 3-6), 0 none

    double theta, g;
    int rk, wd_iter, wd_cells, gcval_eta;
    bool onestep = false;
    bool tl2 = false;          // second-order time levels of the nonlinear terms (one Cartesian grid, no A 549 1)
    bool tp2 = false;          // the estimate pass for the advecting transport without tl2 (curvilinear, multi-block and
                               // refined grids): the passes of tl2 with the other time levels of the step before
    double *accU = nullptr, *accV = nullptr, *accW = nullptr;    // tl2: the implicit acceleration of the step before
    double *vaU = nullptr, *vaV = nullptr, *vaW = nullptr;       // A 510 11 with tl2: its vertical-advection part (col_va)
    double *FSx = nullptr, *FSy = nullptr;    // tl2: the sources of the last predictor stage (wind, Coriolis, ...; the wall sides)
    double thv = 1.0;                          // time level of the vertical advection in the column operator (1: implicit Euler)
    std::vector<double> vau, val, vad;         // its coefficients (unscaled), for the explicit part
    void col_adv_rhs(lexer*, const double*, double*);      // A 510 11: one stage, BCG advection, the nonlinear terms at t^{n+1/2}
    void col_va(lexer*, const double*, const double*, double, double*);   // the vertical advection of a cell column over the step
    vector<double> wv0, wv1, wva;              // its columns (knoz)
    bool qini;
    // the pass of step_single (tl2: 0 the estimate, 1 the final pass); faces_set: face_update of the step done by the
    // coupled non-hydrostatic solve (nh_solve), si_end_faces skips it
    int si_pass=0;
    // wave breaking (A 550 1/3/4): the breaking viscosity vb of the columns from the water level of the step, used in the
    // next step (horizontal diffusion of the predictor, column operator); nbrk_lim: its cap for the explicit diffusion
    nhflow_breaking *pbrk = nullptr;
    double vb_layer(lexer*, fdm_nhf*, int, int, int);
    bool faces_set=false;

    double Qin, Qout;       // discharge over the open boundaries of the step (theta weighted)
    double vrelax;          // volume change of the relaxation zones of the step
    ofstream silog;
};

#endif
