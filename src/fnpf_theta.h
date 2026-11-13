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
Architect: Hans Bihs
--------------------------------------------------------------------*/


#ifndef FNPF_THETA_H_
#define FNPF_THETA_H_

#include"fnpf.h"
#include"fnpf_ini.h"
#include"fnpf_sigma.h"
#include"slice4.h"
#include<memory>
#include<vector>

class fnpf_laplace;
class fnpf_fsf;
class fnpf_body;
class fnpf_amr;
class fnpf_laplace_cds2;
class solver2D;
class matrix2D;
class vec2D;

using namespace std;

// ---------------------------------------------------------------------------
// FNPF theta-Robin time integration (A 310 5): 2D and 3D, MPI (REEFMG), wetting-drying A 343 1 and 2, resolved
// floating bodies (X 10 1, fnpf_6DOF), mesh refinement (G 1, fnpf_amr_theta.cpp).
//
// Method. The linear part of the free-surface conditions, eta_t = w and Fifsf_t = -g eta, is integrated with the
// theta method, the nonlinear rest N (Zakharov terms, wind) at t^n+1/2:
//
//   eta^n+1   = eta^n   + dt [ theta w^n+1 + (1-theta) w^n + N_eta^n+1/2 ]
//   Fifsf^n+1 = Fifsf^n + dt [ -g (theta eta^n+1 + (1-theta) eta^n) + N_Fi^n+1/2 ]
//
// Eliminating eta^n+1 gives one Laplace problem per time step with the Robin condition at the free surface
//
//   Fifsf^n+1 + g theta^2 dt^2 w^n+1 = R,   R = Fifsf^n - g dt eta^n - g theta (1-theta) dt^2 w^n - g theta dt^2 N_eta + dt N_Fi
//
// with w = sigz sz(Fi) the free-surface vertical velocity of the sz stencil of A 311 (Fornberg weights). For
// theta >= 1/2 the linear gravity waves are unconditionally stable; theta = 1/2 is Crank-Nicolson (no damping),
// theta > 1/2 damps the grid-scale waves.
//
// Linear solver (A 306). The Laplace rows are those of fnpf_laplace_cds2 (assembled, not solved).
//   1 (default, REEFMG): the Robin row of each column is eliminated into its top Laplace row, which keeps the
//     7-point structure for the V-cycle; the couplings of the sz stencil beyond k-1 are in the Krylov operator
//     only (matrix_diag::xcol). A 307: inexact solves, residual reduction per pass (default 1e-2; 0: N 44).
//   2 (SparseLU, 2D serial): the surface nodes are extra unknowns, direct solve.
//
// Time levels:
//   - predictor: the sigma grid of the solve is that of eta^n+1 predicted with Adams-Bashforth 2 (Euler in the
//     first step);
//   - pass 0: N at t^n+1/2 extrapolated with Adams-Bashforth 2, the sigma cross-derivative terms of the Laplace
//     rows with Fi^n;
//   - corrector passes (A 302, default 4; A 304: stop when the change of the solution is below A 304 times its
//     size, default 1e-3): the sigma grid of eta^n+1 of the last solution, N^n+1/2 as the mean of N^n and N^n+1
//     (trapezoidal), the cross terms with the last solution. Each pass re-solves the whole problem (SparseLU: one
//     defect-correction step with the factorisation of pass 0; REEFMG: warm start), so the solver error does not
//     accumulate. Where the surface Courant number |Fx| dt/dx exceeds 0.5 the change of the nonlinear terms is
//     preconditioned with the advection (precond_adv).
//
// Relaxation zones (B 98 2, B 99 1/2), targets at t^n+1: the relaxation of Fifsf, f -> r f + T, is part of the
// Robin rows,
//   Fifsf^n+1 + r g theta^2 dt^2 w^n+1 = r R + T,
// so the solution Fi is that of the relaxed surface potential (r = 0: Dirichlet, the target). A 303 sets the
// strength per time: > 0 factors r^(dt/A 303), 0 r per step, -1 < A 303 < 0 per time as the three RK3 stages of
// an RK3 run at N 47 -A 303, <= -1 (default -200) as an RK3 run with -A 303 steps per period of the wave input
// (fnpf_wave_period.h); without a wave period as RK3 at N 47 0.1.
//
// Wave generation and absorption at the edges: B 98 3/4 (and 6-8) and B 99 3-5, 7 (and 6, 8) through
// iowave::inflow_fnpf as in RK3, with the targets at t^n+1, in the predictor and in each corrector pass.
//
// Wetting-drying (A 343 >= 1): the free-surface conditions of fnpf_fsfbc_wd as RK3. The coastline relaxation
// (A 341/A 342/A 347; strength A 333) and Fifsf = 0 on land are part of the affine maps of the zones
// (relax_affine); the dry columns are not in the solve (A 343 1: Fifsf of the land map, A 343 2: held, set by
// wetdry_dynamic after the step). With a coastline N 48 1 uses at most the RK3 step of N 48 2 at N 47 0.75.
// A 379 delta d: theta + delta within d [m] of the still-water coastline, damps the short waves in the swash.
//
// Resolved floating bodies (X 10 1, fnpf_6DOF): the body moves with Heun (sixdof_obj::rk2): stage 0 with the loads
// of the state at t^n, stage 1 with the loads of the first corrector state (kE1, kF1) in corrector pass 0; the later
// passes solve with the body at t^n+1. The footprint columns are Dirichlet nodes with the extended Fifsf; the body
// geometry on the sigma grid of a pass before its assembly (solve_pre), the body band of Fi after the solve
// (solve_post). N 48 1 stays at the RK3 step. Not with externally driven bodies, SparseLU or without corrector
// passes.
//
// Mesh refinement (G 1, fnpf_amr_theta.cpp): the patches take the same step and the Robin rows of all grids are
// solved together in the composite solve of fnpf_amr; after each update the state of all grids is made
// consistent. Breaking A 350 1-3 on the patches; no bodies, wetting-drying (G 30 0) or subcycling.
//
// Damping with the viscosity vb of breaking and of the coastline (A 350, A 357): at the start of the step, as a
// source of the step, so the solve gives Fi of the damped surface. Edits of the surface after the solve (the
// breaking filters A 352, the wetting-drying clamps, the coastline of A 343 0): when an edit changed eta or Fifsf,
// a Dirichlet Laplace solve on the final grid follows (as the last RK3 stage).
//
// A 308 c > 0 (theta 0.5, N 48 2): 4th-order phase. The step is the Newmark scheme (gamma = theta) with
// beta = theta^2 = 1/4 in the Fifsf update,
//   Fifsf^n+1 + g beta dt^2 w^n+1 = Fifsf^n - g dt eta^n - g (theta - beta) dt^2 w^n - g theta dt^2 N_eta + dt N_Fi,
// whose relative frequency error is (beta - 1/12) (omega dt)^2/2, of 4th order for beta = 1/12 (Fox-Goodwin).
// Below 1/4 the scheme is stable for omega dt < 1/sqrt(1/4 - beta); beta = max(1/12, 1/4 - 1/(c omega_max dt)^2)
// per step, omega_max dt = sqrt(3) dt/dt_crit the grid cut-off mode of N 48 2.
//
// Parameters: A 302 corrector passes, A 303 zone strength, A 304 corrector tolerance, A 305 theta (default 0.5),
// A 306 linear solver, A 307 inexact solves, A 308 Newmark beta, A 333 coastline relaxation, A 379 theta near the
// coastline, A 406 n: a failed step is repeated with half the step, at most n times (default 3, 0 off).
// ---------------------------------------------------------------------------

class fnpf_theta final : public fnpf_ini, public fnpf_sigma
{
public:
	fnpf_theta(lexer*, fdm_fnpf*, ghostcell*);
	virtual ~fnpf_theta();
    
    void start(lexer*, fdm_fnpf*, ghostcell*, solver*, convection*, ioflow*, reini*) override final;
    void inidisc_step1(lexer*, fdm_fnpf*, ghostcell*, ioflow*, solver*) override final;
    void inidisc_step2(lexer*, fdm_fnpf*, ghostcell*, ioflow*, solver*) override final;
    void ini_wetdry(lexer*, fdm_fnpf*, ghostcell*) override final;
    fnpf_body* body() override final {return pbody;}
    void attach_amr(fnpf_amr*) override final;
    
private:
    // one theta step; false when a Laplace solve of it failed (diverged, NaN, iteration limit N 46), the largest
    // surface velocity jumped (more than 3 times that at t^n plus 1 m/s) or the state is not finite (A 406: the
    // step is then repeated from t^n with half the step)
    bool step(lexer*, fdm_fnpf*, ghostcell*, solver*, convection*, ioflow*, reini*);
    void solve_check(lexer*);
    void state_save(lexer*, fdm_fnpf*);
    void state_restore(lexer*, fdm_fnpf*, ghostcell*, ioflow*);
    bool step_fail = false;
    int fail_why = 0;              // 1 Laplace solve, 2 surface velocity jump, 3 state not finite
    int retry_steps = 0;           // steps since the last rejected step (the step may grow by 2 per step then)
    long retrysum = 0;
    std::vector<double> sv_fi;     // the state at t^n of a step that may be repeated (A 406)
    slice4 sv_eta, sv_fifsf, sv_vb, sv_neo, sv_nfo, sv_eko, sv_fko;
    double sv_dtprev = 0.0;
    int sv_nstep = 0;

    // the system on the current sigma grid: Laplace rows (cross terms of the current Fi and fs) and Robin rows (R)
    void robin_build(lexer*, fdm_fnpf*, ghostcell*, solver*, slice&, slice&, slice&);
    // pass 0: build, factorise, solve
    void robin_solve(lexer*, fdm_fnpf*, ghostcell*, solver*, slice&, slice&, slice&);
    // corrector pass: build on the current grid, one defect-correction step with the factorisation of pass 0
    double robin_correct(lexer*, fdm_fnpf*, ghostcell*, solver*, slice&, slice&);
    // the solution back to Fi (interior), fout (Fifsf^n+1) and wnew (w^n+1)
    void robin_unpack(lexer*, fdm_fnpf*, ghostcell*, slice&);
    // eta^n+1 of the current solution, N^n+1/2 = (N^n + N^n+1)/2, R, and the sigma grid of eta^n+1
    void corrector_state(lexer*, fdm_fnpf*, ghostcell*, ioflow*);
    // the relaxation of the zones as f -> R f + T (from f = 0 and f = 1), with the factors r^(dt/A 303)
    void relax_affine(lexer*, fdm_fnpf*, ghostcell*, ioflow*, int, slice&, slice&);
    void relax_scale(lexer*, double, slice&, slice&, bool);   // the strength per time of a relaxation (A 303, A 333)
    
    // A 306 1 (REEFMG): the Robin rows eliminated into the top Laplace rows, one solve per pass (N 44, or
    // inexact with A 307); returns the change of the solution
    double robin_solve_mg(lexer*, fdm_fnpf*, ghostcell*, solver*, slice&, slice&, slice&, int);
    void robin_build_mg(lexer*, fdm_fnpf*, ghostcell*, solver*, slice&, slice&);
    void sz_stencil(lexer*, ghostcell*);
    
    int gcval,gcval_u,gcval_v,gcval_w;
    int gcval_eta,gcval_fifsf;
    
    slice4 ek,fk;          // tendencies of the full free-surface conditions at t^n
    slice4 Ne,Nf,Neo,Nfo;  // nonlinear parts at t^n and t^n-1
    slice4 nEm,nFm;        // nonlinear parts at t^n+1/2 of the step
    slice4 eko,fko;        // full tendencies at t^n-1 (predictor)
    slice4 wn;             // w^n
    slice4 erk1,frk1;      // predicted eta and Fifsf (sigma grid of the solve, cross terms of pass 0)
    slice4 erk2;           // eta^n+1 of the current solution (corrector)
    slice4 R0,Rs;          // right-hand side of the Robin row without and with the nonlinear terms
    slice4 relR,relT;      // relaxation of Fifsf in the zones, f -> relR f + relT (relR 1, relT 0 outside)
    slice4 relRe,relTe;    // relaxation of eta in the zones
    slice4 thl;            // A 379: theta of the column (theta + A 379 delta within A 379 d of the coastline)
    slice4 fnew,wnew;      // Fifsf^n+1, w^n+1 of the solve
    slice4 en;
    slice4 colD,colA,colB; // REEFMG: Robin row per column, D Fs + A sum_m>=1 w_m Fi(knoz-m) = B
    slice4 fsold;          // REEFMG: Fifsf^n+1 of the last pass
    slice4 nEo,nFo,pxe,pxf;    // corrector: N of the last pass and the preconditioned changes of N
    std::unique_ptr<solver2D> psolv2d;    // the 2D solve of precond_adv
    std::unique_ptr<matrix2D> Nadv;
    std::unique_ptr<vec2D> xadv,radv;
    void precond_adv(lexer*, fdm_fnpf*, ghostcell*, slice&);
    
    double theta;
    bool thloc = false;    // A 379 on: theta per column (thl), set at the first step
    bool thset = false;
    double beta,betasum;   // Newmark beta of the Fifsf update (theta^2; A 308: 1/12 .. 1/4)
    int lsolv;             // A 306: 1 REEFMG, 2 SparseLU
    double ctol;           // A 307: residual reduction of the inexact solves (REEFMG, 0: N 44)
    long cycsum0, cycsumc; // REEFMG cycles of pass 0 and of the corrector passes
    int npt;               // points of the sz stencil
    const double *wsz;     // its weights
    std::vector<double> mts;   // REEFMG: coupling of the top rows to the surface node
    std::vector<double> fiold; // REEFMG: Fi of the last pass
    int npass;
    double passtol;
    long passsum, stepsum, editsum, padvsum;   // padvsum: passes with the advection preconditioner
    double csmaxsum;                          // largest surface Courant number of the run
    int nstep;
    double dt_prev;
    
    std::vector<double> mt0;   // surface coupling of the top Laplace row
    std::vector<double> rhsA;
    
    fnpf_laplace *plap;
    fnpf_laplace_cds2 *plap0;
    fnpf_fsf *pf;
    fnpf_body *pbody;
    solver *psolv_step=nullptr;   // the Laplace solver of the step (the body loads of the corrector)
    int passq=-1;                 // the corrector pass under way (-1: pass 0)
    slice4 kE1,kF1;               // tendencies of eta and Fifsf of the first corrector state (the body loads at t^n+1)
    slice4 cR,cT;                 // the coastline relaxation as f -> cR f + cT (relax_affine)
    fnpf_amr *pamr=nullptr;       // mesh refinement (G 1): the patches (fnpf_amr_theta.cpp)
    
    double solvetime;
};

#endif
