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
// Wetting-drying (A 343 >= 1): the free-surface conditions of fnpf_fsfbc_wd as RK3. The coastline relaxation
// (A 341/A 342/A 347) and Fifsf = 0 on land are part of the affine maps of the zones (relax_affine), i.e. of the
// Robin rows and the eta update; the dry columns are not in the solve (A 343 1: Fifsf of the land map, A 343 2:
// held, set by wetdry_dynamic after the step); the clamps and the dynamic front after the step are surface edits
// (consistency solve). With a coastline N 51 3 uses at most the RK3 step of N 51 2 at N 47 0.75. The corrector
// passes precondition the change of the nonlinear terms with the advection (precond_adv) where the surface Courant
// number |Fx| dt/dx exceeds 0.5; without that they diverged at the swash at 2 x the RK3 step. With it,
// a static coastline runs at 2-3 x the RK3 step with N 51 2 and N 47 1.5-2.25 (Ting & Kirby: the breaking point as at
// the RK3 step); but in Mase & Kirby the Hm0 of the shallowest gauges grows with the step (no numerical damping),
// and the moving coastline (A 343 2) is unstable at 2 x at the wetted front, so N 51 3 keeps the RK3 step there
// (validation 09).
//
// Resolved floating bodies (X 10 1, fnpf_6DOF): the body moves with Heun (sixdof_obj::rk2): stage 0 with the loads of
// the state at t^n at the start of the step, stage 1 with the loads of the first corrector state (t^n+1, its tendencies
// kE1, kF1) in corrector pass 0; the later passes solve with the body at t^n+1. The footprint columns are Dirichlet
// nodes with the extended Fifsf, not Robin rows; the body geometry on the sigma grid of a pass before its assembly
// (solve_pre), the body band of Fi after the solve (solve_post); the psi solves of the loads with the plain Laplace
// solver. N 51 3 stays at the RK3 step: the natural periods of the body limit it. Not with externally driven bodies,
// SparseLU or without corrector passes (A 302 >= 1).
//
// Mesh refinement (G 1, fnpf_amr_theta.cpp): the patches take the same step (tendencies, predictor, corrector passes)
// and the Robin rows of all grids are solved together in the composite solve of fnpf_amr, with the couplings of the sz
// stencil beyond k-1 in its operator and lumped into k-1 in its V-cycles; after each update the state of all grids is
// made consistent (restriction, the cells around the patches). N 51 3 with the stability step of the finest cells. The
// patches with breaking A 350 1-3 (damping at the start of the step, detection and filters at its end, as level 0),
// without bodies, wetting-drying (G 30 0) or subcycling (refused).
//
// The linear part of the free-surface conditions, eta_t = w and Fifsf_t = -g eta, is integrated with the
// theta method, the nonlinear rest N (Zakharov terms, wind) at t^n+1/2:
//
//   eta^n+1   = eta^n   + dt [ theta w^n+1 + (1-theta) w^n + N_eta^n+1/2 ]
//   Fifsf^n+1 = Fifsf^n + dt [ -g (theta eta^n+1 + (1-theta) eta^n) + N_Fi^n+1/2 ]
//
// Eliminating eta^n+1 gives one Laplace problem per time step with the Robin condition at the free surface
//
//   Fifsf^n+1 + g theta^2 dt^2 w^n+1 = R,   R = Fifsf^n - g dt eta^n - g theta (1-theta) dt^2 w^n - g theta dt^2 N_eta + dt N_Fi
//
// with w = sigz sz(Fi) the free-surface vertical velocity of the sz stencil of A311 (A327 1). For theta >= 1/2
// the linear gravity waves are unconditionally stable; theta = 1/2 is Crank-Nicolson (no damping), theta > 1/2
// damps the grid-scale waves.
//
// The Laplace rows are those of fnpf_laplace_cds2 (assembled, not solved). A 306 1 (default, REEFMG): the Robin
// row of each column is eliminated into its top Laplace row (Fs from the Robin row), which keeps the 7-point
// structure for the V-cycle; the couplings of the sz stencil beyond k-1 (5 and 7 points) are in REEFMG's Krylov
// operator only (matrix_diag::xcol). A 306 2: the surface nodes are added as unknowns and the system is solved
// with a sparse direct solver (Eigen SparseLU, 2D serial, the prototype of steps 1-3).
//
// Time levels (test-plan step 2):
//   - predictor: the sigma grid of the solve is that of eta^n+1 predicted with Adams-Bashforth 2 (second order;
//     Euler in the first step);
//   - pass 0: N at t^n+1/2 extrapolated with Adams-Bashforth 2 from t^n and t^n-1, the sigma cross-derivative
//     terms of the Laplace rows with Fi^n (lagged by one step);
//   - A 302 corrector passes (at most, default 4; A 304: stop when the change of the solution is below A 304
//     times its size, default 1e-3): the sigma grid of eta^n+1 of the last solution, N at t^n+1/2 as the
//     mean of N^n and N^n+1 of the last solution (trapezoidal), and the cross-derivative terms with the last
//     solution. A pass assembles the system on that grid; SparseLU makes one defect-correction step with the
//     factorisation of pass 0, REEFMG solves it from the last solution (warm start). With REEFMG the solves are
//     inexact (A 307, default 1e-2): each reduces its residual by A 307 relative to the warm start (BiCGStab,
//     about one iteration). The passes re-solve the whole problem, so the solver error does not accumulate; the
//     convergence of the passes is set by the nonlinear coupling, not by the linear solves (5 November 2026).
//
// Relaxation zones (B 98 2, B 99 1/2), targets at t^n+1 (wavegen_precalc_fnpf again in start()): the relaxation of Fifsf, f -> r f + T, is part of the Robin rows,
//   Fifsf^n+1 + r g theta^2 dt^2 w^n+1 = r R + T,
// so the solution Fi is that of the relaxed surface potential (r = 0: Dirichlet, the target). Relaxing Fifsf
// after the solve (as the RK3 stages do before their solve) leaves Fi of the unrelaxed potential below the
// relaxed surface: w is wrong in the zones and the step is unstable (wave flume, 4 November 2026).
//
//   Without the corrector the lagged cross terms and the predicted grid feed the even harmonics (cos 2kx) of a
//   standing wave, which grow in time (first order in dt, proportional to the wave slope); one pass removes it.
//
//   One pass removes the growth of the even harmonics of step 2. At steps of 2-4 times the RK3 step one
//   defect-correction step does not converge: the unconverged corrector loses wave energy (JONSWAP Hm0 -16 %
//   at 2x with one pass, -2.5 % converged, A 303 0; step 3, 5 November 2026), hence the passes to convergence.
//
//   Zones once per theta step relax less per time than in the three RK3 stages (weaker beach, 6-8 % higher
//   Stokes waves from reflection); A 303 < 0 (default -0.1) relaxes per time as RK3 at N 47 -A 303 (dt_crit of
//   fnpf_timestep), and per step at least as one RK3 step: at smaller steps as an RK3 run with the same step.
//   The zone strength per time matters: in the 3D short-crested sea (step 5) the gauges after the generation
//   zone change by up to 16 % in Hm0 between zones per time as RK3 at N 47 0.75 and at N 47 0.04; Stokes
//   heights 1.035-1.044 with -0.75, 0.986-0.992 with -0.1 (the default from 6 to 7 November 2026).
//   With N 47 -A 303 the zones depended on the grid through dt_crit (cell size, depth, the velocity margin of the
//   criterion): -0.1 was an RK3 run with 160-180 steps per period in the 3D production sea (step 5) and in the
//   JONSWAP and Stokes benchmarks, 450 in Beji & Battjes, 740 in Ting & Kirby, and 710 and 1420 in a 2D flume at
//   dx 0.1 and 0.05 m; with a fixed step a refined patch did not change it (that flume with a patch at dx/2 was
//   4.6-4.9 % from the fine grid, 1.6-1.7 % with zones per period). A 303 <= -1 (default -200, patch 0032): per
//   time as an RK3 run with -A 303 steps per period of the wave input (wT, wTp for spectra, the peak of the
//   components for reconstructed waves B 92 51-53 since patch 0033, fnpf_wave_period.h), independent of the
//   grid; the 3D sea stays as before (Hm0 bias -2.6 %), the 2D benchmarks are the same or better (validation 13).
//   Without a wave period as RK3 at N 47 0.1, as before. Mase & Kirby (B 92 51, Tp 0.988 s) and Tian tp1 (0.937 s)
//   with zones per period: checker errors 0.281 and 0.0222 (before 0.276 and 0.0129); Mase & Kirby within its
//   run-to-run spread (0.28-0.30), Tian with an energy loss of 0.111 for 0.120 (reference 0.133; validation 14).
//
// Damping with the viscosity vb of breaking and of the coastline (A 350, A 357): at the start of the step, as a
// source of the step (patch 0022), so the solve gives Fi of the damped surface.
// Edits of the surface after the solve (the breaking filters A 352, the wetting-drying clamps, the coastline of
// A 343 0): Fi below would be that of the unedited surface (as with the relaxation zones before patch 0006); when
// an edit changed eta or Fifsf, a Dirichlet Laplace solve on the final grid follows (as the last RK3 stage).
//
// A 306 linear solver (1 REEFMG, default; 2 SparseLU), A 307 residual reduction of the REEFMG solves (default
// 1e-2; 0: N 44).
//
// A 308 c > 0 (theta 0.5, N 51 2): 4th-order phase. The step above is the Newmark scheme (gamma = theta) with
// beta = theta^2 = 1/4 in the Fifsf update,
//   Fifsf^n+1 + g beta dt^2 w^n+1 = Fifsf^n - g dt eta^n - g (theta - beta) dt^2 w^n - g theta dt^2 N_eta + dt N_Fi,
// whose relative frequency error is (beta - 1/12) (omega dt)^2/2: (omega dt)^2/12 for Crank-Nicolson, of 4th
// order for beta = 1/12 (Fox-Goodwin). Below 1/4 the scheme is stable for omega dt < 1/sqrt(1/4 - beta);
// beta = max(1/12, 1/4 - 1/(c omega_max dt)^2) per step, omega_max dt = sqrt(3) dt/dt_crit the grid cut-off
// mode of N 51 2 (Doppler shift included). The eta update and the cost are unchanged.
//   Near the limit the discrete frequencies of the shortest resolved waves approach pi/dt (Crank-Nicolson
//   compresses them below the exact ones), so the gain is in the waves well below the cut-off: at twice the RK3
//   step Beji & Battjes 0.233 instead of 0.351 (rel. rms), Stokes 0.03 instead of 1.26 (A 308 1); JONSWAP Hm0
//   -5.5 % instead of -4.6 %; at four times the RK3 step no gain. Unstable at 4x with c = 0.6, stable with
//   c = 0.8: the cut-off estimate of N 51 2 is conservative by about 1.3-1.6 (6 November 2026).
// A 305 theta (default 0.5: 0.55 damps 9 % per period at 21 steps), A 302 corrector passes (default 4),
// A 303 relaxation zones: > 0 factors r^(dt/A 303), 0 r per step, -1 < A 303 < 0 per time as the three RK3 stages
// of an RK3 run at N 47 -A 303, <= -1 (default -200) as an RK3 run with -A 303 steps per period of the wave input;
// A 333 the coastline relaxation (A 343 >= 1: > 0 as RK3 at N 47 A 333, 0 as the zones); A 304 tolerance of the
// corrector (default 1e-3); A 379 delta d: theta + delta within d [m] of the still-water coastline (default off).
//
// A 379 (patch 0034): in the inner surf and swash zone the theta step at the RK3 step (N 51 3 with a coastline)
// carries more short-wave energy than the converged solution: Mase & Kirby, Hm0 above 0.3 Hz at h = 2.5 cm 9 % above
// RK3 at N 47 0.25, 4 % at half the step, 2 % at a quarter (passes, filters, breaking thresholds and zone strength
// do not change it). Off-centring theta there damps it: A 379 0.1 2.0 puts the inner gauges on the converged RK3
// ones (checker 0.220 for 0.28, RK3 0.236 converged, 0.217 at N 47 0.75 with its own damping). It is a damping of
// order (theta - 1/2) (omega dt)^2 per step and so depends on d, delta and dt: d = 3 m over-damps (0.176). A
// band by depth instead of by distance damps the harmonics over a bar (Beji & Battjes 0.248 -> 0.280).
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
