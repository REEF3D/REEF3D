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
#include<vector>

class fnpf_laplace;
class fnpf_fsf;
class fnpf_body;
class fnpf_amr;
class fnpf_laplace_cds2;

using namespace std;

// ---------------------------------------------------------------------------
// FNPF theta-Robin time integration (A 310 5): 2D and 3D, MPI (REEFMG); no wetting-drying, no bodies, no AMR.
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
// The Laplace rows are those of fnpf_laplace_cds2 (assembled, not solved). A 305 1 (default, REEFMG): the Robin
// row of each column is eliminated into its top Laplace row (Fs from the Robin row), which keeps the 7-point
// structure for the V-cycle; the couplings of the sz stencil beyond k-1 (5 and 7 points) are in REEFMG's Krylov
// operator only (matrix_diag::xcol). A 305 2: the surface nodes are added as unknowns and the system is solved
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
//     inexact (A 306, default 1e-2): each reduces its residual by A 306 relative to the warm start (BiCGStab,
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
//   Stokes waves from reflection); A 303 < 0 (default) relaxes per time as RK3 at N 47 0.75 (dt_crit of
//   fnpf_timestep).
//
// Edits of the surface after the solve (breaking damping and filters A 350/A 352, coastline): Fi below would be
// that of the unedited surface (as with the relaxation zones before patch 0006); when an edit changed eta or
// Fifsf, a Dirichlet Laplace solve on the final grid follows (as the last RK3 stage). Breaking events only.
//
// A 305 linear solver (1 REEFMG, default; 2 SparseLU), A 306 residual reduction of the REEFMG solves (default
// 1e-2; 0: N 44).
// A 300 theta (default 0.5: 0.55 damps 9 % per period at 21 steps), A 302 corrector passes (default 4),
// A 303 relaxation zones: > 0 factors r^(dt/A 303), 0 r per step, < 0 (default) per time as the three RK3 stages
// of an RK3 run at N 47 0.75; A 304 tolerance of the corrector (default 1e-3).
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
    void relax_affine(lexer*, ghostcell*, ioflow*, int, slice&, slice&);
    
    // A 305 1 (REEFMG): the Robin rows eliminated into the top Laplace rows, one solve per pass (N 44, or
    // inexact with A 306); returns the change of the solution
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
    slice4 fnew,wnew;      // Fifsf^n+1, w^n+1 of the solve
    slice4 en;
    slice4 colD,colA,colB; // REEFMG: Robin row per column, D Fs + A sum_m>=1 w_m Fi(knoz-m) = B
    slice4 fsold;          // REEFMG: Fifsf^n+1 of the last pass
    
    double theta;
    int lsolv;             // A 305: 1 REEFMG, 2 SparseLU
    double ctol;           // A 306: residual reduction of the inexact solves (REEFMG, 0: N 44)
    long cycsum0, cycsumc; // REEFMG cycles of pass 0 and of the corrector passes
    int npt;               // points of the sz stencil
    const double *wsz;     // its weights
    std::vector<double> mts;   // REEFMG: coupling of the top rows to the surface node
    std::vector<double> fiold; // REEFMG: Fi of the last pass
    int npass;
    double passtol;
    long passsum, stepsum, editsum;
    int nstep;
    double dt_prev;
    
    std::vector<double> mt0;   // surface coupling of the top Laplace row
    std::vector<double> rhsA;
    
    fnpf_laplace *plap;
    fnpf_laplace_cds2 *plap0;
    fnpf_fsf *pf;
    fnpf_body *pbody;
    
    double solvetime;
};

#endif
