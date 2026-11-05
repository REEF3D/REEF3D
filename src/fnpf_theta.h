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
// FNPF theta-Robin time integration (A 310 5), prototype: 2D, serial, no wetting-drying, no bodies, no AMR.
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
// The Laplace rows are those of fnpf_laplace_cds2 (assembled, not solved); the surface nodes are added as
// unknowns with the Robin row and the system is solved with a sparse direct solver (Eigen SparseLU). For the
// production version the Robin row goes into REEFMG.
//
// Time levels (test-plan step 2):
//   - predictor: the sigma grid of the solve is that of eta^n+1 predicted with Adams-Bashforth 2 (second order;
//     Euler in the first step);
//   - pass 0: N at t^n+1/2 extrapolated with Adams-Bashforth 2 from t^n and t^n-1, the sigma cross-derivative
//     terms of the Laplace rows with Fi^n (lagged by one step);
//   - A 302 corrector passes (at most, default 4; A 304: stop when the change of the solution is below A 304
//     times its size, default 1e-3): the sigma grid of eta^n+1 of the last solution, N at t^n+1/2 as the
//     mean of N^n and N^n+1 of the last solution (trapezoidal), and the cross-derivative terms with the last
//     solution. A pass assembles the system on that grid and makes one defect-correction step with the
//     factorisation of pass 0 (no new factorisation).
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
//   at 2x with one pass, -2.5 % converged; step 3, 5 November 2026), hence the passes to convergence.
//
// A 300 theta (default 0.5: 0.55 damps 9 % per period at 21 steps), A 302 corrector passes (default 4),
// A 303 reference step of the relaxation zones (factors r^(dt/A 303), default 0: r per step), A 304 tolerance
// of the corrector (default 1e-3).
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
    
    double theta;
    int npass;
    double passtol;
    long passsum, stepsum;
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
