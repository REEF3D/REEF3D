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

#ifndef MOMENTUM_RK_H_
#define MOMENTUM_RK_H_

#include"momentum.h"
#include"momentum_forcing.h"
#include"bcmom.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include<vector>

class convection;
class diffusion;
class pressure;
class turbulence;
class solver;
class poisson;
class fluid_update;
class reini;
class picard;
class heat;
class concentration;
class density;
class sixdof;
class fsi;

using namespace std;

/*--------------------------------------------------------------------
momentum_rk: Runge-Kutta time stepping of the CFD momentum equations
(primitive variables), one class for the schemes

  N40   time integration                         level set
  12    SSP-RK2 (Heun)                           outside (free surface module)
  13    SSP-RK3 (Shu-Osher)                      outside
  44    low-storage RK3 (Spalart et al.)         outside
  2     SSP-RK2                                  inside every stage (fully coupled)
  3     SSP-RK3                                  inside every stage (default)
  4     low-storage RK3                          inside every stage
  33    SSP-RK3, conservative form (rho u)       inside every stage
  5     one step, Godunov convection (godunov),  after the projection, with u(n+1)
        one projection
  6     one step with two projections: the      with u(n+1/2), between the
        half step u(n+1/2), p(n+1/2) projected,  two projections
        then the full step with them

N40 6 (step_twoproj), the predictor-corrector of Bell, Colella and Glaz on the staggered grid:
  1. u* = u(n) + dt/2 [C(u(n); u(n)) + S - grad p(n)/rho(n) + g] + dt/2 D(u*)   (D: implicit)
     projected with dt/2 and the face density of phi(n): u(n+1/2), p(n+1/2)
  2. phi(n+1) = phi(n) + dt L(phi(n); u(n+1/2)), reinitialised; density and viscosity of
     phi(n+1/2) = (phi(n) + phi(n+1))/2
  3. u* = u(n) + dt [C(u(n); u(n+1/2)) + S - grad p(n+1/2)/rho(n+1/2) + g]
         + dt/2 [D(u(n)) + D(u*)]                                    (Crank-Nicolson)
     projected with dt and the face density of phi(n+1/2): u(n+1), p(n+1)
  4. density and viscosity of phi(n+1)
C(f; U): the Godunov convection of f (face states at n+1/2) with the advecting velocity U. The
forces of the half step are projected, so that the pressure gradient over the face density
is balanced at the interface (the half-step forces in the Godunov face states of N 40 5,
without a projection, made air jets in the dam break).
A prototype (validation 75): second order in time on smooth interface flows, but in violent flows
(the 2D dam break) stable up to N 47 0.3, not at 0.45 (N 56 2 delays the growth, it does not stop
it): the half-step velocity carries the air accelerations of the half step and convects the air
momentum with them. N 40 5 is the robust one-step scheme.

On the mesh refinement (G 7 0, cfd_amr::step_comp) N 40 6 takes the same steps on all grids, with
the composite projection after each of the two momentum steps and the level set of all levels
between them.

With the compressible air (W 32, pjm_corr_comp.cpp) both projections take p(n) of the start of the
step as the reference of the pressure change: the half step over dt/2, the full step over dt.

N 59 1 (N 40 5, 6 with the level set, also on the mesh refinement, whose patches take the
same Godunov convection): the convection of the momentum as rho u, moved with the mass flux of the
Godunov convection (godunov::start, see godunov.h); the density of the step starts from the level
set of the step (phi(n) also in the full step of N 40 6, kept in frk1). Validation 77:
the impact pressures of the 3D dam break as RK3 without the oscillations of the one-step schemes,
N 40 6 stable at N 47 0.45 in the 2D dam break.

SSP stage s:   u(s+1) = a_s u(n) + b_s [ u(s) + dt (F + D) ]
low storage:   u(s+1) = u(s) + 2 alpha_s dt (S + D) + gamma_s dt C(s) + zeta_s dt C(s-1)
with F the explicit terms (sources, pressure gradient, convection C), D the
diffusion (implicit for D 20 2/3: solved with the weight b_s resp. 2 alpha_s
after the explicit terms are combined), followed by the direct forcing and the
pressure projection of each stage.

Implicit diffusion (D 20 2), D 23 2 (default): second order in time. The
stages solve with their own implicit weight and add explicit corrections
dt sum_j d_sj D(u_j) with the diffusion D(u_j) of earlier stages, evaluated at the projected stage velocities with the operator of the
implicit scheme (diffusion::apply_u/v/w), so that the implicit part is a second-order, L-stable DIRK with the same stage times as
the explicit part (IMEX order conditions):
  SSP-RK3:   rows (u1, u2, u3) = [1], [1/4, 1/4], [-1/2, 1, 1/2]
  SSP-RK2:   rows (u0, u1, u2) = [0, 1], [1/2, -1/2, 1]   (D(u0) from the last step)
  low-storage RK3: stage s weights beta_s D(u_s) + (2 alpha_s - beta_s) D(u_s+1),
             beta = (4/15, 0, 29/150)                    (D(u0) from the last step)
D 23 1: all implicit weight on the new stage value (first order, as before).
SSP-RK2 and low storage use D(u(n)) of the last step and start with the
first-order weights.

Conservative form (N40 33): the convection is advanced for u, rho u and
the face density rho in each stage, and the convected velocity is
reconstructed from rho u / rho (vel_limiter) before the sources and the
diffusion are added.
--------------------------------------------------------------------*/

class pic;

class momentum_rk final : public momentum, public momentum_forcing, public bcmom
{
public:
	momentum_rk(lexer*, fdm*, ghostcell*, convection*, convection*, diffusion*, pressure*, poisson*,
                turbulence*, solver*, solver*, ioflow*, heat*&, concentration*&, reini*, fsi*);
	virtual ~momentum_rk();
	void start(lexer*, fdm*, ghostcell*, vrans*, sixdof*) override final;

    // PIC (F 81 2, pic_f): the convection is left out on the faces in the water
    static pic *ppic;

    // mesh refinement (cfd_amr, G 15): the SSP step in parts, called grid by grid with the fills of
    // the patches and the composite projection in between (step_ssp calls them in this order)
    int amr_stages() const {return stages;}
    bool amr_ssp() const {return scheme==SSP && !conservative;}
    // N 40 6: the weight of the projection of the half step (s 0: dt/2) and of the full step (s 1: dt)
    double amr_alpha(int s) const {return twoproj ? (s==0 ? 0.5 : 1.0) : ssp_b[s];}
    // N 40 6 (step_twoproj) in parts, called grid by grid: the momentum of the half step (half true:
    // u(n) -> u* in urk2) or of the full step (half false: u(n) -> u* in a->u with the advecting
    // velocity u(n+1/2)) with the forcing (1: all, 2: solids and topography only, 0: none), and the
    // level set, density and viscosity of the middle of the step. The stage index of the accessors
    // (amr_vel, amr_velout, amr_alpha) is 0 for the half step, 1 for the full step; the level set
    // of both is ls (amr_phi_in, amr_phi_out)
    bool amr_twoproj() const {return twoproj;}
    void amr_twoproj_momentum(lexer*, fdm*, ghostcell*, vrans*, sixdof*, bool, int);
    void amr_twoproj_mid(lexer*, fdm*, ghostcell*);
    void amr_step_begin(lexer*, fdm*, ghostcell*);
    void amr_ls_transport(lexer*, fdm*, ghostcell*, int);
    void amr_ls_finish(lexer*, fdm*, ghostcell*, int, bool, int);
    void amr_ls_reini(lexer*, fdm*, ghostcell*, int, int);
    int amr_reini_iters(lexer*, int) const;
    void amr_momentum(lexer*, fdm*, ghostcell*, vrans*, sixdof*, int, int);
    // forcing 3 of amr_momentum: the predictor alone, the forcing of the bodies later (X 198: the
    // fluid-momentum coupling fits the bodies to the predicted velocity of all patches first)
    void amr_forcing(lexer*, fdm*, ghostcell*, sixdof*, int);
    void amr_project_after(lexer*, fdm*, ghostcell*, field&, field&, field&);
    void amr_projection(lexer*, fdm*, ghostcell*, int);
    void amr_solid_forcing(lexer*, fdm*, ghostcell*, field&, field&, field&, double);
    void amr_stage_end(lexer*, fdm*, ghostcell*, int);
    field& amr_vel(fdm*, int c, int s);
    field& amr_velout(fdm*, int c, int s);
    field4& amr_phi_in(int s);
    field4& amr_phi_out(int s);

    // N 63 2, 3: the level set of the one-step schemes with the scheme of F 35 in SSP-RK2 or SSP-RK3
    // stages with the velocity of the step (U, V, W fixed); the mesh refinement calls the stages grid
    // by grid with the fills of the patches in between (cfd_amr), level 0 alone takes them in
    // amr_ls_transport.  ls_rk_n: the number of stages (1: one forward-Euler step, N 63 0, 1)
    int ls_rk_n() const {return lsrk_n;}
    void ls_rk_begin(lexer*, fdm*);
    void ls_rk_stage(lexer*, fdm*, int st, field&, field&, field&);
    void ls_rk_ghost(lexer*, ghostcell*);
    field& ls_vel(fdm*, int c);

    // the diffusion and its solver (the mesh refinement: the implicit diffusion of heat and
    // concentration on the grid of this object)
    diffusion* amr_diffusion() { return pdiff; }
    solver* amr_solver() { return psolv; }

private:
    enum {SSP=0, LOWSTORAGE=1};

    void step_ssp(lexer*, fdm*, ghostcell*, vrans*, sixdof*);
    void step_twoproj(lexer*, fdm*, ghostcell*, vrans*, sixdof*);
    void component_onestep(lexer*, fdm*, ghostcell*, vrans*, int, field&, field&, field&, field&, field&, double, bool);
    void twoproj_forcing(lexer*, fdm*, ghostcell*, sixdof*, field&, field&, field&, double, bool, int);

    void step_lowstorage(lexer*, fdm*, ghostcell*, vrans*, sixdof*);

    void levelset_transport(lexer*, fdm*, ghostcell*, int);
    void levelset_finish(lexer*, fdm*, ghostcell*, int, bool, int);
    void levelset_lowstorage(lexer*, fdm*, ghostcell*, int);

    void component_ssp(lexer*, fdm*, ghostcell*, vrans*, int, int);
    void convection_conservative(lexer*, fdm*, ghostcell*, int);
    void component_lowstorage(lexer*, fdm*, ghostcell*, vrans*, int, int);

    void sources(lexer*, fdm*, ghostcell*, vrans*, int, field&);
    void rhs(lexer*, fdm*, int);
    void convection_start(lexer*, fdm*, int, field&, field&, field&, field&);
    void diffusion_start(lexer*, fdm*, ghostcell*, int, field&, field&, field&, field&, double);
    void projection(lexer*, fdm*, ghostcell*, field&, field&, field&, double);

    // second-order implicit diffusion
    void imex_step_setup(lexer*);
    void imex_pre(lexer*, fdm*, int, int, field&);
    void imex_accumulate(lexer*, fdm*, ghostcell*, int, field&, field&, field&);
    field& DC(int c);

    // stage fields
    field& vel(fdm*, int c, int s);         // input velocity component c of SSP stage s
    field& velout(fdm*, int c, int s);      // output velocity component c of SSP stage s
    field& comp(fdm*, field&, field&, field&, int c);
    field& FGH(fdm*, int c);
    void clear_FGH(lexer*, fdm*);

    // conservative form
    void face_density(lexer*, fdm*, ghostcell*, field&, field&, field&);
    double vel_limiter(lexer*, fdm*, field&, field&, field&, field&);
    double ro_filter(lexer*, fdm*, field&);
    field& M(int c, int idx);
    field& RO(int c, int idx);
    field& UR(int c);

    int scheme;
    int stages;
    bool twoproj;       // N 40 6: one step with two projections (step_twoproj)
    bool levelset;
    bool conservative;
    std::vector<double> csave;   // the right-hand side before the convection (PIC, F 81 2)

    // SSP: u(s+1) = ssp_a[s] u(n) + ssp_b[s] (u(s) + dt F)
    double ssp_a[3], ssp_b[3];

    // low storage
    double ls_alpha[3], ls_gamma[3], ls_zeta[3];

    // SSP: stage velocities; low storage: urk1 = stage velocity, urk2 = convection of the previous stage
    field1 urk1, urk2, fx;
	field2 vrk1, vrk2, fy;
	field3 wrk1, wrk2, fz;

    // level set (N40 = 2, 3, 4): ls = phi(n), frk1/frk2 = stages; low storage: frk1 = convection of the previous stage;
    // N 40 5 with N 59: frk1 = the level set before its last transport (phiold)
    field4 *ls, *frk1, *frk2;

    // N 63 2, 3: phi(n) for the stages of the level set of the one-step schemes
    field4 *lsrk = nullptr;
    int lsrk_n = 1;
    void ls_rk_all(lexer*, fdm*, ghostcell*, field&, field&, field&);

    // N 40 5 with N 59: the mass transport of the momentum starts from the level set before its last
    // transport (frk1), kept by amr_ls_transport; phiold_valid: kept in the last step of this grid (1),
    // a double for the checkpoints (ckpt_registry)
    bool phiold;
    double phiold_valid;
    unsigned long long ckid_phiold;

    // second-order implicit diffusion: implicit weight, accumulation weight of the solved D and
    // whether the accumulated correction is added, per stage
    bool imex;
    double imex_valid;          // D(u(n)) of the last step available (1); a double for the checkpoints (ckpt_registry)
    unsigned long long ckid_imex;
    double imex_g[3], imex_w[3];
    bool imex_use[3];
    field1 *Dtmp;
    field1 *Du;
    field2 *Dv;
    field3 *Dw;

    // conservative form (N40 = 33): reconstructed velocity, rho u and face density of the stages
    // (index 0: level n, 1, 2: stages; the last stage writes index 0)
    field1 *ur, *Mx[3], *rox[3];
    field2 *vr, *My[3], *roy[3];
    field3 *wr, *Mz[3], *roz[3];
    density *pd;
    double ro_threshold;
    double val;
    int gcval_ro;

    fluid_update *pupdate;
    picard *ppicard;

	int gcval_u, gcval_v, gcval_w;
    int gcval_phi;
	double starttime;

	convection *pconvec;
    convection *pfsfdisc;
	diffusion *pdiff;
	pressure *ppress;
	poisson *ppois;
	turbulence *pturb;
	solver *psolv;
    solver *ppoissonsolv;
	ioflow *pflow;
    reini *preini;
    sixdof *p6dof;
    fsi *pfsi;
};

#endif
