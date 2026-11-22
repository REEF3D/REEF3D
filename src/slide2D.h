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

#ifndef SLIDE2D_H_
#define SLIDE2D_H_

#include"increment.h"
#include"slice4.h"
#include"rheology_core.h"
#include<fstream>
#include<functional>
#include<iosfwd>

class lexer;
class fdm2D;
class fdm_nhf;
class matrix2D;
class vec2D;
class ghostcell;
class slice;
class solver2D;

// Two-layer landslide model (W 150 1): a depth-averaged slide layer (h_s, h_s u_s) under the SFLOW water.
//
// - slide:  d(h)/dt + div(h u) = 0
//           d(h u)/dt + div(h u u) + g' grad(h^2/2) = - g' h grad(z_b) - r sub g h grad(eta_w) + h g_s
//                                                   - tau_b/rho_s + tau_i/rho_s
//           g' = (1 - r sub) g, r = rho_w/rho_s, sub = min(1, h_w/W 156) (0 subaerial, 1 submerged),
//           eta_w the water surface, g_s the along-slope gravity W 20, W 21 (inclined frame, no water)
// - basal friction from the rheology core (W 120 - W 128) with sigma_n' = (1 - lambda)(rho_s - sub rho_w) g_n h
// - interface drag tau_i = rho_w C_d |u_w - u_s| (u_w - u_s) (W 155), equal and opposite in the water
// - the water sees the bed z_b + h: b->bed, then sflow_eta::depth_update keeps the water column WL and moves the
//   surface with the bed (the water volume is unchanged by the moving bed)
// - W 163 1: the water pressure on the slide surface includes the non-hydrostatic pressure at the bed of the
//   water, - r sub h grad(g eta_w + p_nh,b/rho_w) (NHFLOW: d->P of the bed node; SFLOW A 220 2/3: the bed pressure of the
//   quadratic pressure); 0 (default) the hydrostatic water pressure only
// - W 162 a > 0: the water sees z_b + h_f instead, h_f - div(D grad h_f) = h, D = a d^2 (d the still water depth
//   over z_b, F 60): a conservative low-pass filter of the slide under the water, response 1/(1 + a (k d)^2) against
//   1/cosh(k d) of the potential-flow transfer (Kajiura 1963; a = 0.5 matches for small k d); no filter where dry
//
// Numerics: finite volumes on the SFLOW slice, MUSCL (minmod) with the hydrostatic reconstruction of Audusse et
// al. (2004), HLL fluxes, SSP-RK2 sub-steps with their own CFL (W 157) inside the water step; friction and
// drag as backward-Euler steps after each stage (as the SFLOW adapter of the core); pore pressure as in W 126.

class slide2D final : public increment
{
public:
    slide2D(lexer*, fdm2D*, ghostcell*);
    slide2D(lexer*, fdm_nhf*, ghostcell*);     // NHFLOW (W 150 with A 10 5)
    ~slide2D() = default;

    // the initial slide (W 151 - W 153) on the bed of the grid; sets the bed of the water to z_b + h
    void ini(lexer*, fdm2D*, ghostcell*);
    void ini(lexer*, fdm_nhf*, ghostcell*);

    // one water step: the slide over p->dt in sub-steps, with the water state of the end of the step
    void start(lexer*, fdm2D*, ghostcell*);
    void start(lexer*, fdm_nhf*, ghostcell*);

    // the largest water step: W 159 slide steps (no water: the slide alone sets the step)
    double dt_cap(lexer*, ghostcell*);

    // drag of the slide on the water: the rhs of the water momentum (b->F, b->G), called in the water stages
    void u_source(lexer*, fdm2D*, slice&, slice&, slice&);
    void v_source(lexer*, fdm2D*, slice&, slice&, slice&);
    // NHFLOW: the drag on the bottom water layer (comp 0: U into F, 1: V into G), called with the bed friction
    void drag_nhflow(lexer*, fdm_nhf*, int, double*, double*, slice&);
    // semi-implicit NHFLOW (W 155): the drag rate of the bottom layer of a cell (thickness dz) with the water velocity
    // uw, vw; stores the rate and the slide velocity the water relaxes to (sgc; suc, svc)
    double si_drag(lexer*, int, int, double, double, double);
    void si_drag_reset(lexer*);
    slice4 sgc, suc, svc;

    // output (SFLOW FSF VTP and the LAGOON store, P 18)
    void offset_vtp(lexer*, int*, int&);
    void name_vtp(std::ostream&, int*, int&);
    static void name_pvtp(std::ostream&);
    void print_vtp(lexer*, std::ostream&);

    // checkpoints (P 48)
    void ckpt_fields(const std::function<void(const char*, char, double*)>&, const std::function<void(const char*, char, int*)>&);

    // the volumes of the slide and of the water [m^3] (global)
    double volume_slide(lexer*, ghostcell*);
    double volume_water(lexer*, ghostcell*);

    slice4 h, qx, qy;       // the slide depth and discharges
    slice4 lambda;          // basal pore pressure ratio
    slice4 zb;              // the bed under the slide
    slice4 pnb;             // W 163: the non-hydrostatic pressure of the water at the bed / rho_w (set by the water model)
    slice4 wb;              // the velocity of the bed seen by the water, (bed^{n+1} - bed^n)/dt (non-hydrostatic water)

private:
    void ghosts(lexer*, ghostcell*, slice&, slice&, slice&);
    void rhs(lexer*, slice&, slice&, slice&);
    void flux_x(lexer*, slice&, slice&, slice&);
    void flux_y(lexer*, slice&, slice&, slice&);
    void sources(lexer*, slice&, slice&, slice&);
    void friction(lexer*, slice&, slice&, slice&, double);
    void positivity(lexer*, slice&, slice&, slice&);
    void pore_update(lexer*, ghostcell*);
    void bed_update(lexer*, ghostcell*);
    double timestep(lexer*, ghostcell*);
    void log(lexer*, ghostcell*, int);

    double gn_cell(lexer*);
    double sub_cell(lexer*, int, int);
    double gred(lexer*, int, int);
    double vel(slice&, slice&, int, int);

    static double minmod(double, double);

    rheo::basal_params q;

    // work arrays: stage states, rhs, face values
    slice4 h0, qx0, qy0;
    slice4 h1, qx1, qy1;
    slice4 Rh, Rx, Ry;
    slice4 FhX, FxLX, FxRX, FyX;     // x faces i+1/2: mass, x-momentum seen from the left / right cell, y-momentum
    slice4 FhY, FxY, FyLY, FyRY;     // y faces j+1/2
    slice4 hE, hW, zE, zW, hN, hS, zN, zS;     // reconstructed depth and bed at the faces of each cell
    slice4 lnew;
    slice4 hstep;           // the slide depth at the start of the water step
    slice4 gr;              // the reduced gravity g' of the cells (the water state of the step)
    slice4 hf, dfil;        // W 162: the filtered slide seen by the water, the filter diffusivity
    solver2D *pfilt=nullptr;

    void filter(lexer*, ghostcell*);
    matrix2D *Mf=nullptr;
    vec2D *xf=nullptr, *rf=nullptr;

    // the water: SFLOW or NHFLOW
    void bind(lexer*, fdm2D*);
    void bind(lexer*, fdm_nhf*);
    void ini_common(lexer*, ghostcell*);
    void step(lexer*, ghostcell*);
    slice *wbed=nullptr, *weta=nullptr, *wWL=nullptr, *wU=nullptr, *wV=nullptr;
    slice4 uwb, vwb;        // NHFLOW: the velocity of the bottom layer
    double wdry=0.0;        // dry below this water depth (A 244 SFLOW, A 544 NHFLOW)
    bool nhf=false;

    double g, rhos, rhow, r, Cd, dsub, cfl, hmin;
    int nsub, nsub_max;
    double vol0;
    double p_W162=0.0;     // W 162 of the run (checkpoint fields)
    bool logfresh=true;     // the log starts new (not continued from a checkpoint)

    std::ofstream logout;
};

#endif
