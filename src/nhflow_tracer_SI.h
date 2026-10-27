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

#ifndef NHFLOW_TRACER_SI_H_
#define NHFLOW_TRACER_SI_H_

#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include<vector>
#include<fstream>
#include<functional>

class lexer;
class fdm_nhf;
class ghostcell;
class nhflow_momentum_SI;

using namespace std;

// ---------------------------------------------------------------------------
// Tracers of the hydrostatic semi-implicit NHFLOW (A 534 1, A 600): salinity (A 600 1), and
// temperature (A 600 2), passive in P12a (no density yet).
//
// Finite volumes on the sigma layers with the transport of the step: the layer face transport
// T_k = theta Q_k^{n+1} + (1-theta) Q_k^n that moved the water (nhflow_momentum_SI::face_update), and
// the vertical transport omega from the layer continuity with the same T_k and the water depths WL^n,
// WL^{n+1} of the step:
//
//     omega_{k+1/2} = omega_{k-1/2} - dsig_k ( (WL^{n+1} - WL^n)/dt + div_h T_k ),  omega_{-1/2} = 0
//
// so that a uniform tracer stays uniform and the tracer mass changes only by the open boundaries
// (to round-off, also with wetting and drying, relaxation zones and mass sources: the remainder
// omega_{K+1/2} at the surface carries the concentration of the top layer).
// One step (after the momentum step), in nsub substeps of dt/nsub, the water depth linear between
// WL^n and WL^{n+1}:
//  1. horizontal, explicit: the face values of the upwind side, TVD limited (A 609: 1 upwind, 2 van
//     Leer, 3 MC) along the grid index with the Lax-Wendroff factor (1 - nu)/2 of the face Courant number
//     (Sweby 1984: second order in time and space), first order at the open boundaries; horizontal
//     diffusion along the sigma layers (A 604, through the open faces only)
//  2. vertical, implicit per column: upwind advection with omega, diffusion K_v = EV/sigma_t + K0
//     (A 605; EV of the turbulence model or of A 547), an M-matrix: positive for any omega
// nsub (the same on all grids and ranks, at most A 610) keeps dt (in + 2 out + diffusion) / content of every
// layer below 0.9 (upwind A 609 1: dt (out + diffusion)), so that the horizontal step is bounded; the
// columns shallower than 10 A 544 do not count. A layer beyond the bound in its substep (draining
// shoreline cells, or beyond A 610) is flushed: all its outflow carries its content mixed with its inflow
// of the substep, (V C + dt sum_in T C_in)/(V + dt in), first order, bounded and conservative (pass 1, filled
// to the neighbours before the fluxes). A value outside the bounds of its neighbour columns and inflow
// values (the vertical step of a layer drained below zero) is clipped, the clipped mass in the log.
// Open boundaries (nhflow_momentum_SI::bcf): inflow carries the boundary value, outflow the cell value;
// the discharge boundaries (rivers: B 60, G 54) A 606, the water level and radiating boundaries (sea)
// A 607 with the Thatcher-Harleman return (A 608, per layer: after the flow turns inward the value
// goes from the last outflow value to the sea value as 1/2 (1 - cos(pi t/T_TH))).
// Initial values: A 601 uniform, A 602 / A 603 linear between two points.
// Multi-block grids (G 50 1): nhflow_amr::si_tracer fills the cells around the blocks before every
// substep; a face between two blocks gets the same flux on both sides (same T_k and face depth, same
// cell values).
// Output: the log REEF3D_NHFLOW_Log/REEF3D_NHFLOW_tracer.dat (mass, boundary flux, clipping, range),
// the gauges P 51 (REEF3D_NHFLOW_Log/REEF3D_NHFLOW_tracer_gauges.dat: depth mean, bottom, top of the
// nearest cell), the 3D output (P 10) and the block output (depth mean, bottom, top).
// ---------------------------------------------------------------------------

class nhflow_tracer_SI
{
public:
    nhflow_tracer_SI(lexer*, fdm_nhf*, nhflow_momentum_SI*);
    ~nhflow_tracer_SI();

    // the step: prepare (initial values, omega, boundary values, the ratio of the substeps), nsub
    // substeps (the caller fills the cells around the grid before each), finish
    double prepare(lexer*, fdm_nhf*);
    void flushed(lexer*, fdm_nhf*, int s, int nsub);   // pass 1 of a substep: the outflow value of flushed cells (Cx)
    void substep(lexer*, fdm_nhf*, int s, int nsub);   // pass 2: the fluxes and the column (the caller fills Cx between)
    static int substeps(lexer*, double rmax);

    // the values of the solid and outside cells next to fluid cells, and of the layers below the bed and
    // above the surface: the nearest fluid value (output at the nodes); after the step
    void ghosts(lexer*, fdm_nhf*);

    // the field of tracer m of a grid (0 salinity, 1 temperature)
    static double* field(fdm_nhf*, int m);

    // the sums of the step over the cells own(i,j) and the open faces of this grid:
    // per tracer mass, boundary inflow, clipped mass, min, max; violations of the bounds
    struct stats
    {
        double mass[2], bin[2], clip[2], cmin[2], cmax[2];
        double nviol;
    };
    void sums(lexer*, fdm_nhf*, const std::function<bool(int,int)> &own, stats &s);
    static void reduce(lexer*, stats&);         // over the ranks (MPI)
    static void stats_zero(stats&);
    void write_log(lexer*, const stats&, int nsub, double rmax);

    // gauges P 51: the nearest own cell centre per gauge (distance, cell); values depth mean, bottom, top
    void gauge_find(lexer*, const std::function<bool(int,int)> &own, vector<double> &dist, vector<int> &cell);
    void gauge_values(lexer*, fdm_nhf*, const vector<int> &cell, vector<double> &val);
    static void gauge_combine(lexer*, vector<double> &dist, vector<double> &val);   // over the ranks: the nearest
    static void gauge_write(lexer*, const vector<double> &val);

    // the fluid mask (1 in fluid cells; the caller exchanges / fills it with the tracers)
    double *mask;

    // the outflow value of a flushed cell (its content mixed with its inflow), -1e30 elsewhere; per tracer
    double *Cx[2];

private:
    // pass 2: the face values of the substep (x and y faces of the own cells), the range of every column
    double *fvx[2], *fvy[2];
    vector<double> cmn[2], cmx[2];

public:
    void mask_own(lexer*);

    // checkpoints: the Thatcher-Harleman state, the initialisation flag
    void ckpt_fields(const std::function<void(const char*, char, double*)> &fd);

    int nt;                     // number of tracers

private:
    double bvalue(int m, int type, double thc, double tht) const;
    double face_value(lexer*, int mm, int dir, int ii, int jj, int kk, double T);
    double lim(double a, double b) const;
    int cidx(lexer*, int, int, int) const;
    int sidx(lexer*, int, int) const;
    bool fluid(lexer*, int, int) const;
    void initial(lexer*, fdm_nhf*);

    nhflow_momentum_SI *m;
    int K, N, scheme;
    double dt;
    double riv[2], sea[2], tth, kh, sigt, k0;     // A 606, A 607, A 608, A 604, A 605
    double sdt = 0.0, sfa = 0.0;                  // the substep: its length, its start (fraction of the step)
    bool cx_on = false;                           // pass 2: the face values of flushed upwind cells from Cx

    // geometry: face lengths, normal distance of the centres, cell area
    slice1 lx, dnx;
    slice2 ly, dny;
    slice4 ar;

    // omega of the transport per column (K+1 levels), the water depths of the step
    vector<double> om;
    slice4 wl0, wl1;

    // the values at the start of the substep (all cells), the work column
    double *Co[2];
    vector<double> ca, cb, cc, cr, cx;

    // open faces of the step: face index -> boundary face (x and y faces), boundary values per layer
    vector<int> bfx, bfy;
    vector<double> bval[2];

    // Thatcher-Harleman: per cell and layer of the boundary cells the last outflow value and the time of
    // inflow since; initialisation flag (checkpoint)
    double *thc[2], *tht;
    slice4 trini;
    bool ini;

    // step sums of this grid
    vector<double> clipc[2], binc[2], violc;     // per column: clipped mass, boundary inflow, bound violations (own cells summed)

    ofstream tlog;
};

#endif
