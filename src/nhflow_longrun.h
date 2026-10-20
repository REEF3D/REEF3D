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

#ifndef NHFLOW_LONGRUN_H_
#define NHFLOW_LONGRUN_H_

#include<vector>
#include<string>

class lexer;
class fdm_nhf;
class ghostcell;
class nhflow_amr;

using namespace std;

// ---------------------------------------------------------------------------
// NHFLOW long-run statistics (P 160 t0 dt_out): from t0 on, every time step adds the state of every
// column to running sums, and every dt_out (and at the end of the run) the maps of the record
// [t0, t] are written to REEF3D_NHFLOW_Tidal (.vts per grid and rank, one .vtm):
//
//  - harmonic analysis on the fly: least squares fit x(t) = Z0 + sum_k a_k cos(w_k t) + b_k sin(w_k t)
//    of the water level eta and the depth-averaged velocity (u, v) of every column. The normal
//    equations are integrals over time (weights dt, so a varying time step does not bias the fit);
//    the matrix is the same for all columns, the right-hand sides are summed per column. At output
//    the constituents that the record resolves are chosen (Rayleigh criterion P 161, default 1): in the
//    order M2 S2 K1 O1 M4 N2 K2 P1 Q1 MS4 MN4 M6 2N2, a constituent is taken if |f_k - f_l| T >= R for
//    every constituent l before it (taken or not) and f_k T >= R (separation from the mean).
//    Output per constituent: amplitude and phase of eta (x = A cos(w t - g), g in degrees, relative to
//    t = 0 of the run, no nodal corrections), and the tidal current ellipse (major and minor semi-axis,
//    minor > 0 anticlockwise, inclination of the major axis in degrees from x).
//  - time averages: eta, depth-averaged u, v (Eulerian residual current), the depth-integrated
//    transport qx, qy (residual transport), and min/max of eta (wet times), max of the speed, and the
//    fraction of the time the column is wet (inundation frequency).
//
// Dry columns: eta is their bed level, u = v = 0 (part of the fit as such); use wet_fraction to judge
// intertidal columns. Grids: the native grid of the rank (incl. the ring of a multi-block grid) and the
// blocks 1.. of a multi-block grid (G 50 1). The refined patches of G 1 are not included.
// The sums are not part of a restart.
// ---------------------------------------------------------------------------

class nhflow_longrun
{
public:
    nhflow_longrun(lexer*);
    ~nhflow_longrun();

    // after the momentum step (simtime advanced, dt the step just done)
    void step(lexer*, fdm_nhf*, ghostcell*, nhflow_amr*);

    // a grid: its cells [i0,i1) x [j0,j1) are written
    void add_grid(lexer*, fdm_nhf*, int i0, int i1, int j0, int j1);

    static const int NC = 13;

private:
    struct grid
    {
        lexer *q;
        fdm_nhf *d;
        int i0,i1,j0,j1;
        vector<double> s;           // per column (slice index): NV sums
    };

    void accumulate(grid&, const double *phi, double w);
    void write(lexer*, ghostcell*);
    int select(double T, vector<int> &sel) const;

    vector<grid> G;
    bool ini = false;
    double t0, dtout, tnext;
    double Twin = 0.0;              // record length (s)
    int count = 0;
    int NB;                         // basis size 2 NC + 1
    int NV;                         // sums per column
    vector<double> M;               // NB x NB, integral phi phi^T dt
    vector<double> omega;           // rad/s
    vector<string> name;
    double rayleigh;
};

#endif
