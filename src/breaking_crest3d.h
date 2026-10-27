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

#ifndef BREAKING_CREST3D_H_
#define BREAKING_CREST3D_H_

#include"increment.h"
#include<vector>
#include<functional>
#include<fstream>

class lexer;
class ghostcell;
class slice;

using namespace std;

// Crest tracking for the breaking onset B = u_s/C (Barthelemy et al. 2018) on the free surface of FNPF and NHFLOW,
// without a preferred direction (waves in any direction on the horizontal grid) and without gathering the surface.
//
// 1. crest points: a cell is a crest point when eta is a maximum along the ridge normal n, the eigenvector of the
//    most negative curvature of eta (Hessian; 2D grids: n = x). Sub-grid position and crest elevation from the
//    parabola through eta at -d, 0, +d along n (bilinear values), d the distance across one cell along n.
//    The sign of n follows the crest point of the last step it continues, a new one takes the sign of u_s.n.
// 2. tracking: a crest moves less than one cell per step, so every crest point continues the crest point of the
//    last step nearest to its predicted position (moved by C dt along n) within two cells (same orientation). A crest
//    point that continues no breaking crest point but lies on the front face of one of the last step takes over its
//    breaking (and its C until it has its own): when a breaking crest turns into a bore (NHFLOW), the maximum jumps
//    from the collapsing crest to the top of the bore ahead of it. The crest points within reach of a subdomain border are
//    exchanged between the ranks, so the result does not depend on the decomposition. The crest speed C is the
//    least-squares slope of the positions projected on n (samples every quarter cell of travel, up to 8, C from
//    at least 4), B = u_s.n / C.
// 3. rays: the crest points with B > 0.5 or breaking are gathered. Along each ray x_c + s n every rank samples its
//    own part of the surface (every half cell): the first minimum ahead and behind (troughs; a dry sample ends the
//    ray), the crest-to-trough height H, the water depths at crest and trough, the largest front-face slope, and the
//    caller's dissipation per unit viscosity along the painted interval [-s_back, s_ahead]. Minima, maxima and
//    sums over the ranks; every rank then has the same numbers for every ray.
// 4. state: onset at B > Bon with H > Hmin and C > Cmin c_lin, c_lin the linear phase speed at the crest depth for
//    the period of the generated waves (or k = 2 pi/L): no onset of short waves riding on the waves, end when the largest front-face slope falls below slope_end (not in the
//    step of onset). Strength per unit crest length eps (caller), viscosity nu = eps / (dissipation per nu),
//    at most |C| (s_back + s_ahead).
// 5. painting: cells within 1.5 cells of a ray, with the projection s of their centre on [-s_back, s_ahead],
//    take max(vb, nu S(s)), S the sin^2 ramp over the first and last 10 % of the interval; flags on the interval
//    +- 2 cells.
//
// Cartesian grids (geo_curv 0), uniform or stretched.

struct bc3_point
{
    double x=0.0, y=0.0, nx=1.0, ny=0.0;   // crest position, unit ridge normal (direction of travel)
    double eta=0.0, us=0.0, wl=0.0, d=0.0; // crest elevation, u_s.n, total water depth at the crest, step along n
    double C=0.0, B=0.0, b=0.0, ton=-1.0e20;
    double eps=0.0, nu=0.0, dpn=0.0, H=0.0, ht=0.0, L=0.0, sa=0.0, sb=0.0, sback=0.0, kw=0.0;
    long long id=0;
    int state=0, hasC=0, nh=0, nbh=0, i=0, j=0, rank=0;
    double th[8], xh[8], yh[8];            // position samples (t, x, y), oldest first
    double tb[8], bh[8];                   // B samples (t, B)
};

class breaking_crest3d : public increment
{
public:
    breaking_crest3d(lexer*);
    virtual ~breaking_crest3d();

    // parameters
    double Bon = 0.85;          // onset B
    double Bcand = 0.5;         // rays and log for B above this
    double slope_end = 0.1405;  // tan(8 deg)
    double Hmin = 0.0;          // smallest crest-to-trough height
    double Cmin = 0.5;          // onset only for crests travelling at least at Cmin of the linear phase speed
    double Tdom = 0.0;          // period of the generated waves (0: unknown, the phase speed from k = 2 pi / L)
    int deep = 1;               // deep-water strength b(Gamma) (Derakhti et al. 2018)

    // one tracking step at time t: eta, horizontal surface velocity (us, vs), total water depth WL, wet flags,
    // zone (int per slice cell, 1: no crests, may be nullptr). With the velocity potential at the surface (phi,
    // FNPF) u_s.n is the derivative of the parabola of phi along n at the crest, otherwise the parabola of
    // (us, vs).n through -d, 0, +d at the crest.
    void detect(lexer*, ghostcell*, slice &eta, slice &us, slice &vs, slice &WL, int *zone, double t, slice *phi=nullptr);

    // rays of the candidates, state, strength and painting.
    //   sback(point)              length of the interval behind the crest that is painted (0: crest to trough)
    //   integrand(point, x, y)    dissipation per unit viscosity, crest length and ray length at (x, y)
    //   epsf(point)               strength per unit crest length
    // vb: viscosity (max), brk: flags 1, kw: wavenumber 2 pi / L of the ray of the largest viscosity (may be nullptr)
    void rays(lexer*, ghostcell*, slice &eta, slice &WL, double t,
              const std::function<double(const bc3_point&)> &sback,
              const std::function<double(const bc3_point&, double, double)> &integrand,
              const std::function<double(const bc3_point&)> &epsf,
              slice &vb, int *brk, slice *kw);

    // the rays of the last call (all ranks have them) and the number of breaking crest points
    std::vector<bc3_point> cand;
    int nbreak = 0;

    // still-water depth (eta - WL taken positive down) of the deepest onset so far
    double deepest_onset = -1.0e20;

    // crest log (P 59 1): rank 0 writes the candidates
    void log(lexer*, double t, const char *file);

    // bilinear value of a slice at (x, y); valid within one cell outside the own subdomain
    double ip(lexer*, slice&, double x, double y);
    // cell containing (x, y) (own cells and ghosts), -1e9 when outside
    int cell_i(lexer*, double x);
    int cell_j(lexer*, double y);
    bool own(lexer*, double x, double y);

private:
    std::vector<bc3_point> pts, prev, halo;   // crest points of this step, of the last step, neighbours' of the last step
    long long nextid = 0;
    double tprev = -1.0e20;   // time of the last detect
    int nproc = 1, rank = 0;
    std::ofstream logf;

    bool wetcell(lexer*, int, int);
    double celld(lexer*, int, int, double, double);
    void exchange_halo(lexer*);
};

#endif
