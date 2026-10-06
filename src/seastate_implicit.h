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

#ifndef SEASTATE_IMPLICIT_H_
#define SEASTATE_IMPLICIT_H_

#include"increment.h"
#include<vector>

class lexer;
class ghostcell;
class fdm_seastate;
class seastate_store;
class seastate_exchange;
class seastate_source;
class slice;
class sliceint;

// spectra on the faces of the cells covered by a finer grid (REEFAMR, seastate_amr): face(i,j,s)
// is the spectrum that leaves covered cell (i,j) through its side s (0 x-, 1 x+, 2 y-, 3 y+),
// nullptr: the spectrum of the cell
struct seastate_faces
{
    virtual ~seastate_faces() {}
    virtual const float *face(int i, int j, int s) const = 0;
};

using namespace std;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - implicit transport of the wave action density

  (N - N0)/dt + d(cx N)/dx + d(cy N)/dy + d(c_sig N)/dsig + d(c_theta N)/dtheta = 0

Finite volumes, first-order upwind in x, y, sigma and theta, implicit
(BDF1; stationary: 1/dt = 0). Face velocities: mean of the two cells
(geographic, sigma), analytic at the face angle (theta).

The matrix is an M-matrix: positive diagonal (1/dt + outflow
coefficients), non-positive off-diagonals (inflow coefficients), and
every column sums to 1/dt >= 0 (each outflow of a cell is the inflow of
its neighbour). Its inverse is non-negative, so N >= 0 without limiters
or clipping, for any time step.

Solver: Gauss-Seidel with four sweeps per iteration, one per directional
quadrant (theta in [q pi/2, (q+1) pi/2)), each in the downstream cell
order of its quadrant, so that geographic upwinding is exact within a
sweep. In each cell the directions of the quadrant are solved together
per frequency (tridiagonal in theta, Thomas algorithm); frequencies and
the directions of the other quadrants use the latest values. Halo
exchange after each sweep (block Gauss-Seidel across ranks).

Source terms (Phase 2, seastate_source): dN/dt = P - D N with P >= 0
added to the right-hand side and D >= 0 to the diagonal, computed per
cell from the latest spectrum in every sweep (Gauss-Seidel), so the
matrix stays an M-matrix. With the deep-water physics the change of N
per iteration is limited by the action density limiter of the source
terms (seastate_source::limit), which keeps N >= 0.

Boundaries: land (inactive cells) absorbs (outflow, no inflow). Domain
sides with side[s] = 1 (x-, x+, y-, y+) let the boundary spectrum Nb
enter; sides with side[s] = 2 are zero-gradient (the inflow spectrum
is the spectrum of the boundary cell itself, implicit as long as the
diagonal stays dominant, otherwise with the latest value); other sides
are open without incoming waves. Surfbeat (boundary_rows): the x- side
takes a spectrum per row j instead of Nb; with a boundary series (A 711 3,
boundary_sides) every side with side[s] = 1 takes a spectrum per boundary
cell. Frequency range ends: outflow only.

Surfbeat (second_order, A 775 2, cell_surfbeat): the wave groups travel
far compared with their length, so the transport is Crank-Nicolson in
time (theta = 1/2 with the old values N0 of the step) and the
geographic fluxes get a second-order correction
F2 - F1 = c/2 phi(r) (N_down - N_up) (van Leer limiter, r = (N_up -
N_upup)/(N_down - N_up)) between active cells, with the latest values
(deferred correction); the source terms stay implicit, N is clipped at
zero. The first-order path (cell) is unchanged.
Directions: periodic (full circle).

Mesh refinement (Phase 5b, seastate_amr): the same solver runs on every
grid of the hierarchy. range() limits the sweeps to the interior of a
patch (the cells around it hold the spectra of the neighbour grids and
act as inflow); skip() leaves out the cells covered by a finer grid,
whose spectra are restricted from it, and faces() gives the spectra
that leave the finer grid through the faces of the covered cells (the
mean of the two fine cells next to the face), so that the coarse grid
takes the outflow of the fine grid. Without these calls the solver
works on the whole grid as before.
--------------------------------------------------------------------*/

class seastate_implicit : public increment
{
public:
    seastate_implicit(lexer*, fdm_seastate*);

    void sources(seastate_source *s) {src = s;}

    // surfbeat: spectra of the x- boundary per row j (local, 0..knoy-1), nbin each, at the new
    // and the previous time level; nullptr: Nb
    void boundary_rows(const vector<float> *rows, const vector<float> *rows_old) {Nbs[0] = rows; Nbs0[0] = rows_old;}

    // boundary spectra per side and boundary cell (A 711 3): x-, x+ per local row j, y-, y+ per
    // local column i, nbin each; nullptr: Nb
    void boundary_sides(const vector<float> *const s[4]) {for(int k=0; k<4; ++k) Nbs[k] = s[k];}

    // wind field (A 730 2): U10 [m/s] and direction [rad] per cell for the source terms
    void wind_field(slice *U10, slice *dir) {wU = U10; wD = dir;}

    // surfbeat (A 775 2): Crank-Nicolson, second-order geographic fluxes (cell_surfbeat);
    // needs N0 (2 iterations) and 2 halo layers
    void second_order(bool s) {second = s;}

    void iterate(lexer*, ghostcell*, fdm_seastate*, seastate_exchange*, const seastate_store *N0,
                 double rdt, const vector<float> &Nb, const int side[4], bool refraction, bool fshift);

    // one sweep of the directional quadrant q (iterate: q = 0..3, halo exchange after each)
    void sweep(lexer*, fdm_seastate*, int q, const seastate_store *N0, double rdt,
               const vector<float> &Nb, const int side[4], bool refraction, bool fshift);

    // mesh refinement: sweeps over [i0,i1]x[j0,j1] only; cells with skip(i,j) == 1 not solved;
    // spectra on the faces of the skipped cells
    void range(int i0_, int i1_, int j0_, int j1_) {ri0 = i0_; ri1 = i1_; rj0 = j0_; rj1 = j1_; ranged = true;}
    void unrange() {ranged = false;}
    void skip(sliceint *s) {skp = s;}
    void faces(const seastate_faces *f) {fcs = f;}

private:

    void cell(lexer*, fdm_seastate*, int q, const float *N0, double rdt,
              const vector<float> &Nb, const int side[4], bool refraction, bool fshift);

    void cell_surfbeat(lexer*, fdm_seastate*, int q, const seastate_store *N0, double rdt,
                       const vector<float> &Nb, const int side[4], bool refraction);

    int nsig, ndir;
    seastate_source *src;
    const vector<float> *Nbs[4], *Nbs0[4];
    slice *wU, *wD;
    bool second;
    bool ranged = false;
    int ri0 = 0, ri1 = -1, rj0 = 0, rj1 = -1;
    sliceint *skp = nullptr;
    const seastate_faces *fcs = nullptr;
    vector<double> P, D;                // source terms of the cell (src != nullptr)
    vector<int> m0, m1;                 // first and last direction of each quadrant
    vector<double> csig;                // c_sigma of the cell, per frequency, for the current direction
    vector<double> la, di, up, rhs, sol, cp, dp;   // tridiagonal system of one frequency
};

#endif
