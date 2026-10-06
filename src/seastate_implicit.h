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
takes a spectrum per row j instead of Nb. Frequency range ends: outflow only.

Surfbeat (second_order, A 775 2, cell_surfbeat): the wave groups travel
far compared with their length, so the transport is Crank-Nicolson in
time (theta = 1/2 with the old values N0 of the step) and the
geographic fluxes get a second-order correction
F2 - F1 = c/2 phi(r) (N_down - N_up) (van Leer limiter, r = (N_up -
N_upup)/(N_down - N_up)) between active cells, with the latest values
(deferred correction); the source terms stay implicit, N is clipped at
zero. The first-order path (cell) is unchanged.
Directions: periodic (full circle).
--------------------------------------------------------------------*/

class seastate_implicit : public increment
{
public:
    seastate_implicit(lexer*, fdm_seastate*);

    void sources(seastate_source *s) {src = s;}

    // surfbeat: spectra of the x- boundary per row j (local, 0..knoy-1), nbin each, at the new
    // and the previous time level; nullptr: Nb
    void boundary_rows(const vector<float> *rows, const vector<float> *rows_old) {Nbx = rows; Nbx0 = rows_old;}

    // surfbeat (A 775 2): Crank-Nicolson, second-order geographic fluxes (cell_surfbeat);
    // needs N0 (2 iterations) and 2 halo layers
    void second_order(bool s) {second = s;}

    void iterate(lexer*, ghostcell*, fdm_seastate*, seastate_exchange*, const seastate_store *N0,
                 double rdt, const vector<float> &Nb, const int side[4], bool refraction, bool fshift);

private:
    void sweep(lexer*, fdm_seastate*, int q, const seastate_store *N0, double rdt,
               const vector<float> &Nb, const int side[4], bool refraction, bool fshift);

    void cell(lexer*, fdm_seastate*, int q, const float *N0, double rdt,
              const vector<float> &Nb, const int side[4], bool refraction, bool fshift);

    void cell_surfbeat(lexer*, fdm_seastate*, int q, const seastate_store *N0, double rdt,
                       const vector<float> &Nb, const int side[4], bool refraction);

    int nsig, ndir;
    seastate_source *src;
    const vector<float> *Nbx, *Nbx0;
    bool second;
    vector<double> P, D;                // source terms of the cell (src != nullptr)
    vector<int> m0, m1;                 // first and last direction of each quadrant
    vector<double> csig;                // c_sigma of the cell, per frequency, for the current direction
    vector<double> la, di, up, rhs, sol, cp, dp;   // tridiagonal system of one frequency
};

#endif
