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

Boundaries: land (inactive cells) absorbs (outflow, no inflow). Domain
sides with side[s] = 1 (x-, x+, y-, y+) let the boundary spectrum Nb
enter; other sides are open without incoming waves. Frequency range
ends: outflow only. Directions: periodic (full circle).
--------------------------------------------------------------------*/

class seastate_implicit : public increment
{
public:
    seastate_implicit(lexer*, fdm_seastate*);

    void iterate(lexer*, ghostcell*, fdm_seastate*, seastate_exchange*, const seastate_store *N0,
                 double rdt, const vector<float> &Nb, const int side[4], bool refraction, bool fshift);

private:
    void sweep(lexer*, fdm_seastate*, int q, const seastate_store *N0, double rdt,
               const vector<float> &Nb, const int side[4], bool refraction, bool fshift);

    void cell(lexer*, fdm_seastate*, int q, const float *N0, double rdt,
              const vector<float> &Nb, const int side[4], bool refraction, bool fshift);

    int nsig, ndir;
    vector<int> m0, m1;                 // first and last direction of each quadrant
    vector<double> csig;                // c_sigma of the cell, per frequency, for the current direction
    vector<double> la, di, up, rhs, sol, cp, dp;   // tridiagonal system of one frequency
};

#endif
