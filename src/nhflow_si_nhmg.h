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

#ifndef NHFLOW_SI_NHMG_H_
#define NHFLOW_SI_NHMG_H_

#include<vector>
#include<memory>
#include<mpi.h>

// Multigrid for the coupled eta-phi system of the semi-implicit non-hydrostatic NHFLOW (nhflow_momentum_SI_nh.cpp),
// after the solver of Basilisk's multilayer model (layered/nh.h): the relaxation solves the unknowns of a whole column at
// once (eta and the K interface pressures, B = K+1), the neighbouring columns fixed.
//
// The system has a block 5-point structure over the columns of a rank (blocks B x B: the cell and its neighbours in
// +x, -x, +y, -y). Level 0 is the assembled operator (double: the products of FGMRES); the preconditioner (single
// precision by default, nhflow_si_nhmg_core) has its own copy of level 0 and the coarse levels, which aggregate 2 x 2
// columns (2 x 1 with one cell in y) of a rank and take the Galerkin operator R A P with the piecewise constant
// prolongation P and R = P^T (summed in double), which keeps the block 5-point structure.
// Smoother: block Gauss-Seidel in red-black order (the inverses of the diagonal blocks, LU with pivoting in double).
// Cycle: the K-cycle of Notay (AGMG): below the finest level each coarse correction is two GCR steps preconditioned
// with the cycle of the level; with the piecewise constant aggregation the plain V-cycle needed 14 / 32 / 66 FGMRES
// iterations for 32 / 128 / 512 columns, the K-cycle 13 / 18 / 18. MPI: the halo of one column per level from the
// neighbouring ranks of the Cartesian partition (p->nb1..nb4); all ranks have the same number of levels; the coarsest
// level by sweeps.
template<class T> class nhflow_si_nhmg_core;

class nhflow_si_nhmg
{
public:
    // the level-0 box of the rank: nx x ny columns, B unknowns per column; nb: the ranks in -x, +x, -y, +y (< 0 none)
    void setup(MPI_Comm comm, int B, int nx, int ny, const int nb[4]);

    // block d (0 the cell, 1 +x, 2 -x, 3 +y, 4 -y) of the row of column (ix,iy) on level 0, B x B row-major; set(): an
    // entry and its place in the pattern of the products (entries written through block() only: mark_all())
    double* block(int ix, int iy, int d) { return &A0[((size_t)(ix*ny0 + iy)*5 + d)*BB]; }
    void set(int ix, int iy, int d, int r, int c, double v) { block(ix,iy,d)[r*B+c] = v; nz0[(size_t)d*BB + r*B + c] = 1; }
    void mark_all();
    void zero();

    // the pattern of level 0 (the products skip the entries that are zero in every column); refresh: the preconditioner
    // from the current level-0 blocks (without: the one of an earlier step, FGMRES is flexible)
    void finalize(bool refresh=true);
    void pattern();     // the pattern of level 0 only (the products, no preconditioner)

    // y = A x (double) and z = M^-1 r (one cycle) on level 0, vectors of nx*ny*B (column (ix,iy), unknown t at (ix*ny + iy)*B + t)
    void matvec(const double *x, double *y);
    void vcycle(const double *r, double *z);

    int levels() const;

    // options (read by finalize and vcycle); 2D standing wave 256^2 x 4, T/40, 5.1 iterations per step either way:
    // single precision about 10 % less time for the step than double, the fused residual about 5 %. Tried and not taken
    // (128^2 x 4): the V-cycle 13.6 iterations instead of 5.3, 13.9 s instead of 11.2 (omega_c 1.5: 14.0 iterations);
    // ktol 0.25: 5.7 iterations, no faster; two sweeps: 4.6 iterations, 7 % slower
    bool single = true;         // the preconditioner in single precision (vectors, blocks, the halo messages)
    bool kcycle = true;         // Krylov-accelerated coarse corrections (K-cycle); false: V-cycle
    bool fused = true;          // the residual after the pre-smoothing only on the colour relaxed first (zero on the one relaxed last)
    double omega_c = 1.0;       // over-correction of the coarse correction (V-cycle)
    double ktol = 0.0625;       // K-cycle: the second GCR step when |r|^2 after the first > ktol |b|^2
    int npre = 1, npost = 1, ncoarse = 20;

private:
    MPI_Comm comm;
    int B=0, BB=0, nx0=0, ny0=0;
    int nb[4];
    bool alone = true;      // no neighbouring rank: no halo exchange
    std::vector<double> A0, xh0, sbuf, rbuf;
    std::vector<char> nz0;
    std::vector<std::vector<int>> pat0;
    std::shared_ptr<nhflow_si_nhmg_core<float>> cf;
    std::shared_ptr<nhflow_si_nhmg_core<double>> cd;
};

#endif
