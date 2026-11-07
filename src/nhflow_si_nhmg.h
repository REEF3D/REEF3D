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
#include<mpi.h>

// Multigrid for the coupled eta-phi system of the semi-implicit non-hydrostatic NHFLOW (nhflow_momentum_SI_nh.cpp),
// after the solver of Basilisk's multilayer model (layered/nh.h): the relaxation solves the unknowns of a whole column at
// once (eta and the K interface pressures, B = K+1), the neighbouring columns fixed.
//
// The system has a block 5-point structure over the columns of a rank (blocks B x B: the cell and its neighbours in
// +x, -x, +y, -y). Level 0 is the assembled operator; the coarse levels aggregate 2 x 2 columns (2 x 1 with one cell in y)
// of a rank and take the Galerkin operator R A P with the piecewise constant prolongation P and R = P^T, which keeps the
// block 5-point structure (Basilisk rediscretises on its coarse levels; the operator here is assembled from the face derivatives).
// Smoother: block Gauss-Seidel in red-black order (the column solves with the LU of the diagonal blocks). Cycle: the
// K-cycle of Notay (AGMG): below the finest level each coarse correction is two GCR steps preconditioned with the
// cycle of the level; with the piecewise constant aggregation the plain V-cycle needed 14 / 32 / 66 FGMRES iterations
// for 32 / 128 / 512 columns, the K-cycle 13 / 18 / 18. MPI: the halo of one column per level from the neighbouring
// ranks of the Cartesian partition (p->nb1..nb4); all ranks have the same number of levels; the coarsest level by
// sweeps.
class nhflow_si_nhmg
{
public:
    // the level-0 box of the rank: nx x ny columns, B unknowns per column; nb: the ranks in -x, +x, -y, +y (< 0 none)
    void setup(MPI_Comm comm, int B, int nx, int ny, const int nb[4]);

    // block d (0 the cell, 1 +x, 2 -x, 3 +y, 4 -y) of the row of column (ix,iy) on level 0, B x B row-major
    double* block(int ix, int iy, int d) { return &L[0].A[((size_t)(ix*L[0].ny + iy)*5 + d)*BB]; }
    void zero();

    // the coarse operators and the LU of the diagonal blocks (after the level-0 blocks are set)
    void finalize();

    // y = A x and z = M^-1 r on level 0, vectors of nx*ny*B (column (ix,iy), unknown t at (ix*ny + iy)*B + t)
    void matvec(const double *x, double *y);
    void vcycle(const double *r, double *z);      // one cycle (K-cycle unless kcycle = false)

    int levels() const { return (int)L.size(); }
    double omega_c = 1.0;      // over-correction of the coarse correction
    int npre = 1, npost = 1, ncoarse = 20;      // one pre- and post-sweep: 8.1 iterations instead of 7.6 at a quarter less time
    bool kcycle = true;      // Krylov-accelerated coarse corrections (K-cycle)

private:
    struct level
    {
        int nx=0, ny=0;
        std::vector<double> A;      // 5 blocks per column
        std::vector<double> LU;     // the LU of the diagonal block
        std::vector<int> piv;
        std::vector<double> x, b, r;    // with the halo of one column
        std::vector<double> kb, v1, w1, v2, w2;     // K-cycle
        int nbx[4] = {-1,-1,-1,-1};
        // the columns c of row r of block d that are not zero in any column of the level (finalize), in ascending
        // order: pat[d*B + r]; the products skip the zero entries of the blocks (the same sums, bitwise)
        std::vector<std::vector<int>> pat;
    };

    MPI_Comm comm;
    int B=0, BB=0;
    int nb[4];
    std::vector<level> L;
    std::vector<double> sbuf, rbuf, tmp;
    bool alone = true;      // no neighbouring rank: no halo exchange
    void block_mul(const level &l, int d, const double *a, const double *x, double *y, double sign) const;

    int hix(const level &l, int ix, int iy) const { return (ix+1)*(l.ny+2) + (iy+1); }
    void exchange(level &l, std::vector<double> &v);
    void smooth(int lev, int sweeps, bool zero=false);
    void residual(int lev);
    void cycle(int lev);
    void kcorrect(int lev);
    void apply(int lev, std::vector<double> &x, std::vector<double> &y);
    void dots(int lev, const std::vector<double> *const *u, const std::vector<double> *const *w, int n, double *res);
    void lu(double *a, int *pv);
    void lusolve(const double *a, const int *pv, double *x);
};

#endif
