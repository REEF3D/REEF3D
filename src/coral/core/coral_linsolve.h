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
// REEF3D::CORAL: linear solver of the Newton systems (time integration, statics): the matrix comes as
// triplets; up to dense_max unknowns it is factorized dense (partial pivoting LU), above that sparse
// (SuperLU-type supernodal LU with COLAMD ordering), so frames of thousands of nodes stay cheap.

#ifndef CORAL_LINSOLVE_H_
#define CORAL_LINSOLVE_H_

#include"coral/core/coral_types.h"
#include<Eigen/Sparse>
#include<memory>
#include<vector>

namespace coral
{
    class linear_solver
    {
    public:
        linear_solver();
        ~linear_solver();
        linear_solver(const linear_solver&) = delete;
        linear_solver &operator=(const linear_solver&) = delete;

        // factorize the N x N matrix of the triplets (duplicates are summed); false if singular
        bool compute(int N, const std::vector<Eigen::Triplet<double>> &T);
        // solution; NaN if the sparse factorization found the matrix singular
        VecX solve(const VecX &b) const;
        bool sparse() const {return sparse_;}

        int dense_max = 400;

    private:
        struct impl;
        std::unique_ptr<impl> p_;
        bool sparse_ = false;
        bool ok_ = false;
    };
}

#endif
