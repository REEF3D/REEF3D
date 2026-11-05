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

#include"coral/core/coral_linsolve.h"
#include<Eigen/SparseLU>
#include<limits>

namespace coral
{

struct linear_solver::impl
{
    Eigen::PartialPivLU<MatX> dense;
    Eigen::SparseMatrix<double> A;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> lu;
    bool pattern = false;
};

linear_solver::linear_solver() : p_(new impl) {}
linear_solver::~linear_solver() = default;

bool linear_solver::compute(int N, const std::vector<Eigen::Triplet<double>> &T)
{
    sparse_ = N > dense_max;
    if(!sparse_)
    {
        MatX J = MatX::Zero(N, N);
        for(const auto &t : T)
            J(t.row(), t.col()) += t.value();
        p_->dense.compute(J);
        ok_ = true;
        return true;
    }
    Eigen::SparseMatrix<double> A(N, N);
    A.setFromTriplets(T.begin(), T.end());
    A.makeCompressed();
    // the symbolic analysis is kept while the sparsity pattern does not change
    const bool same = p_->pattern && A.nonZeros()==p_->A.nonZeros()
        && std::equal(A.outerIndexPtr(), A.outerIndexPtr() + N + 1, p_->A.outerIndexPtr())
        && std::equal(A.innerIndexPtr(), A.innerIndexPtr() + A.nonZeros(), p_->A.innerIndexPtr());
    p_->A = A;
    if(!same)
    {
        p_->lu.analyzePattern(p_->A);
        p_->pattern = true;
    }
    p_->lu.factorize(p_->A);
    ok_ = p_->lu.info()==Eigen::Success;
    return ok_;
}

VecX linear_solver::solve(const VecX &b) const
{
    if(!sparse_)
        return p_->dense.solve(b);
    if(!ok_)
        return VecX::Constant(b.size(), std::numeric_limits<double>::quiet_NaN());    // singular: callers see a failed step
    return p_->lu.solve(b);
}

}
