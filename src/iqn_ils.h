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

#ifndef IQN_ILS_H_
#define IQN_ILS_H_

// Interface quasi-Newton with an inverse Jacobian from a least-squares model (IQN-ILS, Degroote et al. 2009)
// for a partitioned fixed point x = H(x) (here: membrane node velocities -> fluid stage -> structure -> velocities).
//
//   r^k = H(x^k) - x^k,   columns  V = [r^{k} - r^{k-1}, ...],  W = [H(x^k) - H(x^{k-1}), ...]
//   c   = argmin |V c + r^k|        (QR by modified Gram-Schmidt, nearly dependent columns dropped)
//   x^{k+1} = H(x^k) + W c
//
// Without columns a relaxed step x^{k+1} = x^k + omega r^k. At most min(maxcols, n/2) columns are used, so that
// the least-squares model does not fit the noise of the fluid solution (iterative pressure solver). The columns of
// the converged sequences of the last 'reuse' calls are kept and used first in the next ones (history per slot,
// e.g. per Runge-Kutta stage, whose fixed points differ by the stage step); restart() drops them.
//
// Multi-vector variant (imvj = true, IQN-IMVJ, Lindner et al. 2015): instead of reusing columns, the inverse Jacobian
// M (dx~ = M dr) of the end of the last sequence is carried over and updated with the columns of the current one,
//   M_k = M_prev + (W - M_prev V)(V^T V)^-1 V^T,     x^{k+1} = H(x^k) - M_k r^k,
// which keeps the information of all earlier sequences at the cost of an n x n matrix per slot (one matrix-vector
// product per iteration, one rank-m update per sequence). Pure numerics, no MPI: every rank calling with identical
// data gets identical results.

#include<Eigen/Dense>
#include<vector>
#include<deque>

class iqn_ils
{
public:
    iqn_ils(int slots=4, double omega=0.5, int reuse=8, int maxcols=100, double eps=1.0e-2, bool imvj=false);

    void begin(int slot);                                                   // start a new fixed-point sequence
    Eigen::VectorXd next(const Eigen::VectorXd &x, const Eigen::VectorXd &xt);  // input x, output H(x): next input
    void end();                                                             // sequence finished: keep its columns
    void restart();                                                         // drop the columns of earlier sequences

    int columns() const {return ncols_;}                                    // columns used in the last update

private:
    struct block { Eigen::MatrixXd V, W; };

    Eigen::VectorXd next_imvj(const Eigen::VectorXd&, const Eigen::VectorXd&);
    void end_imvj();
    int qr(const std::vector<Eigen::VectorXd>&, int, Eigen::MatrixXd&, Eigen::MatrixXd&, std::vector<int>&) const;

    bool imvj_;
    std::vector<Eigen::MatrixXd> M_;            // IMVJ: inverse Jacobian per slot
    std::vector<char> mset_;
    std::vector<Eigen::VectorXd> mr_;           // IMVJ: M_prev r^j of the current sequence

    std::vector<std::deque<block> > hist_;
    std::vector<Eigen::VectorXd> xt_, r_;
    int slot_;
    double omega_, eps_;
    int reuse_, maxcols_, ncols_;
};

#endif
