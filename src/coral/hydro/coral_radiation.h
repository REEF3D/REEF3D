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
// REEF3D::CORAL: radiation memory of the Cummins equation
//   F_rad(t) = -A_inf xi''(t) - int_0^t K(t - tau) xi'(tau) dtau,
//   K(t) = 2/pi int_0^inf B(w) cos(w t) dw
// K is integrated exactly for B piecewise linear in w (B(0) = 0, zero above the last frequency), so
// it has no aliasing at large t. A_inf from the database or, if missing, from Ogilvie's relation
//   A_inf = A(w) + 1/w int_0^inf K(t) sin(w t) dt
// averaged over the upper half of the frequencies.
//
// State-space alternative: each K_ij(t) is realized as sum_q Re(c_q b_q exp(lambda_q t)) (eigensystem
// realization from the Hankel matrix of the samples of K_ij, order from the Hankel singular values,
// in diagonal form), so the memory is carried by linear states instead of the convolution history.
// The states are advanced exactly for a velocity linear over the step.

#ifndef CORAL_RADIATION_H_
#define CORAL_RADIATION_H_

#include"coral/core/coral_types.h"
#include<vector>

namespace coral
{
    class hydro_db;

    class retardation
    {
    public:
        // samples K(k dt), k = 0 .. tmax/dt
        void compute(const hydro_db&, double tmax, double dt);
        MatX at(double t) const;      // linear interpolation, zero beyond tmax
        double dt() const {return dt_;}
        double tmax() const {return tmax_;}
        int nsample() const {return (int)K_.size();}
        const MatX &sample(int k) const {return K_[k];}
        int ndof() const {return K_.empty() ? 0 : (int)K_[0].rows();}

    private:
        double dt_ = 0.0, tmax_ = 0.0;
        std::vector<MatX> K_;
    };

    // infinite-frequency added mass from A(w) and K(t) (Ogilvie)
    MatX ogilvie_ainf(const hydro_db&, const retardation&);

    class radiation_ss
    {
    public:
        struct pair
        {
            int i = 0, j = 0;            // force i from velocity j
            Eigen::VectorXcd lam;        // continuous poles
            Eigen::VectorXcd b, c;       // input and output weights (diagonal form)
            double error = 0.0;          // max |K_ss - K| / max |K_ij| over the samples
        };

        // realization of every K_ij with max |K_ij| > skip * max |K|; order up to rmax, chosen so that
        // the discarded Hankel singular values are below tol of their sum
        void identify(const retardation&, double tol=1.0e-4, int rmax=12, double skip=1.0e-6);

        double kernel(int p, double t) const;    // K of pair p at t

        // states x of pair p over a step h with the velocity going linearly from u0 to u1
        static void advance(const pair&, Eigen::VectorXcd &x, double u0, double u1, double h);
        // radiation force of pair p (the memory part, K * u) for the states x
        static double output(const pair &q, const Eigen::VectorXcd &x) {return (q.c.cwiseProduct(x)).sum().real();}
        int nstates() const;
        int ndof() const {return ndof_;}

        std::vector<pair> pairs;

    private:
        int ndof_ = 0;
    };
}

#endif
