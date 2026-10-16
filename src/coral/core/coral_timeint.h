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
// REEF3D::CORAL: time integration of the constrained multibody system.
//
// Lie group generalized-alpha method (Chung & Hulbert 1993; Arnold & Bruels 2007; Bruels, Cardona &
// Arnold 2012) on R3 x SO(3) per body, index-3 constraints with Lagrange multipliers:
//   M(q) dnu + h(q,nu) = f(q,nu,dnu,t) - G(q)^T lambda,   g(q) = 0      at t_{n+1}
//   q_{n+1}  = q_n o exp(h nu_n + h^2 (1/2 - beta) a_n + h^2 beta a_{n+1})
//   nu_{n+1} = nu_n + h (1 - gamma) a_n + h gamma a_{n+1}
//   (1 - alpha_m) a_{n+1} + alpha_m a_n = (1 - alpha_f) dnu_{n+1} + alpha_f dnu_n
// alpha_m, alpha_f, gamma, beta from the spectral radius rho_inf (1: no numerical damping, second order).
// With constraints use rho_inf < 1 (default 0.8): with rho_inf = 1 the algebraic components (lambda)
// are not damped, and fast gyroscopic motion of constrained bodies can become unstable.
//
// The unknowns dnu_{n+1} and lambda_{n+1} are solved with Newton's method. The Jacobian is formed
// by finite differences of the residual (loads may depend on positions, velocities and
// accelerations, e.g. added mass), the lambda columns exactly (G^T); it is kept over iterations and
// steps while Newton converges fast (modified Newton) and rebuilt otherwise.

#ifndef CORAL_TIMEINT_H_
#define CORAL_TIMEINT_H_

#include"coral/core/coral_types.h"
#include<vector>

namespace coral
{
    class model;

    struct timeint_settings
    {
        double rho_inf = 0.8;       // spectral radius at infinite frequency, [0,1]
        double tol = 1.0e-12;       // Newton: position (velocity) increments <= tol (1 + max |x| (|nu|))
        int max_iter = 30;
        bool full_newton = false;   // new Jacobian in every iteration
        double fd_eps = 1.0e-7;     // relative perturbation of the finite-difference Jacobian
    };

    struct timeint_stats
    {
        long steps = 0;
        long iterations = 0;
        long jacobians = 0;
        long retries = 0;
        long failures = 0;
        double max_constraint = 0.0;    // max |g| after the steps
    };

    class timeint
    {
    public:
        explicit timeint(const timeint_settings &s = timeint_settings());

        // consistent accelerations and multipliers at t0 (initializes the loads)
        void initialize(model&, double t0);

        // advance from t to t+h; on return the model holds the state at t+h.
        // false if Newton did not converge (the model is then back at t).
        bool step(model&, double t, double h);

        // Lagrange multipliers of the last step (constraint forces on the bodies: -G^T lambda)
        const VecX &multipliers() const {return lam_;}

        // generalized constraint forces -G^T lambda of the last step
        VecX reactions(const model&) const;

        const timeint_settings &settings() const {return set_;}
        double alpha_m() const {return am_;}
        double alpha_f() const {return af_;}
        double gamma() const {return gam_;}
        double beta() const {return beta_;}

        timeint_stats stats;

    private:
        bool converged(const model&, double ddnu, double h) const;
        bool stagnated(const model&, double ddnu, double h) const;
        void store(const model&);
        void restore(model&) const;
        void trial(model&, const VecX &dnu, double h) const;
        void residual(model&, double t1, double h, const VecX &u, VecX &r);
        void residual_initial(model&, double t0, const VecX &c, const VecX &u, VecX &r);
        void build_jacobian(model&, double t1, double h, const VecX &u, const VecX &r, bool initial, const VecX &c);

        timeint_settings set_;
        double am_, af_, gam_, beta_;
        double scale_ = 1.0;           // multiplier scaling (typical mass)

        std::vector<Vec3> xn_;
        std::vector<Mat3> Rn_;
        VecX nun_, dnun_, an_, lam_;

        MatX J_;
        Eigen::PartialPivLU<MatX> lu_;
        bool jac_ok_ = false;
        double jac_h_ = 0.0;

        VecX f_, in_, g_;
        MatX G_;
    };
}

#endif
