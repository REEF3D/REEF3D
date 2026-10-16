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
// REEF3D::CORAL: interface of all loads on the multibody system (gravity, springs, hydrodynamics,
// mooring, power take-offs, aerodynamics, ...).
//
// eval() adds the generalized forces of the load for the current state of the model (positions,
// velocities and accelerations of the bodies). The time integrator calls it several times per step
// with trial states (Newton iterations), so eval() must not change the load.
// commit() is called once per converged step with the converged state; loads with a history
// (radiation memory, internal states, line dynamics) update it there.
// coupling() tells the Jacobian which bodies the load couples: global (any body may act on any other;
// the default), local (only the bodies listed by coupled_bodies() and the generalized coordinates
// listed by coupled_coordinates()), diagonal (each body only on itself, e.g. gravity). Local and
// diagonal loads let the integrators perturb several bodies at once.
// has_tangent()/tangent(): a load with an analytic tangent adds cK dF/dq + cC dF/dnu + cM dF/dnu'
// itself (q: configuration increments as the integrators take them: x += dx global, R = R exp(dtheta)
// body frame, coordinates += dq; nu: velocities, nu': accelerations); it is then left out of the
// finite-difference Jacobian and of the colouring.

#ifndef CORAL_LOAD_H_
#define CORAL_LOAD_H_

#include"coral/core/coral_types.h"
#include<Eigen/Sparse>
#include<string>
#include<vector>

namespace coral
{
    class model;

    enum class coupling_kind {global, local, diagonal};
    using triplet = Eigen::Triplet<double>;

    class load
    {
    public:
        virtual ~load() = default;

        virtual std::string name() const = 0;

        // called once before the first step (state of the model at t0)
        virtual void initialize(const model&, double /*t0*/) {}

        // add the generalized forces (size model::ndof()) of the current state at time t
        virtual void eval(const model&, double t, VecX &f) const = 0;

        // converged step from t-dt to t
        virtual void commit(const model&, double /*t*/, double /*dt*/) {}

        // bodies whose state changes the forces of the load on each other
        virtual coupling_kind coupling() const {return coupling_kind::global;}
        virtual void coupled_bodies(std::vector<int>&) const {}
        virtual void coupled_coordinates(std::vector<int>&) const {}

        // analytic tangent (rows and columns: generalized degrees of freedom of the model)
        virtual bool has_tangent() const {return false;}
        virtual void tangent(const model&, double /*t*/, double /*cK*/, double /*cC*/, double /*cM*/,
                             std::vector<triplet>&) const {}

        // potential energy of conservative loads (energy checks), 0 otherwise
        virtual double potential(const model&) const {return 0.0;}

        // output channels: names once, values per output
        virtual void channel_names(std::vector<std::string>&) const {}
        virtual void channel_values(const model&, double /*t*/, std::vector<double>&) const {}
    };
}

#endif
