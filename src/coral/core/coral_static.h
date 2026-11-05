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
// REEF3D::CORAL: static equilibrium of the model at time t: the loads balance with the constraint
// forces, velocities and accelerations zero:
//   f(q, 0, 0, t) - G(q)^T lambda = 0,   g(q) = 0
// Newton on the configuration (increments dx global, dtheta body frame per body: x += dx,
// R = R exp(dtheta)) with the coloured finite-difference Jacobian (coral_fdjac.h). Full Newton steps
// first; if they diverge, Newton again from the start configuration with steps damped by the natural
// monotonicity test (Deuflhard): a step is accepted when the simplified Newton increment shrinks.
// Converged when the increments are below the tolerance. Strongly nonlinear problems (large rotations)
// are solved in load steps by the caller (change the load, solve again).

#ifndef CORAL_STATIC_H_
#define CORAL_STATIC_H_

#include"coral/core/coral_types.h"

namespace coral
{
    class model;

    struct static_settings
    {
        int max_iter = 60;
        double tol = 1.0e-9;        // increments of positions and rotations <= tol (1 + max |x|)
        double fd_eps = 1.0e-7;
        int dense_max = 400;        // dense LU up to this many unknowns, sparse LU above
    };

    struct static_result
    {
        bool converged = false;
        int iterations = 0;
        double residual = 0.0;      // max |r| of the force rows at the end
    };

    static_result static_solve(model&, double t, const static_settings& = static_settings());
}

#endif
