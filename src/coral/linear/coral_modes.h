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
// REEF3D::CORAL: linear analysis of the model about its current configuration (velocities zero):
//   M dq'' + C dq' + K dq = 0  on the generalized coordinates (6 per body: position global, rotation
//   body frame), restricted to the constraint manifold (null space of G).
// M: mass of the bodies and acceleration-dependent loads (added mass A_inf), C: velocity derivatives
// of the loads (damping), K: derivatives of the loads and of the constraint forces G^T lambda (lambda
// of the static balance in the current configuration) with respect to the configuration; all by
// finite differences. Radiation memory is not part of the linear model (A_inf and the instantaneous
// damping only). Undamped modes from K phi = w^2 M phi, damped ones from the first-order form.

#ifndef CORAL_MODES_H_
#define CORAL_MODES_H_

#include"coral/core/coral_types.h"

namespace coral
{
    class model;

    struct linear_model
    {
        MatX M, C, K;      // full, 6 per body
        MatX N;            // basis of the constraint null space (columns)
        MatX Mr, Cr, Kr;   // reduced: N^T X N
    };

    struct mode_result
    {
        VecX frequency;    // [Hz], ascending
        VecX damping;      // damping ratio (damped modes), 0 for undamped
        MatX shape;        // columns: modes in the full coordinates (6 per body), largest entry 1
        bool ok = false;
    };

    linear_model linearize(model&, double t);

    // the first n modes (n <= number of free coordinates); damped: eigenvalues of the first-order form
    mode_result modes(model&, double t, int n, bool damped=false);
}

#endif
