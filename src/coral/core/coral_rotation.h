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
// REEF3D::CORAL: rotations on SO(3).
// A rotation matrix R maps body-frame vectors to the global frame (x_global = R x_body).
// Rotation vectors psi: R = exp(psi~), |psi| = angle, psi/|psi| = axis.
// Variations are taken in the body frame: dR = R dtheta~, so d(R s) = -R s~ dtheta.

#ifndef CORAL_ROTATION_H_
#define CORAL_ROTATION_H_

#include"coral/core/coral_types.h"

namespace coral
{
    // skew-symmetric matrix: skew(a) b = a x b
    Mat3 skew(const Vec3 &a);

    // exponential map: rotation vector -> rotation matrix (Rodrigues)
    Mat3 exp_so3(const Vec3 &psi);

    // logarithm: rotation matrix -> rotation vector, |psi| in [0, pi]
    Vec3 log_so3(const Mat3 &R);

    // tangent operator T(psi): for R = exp(psi), dR = R (T(psi) dpsi)~
    Mat3 tangent_so3(const Vec3 &psi);

    // Euler angles (roll phi, pitch theta, yaw psi; R = Rz(psi) Ry(theta) Rx(phi)), for output and input
    Vec3 euler_from_matrix(const Mat3 &R);
    Mat3 matrix_from_euler(const Vec3 &rpy);

    // nearest rotation matrix (polar decomposition), against round-off drift
    Mat3 orthonormalize(const Mat3 &R);
}

#endif
