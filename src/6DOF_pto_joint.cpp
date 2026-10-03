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

#include"6DOF_pto_joint.h"

void pto_joint_prismatic::initialize(const Eigen::Vector3d &axis, const Eigen::Vector3d &x_att,
                                     const Eigen::Vector3d &c, const Eigen::Matrix3d &R)
{
    // axis: inertial, x_att: inertial position of the attachment at t=0, c/R: initial body pose
    a  = axis.normalized();
    rb = R.transpose()*(x_att - c);
    x0 = x_att;

    J.head<3>() = a;
    J.tail<3>() = (R*rb).cross(a);
}

pto_state pto_joint_prismatic::state(const Eigen::Vector3d &c, const Eigen::Matrix3d &R,
                                     const Eigen::Vector3d &v, const Eigen::Vector3d &w, double t)
{
    const Eigen::Vector3d r = R*rb;

    J.head<3>() = a;
    J.tail<3>() = r.cross(a);

    pto_state s;

    s.q  = a.dot(c + r - x0);
    s.qd = a.dot(v) + J.tail<3>().dot(w);
    s.t  = t;

    return s;
}
