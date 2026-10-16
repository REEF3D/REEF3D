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
// REEF3D::CORAL: the types used across CORAL (core layer: no REEF3D headers, no MPI).

#ifndef CORAL_TYPES_H_
#define CORAL_TYPES_H_

#include<Eigen/Dense>
#include<cstddef>

namespace coral
{
    using index = std::ptrdiff_t;

    using Vec3 = Eigen::Vector3d;
    using Vec6 = Eigen::Matrix<double,6,1>;
    using VecX = Eigen::VectorXd;
    using Mat3 = Eigen::Matrix3d;
    using Mat6 = Eigen::Matrix<double,6,6>;
    using MatX = Eigen::MatrixXd;
    using Quat = Eigen::Quaterniond;
}

#endif
