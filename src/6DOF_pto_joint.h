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

#ifndef SIXDOF_PTO_JOINT_H_
#define SIXDOF_PTO_JOINT_H_

#include"6DOF_pto.h"
#include<Eigen/Dense>

// Prismatic PTO joint between the body and the ground (X 501).
// Axis a fixed in the inertial frame, attachment point x_att on the body (body-fixed offset r_b
// from the CoG), reference x0 = initial attachment point:
//
//   q  = a . (c + R r_b - x0)
//   qd = a . (v + w x r) = J^T [v; w],   J = [a; r x a],   r = R r_b
//
// The PTO force F acts along a at the attachment point: load on the body F J (force F a,
// moment about the CoG F r x a). Heaving point absorber: a = (0,0,1), x_att = CoG.

class pto_joint_prismatic
{
public:
    void initialize(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Matrix3d&);

    pto_state state(const Eigen::Vector3d&, const Eigen::Matrix3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, double);

    // generalised direction of the last state() call
    Eigen::Matrix<double,6,1> J = Eigen::Matrix<double,6,1>::Zero();

    Eigen::Vector3d a  = Eigen::Vector3d(0.0,0.0,1.0);
    Eigen::Vector3d rb = Eigen::Vector3d::Zero();
    Eigen::Vector3d x0 = Eigen::Vector3d::Zero();
};

#endif
