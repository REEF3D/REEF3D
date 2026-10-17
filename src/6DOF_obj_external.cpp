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

// External dynamics of a 6DOF body (6DOF_external.h): the state comes from the external solver, the
// loads of the fluid go to it.

#include"6DOF_obj.h"
#include"lexer.h"

void sixdof_obj::external_set_state(const Eigen::Vector3d &c, const Eigen::Matrix3d &R, const Eigen::Vector3d &v, const Eigen::Vector3d &w)
{
    // CoG, quaternion (R = E G^T, e = (e0, e1, e2, e3), sign kept continuous), momentum and the body
    // frame angular momentum; the velocities themselves are kept for update_fbvel
    Eigen::Quaterniond q(R);
    q.normalize();
    Eigen::Vector4d e(q.w(), q.x(), q.y(), q.z());
    
    if(e.dot(rb.e)<0.0)
    e = -e;
    
    rb.c = c;
    rb.e = e;
    rb.quat_matrices();
    rb.p = rb.mass*v;
    rb.h = rb.I*(rb.R.transpose()*w);
    rb.dc = v;
    
    ext_v_ = v;
    ext_w_ = w;
    ext_on_ = true;
}

void sixdof_obj::external_fluid_load(lexer *p, Eigen::Matrix<double,6,1> &F) const
{
    // the fluid load of the last force evaluation without the weight (Xe .. Ne include gravity)
    F(0) = Xe - p->W20*Mass_fb;
    F(1) = Ye - p->W21*Mass_fb;
    F(2) = Ze - p->W22*Mass_fb;
    F(3) = Ke;
    F(4) = Me;
    F(5) = Ne;
}

void sixdof_obj::external_properties(double &m, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Matrix3d &I) const
{
    m = Mass_fb;
    c = c_;
    R = R_;
    I = I_;
}
