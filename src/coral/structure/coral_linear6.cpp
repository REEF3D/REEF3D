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

#include"coral/structure/coral_linear6.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_rotation.h"

namespace coral
{

linear6::linear6(const model &m, int body, const Vec3 &P, const Mat6 &K, const Mat6 &C, const Vec6 &F0, const std::string &name)
    : name_(name), b_(body), K_(K), C_(C), F0_(F0)
{
    const rigid_body &b = m.body(body);
    s_ = b.to_body(P);
    p0_ = P;
    R0_ = b.R;
}

Vec6 linear6::displacement(const model &m) const
{
    const rigid_body &b = m.body(b_);
    Vec6 xi;
    xi.head<3>() = b.point(s_) - p0_;
    xi.tail<3>() = log_so3(b.R*R0_.transpose());
    return xi;
}

Vec6 linear6::velocity(const model &m) const
{
    const rigid_body &b = m.body(b_);
    Vec6 v;
    v.head<3>() = b.point_velocity(s_);
    v.tail<3>() = b.omega_global();
    return v;
}

void linear6::eval(const model &m, double, VecX &f) const
{
    const Vec6 F = F0_ - K_*displacement(m) - C_*velocity(m);
    const rigid_body &b = m.body(b_);
    m.add_force(f, b_, b.point(s_), F.head<3>());
    m.add_moment(f, b_, F.tail<3>());
}

double linear6::potential(const model &m) const
{
    const Vec6 xi = displacement(m);
    return 0.5*xi.dot(0.5*(K_ + K_.transpose())*xi) - F0_.dot(xi);
}

}
