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

#include"coral/hydro/coral_body_drag.h"
#include"coral/core/coral_model.h"
#include<cmath>

namespace coral
{

body_drag::body_drag(const model &m, int body, const Vec3 &P, const Vec6 &Q, const std::string &name)
    : name_(name), b_(body), Q_(Q)
{
    s_ = m.body(body).to_body(P);
}

body_drag::body_drag(const model &m, int body, const Vec3 &P, const Vec3 &e, double Qe, const std::string &name)
    : name_(name), b_(body), e_(e.normalized()), Qe_(Qe), axis_(true)
{
    s_ = m.body(body).to_body(P);
}

Vec6 body_drag::force(const model &m) const
{
    const rigid_body &b = m.body(b_);
    Vec6 F = Vec6::Zero();
    if(axis_)
    {
        const double u = b.point_velocity(s_).dot(e_);
        F.head<3>() = -Qe_*std::fabs(u)*u*e_;
        return F;
    }
    Vec6 u;
    u.head<3>() = b.point_velocity(s_);
    u.tail<3>() = b.omega_global();
    for(int i=0; i<6; ++i)
        F(i) = -Q_(i)*std::fabs(u(i))*u(i);
    return F;
}

double body_drag::power(const model &m) const
{
    const rigid_body &b = m.body(b_);
    const Vec6 F = force(m);
    return F.head<3>().dot(b.point_velocity(s_)) + F.tail<3>().dot(b.omega_global());
}

void body_drag::eval(const model &m, double, VecX &f) const
{
    const Vec6 F = force(m);
    const rigid_body &b = m.body(b_);
    m.add_force(f, b_, b.point(s_), F.head<3>());
    m.add_moment(f, b_, F.tail<3>());
}

void body_drag::channel_names(std::vector<std::string> &c) const
{
    for(const char *x : {"Fx", "Fy", "Fz", "Mx", "My", "Mz", "P"})
        c.push_back(name_ + "_" + x);
}

void body_drag::channel_values(const model &m, double, std::vector<double> &v) const
{
    const Vec6 F = force(m);
    for(int i=0; i<6; ++i)
        v.push_back(F(i));
    v.push_back(power(m));
}

}
