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

#include"coral/structure/coral_spring.h"
#include"coral/core/coral_model.h"
#include<algorithm>

namespace coral
{

spring_damper::spring_damper(const model &m, int a, const Vec3 &pa, int b, const Vec3 &pb,
                             double k, double c, double L0, double T0, bool tension_only, const std::string &name)
    : name_(name), a_(a), b_(b), k_(k), c_(c), L0_(L0), T0_(T0), tension_only_(tension_only)
{
    sa_ = m.body(a).to_body(pa);
    sb_ = m.body(b).to_body(pb);
    if(L0_ < 0.0)
        L0_ = (pb - pa).norm();
}

double spring_damper::length(const model &m) const
{
    return (m.body(b_).point(sb_) - m.body(a_).point(sa_)).norm();
}

double spring_damper::tension(const model &m) const
{
    const Vec3 pa = m.body(a_).point(sa_), pb = m.body(b_).point(sb_);
    const Vec3 d = pb - pa;
    const double L = d.norm();
    if(L < 1.0e-14)
        return 0.0;
    const Vec3 e = d/L;
    const double Ld = (m.body(b_).point_velocity(sb_) - m.body(a_).point_velocity(sa_)).dot(e);
    double T = T0_ + k_*(L - L0_) + c_*Ld;
    if(tension_only_)
        T = std::max(T, 0.0);
    return T;
}

void spring_damper::eval(const model &m, double, VecX &f) const
{
    const Vec3 pa = m.body(a_).point(sa_), pb = m.body(b_).point(sb_);
    const Vec3 d = pb - pa;
    const double L = d.norm();
    if(L < 1.0e-14)
        return;
    const Vec3 e = d/L;
    const double T = tension(m);
    m.add_force(f, a_, pa, T*e);
    m.add_force(f, b_, pb, -T*e);
}

double spring_damper::potential(const model &m) const
{
    // conservative part (no damper); for a slack tension-only line the energy of the slack range is 0
    const double s = length(m) - L0_;
    if(tension_only_)
    {
        if(k_ <= 0.0)
            return T0_ > 0.0 ? T0_*s : 0.0;
        const double s0 = -T0_/k_;    // T = 0 at s0
        if(s <= s0)
            return 0.0;
        return 0.5*k_*(s - s0)*(s - s0);
    }
    return T0_*s + 0.5*k_*s*s;
}

void spring_damper::channel_names(std::vector<std::string> &c) const
{
    c.push_back(name_ + "_length");
    c.push_back(name_ + "_tension");
}

void spring_damper::channel_values(const model &m, double, std::vector<double> &v) const
{
    v.push_back(length(m));
    v.push_back(tension(m));
}

}
