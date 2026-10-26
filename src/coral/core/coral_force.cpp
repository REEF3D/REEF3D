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

#include"coral/core/coral_force.h"
#include"coral/core/coral_model.h"
#include<algorithm>

namespace coral
{

point_force::point_force(const model &m, int b, const Vec3 &p, const Vec3 &F_, const Vec3 &M_, double ramp, const std::string &name)
    : F(F_), M(M_), name_(name), b_(b), s_(b >= 0 ? m.body(b).to_body(p) : p), ramp_(ramp)
{
}

void point_force::eval(const model &m, double t, VecX &f) const
{
    if(b_ < 0 || t >= until)
        return;
    const double r = ramp_ > 0.0 ? std::min(1.0, std::max(0.0, t/ramp_)) : 1.0;
    m.add_force(f, b_, m.body(b_).point(s_), r*F);
    m.add_moment(f, b_, r*M);
}

}
