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

#include"coral/core/coral_points.h"
#include"coral/core/coral_model.h"

namespace coral
{

int body_points::add(const model &m, int b, const Vec3 &p)
{
    body_.push_back(b);
    s_.push_back(b >= 0 ? m.body(b).to_body(p) : p);
    return (int)body_.size() - 1;
}

void body_points::point_state(const model &m, int k, Vec3 &p, Vec3 &v, Vec3 &a) const
{
    const int b = body_[k];
    if(b < 0)
    {
        p = s_[k];
        v.setZero();
        a.setZero();
        return;
    }
    const rigid_body &B = m.body(b);
    p = B.point(s_[k]);
    v = B.point_velocity(s_[k]);
    a = B.point_acceleration(s_[k]);
}

void body_points::add_point_force(const model &m, int k, const Vec3 &F, VecX &f) const
{
    const int b = body_[k];
    if(b < 0)
        return;
    m.add_force(f, b, m.body(b).point(s_[k]), F);
}

}
