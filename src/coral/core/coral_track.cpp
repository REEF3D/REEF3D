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

#include"coral/core/coral_track.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_rotation.h"
#include<algorithm>

namespace coral
{

void body_track::record(const model &m, double t)
{
    std::vector<body_snapshot> s;
    for(int b : bodies_)
    {
        const rigid_body &rb = m.body(b);
        body_snapshot q;
        q.x = rb.x;
        q.v = rb.v;
        q.a = rb.vd;
        q.R = rb.R;
        q.w = rb.R*rb.w;
        q.alpha = rb.R*rb.wd;
        s.push_back(q);
    }
    if(!t_.empty() && t <= t_.back())
    {
        // the same time again: replace
        s_.back() = s;
        return;
    }
    t_.push_back(t);
    s_.push_back(s);
}

body_snapshot body_track::state(int k, double t) const
{
    if(t_.empty())
        return body_snapshot{};
    if(t <= t_.front())
        return s_.front()[k];
    if(t >= t_.back())
        return s_.back()[k];
    const size_t j = std::upper_bound(t_.begin(), t_.end(), t) - t_.begin();
    const double ta = t_[j-1], tb = t_[j], h = tb - ta;
    const body_snapshot &A = s_[j-1][k], &B = s_[j][k];
    const double s = (t - ta)/h, s2 = s*s, s3 = s2*s;
    // Hermite basis and its derivative (d/dt = d/ds / h)
    const double h00 = 2.0*s3 - 3.0*s2 + 1.0, h10 = s3 - 2.0*s2 + s, h01 = -2.0*s3 + 3.0*s2, h11 = s3 - s2;
    const double d00 = (6.0*s2 - 6.0*s)/h, d10 = (3.0*s2 - 4.0*s + 1.0)/h, d01 = (-6.0*s2 + 6.0*s)/h, d11 = (3.0*s2 - 2.0*s)/h;
    body_snapshot q;
    q.x = h00*A.x + h10*h*A.v + h01*B.x + h11*h*B.v;
    q.v = h00*A.v + h10*h*A.a + h01*B.v + h11*h*B.a;
    q.a = d00*A.v + d10*h*A.a + d01*B.v + d11*h*B.a;
    const Vec3 thb = log_so3(B.R*A.R.transpose());
    const Vec3 th = h10*h*A.w + h01*thb + h11*h*B.w;
    q.R = exp_so3(th)*A.R;
    q.w = h00*A.w + h10*h*A.alpha + h01*B.w + h11*h*B.alpha;
    q.alpha = d00*A.w + d10*h*A.alpha + d01*B.w + d11*h*B.alpha;
    return q;
}

}
