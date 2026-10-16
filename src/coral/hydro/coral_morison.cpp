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

#include"coral/hydro/coral_morison.h"
#include"coral/core/coral_model.h"
#include<cmath>

namespace coral
{

morison_load::morison_load(const model &m, int body, const std::vector<morison_member> &members,
                           std::shared_ptr<const wave_field> waves, double rho, const Vec3 &current, const std::string &name)
    : name_(name), b_(body), waves_(std::move(waves)), rho_(rho), current_(current)
{
    const rigid_body &b = m.body(body);
    for(const auto &mb : members)
    {
        const Vec3 d = mb.p2 - mb.p1;
        const double L = d.norm();
        if(L <= 0.0 || mb.n < 1)
            continue;
        for(int i=0; i<mb.n; ++i)
        {
            strip s;
            s.s = b.to_body(mb.p1 + (i + 0.5)/mb.n*d);
            s.e = b.R.transpose()*(d/L);
            s.L = L/mb.n;
            s.D = mb.D;
            s.Cd = mb.Cd;
            s.Ca = mb.Ca;
            s.Cm = mb.Cm;
            strips_.push_back(s);
        }
    }
}

// force of one strip at its current position p; false if the strip is out of the water
bool morison_load::strip_force(const model &m, const strip &s, double t, Vec3 &p, Vec3 &f) const
{
    const rigid_body &b = m.body(b_);
    p = b.point(s.s);
    if(p(2) > (waves_ ? waves_->elevation(p(0), p(1), t) : 0.0))
        return false;
    const Vec3 e = b.R*s.e;
    Vec3 uw = Vec3::Zero(), aw = Vec3::Zero();
    double pr = 0.0;
    if(waves_)
        waves_->kinematics(p, t, uw, aw, pr);
    const Vec3 ub = b.point_velocity(s.s), ab = b.point_acceleration(s.s);
    auto normal = [&](const Vec3 &v) {return Vec3(v - v.dot(e)*e);};
    const Vec3 ur = normal(uw + current_ - ub);
    const double A = 0.25*M_PI*s.D*s.D;
    f = s.L*(rho_*A*s.Cm*normal(aw) - rho_*A*s.Ca*normal(ab) + 0.5*rho_*s.Cd*s.D*ur.norm()*ur);
    return true;
}

void morison_load::total(const model &m, double t, Vec3 &F, Vec3 &M) const
{
    F.setZero();
    M.setZero();
    Vec3 p, f;
    for(const strip &s : strips_)
        if(strip_force(m, s, t, p, f))
        {
            F += f;
            M += p.cross(f);
        }
}

void morison_load::eval(const model &m, double t, VecX &fg) const
{
    if(b_ < 0)
        return;
    Vec3 p, f;
    for(const strip &s : strips_)
        if(strip_force(m, s, t, p, f))
            m.add_force(fg, b_, p, f);
}

void morison_load::channel_names(std::vector<std::string> &c) const
{
    for(const char *d : {"Fx", "Fy", "Fz"})
        c.push_back(name_ + "_" + d);
}

void morison_load::channel_values(const model &m, double t, std::vector<double> &v) const
{
    Vec3 F, M;
    total(m, t, F, M);
    for(int i=0; i<3; ++i)
        v.push_back(F(i));
}

}
