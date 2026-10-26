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

#include"coral/hydro/coral_strip_hydro.h"
#include"coral/core/coral_model.h"
#include<algorithm>
#include<cmath>

namespace coral
{

strip_hydro::strip_hydro(const model &m, int body, const std::vector<strip_member> &mem, std::shared_ptr<const wave_kinematics> w,
                         double rho, double g, const Vec3 &current, const std::string &name)
    : name_(name), b_(body), waves_(std::move(w)), rho_(rho), g_(g), current_(current)
{
    const rigid_body &b = m.body(body);
    for(const strip_member &mm : mem)
    {
        const Vec3 d = mm.p2 - mm.p1;
        const double L = d.norm();
        const Vec3 eg = d/L, e = b.R.transpose()*eg;
        const int n = std::max(1, mm.n);
        for(int i=0; i<n; ++i)
            strips_.push_back(strip{b.to_body(mm.p1 + (i + 0.5)/n*d), e, L/n, mm.D, mm.Cd, mm.Ca, mm.Cm, mm.buoyancy});
        if(mm.end1) ends_.push_back(end{b.to_body(mm.p1), -e, mm.D, mm.Ca_end, mm.Cd_end});
        if(mm.end2) ends_.push_back(end{b.to_body(mm.p2), e, mm.D, mm.Ca_end, mm.Cd_end});
    }
    for(const strip &st : strips_) p0_.push_back(b.point(st.s));
    for(const end &q : ends_) p0_.push_back(b.point(q.s));
}

void strip_hydro::kinematics(int idx, const Vec3 &p, double t, kin &k) const
{
    k.u.setZero(); k.du.setZero(); k.p = 0.0; k.eta = 0.0;
    if(!waves_)
        return;
    if(exact_kinematics)
    {
        k.eta = waves_->elevation_at(this, idx, p, t);
        waves_->kinematics_at(this, idx, p, t, k.u, k.du, k.p);
        return;
    }
    if(t != t_cache_)
    {
        cache_.resize(p0_.size());
        for(size_t i=0; i<p0_.size(); ++i)
        {
            cache_[i].eta = waves_->elevation_at(this, (int)i, p0_[i], t);
            waves_->kinematics_at(this, (int)i, p0_[i], t, cache_[i].u, cache_[i].du, cache_[i].p);
        }
        t_cache_ = t;
    }
    k = cache_[idx];
}

template<class Add> void strip_hydro::forces(const model &m, double t, Add &&add, double &V) const
{
    V = 0.0;
    const rigid_body &b = m.body(b_);
    for(int i=0; i<(int)strips_.size(); ++i)
    {
        const strip &s = strips_[i];
        const Vec3 c = b.point(s.s), e = b.R*s.e;
        kin k;
        kinematics(i, c, t, k);
        const double eta = k.eta;
        // wet fraction from the vertical extent of the strip, the wet part from its lower end
        const Vec3 h = 0.5*s.L*e;
        const Vec3 lo = h(2) > 0.0 ? Vec3(c - h) : Vec3(c + h), hi = h(2) > 0.0 ? Vec3(c + h) : Vec3(c - h);
        double f;
        if(hi(2) - lo(2) > 1.0e-12*s.L) f = std::clamp((eta - lo(2))/(hi(2) - lo(2)), 0.0, 1.0);
        else f = c(2) < eta ? 1.0 : 0.0;
        if(f <= 0.0)
            continue;
        const Vec3 p = lo + 0.5*f*(hi - lo);
        const double A = 0.25*M_PI*s.D*s.D;
        V += A*f*s.L;
        Vec3 F = Vec3::Zero();
        if(s.buoy)
            F(2) += rho_*g_*A*f*s.L;
        const Vec3 uw = k.u, aw = k.du;
        const Vec3 sp = b.to_body(p);
        const Vec3 ub = b.point_velocity(sp), ab = b.point_acceleration(sp);
        auto normal = [&](const Vec3 &v) {return Vec3(v - v.dot(e)*e);};
        const Vec3 ur = normal(uw + current_ - ub);
        F += f*s.L*(rho_*A*s.Cm*normal(aw) - rho_*A*s.Ca*normal(ab) + 0.5*rho_*s.Cd*s.D*ur.norm()*ur);
        add(p, F);
    }
    for(int k=0; k<(int)ends_.size(); ++k)
    {
        const end &q = ends_[k];
        const Vec3 p = b.point(q.s), nrm = b.R*q.n;
        const int idx = (int)strips_.size() + k;
        kin kk;
        kinematics(idx, p, t, kk);
        if(p(2) >= kk.eta)
            continue;
        const Vec3 uw = kk.u, aw = kk.du;
        const double pr = kk.p;
        const Vec3 ub = b.point_velocity(q.s), ab = b.point_acceleration(q.s);
        const double A = 0.25*M_PI*q.D*q.D, ma = q.Ca*rho_*q.D*q.D*q.D/6.0;
        const double ur = (uw + current_ - ub).dot(nrm);
        const Vec3 F = -rho_*pr*A*nrm + ma*(aw - ab).dot(nrm)*nrm + 0.5*rho_*q.Cd*A*std::fabs(ur)*ur*nrm;
        add(p, F);
    }
}

void strip_hydro::eval(const model &m, double t, VecX &fg) const
{
    if(b_ < 0)
        return;
    double V;
    forces(m, t, [&](const Vec3 &p, const Vec3 &F) {m.add_force(fg, b_, p, F);}, V);
}

void strip_hydro::total(const model &m, double t, Vec3 &F, Vec3 &M, double &V) const
{
    F.setZero();
    M.setZero();
    const Vec3 x = m.body(b_).x;
    forces(m, t, [&](const Vec3 &p, const Vec3 &f) {F += f; M += (p - x).cross(f);}, V);
}

void strip_hydro::kinematics_points(const model &m, const void *source, std::vector<Vec3> &pts) const
{
    if(b_ < 0 || waves_.get() != source)
        return;
    if(!exact_kinematics)
    {
        pts.insert(pts.end(), p0_.begin(), p0_.end());
        return;
    }
    const rigid_body &b = m.body(b_);
    for(const strip &s : strips_) pts.push_back(b.point(s.s));
    for(const end &q : ends_) pts.push_back(b.point(q.s));
}

void strip_hydro::channel_names(std::vector<std::string> &c) const
{
    for(const char *d : {"Fx", "Fy", "Fz", "Mx", "My", "Mz", "V"})
        c.push_back(name_ + "_" + d);
}

void strip_hydro::channel_values(const model &m, double t, std::vector<double> &v) const
{
    Vec3 F, M;
    double V;
    total(m, t, F, M, V);
    for(int i=0; i<3; ++i) v.push_back(F(i));
    for(int i=0; i<3; ++i) v.push_back(M(i));
    v.push_back(V);
}

}
