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

#include"coral/hydro/coral_member_morison.h"
#include"coral/core/coral_model.h"
#include<algorithm>
#include<cmath>

namespace coral
{

member_morison::member_morison(std::shared_ptr<const point_set> pts, const std::vector<std::pair<int,int>> &elements,
                               const member_morison_coef &c, std::shared_ptr<const wave_field> waves, double rho, double g,
                               const Vec3 &current, const std::string &name)
    : name_(name), pts_(std::move(pts)), elem_(elements), c_(c), waves_(std::move(waves)), rho_(rho), g_(g), current_(current)
{
    if(c_.Cm < 0.0)
        c_.Cm = 1.0 + c_.Ca;
    c_.n = std::max(1, c_.n);
    pts_->point_coupling(bodies_, coords_);
    for(const auto &e : elem_) {pidx_.push_back(e.first); pidx_.push_back(e.second);}
    std::sort(bodies_.begin(), bodies_.end());
    bodies_.erase(std::unique(bodies_.begin(), bodies_.end()), bodies_.end());
    bodies_.erase(std::remove(bodies_.begin(), bodies_.end(), -1), bodies_.end());
    std::sort(coords_.begin(), coords_.end());
    coords_.erase(std::unique(coords_.begin(), coords_.end()), coords_.end());
}

// sink(point a, point b, weight b, strip centre, strip force)
template<class Sink>
void member_morison::strips(const model &m, double t, Sink sink) const
{
    const double A = 0.25*M_PI*c_.D*c_.D;
    std::vector<Vec3> P, V, Ac;
    pts_->point_states(m, pidx_, P, V, Ac);
    for(size_t ie=0; ie<elem_.size(); ++ie)
    {
        const auto &el = elem_[ie];
        const Vec3 &pa = P[2*ie], &va = V[2*ie], &aa = Ac[2*ie], &pb = P[2*ie+1], &vb = V[2*ie+1], &ab = Ac[2*ie+1];
        const Vec3 d = pb - pa;
        const double Le = d.norm();
        if(Le <= 0.0)
            continue;
        const Vec3 e = d/Le;
        auto normal = [&](const Vec3 &v) {return Vec3(v - v.dot(e)*e);};
        for(int i=0; i<c_.n; ++i)
        {
            // the strip from s0 to s1; at the free surface only its wet part (elevation at the strip centre),
            // so that the load changes continuously as the surface moves along the member
            double s0 = (double)i/c_.n, s1 = (i + 1.0)/c_.n;
            const Vec3 pm = (1.0 - 0.5*(s0 + s1))*pa + 0.5*(s0 + s1)*pb;
            const double eta = waves_ ? waves_->elevation(pm(0), pm(1), t) : 0.0;
            const double z0 = pa(2) + s0*d(2), z1 = pa(2) + s1*d(2);
            if(z0 > eta && z1 > eta)
                continue;
            if(z0 > eta || z1 > eta)
            {
                const double sw = s0 + (eta - z0)/(z1 - z0)*(s1 - s0);   // where the surface cuts the strip
                if(z0 > eta) s0 = sw; else s1 = sw;
            }
            const double s = 0.5*(s0 + s1), L = Le*(s1 - s0);
            if(L <= 0.0)
                continue;
            const Vec3 p = (1.0 - s)*pa + s*pb;
            const Vec3 us = (1.0 - s)*va + s*vb, as = (1.0 - s)*aa + s*ab;
            Vec3 uw = Vec3::Zero(), aw = Vec3::Zero();
            double pr = 0.0;
            if(waves_)
                waves_->kinematics(p, t, uw, aw, pr);
            const Vec3 ur = normal(uw + current_ - us);
            Vec3 f = L*(rho_*A*c_.Cm*normal(aw) - rho_*A*c_.Ca*normal(as) + 0.5*rho_*c_.Cd*c_.D*ur.norm()*ur);
            if(c_.buoyancy)
                f(2) += rho_*g_*A*L;
            sink(el.first, el.second, s, p, f);
        }
    }
}

void member_morison::eval(const model &m, double t, VecX &f) const
{
    strips(m, t, [&](int a, int b, double s, const Vec3&, const Vec3 &F)
    {
        pts_->add_point_force(m, a, (1.0 - s)*F, f);
        pts_->add_point_force(m, b, s*F, f);
    });
}

void member_morison::total(const model &m, double t, Vec3 &F, Vec3 &M) const
{
    F.setZero();
    M.setZero();
    strips(m, t, [&](int, int, double, const Vec3 &p, const Vec3 &f) {F += f; M += p.cross(f);});
}

void member_morison::channel_names(std::vector<std::string> &c) const
{
    for(const char *d : {"Fx", "Fy", "Fz", "Mx", "My", "Mz"})
        c.push_back(name_ + "_" + d);
}

void member_morison::channel_values(const model &m, double t, std::vector<double> &v) const
{
    Vec3 F, M;
    total(m, t, F, M);
    for(int i=0; i<3; ++i) v.push_back(F(i));
    for(int i=0; i<3; ++i) v.push_back(M(i));
}

}
