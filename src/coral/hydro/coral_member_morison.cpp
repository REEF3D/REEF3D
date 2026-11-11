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
                               const member_morison_coef &c, std::shared_ptr<const wave_kinematics> waves, double rho, double g,
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

void member_morison::init_waterline(const std::vector<Vec3> &P) const
{
    if(wl_init_)
        return;
    wl_init_ = true;
    if(!surface_point())
        return;
    for(size_t ie=0; ie<elem_.size(); ++ie)
    {
        const double za = P[2*ie](2), zb = P[2*ie+1](2);
        if((za <= 0.0 && zb > 0.0) || (zb <= 0.0 && za > 0.0))
            wl_.push_back((int)ie);
    }
}

void member_morison::stencil_gradient(int i6, const Vec3 &pc, double t, const Vec3 &uc, Mat3 &G, Vec3 &slope) const
{
    G.setZero();
    slope.setZero();
    const double h = step();
    for(int q=0; q<3; ++q)
    {
        if(q==1 && c_.plane)
            continue;
        Vec3 up, um, du;
        double pr;
        const Vec3 pp = pc + h*Vec3::Unit(q), pm = pc - h*Vec3::Unit(q);
        waves_->kinematics_at(this, i6 + 2*q, pp, t, up, du, pr);
        waves_->kinematics_at(this, i6 + 2*q + 1, pm, t, um, du, pr);
        const double ep = waves_->elevation_at(this, i6 + 2*q, pp, t), em = waves_->elevation_at(this, i6 + 2*q + 1, pm, t);
        const double bp = waves_->bed_at(this, i6 + 2*q, pp), bm = waves_->bed_at(this, i6 + 2*q + 1, pm);
        const bool wp = pp(2) < ep && !(pp(2) < bp), wm = pm(2) < em && !(pm(2) < bm);
        if(wp && wm) G.col(q) = (up - um)/(2.0*h);
        else if(wp) G.col(q) = (up - uc)/h;
        else if(wm) G.col(q) = (uc - um)/h;
        if(q < 2)
            slope(q) = (ep - em)/(2.0*h);
    }
}

// sink(point a, point b, weight b, point, force)
template<class Sink>
void member_morison::strips(const model &m, double t, Sink sink, bool points) const
{
    const double A = 0.25*M_PI*c_.D*c_.D;
    std::vector<Vec3> P, V, Ac;
    pts_->point_states(m, pidx_, P, V, Ac);
    init_waterline(P);
    const int ns = (int)elem_.size()*c_.n;
    const bool st = stencil();
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
            const int ik = (int)ie*c_.n + i;       // index of the strip for a sampled wave source
            const double eta = waves_ ? waves_->elevation_at(this, ik, pm, t) : 0.0;
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
                waves_->kinematics_at(this, ik, p, t, uw, aw, pr);
            Vec3 fx = Vec3::Zero();
            if(higher() && waves_)
            {
                Mat3 G;
                Vec3 sl;
                if(st) stencil_gradient(ns + 6*ik, pm, t, uw, G, sl);
                else waves_->gradients(p, t, G, sl);
                const Vec3 uf = uw + current_;
                aw += G*uf;                                            // the particle acceleration
                if(c_.model != slender_model::lagrangian)
                    fx = rho_*A*c_.Ca*e.dot(G*e)*normal(uf);           // axial divergence
            }
            const Vec3 ur = normal(uw + current_ - us);
            Vec3 f = L*(rho_*A*c_.Cm*normal(aw) - rho_*A*c_.Ca*normal(as) + 0.5*rho_*c_.Cd*c_.D*ur.norm()*ur + fx);
            if(c_.buoyancy)
                f(2) += rho_*g_*A*L;
            sink(el.first, el.second, s, p, f);
        }
    }
    if(!points || !surface_point() || !waves_)
        return;
    // the point force at the free surface (rainey) or the still water level (kf) of each member crossing z = 0
    const int w0 = ns + (st ? 6*ns : 0);
    for(size_t iw=0; iw<wl_.size(); ++iw)
    {
        const int ie = wl_[iw];
        const auto &el = elem_[ie];
        const Vec3 &pa = P[2*ie], &pb = P[2*ie+1];
        const Vec3 d = pb - pa;
        const double Le = d.norm();
        if(Le <= 0.0 || std::fabs(d(2)) < 1.0e-9*Le)
            continue;
        const Vec3 e = d/Le;
        auto normal = [&](const Vec3 &v) {return Vec3(v - v.dot(e)*e);};
        const double sw0 = -pa(2)/d(2);
        if(sw0 < 0.0 || sw0 > 1.0)
            continue;
        const Vec3 p0 = pa + sw0*d;
        const int iwc = w0 + 7*(int)iw;
        Vec3 u0, du0, slope;
        double pr, eta;
        Mat3 G;
        if(st)
        {
            // sampled: z = 0 in a crest, the surface in a trough; the first-order Taylor step to the other level
            Vec3 us, dus;
            waves_->kinematics_at(this, iwc, p0, t, us, dus, pr);
            eta = waves_->elevation_at(this, iwc, p0, t);
            stencil_gradient(iwc + 1, p0, t, us, G, slope);
            u0 = eta < 0.0 ? Vec3(us - eta*G.col(2)) : us;
            du0 = dus;
            const Vec3 ue = eta < 0.0 ? us : Vec3(us + eta*G.col(2));
            if(c_.model == slender_model::rainey)
            {
                const double se = std::clamp((eta - pa(2))/d(2), 0.0, 1.0);
                const Vec3 un = normal(ue + current_);
                sink(el.first, el.second, se, pa + se*d, Vec3(-0.5*rho_*A*c_.Ca*un.squaredNorm()*normal(slope)));
                continue;
            }
        }
        else
        {
            // a wave source with gradients (linear waves): z = 0 kinematics, at z = eta in a trough
            waves_->kinematics(p0, t, u0, du0, pr);
            waves_->gradients(p0, t, G, slope);
            eta = waves_->elevation(p0(0), p0(1), t);
            if(c_.model == slender_model::rainey)
            {
                Vec3 ue = u0, due;
                if(eta < 0.0)
                {
                    Vec3 pe = p0;
                    pe(2) = eta;
                    waves_->kinematics(pe, t, ue, due, pr);
                }
                const double se = std::clamp((eta - pa(2))/d(2), 0.0, 1.0);
                const Vec3 un = normal(ue + current_);
                sink(el.first, el.second, se, pa + se*d, Vec3(-0.5*rho_*A*c_.Ca*un.squaredNorm()*normal(slope)));
                continue;
            }
        }
        // kf: F_psi = rho A Ca (4/g) |u_n|^2 (du/dt)_n at z = 0
        const Vec3 un = normal(u0 + current_);
        sink(el.first, el.second, sw0, p0, Vec3(rho_*A*c_.Ca*4.0/g_*un.squaredNorm()*normal(du0)));
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

void member_morison::kinematics_points(const model &m, const void *source, std::vector<Vec3> &pts) const
{
    // the centres of the whole strips; with a stencil six neighbours each; the waterline points (and their stencil)
    if(waves_.get() != source)
        return;
    std::vector<Vec3> P, V, Ac;
    pts_->point_states(m, pidx_, P, V, Ac);
    init_waterline(P);
    const bool st = stencil();
    const double h = step();
    auto six = [&](const Vec3 &c)
    {
        for(int q=0; q<3; ++q)
        {
            pts.push_back(c + h*Vec3::Unit(q));
            pts.push_back(c - h*Vec3::Unit(q));
        }
    };
    std::vector<Vec3> centres;
    for(size_t ie=0; ie<elem_.size(); ++ie)
        for(int i=0; i<c_.n; ++i)
        {
            const double s = (i + 0.5)/c_.n;
            centres.push_back((1.0 - s)*P[2*ie] + s*P[2*ie+1]);
        }
    pts.insert(pts.end(), centres.begin(), centres.end());
    if(st)
        for(const Vec3 &c : centres)
            six(c);
    if(surface_point())
        for(int ie : wl_)
        {
            const Vec3 &pa = P[2*ie], &pb = P[2*ie+1];
            const double sw0 = std::clamp(-pa(2)/(pb(2) - pa(2)), 0.0, 1.0);
            const Vec3 p0 = pa + sw0*(pb - pa);
            pts.push_back(p0);
            if(st) six(p0);
            else for(int j=0; j<6; ++j) pts.push_back(p0);
        }
}

void member_morison::total(const model &m, double t, Vec3 &F, Vec3 &M) const
{
    F.setZero();
    M.setZero();
    strips(m, t, [&](int, int, double, const Vec3 &p, const Vec3 &f) {F += f; M += p.cross(f);});
}

Vec3 member_morison::point_force(const model &m, double t) const
{
    Vec3 Ft = Vec3::Zero(), Fs = Vec3::Zero();
    strips(m, t, [&](int, int, double, const Vec3&, const Vec3 &f) {Ft += f;});
    strips(m, t, [&](int, int, double, const Vec3&, const Vec3 &f) {Fs += f;}, false);
    return Ft - Fs;
}

void member_morison::channel_names(std::vector<std::string> &c) const
{
    for(const char *d : {"Fx", "Fy", "Fz", "Mx", "My", "Mz"})
        c.push_back(name_ + "_" + d);
    if(surface_point())
        for(const char *d : {"Fpx", "Fpy", "Fpz"})
            c.push_back(name_ + "_" + d);
}

void member_morison::channel_values(const model &m, double t, std::vector<double> &v) const
{
    Vec3 F, M;
    total(m, t, F, M);
    for(int i=0; i<3; ++i) v.push_back(F(i));
    for(int i=0; i<3; ++i) v.push_back(M(i));
    if(surface_point())
    {
        const Vec3 Fp = point_force(m, t);
        for(int i=0; i<3; ++i) v.push_back(Fp(i));
    }
}

}
