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

#include"coral/hydro/coral_member_slam.h"
#include"coral/core/coral_model.h"
#include"coral/hydro/coral_sampled_kinematics.h"
#include<algorithm>
#include<cmath>

namespace coral
{

double slam_duration(const member_slam_coef &c)
{
    switch(c.model)
    {
        case slam_model::wo: return 13.0/32.0;
        case slam_model::goda: return 1.0;
        case slam_model::cw: return 2.0;
        case slam_model::custom: return c.tau_end;
    }
    return 0.0;
}

double slam_coefficient(const member_slam_coef &c, double tau)
{
    if(tau < 0.0 || tau > slam_duration(c))
        return 0.0;
    switch(c.model)
    {
        case slam_model::wo:
            if(tau <= 0.0)
                return 2.0*M_PI;
            if(tau <= 0.125)
                return 2.0*M_PI - 2.0*std::sqrt(tau)*std::atanh(std::sqrt(1.0 - 0.25*tau));
            else
            {
                const double tp = tau - 1.0/32.0, r = std::sqrt(6.0*tp);
                return M_PI/r - std::pow(8.0*tp/3.0, 0.25)*std::atanh(std::sqrt(std::max(0.0, 1.0 - tp*r)));
            }
        case slam_model::goda: return M_PI*(1.0 - tau);
        case slam_model::cw:
        {
            // s/D = tau/2
            const double sD = 0.5*tau;
            return 5.15*(1.0/(1.0 + 19.0*sD) + 0.107*sD);
        }
        case slam_model::custom: return c.constant ? c.cs0 : c.cs0*(1.0 - tau/c.tau_end);
    }
    return 0.0;
}

member_slam::member_slam(std::shared_ptr<const point_set> pts, const std::vector<std::pair<int,int>> &elements,
                         const member_slam_coef &c, std::shared_ptr<const wave_kinematics> waves, double rho, double g,
                         const std::string &name)
    : name_(name), pts_(std::move(pts)), elem_(elements), c_(c), waves_(std::move(waves)), rho_(rho), g_(g)
{
    c_.n = std::max(1, c_.n);
    std::sort(c_.times.begin(), c_.times.end());
    pts_->point_coupling(bodies_, coords_);
    for(const auto &e : elem_) {pidx_.push_back(e.first); pidx_.push_back(e.second);}
    std::sort(bodies_.begin(), bodies_.end());
    bodies_.erase(std::unique(bodies_.begin(), bodies_.end()), bodies_.end());
    bodies_.erase(std::remove(bodies_.begin(), bodies_.end(), -1), bodies_.end());
    std::sort(coords_.begin(), coords_.end());
    coords_.erase(std::unique(coords_.begin(), coords_.end()), coords_.end());
}

void member_slam::initialize(const model &m, double)
{
    find_waterline(m);
}

void member_slam::find_waterline(const model &m) const
{
    if(iwl_ >= 0)
        return;
    std::vector<Vec3> P, V, A;
    pts_->point_states(m, pidx_, P, V, A);
    for(size_t ie=0; ie<elem_.size() && iwl_ < 0; ++ie)
    {
        const double za = P[2*ie](2), zb = P[2*ie+1](2);
        if((za <= 0.0 && zb > 0.0) || (zb <= 0.0 && za > 0.0))
            iwl_ = (int)ie;
    }
}

Vec3 member_slam::waterline(const model &m) const
{
    // the member at z = 0 (horizontal position of the impact), z far above: a sampled fluid gives the surface velocity
    std::vector<Vec3> P, V, A;
    const int ie = std::max(iwl_, 0);
    pts_->point_states(m, {elem_[ie].first, elem_[ie].second}, P, V, A);
    const Vec3 d = P[1] - P[0];
    const double s = std::fabs(d(2)) > 0.0 ? std::clamp(-P[0](2)/d(2), 0.0, 1.0) : 0.0;
    Vec3 p = P[0] + s*d;
    p(2) = 1.0e3;
    return p;
}

void member_slam::kinematics_points(const model &m, const void *source, std::vector<Vec3> &pts) const
{
    if(waves_.get() != source)
        return;
    find_waterline(m);
    if(iwl_ < 0)
        return;
    pts.push_back(waterline(m));
}

double member_slam::celerity(const Vec3 &p) const
{
    if(c_.period <= 0.0)
        return 0.0;
    double h = c_.depth;
    if(h <= 0.0 && waves_)
        h = -waves_->bed_at(this, 0, p);
    if(!(h > 0.0) || !std::isfinite(h))
    {
        const wave_field *wf = dynamic_cast<const wave_field*>(waves_.get());
        if(const auto *sk = dynamic_cast<const sampled_kinematics*>(waves_.get())) wf = sk->background().get();
        if(wf) h = wf->depth();
    }
    const double w = 2.0*M_PI/c_.period;
    return (h > 0.0 && std::isfinite(h)) ? w/wave_field::wavenumber(w, h, g_) : g_/w;
}

void member_slam::fire(double t0, double eta_b, const Vec3 &uh, double cel)
{
    const double un = uh.norm();
    const double v = c_.V > 0.0 ? c_.V : c_.vfac*(c_.crest_velocity ? un : cel);
    if(v <= 0.0 || eta_b <= 0.0)
        return;
    ev_.push_back(event{t0, v, eta_b, un > 0.0 ? Vec3(uh/un) : Vec3::UnitX()});
    ++nev_;
}

void member_slam::commit(const model &m, double t, double dt)
{
    find_waterline(m);
    if(iwl_ < 0 || (!waves_ && c_.trigger != slam_trigger::time))
        return;
    // the impacts that are over
    const double R = 0.5*c_.D, T = slam_duration(c_);
    ev_.erase(std::remove_if(ev_.begin(), ev_.end(), [&](const event &e) {return t > e.t0 + T*R/e.V;}), ev_.end());
    const Vec3 p = waterline(m);
    Vec3 u = Vec3::Zero(), du;
    double pr;
    if(waves_)
    {
        waves_->kinematics_at(this, 0, p, t, u, du, pr);
        eta_ = waves_->elevation_at(this, 0, p, t);
    }
    // without waves: along x
    const Vec3 uh = waves_ ? Vec3(u(0), u(1), 0.0) : Vec3::UnitX();
    if(c_.trigger == slam_trigger::time)
    {
        // the impacts of the next step start at their times
        while(nt_ < (int)c_.times.size() && c_.times[nt_] < t + dt - 1.0e-9)
        {
            fire(std::max(c_.times[nt_], t), c_.eta_b > 0.0 ? c_.eta_b : eta_, uh, celerity(p));
            ++nt_;
        }
        return;
    }
    // breaking index: the fluid's flag or the kinematic criterion
    if(c_.trigger == slam_trigger::breaking)
        B_ = waves_->breaking_at(this, 0) ? 1.0 : 0.0;
    else
        B_ = uh.norm()/celerity(p);
    if(eta_ <= 0.0)
    {
        armed_ = true;
        seen_ = false;
    }
    else if(armed_)
    {
        if(c_.trigger == slam_trigger::breaking ? B_ > 0.5 : B_ >= c_.Bc)
            seen_ = true;
        // the crest passed the member at the last step
        if(nstep_ >= 2 && eta_ < eta1_ && eta1_ >= eta2_ && eta1_ > 0.0 && seen_)
        {
            fire(t, eta1_, uh1_, celerity(p));
            armed_ = false;
        }
    }
    eta2_ = eta1_;
    eta1_ = eta_;
    uh1_ = uh;
    ++nstep_;
}

// sink(point a, point b, weight b, point, force)
template<class Sink>
void member_slam::impacts(const model &m, double t, Sink sink) const
{
    if(ev_.empty())
        return;
    const double R = 0.5*c_.D;
    std::vector<Vec3> P, V, A;
    pts_->point_states(m, pidx_, P, V, A);
    for(const event &e : ev_)
    {
        const double cs = slam_coefficient(c_, e.V*(t - e.t0)/R);
        if(cs == 0.0)
            continue;
        const double q = rho_*R*e.V*e.V*cs;          // line force
        const double zl = (1.0 - c_.lambda)*e.eta_b, zu = e.eta_b;
        for(size_t ie=0; ie<elem_.size(); ++ie)
        {
            const Vec3 &pa = P[2*ie], &pb = P[2*ie+1];
            const Vec3 d = pb - pa;
            const double Le = d.norm();
            if(Le <= 0.0 || std::fabs(d(2)) < 1.0e-9*Le)
                continue;
            // the part of the element inside the impact zone
            double s0 = (zl - pa(2))/d(2), s1 = (zu - pa(2))/d(2);
            if(s0 > s1) std::swap(s0, s1);
            s0 = std::max(s0, 0.0);
            s1 = std::min(s1, 1.0);
            if(s1 <= s0)
                continue;
            const Vec3 ax = d/Le;
            Vec3 fn = e.dir - e.dir.dot(ax)*ax;
            const double nn = fn.norm();
            if(nn <= 0.0)
                continue;
            fn /= nn;
            for(int i=0; i<c_.n; ++i)
            {
                const double a = s0 + (s1 - s0)*i/c_.n, b = s0 + (s1 - s0)*(i + 1.0)/c_.n, s = 0.5*(a + b);
                // the line force per unit height of the zone, over the height of the sub-strip
                const double dz = std::fabs(d(2))*(b - a);
                sink(elem_[ie].first, elem_[ie].second, s, Vec3(pa + s*d), Vec3(q*dz*fn));
            }
        }
    }
}

void member_slam::eval(const model &m, double t, VecX &f) const
{
    impacts(m, t, [&](int a, int b, double s, const Vec3&, const Vec3 &F)
    {
        pts_->add_point_force(m, a, (1.0 - s)*F, f);
        pts_->add_point_force(m, b, s*F, f);
    });
}

void member_slam::total(const model &m, double t, Vec3 &F, Vec3 &M) const
{
    F.setZero();
    M.setZero();
    impacts(m, t, [&](int, int, double, const Vec3 &p, const Vec3 &f) {F += f; M += p.cross(f);});
}

void member_slam::channel_names(std::vector<std::string> &c) const
{
    for(const char *d : {"Fx", "Fy", "Fz", "Mx", "My", "Mz", "eta", "B", "n"})
        c.push_back(name_ + "_" + d);
}

void member_slam::channel_values(const model &m, double t, std::vector<double> &v) const
{
    Vec3 F, M;
    total(m, t, F, M);
    for(int i=0; i<3; ++i) v.push_back(F(i));
    for(int i=0; i<3; ++i) v.push_back(M(i));
    v.push_back(eta_);
    v.push_back(B_);
    v.push_back(nev_);
}

void member_slam::save_state(std::vector<double> &v) const
{
    v.push_back((double)ev_.size());
    for(const event &e : ev_)
    {
        v.push_back(e.t0); v.push_back(e.V); v.push_back(e.eta_b);
        for(int i=0; i<3; ++i) v.push_back(e.dir(i));
    }
    for(double x : {(double)nev_, armed_ ? 1.0 : 0.0, seen_ ? 1.0 : 0.0, (double)nt_, eta1_, eta2_, B_, eta_,
                    uh1_(0), uh1_(1), uh1_(2), (double)nstep_})
        v.push_back(x);
}

bool member_slam::load_state(const std::vector<double> &v, size_t &k)
{
    if(k >= v.size())
        return false;
    const int ne = (int)v[k++];
    if(ne < 0 || k + 6*(size_t)ne + 12 > v.size())
        return false;
    ev_.clear();
    for(int i=0; i<ne; ++i, k+=6)
        ev_.push_back(event{v[k], v[k+1], v[k+2], Vec3(v[k+3], v[k+4], v[k+5])});
    nev_ = (int)v[k++];
    armed_ = v[k++] > 0.5;
    seen_ = v[k++] > 0.5;
    nt_ = (int)v[k++];
    eta1_ = v[k++]; eta2_ = v[k++]; B_ = v[k++]; eta_ = v[k++];
    uh1_ = Vec3(v[k], v[k+1], v[k+2]);
    k += 3;
    nstep_ = (int)v[k++];
    return true;
}

}
