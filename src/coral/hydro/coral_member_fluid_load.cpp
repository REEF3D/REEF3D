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

#include"coral/hydro/coral_member_fluid_load.h"
#include"coral/core/coral_model.h"
#include<algorithm>
#include<cmath>

namespace coral
{

member_fluid_load::member_fluid_load(const model &m, std::shared_ptr<const point_set> pts, const std::vector<std::pair<int,int>> &elements,
                                     const Vec3 &p1, const Vec3 &p2, int nstrips, const std::string &name)
    : name_(name), pts_(std::move(pts)), a_(p1), ns_(std::max(1, nstrips))
{
    L_ = (p2 - p1).norm();
    e_ = L_ > 0.0 ? Vec3((p2 - p1)/L_) : Vec3::UnitZ();
    pts_->point_coupling(bodies_, coords_);
    std::sort(bodies_.begin(), bodies_.end());
    bodies_.erase(std::unique(bodies_.begin(), bodies_.end()), bodies_.end());
    bodies_.erase(std::remove(bodies_.begin(), bodies_.end(), -1), bodies_.end());
    std::sort(coords_.begin(), coords_.end());
    coords_.erase(std::unique(coords_.begin(), coords_.end()), coords_.end());
    // the element at each strip centre (positions along the axis at the start)
    std::vector<int> idx;
    for(const auto &el : elements) {idx.push_back(el.first); idx.push_back(el.second);}
    std::vector<Vec3> P, V, A;
    pts_->point_states(m, idx, P, V, A);
    for(int k=0; k<ns_; ++k)
    {
        const double s = (k + 0.5)/ns_*L_;
        int best = 0;
        double wbest = 0.0, dbest = 1.0e300;
        for(size_t ie=0; ie<elements.size(); ++ie)
        {
            const double sa = (P[2*ie] - a_).dot(e_), sb = (P[2*ie+1] - a_).dot(e_);
            const double lo = std::min(sa, sb), hi = std::max(sa, sb);
            const double d = s < lo ? lo - s : (s > hi ? s - hi : 0.0);
            if(d < dbest)
            {
                dbest = d;
                best = (int)ie;
                wbest = std::fabs(sb - sa) > 0.0 ? std::min(1.0, std::max(0.0, (s - sa)/(sb - sa))) : 0.5;
            }
        }
        pa_.push_back(elements[best].first);
        pb_.push_back(elements[best].second);
        w_.push_back(wbest);
    }
}

void member_fluid_load::set_sample(double t, const VecX &F)
{
    if(F.size() != 3*ns_)
        return;
    VecX G = F;
    if(excitation)
    {
        if(F0_.size() == 0)
            F0_ = F;
        G -= F0_;
    }
    while(!ts_.empty() && ts_.back() >= t)
    {
        ts_.pop_back();
        Fs_.pop_back();
    }
    ts_.push_back(t);
    Fs_.push_back(G);
    while(ts_.size() > 3)
    {
        ts_.erase(ts_.begin());
        Fs_.erase(Fs_.begin());
    }
}

VecX member_fluid_load::forces(double t) const
{
    const int ns = (int)ts_.size();
    VecX F = VecX::Zero(3*ns_);
    if(ns == 0)
        return F;
    const int k = std::min(ns, std::max(0, order) + 1);
    for(int i=ns-k; i<ns; ++i)
    {
        double l = 1.0;
        for(int j=ns-k; j<ns; ++j)
            if(j != i)
                l *= (t - ts_[j])/(ts_[i] - ts_[j]);
        F += l*Fs_[i];
    }
    return F;
}

void member_fluid_load::eval(const model &m, double t, VecX &f) const
{
    if(ts_.empty())
        return;
    const VecX F = forces(t);
    for(int k=0; k<ns_; ++k)
    {
        const Vec3 Fk = F.segment<3>(3*k);
        pts_->add_point_force(m, pa_[k], (1.0 - w_[k])*Fk, f);
        pts_->add_point_force(m, pb_[k], w_[k]*Fk, f);
    }
}

void member_fluid_load::total(const model &m, double t, Vec3 &Ft, Vec3 &Mt) const
{
    Ft.setZero();
    Mt.setZero();
    if(ts_.empty())
        return;
    const VecX F = forces(t);
    std::vector<int> idx;
    for(int k=0; k<ns_; ++k) {idx.push_back(pa_[k]); idx.push_back(pb_[k]);}
    std::vector<Vec3> P, V, A;
    pts_->point_states(m, idx, P, V, A);
    for(int k=0; k<ns_; ++k)
    {
        const Vec3 Fk = F.segment<3>(3*k);
        const Vec3 p = (1.0 - w_[k])*P[2*k] + w_[k]*P[2*k+1];
        Ft += Fk;
        Mt += p.cross(Fk);
    }
}

void member_fluid_load::channel_names(std::vector<std::string> &c) const
{
    for(const char *d : {"Fx", "Fy", "Fz", "Mx", "My", "Mz"})
        c.push_back(name_ + "_" + d);
}

void member_fluid_load::channel_values(const model &m, double t, std::vector<double> &v) const
{
    Vec3 F, M;
    total(m, t, F, M);
    for(int i=0; i<3; ++i) v.push_back(F(i));
    for(int i=0; i<3; ++i) v.push_back(M(i));
}

// restart from a checkpoint: the samples of the strip loads
void member_fluid_load::save_state(std::vector<double> &v) const
{
    state_out o{v};
    o.vd(ts_);
    o.list(Fs_, [&](const VecX &x){o.m(x);});
    o.m(F0_);
}

bool member_fluid_load::load_state(const std::vector<double> &v, size_t &k)
{
    state_in in{v, k};
    in.vd(ts_);
    in.list(Fs_, [&](VecX &x){in.m(x);});
    in.m(F0_);
    k = in.k;
    return in.ok;
}

}
