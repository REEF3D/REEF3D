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

#include"coral/hydro/coral_fluid_load.h"
#include"coral/core/coral_model.h"
#include<algorithm>
#include<cmath>

namespace coral
{

fluid_load::fluid_load(const model &m, const std::vector<int> &bodies, const std::string &name)
    : name_(name), map_(bodies)
{
    for(int b : map_)
        bname_.push_back(m.body(b).name);
    s_.assign(map_.size(), Vec3::Zero());
    A_ = MatX::Zero(6*nbody(), 6*nbody());
    As_ = MatX::Zero(6*nbody(), 6*nbody());
}

void fluid_load::set_artificial_added_mass(int k, const Mat6 &As)
{
    if(k >= 0 && k < nbody())
        As_.block<6,6>(6*k, 6*k) = As;
}

void fluid_load::set_sample(double t, const VecX &F, const MatX &A, bool provisional, const VecX &a)
{
    const int n = 6*nbody();
    if(F.size() != n)
        return;
    // with the artificial added mass: the sample carries its compensation A_s a_k
    VecX G = F;
    if(a.size() == n)
        G += As_*a;
    // a new sample replaces the provisional one and any sample at or after its time
    if(prov_)
    {
        ts_.pop_back();
        Fs_.pop_back();
        prov_ = false;
    }
    while(!ts_.empty() && ts_.back() >= t)
    {
        ts_.pop_back();
        Fs_.pop_back();
    }
    ts_.push_back(t);
    Fs_.push_back(G);
    prov_ = provisional;
    // keep the samples the polynomials need (3 committed + 1 provisional)
    while(ts_.size() > 4)
    {
        ts_.erase(ts_.begin());
        Fs_.erase(Fs_.begin());
    }
    if(A.rows() == n && A.cols() == n)
        A_ = A;
}

VecX fluid_load::force(double t) const
{
    const int ns = (int)ts_.size();
    if(ns == 0)
        return VecX::Zero(6*nbody());
    // Lagrange polynomial through the last k samples; beyond a provisional sample inside the step (CFD's
    // low-storage RK3: its second stage starts at 8/15 dt) at most linear: a parabola through samples
    // close in time amplifies their noise far beyond them
    int o = std::max(0, order);
    if(prov_ && t > ts_.back() + 1.0e-12*std::max(1.0, std::fabs(t)))
        o = std::min(o, 1);
    const int k = std::min(ns, o + 1);
    VecX F = VecX::Zero(6*nbody());
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

void fluid_load::set_attachment(const model &m, int k, const Vec3 &p)
{
    if(k >= 0 && k < nbody())
        s_[k] = m.body(map_[k]).to_body(p);
}

VecX fluid_load::acceleration(const model &m) const
{
    VecX a(6*nbody());
    for(int i=0; i<nbody(); ++i)
    {
        const rigid_body &b = m.body(map_[i]);
        a.segment<3>(6*i) = b.vd;
        a.segment<3>(6*i+3) = b.R*b.wd;
    }
    return a;
}

VecX fluid_load::velocity(const model &m) const
{
    VecX v(6*nbody());
    for(int i=0; i<nbody(); ++i)
    {
        const rigid_body &b = m.body(map_[i]);
        v.segment<3>(6*i) = b.v;
        v.segment<3>(6*i+3) = b.R*b.w;
    }
    return v;
}

VecX fluid_load::added_mass_force(const model &m) const
{
    return -(A_ + As_)*acceleration(m);
}

void fluid_load::eval(const model &m, double t, VecX &f) const
{
    const VecX F = force(t) + added_mass_force(m);
    for(int i=0; i<nbody(); ++i)
    {
        const int b = map_[i];
        m.add_force(f, b, m.body(b).point(s_[i]), F.segment<3>(6*i));
        m.add_moment(f, b, F.segment<3>(6*i+3));
    }
}

void fluid_load::channel_names(std::vector<std::string> &c) const
{
    for(int i=0; i<nbody(); ++i)
        for(const char *x : {"Fx", "Fy", "Fz", "Mx", "My", "Mz"})
            c.push_back(name_ + "_" + bname_[i] + "_" + x);
}

void fluid_load::channel_values(const model &m, double t, std::vector<double> &v) const
{
    // the total fluid load, added-mass part included
    const VecX F = force(t) + added_mass_force(m);
    for(int k=0; k<F.size(); ++k)
        v.push_back(F(k));
}

// restart from a checkpoint: the points of application and the samples of the fluid loads
void fluid_load::save_state(std::vector<double> &v) const
{
    state_out o{v};
    o.list(s_, [&](const Vec3 &x){o.m(x);});
    o.vd(ts_);
    o.list(Fs_, [&](const VecX &x){o.m(x);});
    o.b(prov_);
    o.m(A_);
    o.m(As_);
}

bool fluid_load::load_state(const std::vector<double> &v, size_t &k)
{
    state_in in{v, k};
    in.list(s_, [&](Vec3 &x){in.m(x);});
    in.vd(ts_);
    in.list(Fs_, [&](VecX &x){in.m(x);});
    in.b(prov_);
    in.m(A_);
    in.m(As_);
    k = in.k;
    return in.ok;
}

}
