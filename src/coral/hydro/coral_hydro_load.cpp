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

#include"coral/hydro/coral_hydro_load.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_rotation.h"
#include<cmath>
#include<complex>
#include<stdexcept>

namespace coral
{

hydro_load::hydro_load(const model &m, std::shared_ptr<const hydro_db> db, const std::vector<int> &bodies,
                       std::shared_ptr<const wave_field> waves, const hydro_settings &s, const std::string &name)
    : name_(name), db_(std::move(db)), waves_(std::move(waves)), set_(s), map_(bodies)
{
    db_->check();
    const int N = db_->nbody();
    if((int)map_.size() != N)
        throw std::runtime_error("coral::hydro_load: the database has " + std::to_string(N) + " bodies");
    s_.resize(N);
    p0_.resize(N);
    R0_.resize(N);
    buoy_.resize(N);
    for(int i=0; i<N; ++i)
    {
        const rigid_body &b = m.body(map_[i]);
        p0_[i] = db_->body[i].ref;
        s_[i] = b.to_body(p0_[i]);
        R0_[i] = b.R;
        buoy_[i] = db_->body[i].volume >= 0.0 ? db_->rho*db_->g*db_->body[i].volume : b.mass*std::fabs(m.gravity(2));
    }
    // radiation
    if(set_.radiation)
    {
        K_.compute(*db_, set_.t_memory, set_.dt_memory);
        Ainf_ = db_->has_ainf ? db_->Ainf : ogilvie_ainf(*db_, K_);
        if(set_.state_space)
            ss_.identify(K_, set_.ss_tol, set_.ss_order);
    }
    else
        Ainf_ = MatX::Zero(6*N, 6*N);
    // excitation per wave component
    if(set_.excitation && waves_)
        for(const auto &c : waves_->components())
            Xc_.push_back(c.a*std::exp(std::complex<double>(0.0, c.phase))*db_->excitation(c.w, c.beta));
}

VecX hydro_load::displacement(const model &m) const
{
    const int N = db_->nbody();
    VecX xi(6*N);
    for(int i=0; i<N; ++i)
    {
        const rigid_body &b = m.body(map_[i]);
        xi.segment<3>(6*i) = b.point(s_[i]) - p0_[i];
        xi.segment<3>(6*i+3) = log_so3(b.R*R0_[i].transpose());
    }
    return xi;
}

VecX hydro_load::velocity(const model &m) const
{
    const int N = db_->nbody();
    VecX v(6*N);
    for(int i=0; i<N; ++i)
    {
        const rigid_body &b = m.body(map_[i]);
        v.segment<3>(6*i) = b.point_velocity(s_[i]);
        v.segment<3>(6*i+3) = b.omega_global();
    }
    return v;
}

VecX hydro_load::acceleration(const model &m) const
{
    const int N = db_->nbody();
    VecX a(6*N);
    for(int i=0; i<N; ++i)
    {
        const rigid_body &b = m.body(map_[i]);
        a.segment<3>(6*i) = b.point_acceleration(s_[i]);
        a.segment<3>(6*i+3) = b.R*b.wd;
    }
    return a;
}

VecX hydro_load::hydrostatic_force(const model &m) const
{
    const int N = db_->nbody();
    VecX F = VecX::Zero(6*N);
    if(set_.hydrostatics)
        F -= db_->C*displacement(m);
    if(set_.buoyancy)
        for(int i=0; i<N; ++i)
            F(6*i+2) += buoy_[i];
    return F;
}

VecX hydro_load::excitation_force(double t) const
{
    const int nd = db_->ndof();
    if(t == exc_t_)
        return exc_f_;
    VecX F = VecX::Zero(nd);
    if(set_.excitation && waves_)
    {
        const auto &comp = waves_->components();
        for(size_t c=0; c<comp.size(); ++c)
        {
            const std::complex<double> e = std::exp(std::complex<double>(0.0, comp[c].w*t));
            F += (Xc_[c]*e).real();
        }
        F *= waves_->ramp(t);
    }
    exc_t_ = t;
    exc_f_ = F;
    return F;
}

// int K(t - tau) v(tau) dtau over the memory: committed history (trapezoid on its times) plus the
// current time t with the trial velocity v1
VecX hydro_load::memory_force(const VecX &v1, double t) const
{
    const int nd = db_->ndof();
    if(!set_.radiation || !started_)
        return VecX::Zero(nd);
    if(set_.state_space)
    {
        VecX F = VecX::Zero(nd);
        const double h = t - t_last_;
        for(size_t p=0; p<ss_.pairs.size(); ++p)
        {
            const auto &q = ss_.pairs[p];
            Eigen::VectorXcd x = xs_[p];
            if(h > 0.0)
                radiation_ss::advance(q, x, v_last_(q.j), v1(q.j), h);
            F(q.i) += radiation_ss::output(q, x);
        }
        return F;
    }
    const size_t n = th_.size();
    if(t != cache_t_ || n != cache_n_)
    {
        // weights of the history points for the trapezoid over [.., t_n, t]
        VecX Fh = VecX::Zero(nd);
        const double T = K_.tmax();
        MatX Kt(nd, nd);
        for(size_t j=n; j-- > 0;)
        {
            const double age = t - th_[j];
            if(age > T)
                break;
            const double left = j > 0 ? th_[j] - th_[j-1] : 0.0;
            const double right = (j+1 < n ? th_[j+1] : t) - th_[j];
            // K(age) by linear interpolation of the samples
            const double x = age/K_.dt();
            int k = (int)x;
            if(k >= K_.nsample()-1) {k = K_.nsample()-2;}
            const double f = x - k;
            Kt = (1.0 - f)*K_.sample(k) + f*K_.sample(k+1);
            Fh += 0.5*(left + right)*(Kt*vh_[j]);
        }
        cache_hist_ = Fh;
        cache_t_ = t;
        cache_n_ = n;
    }
    return cache_hist_ + 0.5*(t - th_.back())*(K_.sample(0)*v1);
}

VecX hydro_load::radiation_force(const model &m, double t) const
{
    const int nd = db_->ndof();
    if(!set_.radiation)
        return VecX::Zero(nd);
    return -Ainf_*acceleration(m) - memory_force(velocity(m), t);
}

VecX hydro_load::drag_force(const model &m) const
{
    const int nd = db_->ndof();
    VecX F = VecX::Zero(nd);
    if(set_.linear_damping.size() || set_.quadratic_drag.size())
    {
        const VecX v = velocity(m);
        if(set_.linear_damping.size())
            F -= set_.linear_damping*v;
        if(set_.quadratic_drag.size())
            F -= set_.quadratic_drag*(v.cwiseAbs().cwiseProduct(v));
    }
    return F;
}

void hydro_load::eval(const model &m, double t, VecX &f) const
{
    const VecX F = hydrostatic_force(m) + excitation_force(t) + radiation_force(m, t) + drag_force(m);
    for(int i=0; i<db_->nbody(); ++i)
    {
        const rigid_body &b = m.body(map_[i]);
        m.add_force(f, map_[i], b.point(s_[i]), F.segment<3>(6*i));
        m.add_moment(f, map_[i], F.segment<3>(6*i+3));
    }
}

void hydro_load::initialize(const model &m, double t0)
{
    th_.clear();
    vh_.clear();
    const VecX v = velocity(m);
    if(!set_.state_space)
    {
        th_.push_back(t0);
        vh_.push_back(v);
    }
    started_ = true;
    v_last_ = v;
    t_last_ = t0;
    xs_.clear();
    for(const auto &q : ss_.pairs)
        xs_.push_back(Eigen::VectorXcd::Zero(q.lam.size()));
    cache_t_ = -1.0e300;
    exc_t_ = -1.0e300;
}

void hydro_load::commit(const model &m, double t, double)
{
    const VecX v = velocity(m);
    if(set_.state_space)
    {
        const double h = t - t_last_;
        for(size_t p=0; p<ss_.pairs.size(); ++p)
            radiation_ss::advance(ss_.pairs[p], xs_[p], v_last_(ss_.pairs[p].j), v(ss_.pairs[p].j), h);
    }
    else
    {
        th_.push_back(t);
        vh_.push_back(v);
        // drop history older than the memory (keep one point beyond it)
        size_t drop = 0;
        while(drop + 1 < th_.size() && t - th_[drop+1] > K_.tmax())
            ++drop;
        if(drop > 1000)
        {
            th_.erase(th_.begin(), th_.begin() + drop);
            vh_.erase(vh_.begin(), vh_.begin() + drop);
        }
    }
    v_last_ = v;
    t_last_ = t;
}

void hydro_load::channel_names(std::vector<std::string> &c) const
{
    const char *const dof[6] = {"x", "y", "z", "rx", "ry", "rz"};
    for(int i=0; i<db_->nbody(); ++i)
    {
        const std::string b = name_ + "_" + (db_->body[i].name.empty() ? std::to_string(i+1) : db_->body[i].name);
        for(const char *part : {"exc", "rad", "hs"})
            for(int d=0; d<6; ++d)
                c.push_back(b + "_" + part + "_" + dof[d]);
    }
}

void hydro_load::channel_values(const model &m, double t, std::vector<double> &v) const
{
    const VecX fe = excitation_force(t), fr = radiation_force(m, t), fh = hydrostatic_force(m);
    for(int i=0; i<db_->nbody(); ++i)
        for(const VecX *F : {&fe, &fr, &fh})
            for(int d=0; d<6; ++d)
                v.push_back((*F)(6*i+d));
}

}
