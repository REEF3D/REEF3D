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

#include"coral/pto/coral_pto.h"
#include"coral/pto/coral_pto_elements.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_constraint.h"
#include<algorithm>
#include<cmath>

namespace coral
{

pto_load::pto_load(const joint &j, bool rotational, const pto_settings &s, const std::string &name)
    : settings(s), name_(name), j_(&j), rot_(rotational)
{
    build();
}

pto_load::pto_load(const model &m, int a, const Vec3 &pa, int b, const Vec3 &pb, const pto_settings &s,
                   const std::string &name)
    : settings(s), name_(name), j_(nullptr), rot_(false), a_(a), b_(b)
{
    sa_ = m.body(a).to_body(pa);
    sb_ = m.body(b).to_body(pb);
    L0_ = (pb - pa).norm();
    build();
}

void pto_load::build()
{
    const pto_settings &s = settings;
    pto_ = pto_composite();
    ctrl_.clear();
    hyd_ = nullptr;
    // generator: spring-damper with preload, Coulomb friction; limits on their sum
    if(s.k != 0.0 || s.c != 0.0 || s.F0 != 0.0 || s.commanded)
        pto_.add(std::make_unique<pto_springdamper>(s.k, s.c, s.q0, s.F0));
    if(s.coulomb > 0.0)
        pto_.add(std::make_unique<pto_coulomb>(s.coulomb, s.v_eps));
    // drive-train losses: act on the bodies, not scaled or limited with the generator
    if(s.c_loss > 0.0)
        pto_.add(std::make_unique<pto_damper>(s.c_loss, true));
    if(s.coulomb_loss > 0.0)
        pto_.add(std::make_unique<pto_coulomb>(s.coulomb_loss, s.v_eps, true));
    pto_.Fmax = s.F_max;
    pto_.Pmax = s.P_max;
    if(s.k_stop > 0.0 || s.c_stop > 0.0)
        pto_.add(std::make_unique<pto_endstop>(s.q_min, s.q_max, s.k_stop, s.c_stop));
    if(s.control == pto_control::latching && s.t_hold > 0.0)
    {
        pto_.add(std::make_unique<pto_latch>(s.k_latch, s.c_latch));
        ctrl_.push_back(std::make_unique<pto_ctrl_latching>(s.t_hold));
    }
    if(s.control == pto_control::declutching && s.t_hold > 0.0)
        ctrl_.push_back(std::make_unique<pto_ctrl_declutching>(s.t_hold));
    // hydraulic PTO: its own element after the others (the checkpoints of the others stay as they were)
    if(s.has_hydraulic)
    {
        auto h = std::make_unique<pto_hydraulic>(s.hydraulic);
        hyd_ = h.get();
        pto_.add(std::move(h));
    }
    gen_ = generator_model(s.generator);
}

void pto_load::coupled_bodies(std::vector<int> &b) const
{
    if(j_)
    {
        b.push_back(j_->body_a);
        b.push_back(j_->body_b);
    }
    else
    {
        b.push_back(a_);
        b.push_back(b_);
    }
}

double pto_load::coordinate(const model &m) const
{
    if(!j_)
        return (m.body(b_).point(sb_) - m.body(a_).point(sa_)).norm() - L0_;
    return rot_ ? j_->rotation_coordinate(m) : j_->coordinate(m);
}

double pto_load::rate(const model &m) const
{
    if(!j_)
    {
        const Vec3 d = m.body(b_).point(sb_) - m.body(a_).point(sa_);
        const double L = d.norm();
        if(L < 1.0e-14)
            return 0.0;
        return (m.body(b_).point_velocity(sb_) - m.body(a_).point_velocity(sa_)).dot(d/L);
    }
    return rot_ ? j_->rotation_rate(m) : j_->rate(m);
}

void pto_load::apply(const model &m, double F, VecX &f) const
{
    if(!j_)
    {
        const Vec3 pa = m.body(a_).point(sa_), pb = m.body(b_).point(sb_);
        const Vec3 d = pb - pa;
        const double L = d.norm();
        if(L < 1.0e-14)
            return;
        const Vec3 e = d/L;
        m.add_force(f, b_, pb, F*e);
        m.add_force(f, a_, pa, -F*e);
    }
    else if(rot_)
        j_->add_axial_torque(m, F, f);
    else
        j_->add_axial_force(m, F, f);
}

generator_power pto_load::electrical(const model &m) const
{
    if(!has_generator())
        return generator_power();
    const pto_output o = evaluate(m, t_);
    return gen_.convert(o.F_gen, rate(m));
}

double pto_load::electrical_energy() const
{
    if(hyd_ && settings.hydraulic.has_gen)
        return hyd_->E_el;
    return E_el_;
}

pto_output pto_load::evaluate(const model &m, double t) const
{
    pto_state st;
    st.q = coordinate(m);
    st.qd = rate(m);
    st.t = t;
    return pto_.force(st);
}

double pto_load::force(const model &m, double t) const
{
    return evaluate(m, t).F;
}

double pto_load::force(const model &m) const
{
    return evaluate(m, t_).F;
}

double pto_load::power(const model &m) const
{
    return evaluate(m, t_).P;
}

double pto_load::loss_power(const model &m) const
{
    return evaluate(m, t_).P_loss;
}

bool pto_load::held() const
{
    for(const auto &c : ctrl_)
        if(c->state() > 0.5)
            return true;
    return false;
}

void pto_load::eval(const model &m, double t, VecX &f) const
{
    apply(m, force(m, t), f);
}

double pto_load::potential(const model &m) const
{
    const pto_settings &s = settings;
    const double q = coordinate(m);
    double e = 0.5*s.k*(q - s.q0)*(q - s.q0) - s.F0*q;
    if(q > s.q_max) e += 0.5*s.k_stop*(q - s.q_max)*(q - s.q_max);
    if(q < s.q_min) e += 0.5*s.k_stop*(q - s.q_min)*(q - s.q_min);
    return e;
}

void pto_load::initialize(const model &m, double t0)
{
    build();
    t_ = t0;
    t0_ = t0;
    energy_ = 0.0;
    energy_loss_ = 0.0;
    E_el_ = 0.0;
    // elements with internal states start at the initial state
    pto_state st;
    st.q = coordinate(m);
    st.qd = rate(m);
    st.t = t0;
    pto_.start(st);
    // the controllers see the initial state first (sign of qd for the zero crossings)
    update_controllers(m, t0, st.q, st.qd);
    P_last_ = power(m);
    Ploss_last_ = loss_power(m);
    Pe_last_ = has_generator() ? electrical(m).P_grid : 0.0;
}

void pto_load::update_controllers(const model &m, double t, double q, double qd)
{
    pto_ctrl_in in;
    in.t = t;
    in.q = q;
    in.qd = qd;
    if(!ext_.empty())
    {
        const pto_output o = evaluate(m, t);
        in.F = o.F;
        in.P = o.P;
    }
    for(auto &c : ctrl_)
        c->update(in, pto_);
    for(auto &c : ext_)
        c->update(in, pto_);
}

void pto_load::commit(const model &m, double t, double dt)
{
    // absorbed energy: trapezoidal rule over the step with the forces of the step (controller
    // state before the update)
    t_ = t;
    const double P = power(m);
    energy_ += 0.5*dt*(P + P_last_);
    const double Pl = loss_power(m);
    energy_loss_ += 0.5*dt*(Pl + Ploss_last_);
    pto_state st;
    st.q = coordinate(m);
    st.qd = rate(m);
    st.t = t;
    const double Pe = has_generator() ? electrical(m).P_grid : 0.0;
    E_el_ += 0.5*dt*(Pe + Pe_last_);
    pto_.accept(st, dt);
    // the controllers act on the converged state at t, for the next step
    update_controllers(m, t, st.q, st.qd);
    if(ext_.empty())
    {
        P_last_ = power(m);
        Ploss_last_ = loss_power(m);
        Pe_last_ = has_generator() ? electrical(m).P_grid : 0.0;
    }
    else
    {
        // a controller that sets the force once per step makes it jump at t; the time integration
        // (generalized-alpha) takes the force of the step's start from the end of the last step, so the
        // force of the step is the mean of the old and the new one: the energies follow it with the
        // powers of the old force at the step's start
        P_last_ = P;
        Ploss_last_ = Pl;
        Pe_last_ = Pe;
    }
}

void pto_load::channel_names(std::vector<std::string> &c) const
{
    for(const char *x : {"q", "qd", "F", "P", "E"})
        c.push_back(name_ + "_" + x);
    if(has_loss())
        for(const char *x : {"Ploss", "Eloss"})
            c.push_back(name_ + "_" + x);
    if(has_generator())
    {
        for(const char *x : {"Pe", "Ee"})
            c.push_back(name_ + "_" + x);
        if(settings.generator.kind == generator_kind::pmsg)
            c.push_back(name_ + "_I");
    }
    pto_.channel_names(name_, c);
    for(const auto &x : ext_)
        x->channel_names(name_, c);
}

void pto_load::channel_values(const model &m, double t, std::vector<double> &v) const
{
    const pto_output o = evaluate(m, t);
    v.push_back(coordinate(m));
    v.push_back(rate(m));
    v.push_back(o.F);
    v.push_back(o.P);
    v.push_back(energy_);
    if(has_loss())
    {
        v.push_back(o.P_loss);
        v.push_back(energy_loss_);
    }
    if(has_generator())
    {
        const generator_power g = gen_.convert(o.F_gen, rate(m));
        v.push_back(g.P_grid);
        v.push_back(E_el_);
        if(settings.generator.kind == generator_kind::pmsg)
            v.push_back(g.I);
    }
    pto_.channel_values(v);
    for(const auto &x : ext_)
        x->channel_values(v);
}

// restart from a checkpoint: the time and power of the last converged state, the absorbed energy, the
// elements (their parameters as the controllers set them) and the controllers
void pto_load::save_state(std::vector<double> &v) const
{
    for(double x : {t_, P_last_, energy_, t0_})
        v.push_back(x);
    pto_.save_state(v);
    v.push_back(double(ctrl_.size()));
    for(const auto &c : ctrl_)
        c->save_state(v);
    // the losses after everything else, so the checkpoints without them stay as they were
    if(has_loss())
        for(double x : {Ploss_last_, energy_loss_})
            v.push_back(x);
    // the generator and the further controllers after that
    if(has_generator())
        for(double x : {Pe_last_, E_el_})
            v.push_back(x);
    if(!ext_.empty())
    {
        v.push_back(double(ext_.size()));
        for(const auto &c : ext_)
            c->save_state(v);
    }
}

bool pto_load::load_state(const std::vector<double> &v, size_t &k)
{
    if(k + 4 > v.size())
        return false;
    t_ = v[k++]; P_last_ = v[k++]; energy_ = v[k++]; t0_ = v[k++];
    pto_.load_state(v, k);
    if(k >= v.size())
        return false;
    const size_t n = size_t(v[k++]);
    for(size_t q=0; q<n && q<ctrl_.size(); ++q)
        ctrl_[q]->load_state(v, k);
    if(has_loss())
    {
        if(k + 2 > v.size())
            return false;
        Ploss_last_ = v[k++]; energy_loss_ = v[k++];
    }
    if(has_generator())
    {
        if(k + 2 > v.size())
            return false;
        Pe_last_ = v[k++]; E_el_ = v[k++];
    }
    if(!ext_.empty())
    {
        if(k >= v.size())
            return false;
        const size_t n = size_t(v[k++]);
        for(size_t q=0; q<n && q<ext_.size(); ++q)
            ext_[q]->load_state(v, k);
    }
    return k <= v.size();
}

}
