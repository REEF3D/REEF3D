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

void pto_load::build()
{
    const pto_settings &s = settings;
    pto_ = pto_composite();
    ctrl_.clear();
    // generator: spring-damper with preload, Coulomb friction; limits on their sum
    if(s.k != 0.0 || s.c != 0.0 || s.F0 != 0.0)
        pto_.add(std::make_unique<pto_springdamper>(s.k, s.c, s.q0, s.F0));
    if(s.coulomb > 0.0)
        pto_.add(std::make_unique<pto_coulomb>(s.coulomb, s.v_eps));
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
}

void pto_load::coupled_bodies(std::vector<int> &b) const
{
    b.push_back(j_->body_a);
    b.push_back(j_->body_b);
}

double pto_load::coordinate(const model &m) const
{
    return rot_ ? j_->rotation_coordinate(m) : j_->coordinate(m);
}

double pto_load::rate(const model &m) const
{
    return rot_ ? j_->rotation_rate(m) : j_->rate(m);
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

bool pto_load::held() const
{
    for(const auto &c : ctrl_)
        if(c->state() > 0.5)
            return true;
    return false;
}

void pto_load::eval(const model &m, double t, VecX &f) const
{
    const double F = force(m, t);
    if(rot_)
        j_->add_axial_torque(m, F, f);
    else
        j_->add_axial_force(m, F, f);
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
    // the controllers see the initial state first (sign of qd for the zero crossings)
    pto_ctrl_in in;
    in.t = t0;
    in.q = coordinate(m);
    in.qd = rate(m);
    for(auto &c : ctrl_)
        c->update(in, pto_);
    P_last_ = power(m);
}

void pto_load::commit(const model &m, double t, double dt)
{
    // absorbed energy: trapezoidal rule over the step with the forces of the step (controller
    // state before the update)
    t_ = t;
    const double P = power(m);
    energy_ += 0.5*dt*(P + P_last_);
    pto_state st;
    st.q = coordinate(m);
    st.qd = rate(m);
    st.t = t;
    pto_.accept(st, dt);
    // the controllers act on the converged state at t, for the next step
    pto_ctrl_in in;
    in.t = t;
    in.q = st.q;
    in.qd = st.qd;
    for(auto &c : ctrl_)
        c->update(in, pto_);
    P_last_ = power(m);
}

void pto_load::channel_names(std::vector<std::string> &c) const
{
    for(const char *x : {"q", "qd", "F", "P", "E"})
        c.push_back(name_ + "_" + x);
}

void pto_load::channel_values(const model &m, double t, std::vector<double> &v) const
{
    const pto_output o = evaluate(m, t);
    v.push_back(coordinate(m));
    v.push_back(rate(m));
    v.push_back(o.F);
    v.push_back(o.P);
    v.push_back(energy_);
}

}
