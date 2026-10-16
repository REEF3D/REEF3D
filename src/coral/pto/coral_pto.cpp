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
#include"coral/core/coral_model.h"
#include"coral/core/coral_constraint.h"
#include<algorithm>
#include<cmath>

namespace coral
{

pto_load::pto_load(const joint &j, bool rotational, const pto_settings &s, const std::string &name)
    : settings(s), name_(name), j_(&j), rot_(rotational)
{}

double pto_load::coordinate(const model &m) const
{
    return rot_ ? j_->rotation_coordinate(m) : j_->coordinate(m);
}

double pto_load::rate(const model &m) const
{
    return rot_ ? j_->rotation_rate(m) : j_->rate(m);
}

double pto_load::force(const model &m) const
{
    const pto_settings &s = settings;
    const double q = coordinate(m), qd = rate(m);
    double F = 0.0;
    if(held_ && s.control==pto_control::latching)
        F = -s.c_latch*qd;
    else if(!(held_ && s.control==pto_control::declutching))
    {
        F = s.F0 - s.k*(q - s.q0) - s.c*qd;
        if(s.coulomb > 0.0)
            F -= s.coulomb*std::tanh(qd/s.v_eps);
        if(s.F_max > 0.0)
            F = std::max(-s.F_max, std::min(s.F_max, F));
        if(s.P_max > 0.0 && std::fabs(F*qd) > s.P_max)
            F = (F > 0.0 ? 1.0 : -1.0)*s.P_max/std::fabs(qd);
    }
    // end stops
    if(q > s.q_max)
        F -= s.k_stop*(q - s.q_max) + (qd > 0.0 ? s.c_stop*qd : 0.0);
    if(q < s.q_min)
        F -= s.k_stop*(q - s.q_min) + (qd < 0.0 ? s.c_stop*qd : 0.0);
    return F;
}

void pto_load::eval(const model &m, double, VecX &f) const
{
    const double F = force(m);
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
    held_ = false;
    armed_ = true;
    t_release_ = t0;
    qd_last_ = rate(m);
    P_last_ = power(m);
    energy_ = 0.0;
    t0_ = t0;
}

void pto_load::commit(const model &m, double t, double dt)
{
    const double qd = rate(m);
    const double P = power(m);
    energy_ += 0.5*dt*(P + P_last_);
    P_last_ = P;
    if(settings.control != pto_control::passive)
    {
        // a zero crossing of qd starts a hold; after a release the PTO is armed again once the body
        // moves (|qd| > 10 v_eps), so the small velocities at the end of a hold do not start the next
        if(!armed_ && std::fabs(qd) > 10.0*settings.v_eps)
            armed_ = true;
        if(held_ && t >= t_release_)
        {
            held_ = false;
            armed_ = false;
        }
        else if(!held_ && armed_ && qd_last_*qd < 0.0 && settings.t_hold > 0.0)
        {
            held_ = true;
            t_release_ = t + settings.t_hold;
        }
    }
    qd_last_ = qd;
}

void pto_load::channel_names(std::vector<std::string> &c) const
{
    for(const char *x : {"q", "qd", "F", "P", "E"})
        c.push_back(name_ + "_" + x);
}

void pto_load::channel_values(const model &m, double, std::vector<double> &v) const
{
    v.push_back(coordinate(m));
    v.push_back(rate(m));
    v.push_back(force(m));
    v.push_back(power(m));
    v.push_back(energy_);
}

}
