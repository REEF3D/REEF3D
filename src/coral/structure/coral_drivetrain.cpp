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

#include"coral/structure/coral_drivetrain.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_model.h"
#include<cmath>

namespace coral
{

drivetrain::drivetrain(const joint &shaft, const joint *gen, const drivetrain_settings &s, const std::string &name)
    : settings(s), name_(name), shaft_(&shaft), gen_(gen) {}

void drivetrain::coupled_bodies(std::vector<int> &b) const
{
    b.push_back(shaft_->body_a);
    b.push_back(shaft_->body_b);
    if(gen_) b.push_back(gen_->body_b);
}

double drivetrain::omega_rotor(const model &m) const {return shaft_->rotation_rate(m);}

double drivetrain::omega_hss(const model &m) const
{
    return settings.N*(gen_ ? gen_->rotation_rate(m) : shaft_->rotation_rate(m));
}

double drivetrain::generator_torque(const model &m) const
{
    const double w = omega_hss(m);
    switch(settings.law)
    {
        case generator_law::constant: return settings.T;
        case generator_law::slip: return settings.slope*(w - settings.omega_sync);
        case generator_law::k_omega2:
        {
            double T = w > 0.0 ? settings.K*w*w : 0.0;
            if(settings.T > 0.0) T = std::min(T, settings.T);
            return T;
        }
        case generator_law::command: return cmd_;
        default: return 0.0;
    }
}

double drivetrain::shaft_torque(const model &m) const
{
    if(!gen_)
        return 0.0;
    double d = shaft_->rotation_coordinate(m) - gen_->rotation_coordinate(m);
    d = std::remainder(d, 2.0*M_PI);
    return settings.k*d + settings.c*(shaft_->rotation_rate(m) - gen_->rotation_rate(m));
}

void drivetrain::eval(const model &m, double, VecX &f) const
{
    const double Tl = -settings.N*generator_torque(m);
    if(gen_)
    {
        const double Ts = shaft_torque(m);
        shaft_->add_axial_torque(m, -Ts, f);
        gen_->add_axial_torque(m, Ts + Tl, f);
    }
    else
        shaft_->add_axial_torque(m, Tl, f);
}

void drivetrain::channel_names(std::vector<std::string> &c) const
{
    for(const char *x : {"rpm_rotor", "rpm_gen", "T_gen", "P_gen", "T_shaft"})
        c.push_back(name_ + "_" + x);
}

void drivetrain::channel_values(const model &m, double, std::vector<double> &v) const
{
    const double w = omega_hss(m), T = generator_torque(m);
    v.push_back(omega_rotor(m)*30.0/M_PI);
    v.push_back(w*30.0/M_PI);
    v.push_back(T);
    v.push_back(settings.eta*T*w);
    v.push_back(gen_ ? shaft_torque(m) : settings.N*T);
}

}
