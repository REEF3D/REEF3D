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

#include"coral/pto/coral_joint_drive.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_model.h"
#include<cmath>

namespace coral
{

joint_drive::joint_drive(const joint &j, bool rotational, const drive_settings &s, const std::string &name)
    : settings(s), name_(name), j_(&j), rot_(rotational)
{
}

void joint_drive::coupled_bodies(std::vector<int> &b) const
{
    b.push_back(j_->body_a);
    b.push_back(j_->body_b);
}

void joint_drive::prescribed(double t, double &q, double &qd, double &qdd) const
{
    const drive_settings &s = settings;
    const double w = 2.0*M_PI/s.period, a = w*t + s.phase*M_PI/180.0;
    double r = 1.0, rd = 0.0, rdd = 0.0;
    if(s.ramp > 0.0 && t < s.ramp)
    {
        const double v = M_PI/s.ramp;
        r = 0.5*(1.0 - std::cos(v*t));
        rd = 0.5*v*std::sin(v*t);
        rdd = 0.5*v*v*std::cos(v*t);
    }
    if(t < 0.0)
        r = rd = rdd = 0.0;
    const double sn = std::sin(a), cs = std::cos(a);
    q = s.q0 + s.amp*r*sn;
    qd = s.amp*(rd*sn + r*w*cs);
    qdd = s.amp*(rdd*sn + 2.0*rd*w*cs - r*w*w*sn);
}

double joint_drive::force(const model &m, double t) const
{
    double qp, qdp, qddp;
    prescribed(t, qp, qdp, qddp);
    const double q = rot_ ? j_->rotation_coordinate(m) : j_->coordinate(m);
    const double qd = rot_ ? j_->rotation_rate(m) : j_->rate(m);
    return -settings.k*(q - qp) - settings.c*(qd - qdp);
}

void joint_drive::eval(const model &m, double t, VecX &f) const
{
    const double F = force(m, t);
    if(rot_)
        j_->add_axial_torque(m, F, f);
    else
        j_->add_axial_force(m, F, f);
}

void joint_drive::channel_names(std::vector<std::string> &c) const
{
    for(const char *x : {"F", "q", "q_p"})
        c.push_back(name_ + "_" + x);
}

void joint_drive::channel_values(const model &m, double t, std::vector<double> &v) const
{
    double qp, qdp, qddp;
    prescribed(t, qp, qdp, qddp);
    v.push_back(force(m, t));
    v.push_back(rot_ ? j_->rotation_coordinate(m) : j_->coordinate(m));
    v.push_back(qp);
}

}
