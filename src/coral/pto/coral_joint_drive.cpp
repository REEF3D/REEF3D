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
#include"coral/pto/coral_joint_force.h"
#include<algorithm>
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
    if(!s.t.empty())
    {
        double Q, Qd, Qdd;
        table(t, Q, Qd, Qdd);
        q = s.q0 + s.scale*r*Q;
        qd = s.scale*(rd*Q + r*Qd);
        qdd = s.scale*(rdd*Q + 2.0*rd*Qd + r*Qdd);
        return;
    }
    const double sn = std::sin(a), cs = std::cos(a);
    q = s.q0 + s.amp*r*sn;
    qd = s.amp*(rd*sn + r*w*cs);
    qdd = s.amp*(rdd*sn + 2.0*rd*w*cs - r*w*w*sn);
}

void joint_drive::table(double t, double &Q, double &Qd, double &Qdd) const
{
    const std::vector<double> &T = settings.t, &V = settings.q;
    const int n = (int)T.size();
    Q = Qd = Qdd = 0.0;
    if(n == 0)
        return;
    if(n == 1)
    {
        Q = V[0];
        return;
    }
    const double P = T[n-1] - T[0];
    if(settings.periodic)
        t = T[0] + (t - T[0]) - P*std::floor((t - T[0])/P);
    else if(t <= T[0] || t >= T[n-1])
    {
        Q = (t <= T[0]) ? V[0] : V[n-1];
        return;
    }
    // slope at point i (Catmull-Rom on the uneven grid; periodic: the table wraps, the last point is the first)
    auto slope = [&](int i) -> double
    {
        if(settings.periodic && (i == 0 || i == n-1))
            return (V[1] - V[n-2])/(T[1] - T[0] + T[n-1] - T[n-2]);
        if(i == 0) return (V[1] - V[0])/(T[1] - T[0]);
        if(i == n-1) return (V[n-1] - V[n-2])/(T[n-1] - T[n-2]);
        return (V[i+1] - V[i-1])/(T[i+1] - T[i-1]);
    };
    int i = (int)(std::upper_bound(T.begin(), T.end(), t) - T.begin()) - 1;
    i = std::clamp(i, 0, n-2);
    const double h = T[i+1] - T[i], u = (t - T[i])/h;
    const double m0 = slope(i)*h, m1 = slope(i+1)*h;
    const double u2 = u*u, u3 = u2*u;
    Q = (2.0*u3 - 3.0*u2 + 1.0)*V[i] + (u3 - 2.0*u2 + u)*m0 + (-2.0*u3 + 3.0*u2)*V[i+1] + (u3 - u2)*m1;
    Qd = ((6.0*u2 - 6.0*u)*V[i] + (3.0*u2 - 4.0*u + 1.0)*m0 + (-6.0*u2 + 6.0*u)*V[i+1] + (3.0*u2 - 2.0*u)*m1)/h;
    Qdd = ((12.0*u - 6.0)*V[i] + (6.0*u - 4.0)*m0 + (-12.0*u + 6.0)*V[i+1] + (6.0*u - 2.0)*m1)/(h*h);
}

bool joint_drive::read_table(const std::string &file, std::vector<double> &t, std::vector<double> &q, std::string &msg)
{
    // the same format as the force table: two numbers per line, times increasing
    if(!joint_force::read_table(file, t, q, msg))
    {
        const std::string tf = "(t F)";
        const size_t p = msg.find(tf);
        if(p != std::string::npos) msg.replace(p, tf.size(), "(t q)");
        return false;
    }
    return true;
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
