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

#include"coral/pto/coral_joint_force.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_model.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<sstream>

namespace coral
{

joint_force::joint_force(const joint &j, bool rotational, const jforce_settings &s, const std::string &name)
    : settings(s), name_(name), j_(&j), rot_(rotational)
{
}

void joint_force::coupled_bodies(std::vector<int> &b) const
{
    b.push_back(j_->body_a);
    b.push_back(j_->body_b);
}

double joint_force::force(double t) const
{
    const jforce_settings &s = settings;
    const std::vector<double> &T = s.t, &F = s.F;
    if(T.empty() || t < 0.0)
        return 0.0;
    double tt = t;
    if(s.periodic && T.size() > 1)
    {
        const double P = T.back() - T.front();
        tt = T.front() + std::fmod(t - T.front(), P);
        if(tt < T.front()) tt += P;
    }
    double v = 0.0;
    if(tt >= T.front() && tt <= T.back())
    {
        const std::size_t i = std::min<std::size_t>(std::upper_bound(T.begin(), T.end(), tt) - T.begin(), T.size() - 1);
        const std::size_t i0 = i == 0 ? 0 : i - 1;
        const double dt = T[i] - T[i0];
        v = dt > 0.0 ? F[i0] + (F[i] - F[i0])*(tt - T[i0])/dt : F[i];
    }
    double r = 1.0;
    if(s.ramp > 0.0 && t < s.ramp)
        r = 0.5*(1.0 - std::cos(M_PI*t/s.ramp));
    return s.scale*r*v;
}

void joint_force::eval(const model &m, double t, VecX &f) const
{
    const double F = force(t);
    if(rot_)
        j_->add_axial_torque(m, F, f);
    else
        j_->add_axial_force(m, F, f);
}

void joint_force::channel_names(std::vector<std::string> &c) const
{
    for(const char *x : {"F", "q", "qd"})
        c.push_back(name_ + "_" + x);
}

void joint_force::channel_values(const model &m, double t, std::vector<double> &v) const
{
    v.push_back(force(t));
    v.push_back(rot_ ? j_->rotation_coordinate(m) : j_->coordinate(m));
    v.push_back(rot_ ? j_->rotation_rate(m) : j_->rate(m));
}

bool joint_force::read_table(const std::string &file, std::vector<double> &t, std::vector<double> &F, std::string &msg)
{
    std::ifstream in(file);
    if(!in)
    {
        msg = "cannot open '" + file + "'";
        return false;
    }
    t.clear(); F.clear();
    std::string line;
    int n = 0;
    while(std::getline(in, line))
    {
        ++n;
        const std::size_t h = line.find('#');
        if(h != std::string::npos) line.erase(h);
        std::istringstream ls(line);
        double a, b;
        if(!(ls>>a))
            continue;
        if(!(ls>>b))
        {
            msg = file + " line " + std::to_string(n) + ": two numbers (t F) expected";
            return false;
        }
        if(!t.empty() && a <= t.back())
        {
            msg = file + " line " + std::to_string(n) + ": the times must increase";
            return false;
        }
        t.push_back(a);
        F.push_back(b);
    }
    if(t.size() < 2)
    {
        msg = file + ": at least two lines t F needed";
        return false;
    }
    return true;
}

}
