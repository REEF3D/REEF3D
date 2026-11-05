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

#include"coral/structure/coral_pitch_actuator.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_model.h"
#include<algorithm>
#include<cmath>

namespace coral
{

pitch_actuator::pitch_actuator(std::vector<const joint*> joints, double kk, double cc, double theta0, const std::string &name)
    : k(kk), c(cc), j_(std::move(joints)), theta0_(theta0), cmd_(j_.size(), theta0), name_(name) {}

void pitch_actuator::coupled_bodies(std::vector<int> &b) const
{
    for(const joint *j : j_)
    {
        if(std::find(b.begin(), b.end(), j->body_a)==b.end()) b.push_back(j->body_a);
        b.push_back(j->body_b);
    }
}

void pitch_actuator::command(int blade, double t, double &th, double &rt) const
{
    rt = 0.0;
    if(!table_) {th = cmd_[blade]; return;}
    if(t <= tt_.front()) {th = tv_.front(); return;}
    if(t >= tt_.back()) {th = tv_.back(); return;}
    const size_t q = std::upper_bound(tt_.begin(), tt_.end(), t) - tt_.begin() - 1;
    rt = (tv_[q+1] - tv_[q])/(tt_[q+1] - tt_[q]);
    th = tv_[q] + (t - tt_[q])*rt;
}

double pitch_actuator::angle(const model &m, int blade) const {return theta0_ + j_[blade]->rotation_coordinate(m);}
double pitch_actuator::rate(const model &m, int blade) const {return j_[blade]->rotation_rate(m);}

double pitch_actuator::torque(const model &m, int blade, double t) const
{
    double th, rt;
    command(blade, t, th, rt);
    return -k*(angle(m, blade) - th) - c*(rate(m, blade) - rt);
}

void pitch_actuator::eval(const model &m, double t, VecX &f) const
{
    for(int q=0; q<(int)j_.size(); ++q)
        j_[q]->add_axial_torque(m, torque(m, q, t), f);
}

void pitch_actuator::channel_names(std::vector<std::string> &n) const
{
    n.push_back(name_ + "_cmd");
    for(size_t q=0; q<j_.size(); ++q)
    {
        n.push_back(name_ + "_b" + std::to_string(q + 1));
        n.push_back(name_ + "_b" + std::to_string(q + 1) + "_T");
    }
}

void pitch_actuator::channel_values(const model &m, double t, std::vector<double> &v) const
{
    const double deg = 180.0/M_PI;
    double th, rt;
    command(0, t, th, rt);
    v.push_back(th*deg);
    for(int q=0; q<(int)j_.size(); ++q)
    {
        v.push_back(angle(m, q)*deg);
        v.push_back(torque(m, q, t));
    }
}

}
