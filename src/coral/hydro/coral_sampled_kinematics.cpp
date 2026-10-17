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

#include"coral/hydro/coral_sampled_kinematics.h"
#include"coral/core/coral_model.h"
#include<algorithm>

namespace coral
{

double sampled_kinematics::elevation(double x, double y, double t) const
{
    return bg_ ? bg_->elevation(x, y, t) : 0.0;
}

void sampled_kinematics::kinematics(const Vec3 &p, double t, Vec3 &u, Vec3 &du, double &p_rho) const
{
    if(bg_)
        bg_->kinematics(p, t, u, du, p_rho);
    else
    {
        u.setZero();
        du.setZero();
        p_rho = 0.0;
    }
}

void sampled_kinematics::collect(const model &m, std::vector<Vec3> &points)
{
    points.clear();
    owners_.clear();
    for(const auto &l : m.load_list)
    {
        const size_t n0 = points.size();
        l->kinematics_points(m, this, points);
        if(points.size() > n0)
            owners_[l.get()] = {(int)n0, (int)(points.size() - n0)};
    }
    if((int)points.size() != n_)
        lev_.clear();
    n_ = (int)points.size();
}

void sampled_kinematics::add_sample(double t, const std::vector<Vec3> &u, const std::vector<double> &eta, const std::vector<int> &inside)
{
    if((int)u.size() != n_ || (int)eta.size() != n_ || (int)inside.size() != n_)
        return;
    while(!lev_.empty() && lev_.back().t >= t)
        lev_.pop_back();
    lev_.push_back(level{t, u, eta, inside});
    while(lev_.size() > 3)
        lev_.pop_front();
}

int sampled_kinematics::find(const void *owner, int index) const
{
    if(lev_.empty())
        return -1;
    const auto it = owners_.find(owner);
    if(it == owners_.end() || index < 0 || index >= it->second.second)
        return -1;
    const int k = it->second.first + index;
    return lev_.back().in[k] ? k : -1;
}

int sampled_kinematics::weights(double t, double *w, double *dw) const
{
    const int ns = (int)lev_.size();
    const int k = std::min(ns, std::max(0, order) + 1);
    const int i0 = ns - k;
    for(int i=0; i<k; ++i)
    {
        const double ti = lev_[i0+i].t;
        double l = 1.0, dl = 0.0;
        for(int j=0; j<k; ++j)
        {
            if(j == i)
                continue;
            const double tj = lev_[i0+j].t, c = 1.0/(ti - tj);
            // derivative of the product: d/dt prod (t - tj) c
            dl = dl*(t - tj)*c + l*c;
            l *= (t - tj)*c;
        }
        w[i] = l;
        dw[i] = dl;
    }
    return k;
}

double sampled_kinematics::elevation_at(const void *owner, int index, const Vec3 &p, double t) const
{
    const int q = find(owner, index);
    if(q < 0)
        return elevation(p(0), p(1), t);
    double w[3], dw[3];
    const int k = weights(t, w, dw);
    const int i0 = (int)lev_.size() - k;
    double e = 0.0;
    for(int i=0; i<k; ++i)
        e += w[i]*lev_[i0+i].eta[q];
    return e;
}

void sampled_kinematics::kinematics_at(const void *owner, int index, const Vec3 &p, double t, Vec3 &u, Vec3 &du, double &p_rho) const
{
    const int q = find(owner, index);
    if(q < 0)
    {
        kinematics(p, t, u, du, p_rho);
        return;
    }
    double w[3], dw[3];
    const int k = weights(t, w, dw);
    const int i0 = (int)lev_.size() - k;
    u.setZero();
    du.setZero();
    for(int i=0; i<k; ++i)
    {
        u += w[i]*lev_[i0+i].u[q];
        du += dw[i]*lev_[i0+i].u[q];
    }
    p_rho = 0.0;
}

}
