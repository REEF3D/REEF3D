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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_pto_controller.h"
#include<algorithm>
#include<cctype>
#include<cmath>
#include<fstream>
#include<sstream>

// --- zero crossing ---

bool pto_zero_crossing::detect(double t, double qd)
{
    const int s = (qd>0.0) - (qd<0.0);

    if(s==0)
    return false;

    bool hit = false;

    if(sref!=0 && s!=sref)
    {
    hit = true;
    t_cross = (qd_last!=qd) ? t_last + (t - t_last)*qd_last/(qd_last - qd) : t;
    }

    sref = s;
    t_last = t;
    qd_last = qd;

    return hit;
}

// --- latching ---

pto_ctrl_latching::pto_ctrl_latching(double t_hold_) : t_hold(t_hold_)
{
}

void pto_ctrl_latching::update(const pto_ctrl_in &in, pto_composite &pto)
{
    if(latched)
    {
        // the latch element itself releases at the stage time t_release
        if(in.t >= t_release - 1.0e-12)
        {
        pto.set_param("latch",0.0);
        latched = false;
        zc.reset();        // the next zero is the one after the motion has resumed
        }

        return;
    }

    if(zc.detect(in.t,in.qd))
    {
    t_release = zc.t_cross + t_hold;
    pto.set_param("q0",in.q);
    pto.set_param("t_release",t_release);
    pto.set_param("latch",1.0);
    latched = true;
    }
}

// --- declutching ---

pto_ctrl_declutching::pto_ctrl_declutching(double t_off_) : t_off(t_off_)
{
}

void pto_ctrl_declutching::update(const pto_ctrl_in &in, pto_composite &pto)
{
    const bool zero = zc.detect(in.t,in.qd);     // keep tracking the sign while off

    if(off && in.t >= t_on - 1.0e-12)
    off = false;

    if(!off && zero && t_off>0.0)
    {
    off = true;
    t_on = zc.t_cross + t_off;
    pto.t_off0 = in.t;          // from this step on (the zero lies in the last step)
    pto.t_off1 = t_on;          // to the stage time t_on
    }
}

// --- hydrodynamic table ---

bool pto_hydro_table::read(const std::string &file)
{
    std::ifstream in(file);

    if(!in.good())
    return false;

    w_.clear(); A_.clear(); B_.clear();

    std::string line;

    while(std::getline(in,line))
    {
        std::replace(line.begin(),line.end(),',',' ');

        size_t i = line.find_first_not_of(" \t\r");

        if(i==std::string::npos)
        continue;

        const char c0 = line[i];

        if(!(std::isdigit(static_cast<unsigned char>(c0)) || c0=='.' || c0=='-' || c0=='+'))
        continue;

        std::istringstream ss(line);
        double w, A, B;

        if(ss>>w>>A>>B && std::isfinite(w) && std::isfinite(A) && std::isfinite(B))
        {
        w_.push_back(w);
        A_.push_back(A);
        B_.push_back(B);
        }
    }

    // ascending in w
    std::vector<size_t> idx(w_.size());

    for(size_t n=0; n<idx.size(); ++n)
    idx[n]=n;

    std::sort(idx.begin(),idx.end(),[&](size_t a, size_t b){return w_[a]<w_[b];});

    std::vector<double> w2, A2, B2;

    for(size_t n : idx)
    {
    w2.push_back(w_[n]);
    A2.push_back(A_[n]);
    B2.push_back(B_[n]);
    }

    w_.swap(w2); A_.swap(A2); B_.swap(B2);

    return w_.size()>=2;
}

double pto_hydro_table::interp(const std::vector<double> &f, double w) const
{
    if(w<=w_.front())
    return f.front();

    if(w>=w_.back())
    return f.back();

    const size_t n = std::upper_bound(w_.begin(),w_.end(),w) - w_.begin();
    const double s = (w - w_[n-1])/(w_[n] - w_[n-1]);

    return (1.0-s)*f[n-1] + s*f[n];
}

double pto_hydro_table::A(double w) const
{
    return interp(A_,w);
}

double pto_hydro_table::B(double w) const
{
    return interp(B_,w);
}

double pto_bopt(double w, double m, double A, double Brad, double C)
{
    const double X = w*(m + A) - C/w;

    return sqrt(Brad*Brad + X*X);
}
