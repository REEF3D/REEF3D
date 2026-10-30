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

#include"coral/control/coral_turbine_controller.h"
#include<algorithm>
#include<cmath>

namespace coral
{

double turbine_controller_settings::slope25() const
{
    const double ws = w_sync();
    return w_rated > ws ? (P_rated/w_rated)/(w_rated - ws) : 0.0;
}

double turbine_controller_settings::w_trans() const
{
    const double s = slope25(), ws = w_sync();
    if(K <= 0.0)
        return ws;
    // where K w^2 meets the line slope25 (w - w_sync)
    return (s - std::sqrt(std::max(s*(s - 4.0*K*ws), 0.0)))/(2.0*K);
}

void turbine_controller_settings::pi_gains(double J, double O, double N, double dP, double wn, double z, double &KP, double &KI)
{
    KP = 2.0*J*O*z*wn/(N*(-dP));
    KI = J*O*wn*wn/(N*(-dP));
}

double turbine_controller::torque_law(double w, double pitch) const
{
    const turbine_controller_settings &s = set;
    if(w >= s.w_rated || pitch >= s.pitch_rgn3)
        return w > 0.0 ? s.P_rated/w : 0.0;
    if(w <= s.w_cutin)
        return 0.0;
    if(w < s.w_rgn2)
        return s.slope15()*(w - s.w_cutin);
    if(w < s.w_trans())
        return s.K*w*w;
    return s.slope25()*(w - s.w_sync());
}

void turbine_controller::initialize(double t, double w, double pitch, double torque)
{
    t_last_ = t;
    wf_ = w;
    th_ = std::clamp(pitch, set.theta_min, set.theta_max);
    const double GK = 1.0/(1.0 + th_/set.theta_K);
    I_ = set.KI > 0.0 ? th_/(GK*set.KI) : 0.0;
    T_ = set.torque_control ? std::min(torque_law(wf_, th_), set.T_max) : torque;
    wp_ = 0.0;
    wfc_ = w;
    started_ = true;
}

void turbine_controller::update(double t, double w, double v)
{
    const turbine_controller_settings &s = set;
    const double dt = t - t_last_;
    if(!started_ || dt <= 0.0 || (s.dt > 0.0 && dt < s.dt*(1.0 - 1.0e-9)))
        return;
    t_last_ = t;
    const double a = std::exp(-dt*s.omega_c);
    wf_ = (1.0 - a)*w + a*wf_;
    if(s.floating())
    {
        // the platform's part of the speed, exact for v held over the step; the compensated speed through the same filter
        const double eA = std::exp(s.float_A*dt);
        wp_ = eA*wp_ + (std::fabs(s.float_A*dt) > 1.0e-12 ? (eA - 1.0)/s.float_A : dt)*s.float_B*v;
        wfc_ = (1.0 - a)*(w - s.float_gain*wp_) + a*wfc_;
    }
    else
        wfc_ = wf_;
    // torque
    if(s.torque_control)
    {
        const double Tl = std::min(torque_law(wf_, th_), s.T_max);
        T_ += std::clamp(Tl - T_, -s.T_rate*dt, s.T_rate*dt);
        region_ = (wf_ >= s.w_rated || th_ >= s.pitch_rgn3) ? 3 : wf_ <= s.w_cutin ? 1 : wf_ < s.w_rgn2 ? 15
                  : wf_ < s.w_trans() ? 2 : 25;
    }
    // pitch
    if(s.pitch_control)
    {
        const double GK = 1.0/(1.0 + th_/s.theta_K);
        const double e = wfc_ - s.w_ref;
        I_ += e*dt;
        if(s.KI > 0.0)
            I_ = std::clamp(I_, s.theta_min/(GK*s.KI), s.theta_max/(GK*s.KI));
        const double cmd = std::clamp(GK*(s.KP*e + s.KI*I_), s.theta_min, s.theta_max);
        th_ += std::clamp(cmd - th_, -s.theta_rate*dt, s.theta_rate*dt);
        th_ = std::clamp(th_, s.theta_min, s.theta_max);
    }
}

void turbine_controller::save_state(std::vector<double> &v) const
{
    v.push_back(t_last_); v.push_back(wf_); v.push_back(I_); v.push_back(T_); v.push_back(th_);
    v.push_back(region_); v.push_back(started_ ? 1.0 : 0.0); v.push_back(wp_); v.push_back(wfc_);
}

bool turbine_controller::load_state(const std::vector<double> &v, size_t &k)
{
    if(k + 9 > v.size()) return false;
    t_last_ = v[k++]; wf_ = v[k++]; I_ = v[k++]; T_ = v[k++]; th_ = v[k++];
    region_ = (int)v[k++]; started_ = v[k++] > 0.5; wp_ = v[k++]; wfc_ = v[k++];
    return true;
}

controller_load::controller_load(const turbine_controller_settings &s, io f, const std::string &name)
    : ctrl_(s), io_(std::move(f)), name_(name) {}

void controller_load::apply()
{
    if(ctrl_.set.torque_control && io_.set_torque) io_.set_torque(ctrl_.torque());
    if(ctrl_.set.pitch_control && io_.set_pitch) io_.set_pitch(ctrl_.pitch());
}

void controller_load::initialize(const model &m, double t0)
{
    ctrl_.initialize(t0, io_.speed(m), io_.pitch ? io_.pitch(m) : 0.0, io_.torque ? io_.torque(m) : 0.0);
    apply();
}

void controller_load::commit(const model &m, double t, double)
{
    ctrl_.update(t, io_.speed(m), io_.nacelle_velocity ? io_.nacelle_velocity(m) : 0.0);
    apply();
}

bool controller_load::load_state(const std::vector<double> &v, size_t &k)
{
    if(!ctrl_.load_state(v, k)) return false;
    apply();
    return true;
}

void controller_load::channel_names(std::vector<std::string> &c) const
{
    for(const char *x : {"wf", "T", "pitch", "region"})
        c.push_back(name_ + "_" + x);
    if(ctrl_.set.floating())
        for(const char *x : {"wplat", "vnac"})
            c.push_back(name_ + "_" + x);
}

void controller_load::channel_values(const model &m, double, std::vector<double> &v) const
{
    v.push_back(ctrl_.filtered_speed()*30.0/M_PI);
    v.push_back(ctrl_.torque());
    v.push_back(ctrl_.pitch()*180.0/M_PI);
    v.push_back(ctrl_.region());
    if(ctrl_.set.floating())
    {
        v.push_back(ctrl_.platform_speed()*30.0/M_PI);
        v.push_back(io_.nacelle_velocity ? io_.nacelle_velocity(m) : 0.0);
    }
}

}
