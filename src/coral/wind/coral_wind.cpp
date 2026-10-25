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

#include"coral/wind/coral_wind.h"
#include"coral/wind/coral_wind_grid.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<sstream>

namespace coral
{

bool iec_turbine::set_class(const std::string &cls, const std::string &cat)
{
    if(cls=="I") V_ref = 50.0;
    else if(cls=="II") V_ref = 42.5;
    else if(cls=="III") V_ref = 37.5;
    else return false;
    if(cat=="A+") I_ref = 0.18;
    else if(cat=="A") I_ref = 0.16;
    else if(cat=="B") I_ref = 0.14;
    else if(cat=="C") I_ref = 0.12;
    else return false;
    return true;
}

double iec_turbine::sigma_etm(double V_hub) const
{
    const double c = 2.0, V_ave = 0.2*V_ref;
    return c*I_ref*(0.072*(V_ave/c + 3.0)*(V_hub/c - 4.0) + 10.0);
}

double iec_turbine::V_e50(double z) const
{
    return 1.4*V_ref*std::pow(std::max(z, 0.0)/z_hub, 0.11);
}

iec_event iec_event_from(const std::string &s)
{
    if(s=="eog") return iec_event::eog;
    if(s=="edc+") return iec_event::edc_pos;
    if(s=="edc-") return iec_event::edc_neg;
    if(s=="ecd+") return iec_event::ecd_pos;
    if(s=="ecd-") return iec_event::ecd_neg;
    if(s=="ews_v+") return iec_event::ews_v_pos;
    if(s=="ews_v-") return iec_event::ews_v_neg;
    if(s=="ews_h+") return iec_event::ews_h_pos;
    if(s=="ews_h-") return iec_event::ews_h_neg;
    return iec_event::none;
}

wind_field::wind_field(double U_ref, double z_ref, double direction) : U_ref_(U_ref), z_ref_(z_ref), dir_(direction) {}

void wind_field::set_series(const std::vector<series_row> &rows, double L_ref)
{
    series_ = rows;
    L_ref_ = L_ref > 0.0 ? L_ref : 1.0;
}

wind_field::series_row wind_field::row(double t) const
{
    if(series_.empty())
        return series_row{t, U_ref_, dir_, 0.0, alpha_, 0.0, 0.0, 0.0};
    if(t <= series_.front().t) return series_.front();
    if(t >= series_.back().t) return series_.back();
    const auto it = std::upper_bound(series_.begin(), series_.end(), t, [](double x, const series_row &r) {return x < r.t;});
    const series_row &a = *(it - 1), &b = *it;
    const double s = (t - a.t)/(b.t - a.t);
    auto l = [s](double x, double y) {return x + s*(y - x);};
    return series_row{t, l(a.U, b.U), l(a.dir, b.dir), l(a.w, b.w), l(a.alpha, b.alpha), l(a.hshear, b.hshear),
                      l(a.vshear, b.vshear), l(a.gust, b.gust)};
}

double wind_field::shear(double z, double alpha) const
{
    if(z0_ > 0.0)
        return z <= z0_ ? 0.0 : std::log(z/z0_)/std::log(z_ref_/z0_);
    if(z <= 0.0)
        return 0.0;
    return alpha==0.0 ? 1.0 : std::pow(z/z_ref_, alpha);
}

double wind_field::ramp(double t) const
{
    if(t_ramp_ <= 0.0 || t >= t_ramp_) return 1.0;
    if(t <= 0.0) return 0.0;
    return 0.5*(1.0 - std::cos(M_PI*t/t_ramp_));
}

namespace
{
    // the IEC hub speed before the event and sigma_1 of the normal turbulence model there
    struct event_ref {double V_hub, sigma1, Lambda1;};
}

double wind_field::speed(double eta, double z, double t) const
{
    const series_row r = row(t);
    double U = r.U*(shear(z, r.alpha) + r.hshear*eta/L_ref_ + r.vshear*(z - z_ref_)/L_ref_) + r.gust;
    if(event_==iec_event::none || t < t_event_)
        return U;
    const series_row r0 = row(t_event_);
    const event_ref e{r0.U*shear(iec_.z_hub, r0.alpha), iec_.sigma_ntm(r0.U*shear(iec_.z_hub, r0.alpha)), iec_.Lambda1()};
    const double tau = t - t_event_;
    switch(event_)
    {
        case iec_event::eog:
        {
            const double T = 10.5;
            if(tau <= T)
            {
                const double Vg = std::min(1.35*(iec_.V_e1(iec_.z_hub) - e.V_hub), 3.3*e.sigma1/(1.0 + 0.1*iec_.D/e.Lambda1));
                U -= 0.37*Vg*std::sin(3.0*M_PI*tau/T)*(1.0 - std::cos(2.0*M_PI*tau/T));
            }
            break;
        }
        case iec_event::ecd_pos: case iec_event::ecd_neg:
        {
            const double T = 10.0, Vcg = 15.0;
            U += tau <= T ? 0.5*Vcg*(1.0 - std::cos(M_PI*tau/T)) : Vcg;
            break;
        }
        case iec_event::ews_v_pos: case iec_event::ews_v_neg: case iec_event::ews_h_pos: case iec_event::ews_h_neg:
        {
            const double T = 12.0, beta = 6.4;
            if(tau <= T)
            {
                const bool vert = event_==iec_event::ews_v_pos || event_==iec_event::ews_v_neg;
                const double sgn = (event_==iec_event::ews_v_pos || event_==iec_event::ews_h_pos) ? 1.0 : -1.0;
                const double s = vert ? z - iec_.z_hub : eta;
                U += sgn*s/iec_.D*(2.5 + 0.2*beta*e.sigma1*std::pow(iec_.D/e.Lambda1, 0.25))*(1.0 - std::cos(2.0*M_PI*tau/T));
            }
            break;
        }
        default:
            break;
    }
    return U;
}

double wind_field::direction(double z, double t) const
{
    const series_row r = row(t);
    double th = r.dir + veer_*(z - z_ref_);
    if(event_==iec_event::none || t < t_event_)
        return th;
    const double tau = t - t_event_;
    const series_row r0 = row(t_event_);
    const double V = r0.U*shear(iec_.z_hub, r0.alpha);
    if(event_==iec_event::edc_pos || event_==iec_event::edc_neg)
    {
        const double T = 6.0;
        const double the = std::min(4.0*std::atan(iec_.sigma_ntm(V)/(V*(1.0 + 0.1*iec_.D/iec_.Lambda1()))), M_PI);
        const double s = event_==iec_event::edc_pos ? 1.0 : -1.0;
        th += s*(tau <= T ? 0.5*the*(1.0 - std::cos(M_PI*tau/T)) : the);
    }
    else if(event_==iec_event::ecd_pos || event_==iec_event::ecd_neg)
    {
        const double T = 10.0;
        const double thcg = V < 4.0 ? M_PI : std::min(M_PI, 4.0*M_PI/V);     // 720 deg / V_hub
        const double s = event_==iec_event::ecd_pos ? 1.0 : -1.0;
        th += s*(tau <= T ? 0.5*thcg*(1.0 - std::cos(M_PI*tau/T)) : thcg);
    }
    return th;
}

double wind_field::hub_speed(double t) const
{
    const series_row r = row(t);
    return ramp(t)*(r.U + r.gust);
}

Vec3 wind_field::velocity(const Vec3 &p, double t) const
{
    const double s = ramp(t);
    if(s==0.0)
        return Vec3::Zero();
    const series_row r = row(t);
    const double c0 = std::cos(r.dir), s0 = std::sin(r.dir);
    const double dx = p(0) - ox_, dy = p(1) - oy_;
    const double xi = c0*dx + s0*dy, eta = -s0*dx + c0*dy;
    const double U = speed(eta, p(2), t);
    const double th = direction(p(2), t);
    const double ct = std::cos(th), st = std::sin(th);
    const double cu = std::cos(upflow_), su = std::sin(upflow_);
    Vec3 v(U*cu*ct, U*cu*st, U*su + r.w);
    if(turb_)
    {
        const double uc = turb_->u_hub > 0.0 ? turb_->u_hub : std::max(r.U, 1.0e-6);
        const Vec3 f = turb_->local(eta, p(2), t + (turb_->xi0 - xi)/uc);
        v += Vec3(f(0)*ct - f(1)*st, f(0)*st + f(1)*ct, f(2));
    }
    return s*v;
}

bool wind_field::read_uniform(const std::string &file, std::vector<series_row> &rows, std::string &msg)
{
    std::ifstream in(file);
    if(!in)
    {
        msg = "cannot open '" + file + "'";
        return false;
    }
    rows.clear();
    std::string line;
    int n = 0;
    while(std::getline(in, line))
    {
        ++n;
        const size_t b = line.find_first_not_of(" \t");
        if(b==std::string::npos || line[b]=='!' || line[b]=='#')
            continue;
        std::istringstream ls(line);
        std::vector<double> v;
        double x;
        while(ls >> x) v.push_back(x);
        if(v.empty())
            continue;
        if(v.size() < 3)
        {
            msg = file + " line " + std::to_string(n) + ": at least t, U and the direction expected";
            return false;
        }
        v.resize(9, 0.0);
        if(!rows.empty() && v[0] <= rows.back().t)
        {
            msg = file + " line " + std::to_string(n) + ": the times must increase";
            return false;
        }
        // InflowWind: direction positive clockwise from above; CORAL counterclockwise
        rows.push_back(series_row{v[0], v[1], -v[2]*M_PI/180.0, v[3], v[5], v[4], v[6], v[7]});
    }
    if(rows.empty())
    {
        msg = file + ": no data";
        return false;
    }
    return true;
}

}
