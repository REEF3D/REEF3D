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

#include"coral/aero/coral_unsteady_aero.h"
#include<algorithm>
#include<cmath>

namespace coral
{

double lag_step(double y0, double u0, double u1, double tau, double dt)
{
    if(dt <= 0.0)
        return y0;
    if(!(tau > 1.0e-12*dt))
        return u1;
    const double x = dt/tau;
    const double om = -std::expm1(-x);          // 1 - E
    return u1 + (y0 - u0)*(1.0 - om) - (u1 - u0)*om/x;
}

double oye_tau1(double a_mean, double R, double U0)
{
    const double a = std::clamp(a_mean, 0.0, 0.5);
    return 1.1/(1.0 - 1.3*a)*R/std::max(U0, 1.0);
}

double oye_tau2(double tau1, double r_over_R)
{
    const double s = std::clamp(r_over_R, 0.0, 1.0);
    return (0.39 - 0.26*s*s)*tau1;
}

inflow_state inflow_steady(const double Wqs[2])
{
    inflow_state s;
    for(int j=0; j<2; ++j)
        s.Wqs[j] = s.Wint[j] = s.W[j] = Wqs[j];
    return s;
}

inflow_state inflow_step(const inflow_state &s0, const double Wqs[2], double tau1, double tau2, double k, double dt)
{
    if(dt <= 0.0)
        return s0;
    inflow_state s;
    for(int j=0; j<2; ++j)
    {
        // W_int + tau1 W_int' = W_qs + k tau1 W_qs' with W_qs linear over the step: the input W_qs + k tau1 slope
        const double sl = (Wqs[j] - s0.Wqs[j])/dt;
        s.Wqs[j] = Wqs[j];
        s.Wint[j] = lag_step(s0.Wint[j], s0.Wqs[j] + k*tau1*sl, Wqs[j] + k*tau1*sl, tau1, dt);
        s.W[j] = lag_step(s0.W[j], s0.Wint[j], s.Wint[j], tau2, dt);
    }
    return s;
}

namespace
{
    const double DEG = M_PI/180.0;

    double param(const airfoil &af, const char *key, double set, double def)
    {
        if(set >= 0.0) return set;
        if(!af.tables.empty())
        {
            auto it = af.tables.front().params.find(key);
            if(it != af.tables.front().params.end()) return it->second;
        }
        return def;
    }
}

stall_polar::stall_polar(const airfoil &af, double re, const dynamic_stall_settings &s) : af_(&af), re_(re)
{
    // the attached flow's indicial function: Jones' approximation of Wagner's function (Hansen, Gaunaa and Madsen
    // 2004), unless the airfoil file or the input gives the constants (AeroDyn's DEFAULT is 0.3, 0.7, 0.14, 0.53)
    A1 = param(af, "A1", s.A1, 0.165);
    A2 = param(af, "A2", s.A2, 0.335);
    b1 = param(af, "b1", s.b1, 0.0455);
    b2 = param(af, "b2", s.b2, 0.3);
    T_f0 = param(af, "T_f0", s.T_f0, s.model==1 ? 6.0 : 3.0);
    T_p = param(af, "T_p", s.T_p, 1.7);
    if(af.tables.empty())
        return;
    a_ = af.tables.front().alpha;
    const int n = (int)a_.size();
    std::vector<double> cl(n);
    for(int i=0; i<n; ++i)
        cl[i] = af.coefficients(a_[i], re).cl;
    // the lift slope: the largest slope between -20 and 20 deg (or C_lalpha of the file)
    double sl = 0.0;
    for(int i=0; i+1<n; ++i)
        if(a_[i] >= -20.0*DEG && a_[i+1] <= 20.0*DEG && a_[i+1] > a_[i])
            sl = std::max(sl, (cl[i+1] - cl[i])/(a_[i+1] - a_[i]));
    cl_alpha = param(af, "C_lalpha", -1.0, sl > 0.0 ? sl : 2.0*M_PI);
    // the zero lift angle: the upward zero crossing of cl nearest to 0 within +-20 deg (or alpha0 of the file)
    double a0 = 0.0, best = 1.0e9;
    for(int i=0; i+1<n; ++i)
        if(a_[i] >= -20.0*DEG && a_[i+1] <= 20.0*DEG && cl[i] <= 0.0 && cl[i+1] > 0.0)
        {
            const double z = a_[i] - cl[i]*(a_[i+1] - a_[i])/(cl[i+1] - cl[i]);
            if(std::fabs(z) < best) {best = std::fabs(z); a0 = z;}
        }
    auto it = af.tables.front().params.find("alpha0");
    alpha0 = it != af.tables.front().params.end() ? it->second*DEG : a0;
    cd0_ = af.coefficients(alpha0, re).cd;
    // f_st: from alpha0 outwards until 2 sqrt(cl/Cl_inv) - 1 first reaches 0, zero beyond
    f_.assign(n, 0.0);
    fs_.assign(n, 0.0);
    auto fraw = [&](int i)
    {
        const double ci = cl_inv(a_[i]);
        if(std::fabs(a_[i] - alpha0) < 1.0e-9) return 1.0;
        const double r = cl[i]/ci;
        if(!(r > 0.0)) return 0.0;
        const double v = 2.0*std::sqrt(r) - 1.0;
        return v <= 0.0 ? 0.0 : std::min(v*v, 1.0);
    };
    const int iu = int(std::lower_bound(a_.begin(), a_.end(), alpha0) - a_.begin());
    for(int i=iu; i<n; ++i)
    {
        f_[i] = fraw(i);
        if(f_[i] <= 0.0) break;
    }
    for(int i=iu-1; i>=0; --i)
    {
        f_[i] = fraw(i);
        if(f_[i] <= 0.0) break;
    }
    for(int i=0; i<n; ++i)
        fs_[i] = 1.0 - f_[i] < 1.0e-3 ? 0.5*cl[i] : (cl[i] - cl_inv(a_[i])*f_[i])/(1.0 - f_[i]);
}

double stall_polar::interp(const std::vector<double> &y, double a) const
{
    const int n = (int)a_.size();
    if(n==0) return 0.0;
    if(a <= a_.front()) return y.front();
    if(a >= a_.back()) return y.back();
    const int k = int(std::upper_bound(a_.begin(), a_.end(), a) - a_.begin()) - 1;
    const double s = (a - a_[k])/(a_[k+1] - a_[k]);
    return y[k] + s*(y[k+1] - y[k]);
}

airfoil_coef stall_polar::stat(double alpha) const {return af_ ? af_->coefficients(alpha, re_) : airfoil_coef{};}
double stall_polar::f_st(double alpha) const {return interp(f_, wrap_angle(alpha));}
double stall_polar::cl_fs(double alpha) const
{
    // from f_st (interpolated) and the static cl at alpha, so that f_st Cl_inv + (1 - f_st) Cl_fs = cl at any angle
    // (interpolating Cl_fs itself would not); near f_st = 1 the table's value; the static lift where f_st = 0
    const double a = wrap_angle(alpha);
    const double f = interp(f_, a);
    if(f <= 0.0) return stat(a).cl;
    if(1.0 - f > 1.0e-3) return (stat(a).cl - cl_inv(a)*f)/(1.0 - f);
    return interp(fs_, a);
}

double stall_Tu(const stall_input &in)
{
    if(!(in.U > 0.0)) return 50.0;
    return std::clamp(0.5*in.chord/in.U, 0.001, 50.0);
}

stall_state stall_steady(const stall_polar &p, int model, const stall_input &in)
{
    stall_state s;
    const double Tu = stall_Tu(in);
    s.a34 = in.alpha + Tu*in.omega;
    if(model==1)
    {
        s.fst = s.x[3] = p.f_st(s.a34);
        return s;
    }
    s.x[0] = p.A1*s.a34;
    s.x[1] = p.A2*s.a34;
    s.clp = p.cl_inv(s.a34) + M_PI*Tu*in.omega;
    s.x[2] = s.clp;
    s.fF = s.x[3] = p.f_st(s.clp/p.cl_alpha + p.alpha0);
    s.fst = p.f_st(s.a34);
    return s;
}

airfoil_coef stall_step(const stall_polar &p, int model, const stall_state &s0, const stall_input &in, double dt, stall_state &s)
{
    s = s0;
    const double Tu = stall_Tu(in);
    const double a34 = in.alpha + Tu*in.omega;
    airfoil_coef c;
    if(model==1)
    {
        const double fst = p.f_st(a34);
        const double f = std::clamp(lag_step(s0.x[3], s0.fst, fst, p.T_f0*Tu, dt), 0.0, 1.0);
        s.x[3] = f; s.fst = fst; s.a34 = a34;
        const airfoil_coef st = p.stat(in.alpha);
        c.cl = f*p.cl_inv(a34) + (1.0 - f)*p.cl_fs(a34);
        c.cd = st.cd;
        c.cm = st.cm;
    }
    else
    {
        s.x[0] = lag_step(s0.x[0], p.A1*s0.a34, p.A1*a34, Tu/p.b1, dt);
        s.x[1] = lag_step(s0.x[1], p.A2*s0.a34, p.A2*a34, Tu/p.b2, dt);
        const double aE = (1.0 - p.A1 - p.A2)*a34 + s.x[0] + s.x[1];
        const double clp = p.cl_alpha*(aE - p.alpha0) + M_PI*Tu*in.omega;
        s.x[2] = lag_step(s0.x[2], s0.clp, clp, p.T_p*Tu, dt);
        const double fF = p.f_st(s.x[2]/p.cl_alpha + p.alpha0);
        const double x4 = std::clamp(lag_step(s0.x[3], s0.fF, fF, p.T_f0*Tu, dt), 0.0, 1.0);
        s.x[3] = x4; s.a34 = a34; s.clp = clp; s.fF = fF; s.fst = p.f_st(a34);
        const airfoil_coef st = p.stat(aE);
        const double fE = p.f_st(aE);
        const double clc = x4*p.cl_inv(aE) + (1.0 - x4)*p.cl_fs(aE);
        const double dCd = 0.5*(std::sqrt(fE) - std::sqrt(x4)) - 0.25*(fE - x4);
        c.cl = clc + M_PI*Tu*in.omega;
        c.cd = st.cd + (a34 - aE)*clc + (st.cd - p.cd0())*dCd;
        c.cm = st.cm - 0.5*M_PI*Tu*in.omega;
    }
    // beyond 40 deg towards the static coefficients, static from 45 deg (AeroDyn's UACutOut)
    const double w = std::clamp((std::fabs(wrap_angle(in.alpha)) - 40.0*DEG)/(5.0*DEG), 0.0, 1.0);
    if(w > 0.0)
    {
        const airfoil_coef st = p.stat(in.alpha);
        c.cl += w*(st.cl - c.cl);
        c.cd += w*(st.cd - c.cd);
        c.cm += w*(st.cm - c.cm);
    }
    return c;
}

}
