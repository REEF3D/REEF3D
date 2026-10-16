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

#include"coral/hydro/coral_waves.h"
#include<cmath>
#include<random>
#include<stdexcept>

namespace coral
{

wave_field::wave_field(double depth, double g) : depth_(depth), g_(g) {}

double wave_field::wavenumber(double w, double d, double g)
{
    const double k0 = w*w/g;
    if(d <= 0.0 || k0*d > 20.0)
        return k0;
    // Newton on f(k) = g k tanh(k d) - w^2, from the Eckart/Fenton start value
    double k = k0/std::sqrt(std::tanh(k0*d));
    for(int i=0; i<50; ++i)
    {
        const double th = std::tanh(k*d);
        const double f = g*k*th - w*w;
        const double df = g*th + g*k*d*(1.0 - th*th);
        const double dk = f/df;
        k -= dk;
        if(std::fabs(dk) < 1.0e-15*k)
            break;
    }
    return k;
}

void wave_field::add_component(double a, double w, double beta, double phase)
{
    wave_component c;
    c.a = a;
    c.w = w;
    c.k = wavenumber(w, depth_, g_);
    c.beta = beta;
    c.phase = phase;
    comp_.push_back(c);
}

void wave_field::add_regular(double H, double T, double beta, double phase)
{
    add_component(0.5*H, 2.0*M_PI/T, beta, phase);
}

double wave_field::spectrum(wave_spectrum type, double w, double Hs, double Tp, double gamma)
{
    if(w <= 0.0)
        return 0.0;
    const double wp = 2.0*M_PI/Tp;
    auto shape = [&](double ww)
    {
        const double r = wp/ww;
        double s = std::pow(ww, -5.0)*std::exp(-1.25*r*r*r*r);
        if(type==wave_spectrum::jonswap && gamma != 1.0)
        {
            const double sig = ww <= wp ? 0.07 : 0.09;
            const double e = (ww - wp)/(sig*wp);
            s *= std::pow(gamma, std::exp(-0.5*e*e));
        }
        return s;
    };
    // normalization to m0 = Hs^2/16 (Simpson over 0.1 wp .. 20 wp; the tail beyond is below 1e-6)
    const int n = 4000;
    const double a = 0.1*wp, b = 20.0*wp, h = (b - a)/n;
    double m0 = shape(a) + shape(b);
    for(int i=1; i<n; ++i)
        m0 += (i%2 ? 4.0 : 2.0)*shape(a + i*h);
    m0 *= h/3.0;
    return Hs*Hs/16.0*shape(w)/m0;
}

void wave_field::add_spectrum(wave_spectrum type, double Hs, double Tp, double gamma, double beta,
                              double wmin, double wmax, int n, unsigned seed)
{
    if(n < 1 || wmax <= wmin)
        throw std::runtime_error("coral::wave_field: spectrum needs n >= 1 and wmax > wmin");
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    const double dw = (wmax - wmin)/n;
    for(int i=0; i<n; ++i)
    {
        const double w = wmin + (i + 0.25 + 0.5*uni(gen))*dw;
        const double a = std::sqrt(2.0*spectrum(type, w, Hs, Tp, gamma)*dw);
        add_component(a, w, beta, 2.0*M_PI*uni(gen));
    }
}

double wave_field::ramp(double t) const
{
    if(t_ramp_ <= 0.0 || t >= t_ramp_)
        return 1.0;
    if(t <= 0.0)
        return 0.0;
    return 0.5*(1.0 - std::cos(M_PI*t/t_ramp_));
}

double wave_field::elevation(double x, double y, double t) const
{
    double eta = 0.0;
    for(const auto &c : comp_)
        eta += c.a*std::cos(c.w*t - c.k*(x*std::cos(c.beta) + y*std::sin(c.beta)) + c.phase);
    return ramp(t)*eta;
}

void wave_field::kinematics(const Vec3 &p, double t, Vec3 &u, Vec3 &du, double &p_rho) const
{
    u.setZero();
    du.setZero();
    p_rho = 0.0;
    const double z = std::min(p(2), 0.0);
    const double d = depth_;
    for(const auto &c : comp_)
    {
        const double cb = std::cos(c.beta), sb = std::sin(c.beta);
        const double th = c.w*t - c.k*(p(0)*cb + p(1)*sb) + c.phase;
        double fh, fv, fp;    // depth functions of the horizontal and vertical velocity, of the pressure
        if(d <= 0.0 || c.k*d > 20.0)
        {
            fh = fv = fp = std::exp(c.k*z);
        }
        else
        {
            const double kz = c.k*(z + d);
            fh = std::cosh(kz)/std::sinh(c.k*d);
            fv = std::sinh(kz)/std::sinh(c.k*d);
            fp = std::cosh(kz)/std::cosh(c.k*d);
        }
        const double cs = std::cos(th), sn = std::sin(th);
        const double uh = c.a*c.w*fh*cs, duh = -c.a*c.w*c.w*fh*sn;
        u(0) += uh*cb;
        u(1) += uh*sb;
        u(2) += -c.a*c.w*fv*sn;
        du(0) += duh*cb;
        du(1) += duh*sb;
        du(2) += -c.a*c.w*c.w*fv*cs;
        p_rho += g_*c.a*fp*cs;
    }
    const double r = ramp(t);
    u *= r;
    du *= r;
    p_rho *= r;
}

double wave_field::hs() const
{
    double m0 = 0.0;
    for(const auto &c : comp_)
        m0 += 0.5*c.a*c.a;
    return 4.0*std::sqrt(m0);
}

}
