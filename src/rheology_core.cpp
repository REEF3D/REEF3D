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

#include"rheology_core.h"
#include<cmath>
#include<algorithm>

namespace rheo
{

// ---------------------------------------------------------------------------------------------
// 1. local laws
// ---------------------------------------------------------------------------------------------

double inertial_number(double gamma, double p_eff, double d, double rho_s)
{
    const double pp = std::max(p_eff,1.0e-12);
    return std::fabs(gamma)*d/std::sqrt(pp/rho_s);
}

double mu_I_plain(const law3D_params &q, double I)
{
    I = std::max(I,0.0);
    return q.mus + (q.mu2 - q.mus)*I/(q.I0 + I);
}

// nu = I mu'(I)/mu(I) of the plain law
static double nu_plain(const law3D_params &q, double I)
{
    const double dmu = (q.mu2 - q.mus)*q.I0/((q.I0 + I)*(q.I0 + I));
    return I*dmu/mu_I_plain(q,I);
}

double mu_I(const law3D_params &q, double I)
{
    I = std::max(I,0.0);

    if(q.I_N>0.0 && I<q.I_N)
    {
        if(I<=0.0)
        return 0.0;

        return std::sqrt(q.alpha/std::log(q.A/I));
    }

    return mu_I_plain(q,I);
}

double well_posedness(double mu, double nu)
{
    const double t = 1.0 - 0.5*nu;
    return 4.0*nu*nu - 4.0*nu + mu*mu*t*t;
}

static double C_plain(const law3D_params &q, double I)
{
    return well_posedness(mu_I_plain(q,I),nu_plain(q,I));
}

static double bisect_C(const law3D_params &q, double a, double b)
{
    // C changes sign in [a,b] (log bisection)
    const double ca = C_plain(q,a);

    for(int it=0; it<200; ++it)
    {
        const double m = std::sqrt(a*b);
        const double cm = C_plain(q,m);

        if((cm>0.0) == (ca>0.0))
        a = m;
        else
        b = m;

        if(b/a < 1.0 + 1.0e-13)
        break;
    }
    return std::sqrt(a*b);
}

double mu_I_setup(law3D_params &q)
{
    // neutral points of the plain law on a log scan, refined by bisection
    const int N = 4000;
    const double lo = -12.0, hi = 2.0;
    double I1 = 0.0, I2 = 0.0;
    double Iprev = std::pow(10.0,lo);
    double cprev = C_plain(q,Iprev);

    for(int n=1; n<=N; ++n)
    {
        const double I = std::pow(10.0,lo + (hi-lo)*double(n)/double(N));
        const double c = C_plain(q,I);

        if(I1==0.0 && cprev>0.0 && c<=0.0)
        I1 = bisect_C(q,Iprev,I);

        else if(I1>0.0 && I2==0.0 && cprev<=0.0 && c>0.0)
        I2 = bisect_C(q,Iprev,I);

        Iprev = I;
        cprev = c;
    }

    q.I_N = I1;

    if(I1>0.0)
    {
        const double muN = mu_I_plain(q,I1);
        q.A = I1*std::exp(q.alpha/(muN*muN));
    }
    else
    q.A = 0.0;

    return I2;
}

double yield_stress(const law3D_params &q, double p_eff)
{
    if(q.law==LAW_HB)
    return q.tauy;

    if(q.law==LAW_HBC)
    return std::max(0.0,q.tauy + q.tanphi*std::max(p_eff,0.0));

    if(q.law==LAW_MUI)
    return q.mus*std::max(p_eff,0.0);

    if(q.law==LAW_COULOMB_VP)
    return q.tanphi/std::sqrt(1.0 + q.tanphi*q.tanphi)*std::max(p_eff,0.0);

    return 0.0;
}

// (1 - exp(-m g))/g with its limit m at g -> 0
static double papanastasiou(double m, double g)
{
    if(m*g<1.0e-12)
    return m;

    return -std::expm1(-m*g)/g;
}

double viscosity(const law3D_params &q, double gamma, double p_eff)
{
    gamma = std::fabs(gamma);
    double mu = q.mu0;

    if(q.law==LAW_HB || q.law==LAW_HBC)
    {
        const double ty = yield_stress(q,p_eff);
        double visc = 0.0;

        if(q.n==1.0)
        visc = q.K;

        else if(gamma>1.0e-300)
        visc = q.K*std::pow(gamma,q.n-1.0);

        else
        visc = q.n<1.0 ? q.mu_max : 0.0;

        mu = ty*papanastasiou(q.m,gamma) + visc;
    }

    if(q.law==LAW_MUI)
    {
        const double pp = std::max(p_eff,0.0);

        if(gamma>1.0e-300)
        {
            const double I = inertial_number(gamma,pp,q.d,q.rho_s);
            mu = q.mu_min + mu_I(q,I)*pp/gamma;
        }
        else
        mu = q.mu_max;
    }

    if(q.law==LAW_COULOMB_VP)
    mu = q.mu_min + yield_stress(q,p_eff)*papanastasiou(q.m,gamma);

    if(!(mu==mu))
    mu = q.mu_max;

    return std::min(std::max(mu,q.mu_min),q.mu_max);
}

// ---------------------------------------------------------------------------------------------
// 2. depth-integrated closures
// ---------------------------------------------------------------------------------------------

static double sigma_eff(const basal_state &s)
{
    return std::max(0.0,1.0 - s.lambda)*std::max(s.sigma_n,0.0);
}

static double pf_mu_dyn(const basal_params &q, double h, double Fr)
{
    if(Fr<=0.0)
    return q.tz2;

    return q.tz1 + (q.tz2 - q.tz1)/(q.beta*h/(q.L*Fr) + 1.0);
}

static double pf_mu_start(const basal_params &q, double h)
{
    return q.tz3 + (q.tz2 - q.tz1)/(h/q.L + 1.0);
}

static double pf_mu(const basal_params &q, double h, double Fr)
{
    if(Fr>=q.beta)
    return pf_mu_dyn(q,h,Fr);

    const double ms = pf_mu_start(q,h);
    return std::pow(Fr/q.beta,q.kappa)*(pf_mu_dyn(q,h,q.beta) - ms) + ms;
}

double hb_mean_speed(double tau_b, double tauy, double K, double n, double h)
{
    if(tau_b<=tauy || h<=0.0)
    return 0.0;

    // simple shear with tau(z) = tau_b (1 - z/h): sheared layer h0, plug above it
    const double h0 = h*(1.0 - tauy/tau_b);
    const double a = std::pow(tau_b/(K*h),1.0/n);

    return a*n/((n+1.0)*h)*std::pow(h0,(n+1.0)/n)*(h - h0*n/(2.0*n+1.0));
}

double hb_basal_stress(double U, double tauy, double K, double n, double h)
{
    if(U<=0.0 || h<=0.0)
    return tauy;

    double lo = tauy;
    double hi = tauy + K*std::pow(std::max(3.0*U/h,1.0e-12),n) + 1.0e-12;

    while(hb_mean_speed(hi,tauy,K,n,h)<U)
    hi = tauy + 2.0*(hi - tauy);

    for(int it=0; it<200; ++it)
    {
        const double m = 0.5*(lo + hi);

        if(hb_mean_speed(m,tauy,K,n,h)<U)
        lo = m;
        else
        hi = m;

        if(hi - lo < 1.0e-13*hi)
        break;
    }
    return 0.5*(lo + hi);
}

static double hb_tauy(const basal_params &q, const basal_state &s)
{
    if(q.law==BASAL_HBC)
    return q.tauy + q.tandelta*sigma_eff(s);

    return q.tauy;
}

static double mui_basal(const basal_params &q, const basal_state &s, double U)
{
    // depth-averaged Bagnold profile: U = (2/5) I sqrt(phi p'/rho) h/d with p' = sigma_n' (basal)
    const double se = sigma_eff(s);

    if(se<=0.0)
    return 0.0;

    const double I = 2.5*U*q.d/(s.h*std::sqrt(q.phis*se/s.rho));
    return (q.mus + (q.mu2 - q.mus)*I/(q.I0 + I))*se;
}

double basal_stress(const basal_params &q, const basal_state &s, double U)
{
    U = std::max(U,0.0);

    switch(q.law)
    {
    case BASAL_COULOMB:
        return q.tandelta*sigma_eff(s);

    case BASAL_VOELLMY:
        return q.tandelta*sigma_eff(s) + s.rho*s.g*U*U/q.xi;

    case BASAL_PF:
    {
        const double Fr = U/std::sqrt(std::max(s.gn*s.h,1.0e-30));
        return pf_mu(q,s.h,Fr)*sigma_eff(s);
    }

    case BASAL_MUI:
        return mui_basal(q,s,U);

    case BASAL_HB:
    case BASAL_HBC:
        return hb_basal_stress(U,hb_tauy(q,s),q.K,q.n,s.h);

    default:
        return 0.0;
    }
}

double basal_stress_static(const basal_params &q, const basal_state &s)
{
    switch(q.law)
    {
    case BASAL_COULOMB:
    case BASAL_VOELLMY:
        return q.tandelta*sigma_eff(s);

    case BASAL_PF:
        return pf_mu_start(q,s.h)*sigma_eff(s);

    case BASAL_MUI:
        return q.mus*sigma_eff(s);

    case BASAL_HB:
    case BASAL_HBC:
        return hb_tauy(q,s);

    default:
        return 0.0;
    }
}

double implicit_speed(const basal_params &q, const basal_state &s, double Ustar, double dt)
{
    if(q.law==BASAL_OFF || Ustar<=0.0)
    return std::max(Ustar,0.0);

    const double c = dt/(s.rho*s.h);
    const double ts = basal_stress_static(q,s);

    // yield criterion: the driving momentum of the step does not overcome the static threshold
    if(Ustar<=c*ts)
    return 0.0;

    if(q.law==BASAL_COULOMB)
    return Ustar - c*ts;

    if(q.law==BASAL_VOELLMY)
    {
        // s + c (ts + b s^2) = Ustar
        const double bb = c*s.rho*s.g/q.xi;
        const double r = Ustar - c*ts;

        if(bb*r<1.0e-12)
        return r;

        return (std::sqrt(1.0 + 4.0*bb*r) - 1.0)/(2.0*bb);
    }

    if(q.law==BASAL_HB || q.law==BASAL_HBC)
    {
        // in the stress: U_hb(tau) + c tau = Ustar, monotone in tau
        const double ty = hb_tauy(q,s);
        double lo = ty, hi = ty + Ustar/c;

        for(int it=0; it<200; ++it)
        {
            const double m = 0.5*(lo + hi);

            if(hb_mean_speed(m,ty,q.K,q.n,s.h) + c*m < Ustar)
            lo = m;
            else
            hi = m;

            if(hi - lo < 1.0e-14*hi)
            break;
        }
        return std::max(0.0,Ustar - c*0.5*(lo + hi));
    }

    // general closure: bisection in the speed, f(s) = s + c tau(s) - Ustar
    double lo = 0.0, hi = Ustar;

    for(int it=0; it<200; ++it)
    {
        const double m = 0.5*(lo + hi);

        if(m + c*basal_stress(q,s,m) < Ustar)
        lo = m;
        else
        hi = m;

        if(hi - lo < 1.0e-14*Ustar)
        break;
    }
    return 0.5*(lo + hi);
}

// ---------------------------------------------------------------------------------------------
// 3. pore pressure
// ---------------------------------------------------------------------------------------------

double pore_relax(double lambda, double lambda_h, double D, double h, double dt)
{
    if(D<=0.0 || h<=0.0)
    return lambda;

    const double r = std::exp(-PI_*PI_*D*dt/(4.0*h*h));
    return lambda_h + (lambda - lambda_h)*r;
}

}
