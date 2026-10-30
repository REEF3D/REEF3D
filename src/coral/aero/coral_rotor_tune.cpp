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

#include"coral/aero/coral_rotor_tune.h"
#include<cmath>
#include<cstdio>
#include<fstream>

namespace coral
{

rotor_tuning tune_rotor(const rotor_def &rd, double rho, double P_r, double O_r, double U_out, const bem_options &o)
{
    rotor_tuning tu;
    const double R = rd.R_tip;
    auto P = [&](double U, double Om, double th) {return rotor_steady(rd, U, Om, th, rho, o);};
    // Cp(lambda) at zero pitch (U 10 m/s), the maximum refined by a parabola
    const double U0 = 10.0;
    double best = -1.0, lb = 0.0;
    const double dl = 0.05;
    for(double l=1.0; l<=16.0; l+=dl)
    {
        const double cp = P(U0, l*U0/R, 0.0).Cp;
        if(cp > best) {best = cp; lb = l;}
    }
    {
        const double c0 = P(U0, (lb - dl)*U0/R, 0.0).Cp, c2 = P(U0, (lb + dl)*U0/R, 0.0).Cp;
        const double den = c0 - 2.0*best + c2;
        const double x = den < 0.0 ? 0.5*(c0 - c2)/den : 0.0;
        tu.lambda_opt = lb + x*dl;
        tu.cp_max = P(U0, tu.lambda_opt*U0/R, 0.0).Cp;
    }
    tu.K = 0.5*rho*M_PI*std::pow(R, 5)*tu.cp_max/std::pow(tu.lambda_opt, 3);
    auto omega = [&](double U) {return std::min(tu.lambda_opt*U/R, O_r);};
    // rated wind speed: P(U, Omega(U), 0) = P_rated
    double a = 2.0, b = U_out;
    if(P(b, omega(b), 0.0).P < P_r)
        tu.U_rated = U_out;
    else
    {
        for(int i=0; i<80; ++i)
        {
            const double c = 0.5*(a + b);
            (P(c, omega(c), 0.0).P < P_r ? a : b) = c;
        }
        tu.U_rated = 0.5*(a + b);
    }
    const double dth = 0.1*M_PI/180.0;
    std::vector<double> th_fit, d_fit;
    for(double U=3.0; U<=U_out + 1.0e-9; U+=0.5)
    {
        rotor_tuning::point p;
        p.U = U;
        p.Omega = omega(U);
        p.theta = 0.0;
        if(U > tu.U_rated)
        {
            // the pitch towards feather that holds the rated power
            double lo = 0.0, hi = 50.0*M_PI/180.0;
            for(int i=0; i<80; ++i)
            {
                const double c = 0.5*(lo + hi);
                (P(U, p.Omega, c).P > P_r ? lo : hi) = c;
            }
            p.theta = 0.5*(lo + hi);
        }
        const rotor_result r = P(U, p.Omega, p.theta);
        p.P = r.P; p.T = r.T; p.Cp = r.Cp;
        p.dPdtheta = (P(U, p.Omega, p.theta + dth).P - P(U, p.Omega, std::max(p.theta - dth, 0.0)).P)/(p.theta + dth - std::max(p.theta - dth, 0.0));
        tu.curve.push_back(p);
        if(U > tu.U_rated && p.theta >= 1.0*M_PI/180.0)
        {
            th_fit.push_back(p.theta);
            d_fit.push_back(p.dPdtheta);
        }
    }
    // dP/dtheta = a + b theta
    const int n = (int)th_fit.size();
    if(n >= 2)
    {
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for(int i=0; i<n; ++i) {sx += th_fit[i]; sy += d_fit[i]; sxx += th_fit[i]*th_fit[i]; sxy += th_fit[i]*d_fit[i];}
        const double bb = (n*sxy - sx*sy)/(n*sxx - sx*sx), aa = (sy - bb*sx)/n;
        if(aa < 0.0 && bb < 0.0)
        {
            tu.dPdtheta0 = aa;
            tu.theta_K = aa/bb;
        }
        else
        {
            // the line does not reach zero pitch with a negative sensitivity (a rotor whose power rises when it starts
            // to pitch at rated): the sensitivity at the first fitted pitch as dP/dtheta(0), theta_K from the slope
            // (|dP/dtheta| over-estimated by |b| theta_1: lower, safe gains near rated)
            tu.dPdtheta0 = aa + bb*th_fit.front();
            tu.theta_K = bb < 0.0 ? tu.dPdtheta0/bb : 1.0e30;
        }
    }
    return tu;
}

rotor_sensitivity rotor_sensitivities(const rotor_def &rd, const rotor_tuning &tu, double rho, double P_r, double O_r, double U,
                                      const bem_options &o)
{
    auto P = [&](double u, double Om, double th) {return rotor_steady(rd, u, Om, th, rho, o);};
    rotor_sensitivity s;
    s.U = U;
    s.Omega = std::min(tu.lambda_opt*U/rd.R_tip, O_r);
    if(U > tu.U_rated)
    {
        double lo = 0.0, hi = 50.0*M_PI/180.0;
        for(int i=0; i<80; ++i)
        {
            const double c = 0.5*(lo + hi);
            (P(U, s.Omega, c).P > P_r ? lo : hi) = c;
        }
        s.theta = 0.5*(lo + hi);
    }
    const rotor_result r = P(U, s.Omega, s.theta);
    s.T = r.T; s.Q = r.Q;
    const double dU = 0.01*U, dth = 1.0e-4, dO = 1.0e-3*s.Omega;
    const rotor_result u1 = P(U + dU, s.Omega, s.theta), u0 = P(U - dU, s.Omega, s.theta);
    const rotor_result t1 = P(U, s.Omega, s.theta + dth), t0 = P(U, s.Omega, s.theta - dth);
    const rotor_result o1 = P(U, s.Omega + dO, s.theta), o0 = P(U, s.Omega - dO, s.theta);
    s.dT_dU = (u1.T - u0.T)/(2.0*dU); s.dQ_dU = (u1.Q - u0.Q)/(2.0*dU);
    s.dT_dtheta = (t1.T - t0.T)/(2.0*dth); s.dQ_dtheta = (t1.Q - t0.Q)/(2.0*dth);
    s.dT_dOmega = (o1.T - o0.T)/(2.0*dO); s.dQ_dOmega = (o1.Q - o0.Q)/(2.0*dO);
    return s;
}

bool rotor_tuning::write(const std::string &file) const
{
    std::ofstream o(file);
    if(!o) return false;
    char s[256];
    o<<"# REEF3D::CORAL rotor tuning (steady BEM, axial inflow): Cp_max "<<cp_max<<" at tip speed ratio "<<lambda_opt
     <<", K "<<K<<" N m s^2 (rotor side), rated wind "<<U_rated<<" m/s, dP/dtheta(0) "<<dPdtheta0<<" W/rad, theta_K "
     <<theta_K*180.0/M_PI<<" deg\n";
    o<<"# U[m/s] Omega[rpm] pitch[deg] P[W] T[N] Cp dP/dtheta[W/rad]\n";
    for(const point &p : curve)
    {
        std::snprintf(s, sizeof(s), "%6.2f %10.5f %10.5f %14.6e %14.6e %9.6f %14.6e\n", p.U, p.Omega*30.0/M_PI, p.theta*180.0/M_PI,
                      p.P, p.T, p.Cp, p.dPdtheta);
        o<<s;
    }
    return (bool)o;
}

}
