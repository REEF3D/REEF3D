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

#include"seastate_structure.h"
#include"seastate_grid.h"
#include"seastate_dispersion.h"
#include<algorithm>
#include<cmath>
#include<complex>
#include<vector>

double seastate_dangremond(double Rc, double Hs, double Tp, double slope, double B)
{
    const double pi = 3.14159265358979323846;
    const double L0p = std::max(1.0e-8,1.5613*Tp*Tp);
    const double xi = std::tan(slope*pi/180.0)/std::sqrt(Hs/L0p);
    const double fvh = Rc/Hs, bvh = B/Hs;
    double t;

    if(bvh==0.0)
    t = std::min(std::max(-0.40*fvh,0.075),0.9);
    else if(bvh<8.0)
    t = std::min(std::max(-0.40*fvh + 0.64*std::pow(bvh,-0.31)*(1.0-std::exp(-0.50*xi)),0.075),0.9);
    else if(bvh>12.0)
    {
    t = -0.35*fvh + 0.51*std::pow(bvh,-0.65)*(1.0-std::exp(-0.41*xi));
    t = std::max(std::min(t,0.93-0.006*bvh),0.05);
    }
    else
    {
    double f1 = -0.40*fvh + 0.64*std::pow(8.0,-0.31)*(1.0-std::exp(-0.50*xi));
    f1 = std::min(std::max(f1,0.075),0.9);
    double f2 = -0.35*fvh + 0.51*std::pow(12.0,-0.65)*(1.0-std::exp(-0.41*xi));
    f2 = std::min(std::max(f2,0.05),0.858);
    t = 3.0*f1 - 2.0*f2 + bvh*(f2-f1)/4.0;
    }

    return t;
}

double seastate_porous(const seastate_grid &g, const double *E, double h, double B, double n, double D50,
                       float *kt2, float *kr2)
{
    typedef std::complex<double> cplx;
    const double grav = seastate_gravity, nu = 1.0e-6, pi = 3.14159265358979323846;
    const cplx I(0.0,1.0);
    const int ns = g.nsig;

    const double fa = 1000.0*(1.0-n)*(1.0-n)*nu/(n*n*n*grav*D50*D50);
    const double fb = 1.1*(1.0-n)/(n*n*n*grav*D50);
    const double s = 1.0 + 0.34*(1.0-n)/n;
    const double cq = std::sqrt(8.0/pi);

    // Gauss points on the width of the structure (mean of |q|^2)
    const double gx[4] = {0.0694318442,0.3300094782,0.6699905218,0.9305681558};
    const double gw[4] = {0.1739274226,0.3260725774,0.3260725774,0.1739274226};

    std::vector<cplx> T(ns), R(ns), kprev(ns);
    std::vector<double> fprev(ns,0.0);
    double qrms = 0.0;

    for(int it=0; it<50; ++it)
    {
    double q2 = 0.0;

        for(int l=0; l<ns; ++l)
        {
        const double sig = g.sig[l];
        const double k = seastate_wavenumber(sig,h);
        const double f = n*grav*(fa + cq*fb*qrms)/sig;

        // complex dispersion k_s tanh(k_s h) = sig^2 (s + i f)/g: the progressive mode, continued from the root
        // without resistance (real, k tanh(kh) = sig^2 s/g) in steps of the resistance (other starting values
        // can converge to a damped evanescent mode), then from the previous fixed-point iteration
        auto newton = [&](cplx z, const cplx &r)
        {
            for(int m=0; m<60; ++m)
            {
            const cplx th = std::tanh(z*h);
            const cplx dz = (z*th - r)/(th + z*h*(1.0-th*th));
            z -= dz;
            if(std::abs(dz)<1.0e-12*std::abs(z))
            break;
            }
            return z;
        };

        const int nst = (it==0) ? 32 : 4;
        const double f0 = (it==0) ? 0.0 : fprev[l];
        cplx ks = (it==0) ? cplx(seastate_wavenumber(sig*std::sqrt(s),h),0.0) : kprev[l];
            for(int st=1; st<=nst; ++st)
            {
            const double fs = f0 + (f-f0)*double(st)/double(nst);
            ks = newton(ks,sig*sig*(s + I*fs)/grav);
            }
        if(ks.imag()<0.0)
        ks = -ks;
        kprev[l] = ks;
        fprev[l] = f;

        const cplx gam = n*k/ks;
        const cplx e1 = std::exp(I*ks*B), e2 = e1*e1;
        const cplx Dn = (1.0+gam)*(1.0+gam) - (1.0-gam)*(1.0-gam)*e2;
        T[l] = 4.0*gam*e1/Dn;
        R[l] = (1.0-gam*gam)*(1.0-e2)/Dn;

        // discharge velocity in the structure, depth-averaged, per unit incident amplitude:
        // q(x) = gamma sig/(k h) (C e^(i ks x) - D e^(-i ks x)), C = 2 (1+gamma)/Dn, D e^(-i ks x) = -2 (1-gamma) e^(i ks (2B-x))/Dn
        const double a2 = 2.0*E[l]*g.dsig[l];
        double m2 = 0.0;
            for(int gp=0; gp<4; ++gp)
            {
            const double x = gx[gp]*B;
            const cplx qx = gam*sig/(k*h)*(2.0*(1.0+gam)*std::exp(I*ks*x) + 2.0*(1.0-gam)*std::exp(I*ks*(2.0*B-x)))/Dn;
            m2 += gw[gp]*std::norm(qx);
            }
        q2 += 0.5*a2*m2;
        }

    const double qn = std::sqrt(q2);
    const double qr = 0.5*(qrms+qn);
    const bool done = std::fabs(qr-qrms)<=1.0e-5*std::max(qr,1.0e-12);
    qrms = (it==0) ? qn : qr;
    if(done && it>0)
    break;
    }

    for(int l=0; l<ns; ++l)
    {
    const double t2 = std::min(std::norm(T[l]),1.0);
    kt2[l] = float(t2);
    kr2[l] = float(std::min(std::norm(R[l]),1.0-t2));
    }

    return qrms;
}

