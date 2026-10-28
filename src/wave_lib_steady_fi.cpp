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

#include"wave_lib_steady_fi.h"
#include<cmath>
#include<algorithm>

void wave_lib_steady_fi::build(double L_, double c_, double h, double g, double t0_,
                               const std::function<double(double,double)> &uh,
                               const std::function<double(double,double)> &w,
                               const std::function<double(double)> &eta)
{
    L = L_;
    c = c_;
    t0 = t0_;

    const double pi = 3.14159265358979323846;

    std::vector<double> e(N), cs(N), sn(N);

    for(int i=0; i<N; ++i)
    {
    e[i] = eta(i*L/double(N));
    cs[i] = cos(2.0*pi*i/double(N));
    sn[i] = sin(2.0*pi*i/double(N));
    }

    const double emax = *std::max_element(e.begin(),e.end());
    const double emin = *std::min_element(e.begin(),e.end());

    zb = -h;
    zt = emax + 0.1*(emax-emin) + 1.0e-6;
    dz = (zt-zb)/double(nz-1);

    // the zero-mean antiderivative of u along s at every level: F_m = U_m / (i k_m), m = 1 .. N/2-1
    F.assign(size_t(nz)*N,0.0);
    std::vector<double> u(N), ar(N/2), ai(N/2);

    for(int kz=0; kz<nz; ++kz)
    {
        const double z = zb + kz*dz;

        for(int i=0; i<N; ++i)
        u[i] = uh(i*L/double(N),z);

        for(int m=1; m<N/2; ++m)
        {
            double re=0.0, im=0.0;

            for(int i=0; i<N; ++i)
            {
            const int q = (m*i)%N;
            re += u[i]*cs[q];
            im -= u[i]*sn[q];
            }

            // U_m / (i k_m) = (im - i re) / k_m
            const double km = 2.0*pi*m/L;
            ar[m] = im/km;
            ai[m] = -re/km;
        }

        for(int i=0; i<N; ++i)
        {
            double f=0.0;

            for(int m=1; m<N/2; ++m)
            {
            const int q = (m*i)%N;
            f += 2.0*(ar[m]*cs[q] - ai[m]*sn[q]);
            }

            F[size_t(kz)*N + i] = f/double(N);
        }
    }

    // Bernoulli constant at the surface
    b = 0.0;

    for(int i=0; i<N; ++i)
    {
    const double s = i*L/double(N);
    const double us = uh(s,e[i]);
    const double ws = w(s,e[i]);

    b += c*us - 0.5*(us*us + ws*ws) - g*e[i];
    }

    b /= double(N);

    ok = true;
}

// periodic Catmull-Rom interpolation of level kz at s
double wave_lib_steady_fi::row(int kz, double s) const
{
    const double ds = L/double(N);
    double xi = s/ds;
    const int i1 = int(floor(xi));
    const double r = xi - i1;

    auto f = [&](int i) { i = ((i%N)+N)%N; return F[size_t(kz)*N + i]; };

    const double f0=f(i1-1), f1=f(i1), f2=f(i1+1), f3=f(i1+2);

    return f1 + 0.5*r*(f2 - f0 + r*(2.0*f0 - 5.0*f1 + 4.0*f2 - f3 + r*(3.0*(f1 - f2) + f3 - f0)));
}

double wave_lib_steady_fi::fi(double s, double z, double t) const
{
    if(!ok)
    return 0.0;

    // position in the wave at t0, in [0, L)
    double xi = fmod(s - c*(t - t0), L);

    if(xi<0.0)
    xi += L;

    // levels: cubic in z, clamped to the table
    z = std::min(std::max(z,zb),zt);

    double zr = (z-zb)/dz;
    int k1 = std::min(int(floor(zr)),nz-2);
    const double r = zr - k1;

    auto g = [&](int k) { k = std::min(std::max(k,0),nz-1); return row(k,xi); };

    const double f0=g(k1-1), f1=g(k1), f2=g(k1+1), f3=g(k1+2);

    const double val = f1 + 0.5*r*(f2 - f0 + r*(2.0*f0 - 5.0*f1 + 4.0*f2 - f3 + r*(3.0*(f1 - f2) + f3 - f0)));

    return val + b*t;
}
