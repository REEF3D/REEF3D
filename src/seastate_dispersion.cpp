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

#include"seastate_dispersion.h"
#include<cmath>

double seastate_wavenumber(double sig, double d)
{
    const double g = seastate_gravity;
    const double k0 = sig*sig/g;

    if(!(d>0.0))
    return 0.0;

    if(k0*d>30.0)
    return k0;

    // Guo (2002) explicit approximation as start value
    const double x = sig*std::sqrt(d/g);
    double k = k0*std::pow(1.0 - std::exp(-std::pow(x,2.4908)),-1.0/2.4908);

    // Newton on f(k) = g k tanh(kd) - sig^2
    for(int it=0; it<50; ++it)
    {
    const double t = std::tanh(k*d);
    const double f = g*k*t - sig*sig;
    const double df = g*t + g*k*d*(1.0 - t*t);
    const double dk = f/df;
    k -= dk;

        if(std::fabs(dk)<1.0e-14*k)
        break;
    }

    return k;
}

double seastate_cg(double sig, double k, double d)
{
    if(!(k>0.0))
    return 0.0;

    const double kd = std::fmin(k*d,30.0);
    const double n = (kd>=30.0) ? 0.5 : 0.5*(1.0 + 2.0*kd/std::sinh(2.0*kd));

    return n*sig/k;
}

double seastate_refraction(double sig, double k, double d)
{
    const double kd = k*d;

    if(!(kd>0.0) || kd>30.0)
    return 0.0;

    return sig/std::sinh(2.0*kd);
}
