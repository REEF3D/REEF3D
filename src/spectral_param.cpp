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

#include"spectral_param.h"
#include"spectral_grid.h"
#include<cmath>

void spectral_param::compute(const spectral_grid &g, const float *N)
{
    const double pi = 3.14159265358979323846;

    *this = spectral_param();

    if(N==nullptr || g.nbin==0)
    return;

    double a=0.0, b=0.0, epeak=0.0;

    for(int l=0; l<g.nsig; ++l)
    {
    const double sig = g.sig[l];
    double e1 = 0.0;   // E(sig_l) integrated over direction
    double ec = 0.0;
    double es = 0.0;

        for(int m=0; m<g.ndir; ++m)
        {
        const double E = sig*double(N[g.bin(l,m)]);
        e1 += E;
        ec += E*g.costh[m];
        es += E*g.sinth[m];
        }

    e1 *= g.dtheta;
    ec *= g.dtheta;
    es *= g.dtheta;

    m0  += e1*g.dsig[l];
    m1  += e1*sig*g.dsig[l];
    mm1 += e1/sig*g.dsig[l];
    a   += ec*g.dsig[l];
    b   += es*g.dsig[l];

        if(e1>epeak)
        {
        epeak = e1;
        lpeak = l;
        }
    }

    if(!(m0>0.0))
    {
    *this = spectral_param();
    return;
    }

    Hs   = 4.0*std::sqrt(m0);
    Tm01 = 2.0*pi*m0/m1;
    Tm10 = 2.0*pi*mm1/m0;
    Tp   = lpeak>=0 ? 2.0*pi/g.sig[lpeak] : 0.0;

    dir = std::atan2(b,a)*180.0/pi;
    if(dir<0.0)
    dir += 360.0;

    const double r = std::sqrt(a*a + b*b)/m0;
    spread = std::sqrt(std::fmax(2.0*(1.0 - r),0.0))*180.0/pi;
}
