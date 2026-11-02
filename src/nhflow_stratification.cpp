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

#include"nhflow_stratification.h"
#include"nhflow_momentum_SI.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include<cmath>

static inline int sc3(lexer *p, int a, int b, int k)
{
    return (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + k-p->kmin;
}

void nhflow_n2(lexer *p, fdm_nhf *d, double *N2)
{
    const int K = p->knoz;
    const int mg = increment::marge;
    const double g = fabs(p->W22), rho0 = p->A621_rho0;
    double r[1024], z[1024];

    for(int a=p->imin; a<p->imin+p->imax; ++a)
    for(int b=p->jmin; b<p->jmin+p->jmax; ++b)
    {
        const double h = d->WL(a,b);

        if(K<2 || K>1024 || !(h>p->A544))
        {
            for(int k=0; k<K; ++k)
            N2[sc3(p,a,b,k)] = 0.0;
            continue;
        }

        for(int k=0; k<K; ++k)
        {
            const int c = sc3(p,a,b,k);
            const double S = d->SALT!=nullptr ? d->SALT[c] : p->A621_s0;
            const double T = d->TEMP!=nullptr ? d->TEMP[c] : p->A621_t0;
            r[k] = nhflow_momentum_SI::eos(p,S,T);
            z[k] = d->bed(a,b) + p->ZP[k+mg]*h;
        }

        for(int k=0; k<K; ++k)
        {
            const int k0 = k>0 ? k-1 : 0, k1 = k<K-1 ? k+1 : K-1;
            N2[sc3(p,a,b,k)] = -(g/rho0)*(r[k1] - r[k0])/(z[k1] - z[k0]);
        }
    }
}

double nhflow_ri(lexer *p, fdm_nhf *d, const double *N2, int a, int b, int k)
{
    const int K = p->knoz;
    const int mg = increment::marge;
    const double h = d->WL(a,b);
    const double n2 = N2[sc3(p,a,b,k)];

    if(K<2 || !(h>p->A544) || n2<=0.0)
    return 0.0;

    const int k0 = k>0 ? k-1 : 0, k1 = k<K-1 ? k+1 : K-1;
    const double dz = (p->ZP[k1+mg] - p->ZP[k0+mg])*h;
    const double du = (d->U[sc3(p,a,b,k1)] - d->U[sc3(p,a,b,k0)])/dz;
    const double dv = p->j_dir==1 ? (d->V[sc3(p,a,b,k1)] - d->V[sc3(p,a,b,k0)])/dz : 0.0;
    const double s2 = du*du + dv*dv;

    return n2/(s2>1.0e-12 ? s2 : 1.0e-12);
}

double nhflow_fm(lexer *p, double Ri)
{
    return Ri>0.0 ? pow(1.0 + p->A625_am*Ri,p->A625_pm) : 1.0;
}

double nhflow_fh(lexer *p, double Ri)
{
    return Ri>0.0 ? pow(1.0 + p->A625_ah*Ri,p->A625_ph) : 1.0;
}
