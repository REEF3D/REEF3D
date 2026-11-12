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

#include"seastate_coupling.h"
#include"seastate_f.h"
#include"seastate_roller.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"fdm_seastate.h"
#include"lexer.h"
#include"slice.h"
#include"seastate_dispersion.h"
#include<cmath>

int seastate_coupling::wave_active(int ii, int jj) const
{
    return pwave->e->wet(ii,jj);
}

void seastate_coupling::bed_waves(lexer *p, slice &Hs, slice &Tm10, slice &Tp, slice &Urms, slice &Tbot, slice &ex, slice &ey, slice &Dr)
{
    fdm_seastate *e = pwave->e;
    const seastate_grid &g = *e->grid;
    const seastate_roller *proll = pwave->roller();
    const double pi = 3.14159265358979323846;

    SLICELOOP4
    {
    Hs(i,j)=Tm10(i,j)=Tp(i,j)=Urms(i,j)=Tbot(i,j)=ex(i,j)=ey(i,j)=0.0;
    Dr(i,j) = (proll!=nullptr) ? 0.0 : -1.0;

        if(e->wet(i,j)==0 || e->N->spec(i,j)==nullptr)
        continue;

    const float *N  = e->N->spec(i,j);
    const float *kc = e->kw->spec(i,j);
    const double d  = e->depth(i,j);

    // the integral parameters of the last wave step (seastate_f::parameters)
    Hs(i,j)   = e->Hs(i,j);
    Tm10(i,j) = e->Tm10(i,j);
    Tp(i,j)   = e->Tp(i,j);

    double m0u=0.0, m2u=0.0, a=0.0, b=0.0;

        for(int l=0; l<g.nsig; ++l)
        {
        const double sig = g.sig[l];
        const double kd  = double(kc[l])*d;

            // deep water: no orbital motion at the bed
            if(!(kd>0.0) || kd>25.0)
            continue;

        const double sh = std::sinh(kd);
        const double t  = sig*sig/(sh*sh);              // transfer function of the orbital velocity squared
        const double w  = g.dsig[l]*g.dtheta;

            for(int m=0; m<g.ndir; ++m)
            {
            const double E  = sig*double(N[g.bin(l,m)])*w*g.wth[m];     // [m^2]
            const double Eu = t*E;                                     // [m^2/s^2]

            m0u += Eu;
            m2u += Eu*sig*sig;
            a   += Eu*g.costh[m];
            b   += Eu*g.sinth[m];
            }
        }

        if(m0u>0.0)
        {
        Urms(i,j) = std::sqrt(m0u);
        Tbot(i,j) = (m2u>0.0) ? 2.0*pi*std::sqrt(m0u/m2u) : 0.0;

        const double r = std::sqrt(a*a + b*b);

            if(r>0.0)
            {
            ex(i,j) = a/r;
            ey(i,j) = b/r;
            }
        }

        // roller dissipation per area and density (surfbeat, A 748 1)
        if(proll!=nullptr)
        Dr(i,j) = seastate_gravity*const_cast<seastate_roller*>(proll)->Dr(i,j);
    }
}
