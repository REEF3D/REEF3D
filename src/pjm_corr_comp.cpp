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
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#include"pjm_corr.h"
#include"lexer.h"
#include"fdm.h"
#include"interface_width.h"

void pjm_corr::comp_ini(lexer *p)
{
    if(p->W32>0 && (p->N40==4 || p->N40==24 || p->N40==44))
    {
    if(p->mpirank==0)
    cout<<"W 32 requires an SSP Runge-Kutta scheme, N 40 2 or 3"<<endl;

    exit(1);
    }

    if(p->W32>0 && p->mpirank==0 && p->amr_patch==0)
    cout<<"air compressibility  gamma: "<<p->W33<<"  band: "<<p->W34<<"  layer: "<<p->W35<<endl;
}

double pjm_corr::comp_weight(double s)
{
    if(s<=comp_a)
    return 0.0;

    if(s<=comp_a+comp_r)
    return comp_W*(s-comp_a)/comp_r;

    if(s<=comp_a+2.0*comp_r)
    return comp_W + (1.0-comp_W)*(s-comp_a-comp_r)/comp_r;

    return 1.0;
}

void pjm_corr::comp_terms(lexer *p, fdm *a, field &u, field &v, field &w, double alpha)
{
    cstage = 1.0;

    if(fabs(alpha-0.25)<1.0e-12)
    cstage = 0.5;

    fac_t = 1.0/(alpha*cstage*p->dt*p->dt);
    fac_c = 1.0/(alpha*p->dt);

    n=0;
    LOOP
    {
    wt = 0.0;

        if(a->phi(i,j,k)<0.0)
        {
        epsi = interface_width(p,a->phi,i,j,k);
        hn = epsi/p->F45;

        comp_a = epsi + hn;
        comp_r = 2.0*hn;

        if(p->W34>0.0)
        comp_a = p->W34;

        if(p->W35>0.0)
        comp_r = p->W35;

        comp_W = 1.5 + comp_a/comp_r;

        wt = comp_weight(-a->phi(i,j,k));
        }

        if(wt>0.0)
        {
        pabs = MAX(p->W36 + a->press(i,j,k), 0.01*p->W36);
        kappa = wt/(p->W33*pabs);

        a->M.p[n] += kappa*fac_t;
        a->rhsvec.V[n] -= kappa*(a->press(i,j,k)-pn(i,j,k))*fac_t;

            if(p->W32==2)
            {
            uc = 0.5*(u(i-1,j,k)+u(i,j,k));
            wc = 0.5*(w(i,j,k-1)+w(i,j,k));

            ugradp = uc*(a->press(i+1,j,k)-a->press(i-1,j,k))/(p->DXP[IP]+p->DXP[IM1])
                   + wc*(a->press(i,j,k+1)-a->press(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]);

                if(p->j_dir==1)
                {
                vc = 0.5*(v(i,j-1,k)+v(i,j,k));
                ugradp += vc*(a->press(i,j+1,k)-a->press(i,j-1,k))/(p->DYP[JP]+p->DYP[JM1]);
                }

            a->rhsvec.V[n] -= kappa*ugradp*fac_c;
            }
        }
    ++n;
    }
}
