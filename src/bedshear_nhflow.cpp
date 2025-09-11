/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"bedshear.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"turbulence.h"
#include"sliceint.h"

// NHFLOW
void bedshear::taubed(lexer *p, fdm_nhf*d, ghostcell *pgc, sediment_fdm *s)
{
    double uabs,cf,manning,tau;
    double density=p->W1;
    double U,V,W;
    
    k=0;
    SEDSLICELOOP
    {
        
        if(p->S16==1)
        {
        U = d->U[IJK];
        V = d->V[IJK];
        W = d->W[IJK];
        
        dist = 0.5*p->DZN[KP]*d->WL(i,j);
        
        uabs = sqrt(U*U + V*V + W*W);
        
        u_plus = (1.0/kappa)*log(30.0*(dist/ks));

        tau=density*(uabs*uabs)/pow((u_plus>0.0?u_plus:1.0e20),2.0);
        }
        
        if(p->S16==2)
        {
        U = d->U[IJK];
        V = d->V[IJK];
        W = d->W[IJK];
        
        dist = 0.5*p->DZN[KP]*d->WL(i,j);
        
        uabs = sqrt(U*U + V*V + W*W);
        
        u_plus = (1.0/kappa)*log(30.0*(dist/ks));

        tau = MAX(density*(uabs*uabs)/pow((u_plus>0.0?u_plus:1.0e20),2.0), density*d->KIN[IJK]*0.3);
        }
        
        if(p->S16==3)
        {
        double v_t,v_d;

        U = d->U[IJK];
        V = d->V[IJK];
        W = d->W[IJK];
        
        v_d = d->VISC[IJK];
        v_t = d->EV[IJK];
        
        dist = 0.5*p->DZN[KP]*d->WL(i,j);
        
        uabs = sqrt(U*U + V*V + W*W);
    
        tau=density*(v_d + v_t)*(uabs/dist);
        }
        
        if(p->S16==4)
        tau=density*d->KIN[IJK]*0.3;
        
        if(p->S16==6)
        {
        double Cval,wh;
        double bedlevel,waterlevel;
        int count=0;
        U=V=wh=0.0;
        KLOOP
        {

            U += d->U[IJK];
            V += d->V[IJK];
            ++count;
            wh+=p->DZN[KP]*d->WL(i,j);
        }

        U=U/double(count);
        U=V/double(count);

        Cval=18.0*log10((12.0*wh)/ks);

        u_abs = sqrt(U*U + V*V);
	
        tau = density*pow(sqrt(9.81)*(u_abs/Cval),2.0);
        }

    s->tau_eff(i,j) = tau;
    s->shearvel_eff(i,j) = sqrt(tau/density);
    s->shields_eff(i,j) = tau/((p->S22-density)*fabs(p->W22)*p->S20);
    }
}

void bedshear::taucritbed(lexer *p, fdm_nhf* d, ghostcell *pgc, sediment_fdm *s)
{
    double density = p->W1;
    
    SEDSLICELOOP
    {
    k = s->bedk(i,j);
    
    tauc = (p->S30*fabs(p->W22)*(p->S22-density))*p->S20*s->reduce(i,j);
  
    s->tau_crit(i,j) = tauc;
    s->shearvel_crit(i,j) = sqrt(tauc/density);
    s->shields_crit(i,j) = tauc/((p->S22-density)*fabs(p->W22)*p->S20);
    
    s->MOB(i,j) = s->shields_eff(i,j)/(fabs(s->shields_crit(i,j))>1.0e-10?s->shields_crit(i,j):1.0e10);
    }
}