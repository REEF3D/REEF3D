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
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"turbulence.h"
#include"sliceint.h"

#define HP (fabs(b->hp(i,j))>1.0e-20?b->hp(i,j):1.0e20)

// SFLOW
void bedshear::taubed(lexer *p, fdm2D *b, ghostcell *pgc, sediment_fdm *s)
{
    double uabs,cf,manning,tau;
    double U,V;
    
    SLICELOOP4
    {
    U = 0.5*(s->P(i,j) + s->P(i-1,j));
    V = 0.5*(s->Q(i,j) + s->Q(i,j-1));

    uabs = sqrt(U*U + V*V);
    
    manning = pow(p->S21*s->ks(i,j),1.0/6.0)/20.0;
    
    if(p->S16==1)
    {
    cf = pow(manning,2.0)/pow(HP,1.0/3.0);
    
    tau = p->W1*9.81*cf*uabs*uabs; 
    }
    
    if(p->S16==2)
    {    
    cf = 2.5*log(12.0*b->hp(i,j)/(p->S21*s->ks(i,j)));
    
    tau = p->W1*9.81*uabs*uabs/(fabs(cf*cf)>1.0e-20?(cf*cf):1.0e20); 
    }
    
    if(p->S16==3)
    {    
    cf = 2.5*log(12.0*b->hp(i,j)/(p->S21*s->ks(i,j)));
    
    tau = p->W1*9.81*uabs*uabs/(fabs(cf*cf)>1.0e-20?(cf*cf):1.0e20); 
    }
    
    if(p->S16==4)
    {
    tau=p->W1*b->kin(i,j)*0.3;
    }
    
    s->tau_eff(i,j) = taueff_loc(i,j) = tau;
    s->shearvel_eff(i,j) = sqrt(tau/p->W1);
    s->shields_eff(i,j) = tau/(p->W1*((p->S22-p->W1)/p->W1)*fabs(p->W22)*p->S20);
    }
}

void bedshear::taucritbed(lexer *p, fdm2D *b, ghostcell *pgc, sediment_fdm *s)
{
	double r;
    
    SLICELOOP4
    {
    tauc = (p->S30*fabs(p->W22)*(p->S22-p->W1))*p->S20*s->reduce(i,j);
  
    s->tau_crit(i,j) = taucrit_loc(i,j) = tauc;
    s->shearvel_crit(i,j) = sqrt(tauc/p->W1);
    s->shields_crit(i,j) = tauc/(p->W1*((p->S22-p->W1)/p->W1)*fabs(p->W22)*p->S20);
    }
}
