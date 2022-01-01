/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_exner.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void sediment_exner::non_equillibrium_solve(lexer* p,fdm* a, ghostcell *pgc, sediment_fdm *s)
{
    //SLICELOOP4
    //a->bedload(i,j) -= p->dtsed*(0.5*(a->P(i,j)+a->P(i+1,j))*dqx0(i,j) + 0.5*(a->Q(i,j)+a->Q(i,j+1))*dqy0(i,j));
    double rhosed=p->S22;
    double rhowat=p->W1;
    double g=9.81;
    double d50=p->S20;
    double visc=p->W2;
    double kappa=0.4;
    double ks=p->S21*d50;
    double Rstar=(rhosed-rhowat)/rhowat;
    double Ds= d50*pow((Rstar*g)/(visc*visc),1.0/3.0);
    double Ti;
    
    SLICELOOP4
    {
   
    Ti=MAX((s->shearvel_eff(i,j)*s->shearvel_eff(i,j)-s->shearvel_crit(i,j)*s->shearvel_crit(i,j))/(s->shearvel_crit(i,j)*s->shearvel_crit(i,j)),0.0);
        
    Ls = 3.0*d50*pow(Ds,0.6)*pow(Ti,0.9);
    
    
    Ls = 4000.0*MAX(shields_eff-shields_crit, 0.0)*d50;
    
    Ls = p->dtsed/p->DXM*sqrt(pow(0.5*(a->P(i,j)+a->P(i+1,j)),2.0) +  pow(0.5*(a->Q(i,j)+a->Q(i,j+1)),2.0));
    
    //cout<<Ls<<endl;

    a->bedload(i,j) += Ls*(dqx0(i,j) + dqy0(i,j));
    }
    
    SLICELOOP4
    {

    Ti=MAX((s->shearvel_eff(i,j)*s->shearvel_eff(i,j)-s->shearvel_crit(i,j)*s->shearvel_crit(i,j))/(s->shearvel_crit(i,j)*s->shearvel_crit(i,j)),0.0);
        
    Ls = 3.0*d50*pow(Ds,0.6)*pow(Ti,0.9);
    
    Ls = 4000.0*MAX(shields_eff-shields_crit, 0.0)*d50;
    
    Ls = p->dtsed/p->DXM*sqrt(pow(0.5*(a->P(i,j)+a->P(i+1,j)),2.0) +  pow(0.5*(a->Q(i,j)+a->Q(i,j+1)),2.0));
    
    
    KLOOP
    a->test(i,j,k) = Ls*(dqx0(i,j) + dqy0(i,j));
    }
    pgc->start4a(p,a->test,1);
    
}