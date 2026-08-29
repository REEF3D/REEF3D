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
    double uabs,cf,manning,tau_i,tau_eff,tau_s;
    double density=p->W1;
    double visc=p->W2;
    double us, usn;
    double U,V,W;
    
    SLICEBASELOOP
    {
    s->tau_eff(i,j) = 0.0;
    s->tau_crit(i,j) = 0.0;
    }
    
    
    SEDSLICELOOP
    {
        k=0;
        
        if(p->S16==1 && p->B22==1)
        {
        U = d->U[IJK];
        V = d->V[IJK];
        W = d->W[IJK];
        
        dist = 0.5*p->DZN[KP]*d->WL(i,j);
        
        uabs = sqrt(U*U + V*V + W*W);
        
        u_plus = (1.0/kappa)*log(30.0*(dist/s->ks_eff(i,j)));

        tau_eff=density*(uabs*uabs)/pow((u_plus>0.0?u_plus:1.0e20),2.0);
        
        us = 1.0e-6;
        
        for(int qn=0; qn<10; ++qn)
        {
        usn = us;
        
        us = uabs/((1.0/kappa) * log(dist*us/visc) + 5.5);
        
        us += 0.5*(usn - us);
        }
        tau_eff = density * us * us;
        }
        
        if(p->S16==1 && p->B22==2)
        {
        U = d->U[IJK];
        V = d->V[IJK];
        W = d->W[IJK];
        
        dist = 0.5*p->DZN[KP]*d->WL(i,j);
        
        uabs = sqrt(U*U + V*V + W*W);
        
        u_plus = (1.0/kappa)*log(30.0*(dist/s->ks_eff(i,j)));

        tau_eff = density*(uabs*uabs)/pow((u_plus>0.0?u_plus:1.0e20),2.0);
        }
        
        if(p->S16==2)
        {
        U = d->U[IJK];
        V = d->V[IJK];
        W = d->W[IJK];
        
        dist = 0.5*p->DZN[KP]*d->WL(i,j);
        
        uabs = sqrt(U*U + V*V + W*W);
        
        u_plus = (1.0/kappa)*log(30.0*(dist/s->ks_eff(i,j)));

        tau_eff = MAX(density*(uabs*uabs)/pow((u_plus>0.0?u_plus:1.0e20),2.0), density*d->KIN[IJK]*0.3);
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
    
        tau_eff = density*(v_d + v_t)*(uabs/dist);
        }
        
        if(p->S16==4)
        tau_eff = density*d->KIN[IJK]*0.3;
        
        if(p->S16==6)
        {
        double wh;
        double bedlevel,waterlevel;
        double cellcount=0.0;
        U=V=wh=0.0;
        KLOOP
        {

            U += d->U[IJK]*p->DZN[KP]*d->WL(i,j);
            V += d->V[IJK]*p->DZN[KP]*d->WL(i,j);
            cellcount += 1.0;
            wh+=p->DZN[KP]*d->WL(i,j);
        }
        
        U=U/wh;
        V=V/wh;
        
        u_abs = sqrt(U*U + V*V);

        //Cval=18.0*log10((12.0*wh)/s->ks_eff(i,j));	
        //tau_eff = density*pow(sqrt(9.81)*(u_abs/Cval),2.0);
        
        manning = pow(s->ks(i,j),1.0/6.0)/20.0;
        cf = pow(manning,2.0)/pow(d->WL(i,j),1.0/3.0);
        
        //cf = sqrt(9.81)/log(12.0*d->WL(i,j)/(s->ks(i,j)));
    
        tau_eff = p->W1*9.81*cf*u_abs*u_abs; 
        
        //cout<<"u_abs: "<<u_abs<<" cf: "<<cf<<" count: "<<cellcount<<" tau_eff: "<<tau_eff<<endl;
        }
        
        // tau_i for Ti
        U = d->U[IJK];
        V = d->V[IJK];
        W = d->W[IJK];
        
        dist = 0.5*p->DZN[KP]*d->WL(i,j);
        
        uabs = sqrt(U*U + V*V + W*W);
        
        u_plus = (1.0/kappa)*log(30.0*(dist/s->ks(i,j)));

        tau_i = density*(uabs*uabs)/pow((u_plus>0.0?u_plus:1.0e20),2.0);
        
        
        // blend with shallow water
        if(p->S38==1 && d->WL(i,j)<p->S39*p->A544)
        {
            double fac,wh;
            double bedlevel,waterlevel;
            U=V=wh=0.0;
            KLOOP
            {

            U += d->U[IJK]*p->DZN[KP]*d->WL(i,j);
            V += d->V[IJK]*p->DZN[KP]*d->WL(i,j);
            wh+=p->DZN[KP]*d->WL(i,j);
            }
        
            U=U/wh;
            V=V/wh;

            
            u_abs = sqrt(U*U + V*V);

            manning = pow(s->ks(i,j),1.0/6.0)/20.0;
            cf = pow(manning,2.0)/pow(d->WL(i,j),1.0/3.0);

            tau_s = p->W1*9.81*cf*u_abs*u_abs; 
            
            if(d->WL(i,j) < p->S39*p->A544 && d->WL(i,j) >= 0.5*p->S39*p->A544)
            fac = (d->WL(i,j) - 0.5*p->S39*p->A544)/(0.5*p->S39*p->A544);
            
            if(d->WL(i,j) < 0.5*p->S39*p->A544)
            fac = 0.0;
            
            tau_eff = fac*tau_eff + (1.0-fac)*tau_s;
            
            tau_i   = fac*tau_i   + (1.0-fac)*tau_s;

        }
        
        
        // add inertia bed shear stress
        if(p->S37==1)
        {
        
            double cm,cv,cn,ci,tau_acc;
            
            ci=1.0;
            cm = 1.5;
            cv = PI/6.0;
            cn = ci/(cm*cm);
            
            tau_acc = p->W1 * cm*cv*cn * p->S20 * p->S21 * d->dudt(i,j);
            
            //cout<<tau_acc<<endl;
            
            tau_eff += tau_acc;
            tau_i += tau_acc;
            
        }


    s->tau_eff(i,j) = tau_eff;
    s->tau_i(i,j) = tau_i;
    s->shearvel_eff(i,j) = sqrt(tau_eff/density);
    s->shields_eff(i,j) = tau_eff/((p->S22-density)*fabs(p->W22)*p->S20);
    }
}

void bedshear::taucritbed(lexer *p, fdm_nhf* d, ghostcell *pgc, sediment_fdm *s)
{
    double density = p->W1;
    
    SEDSLICELOOP
    {
    tauc = s->reduce(i,j) * (p->S30*fabs(p->W22)*(p->S22-density))*p->S20;
  
    s->tau_crit(i,j) = tauc;
    s->shearvel_crit(i,j) = sqrt(tauc/density);
    s->shields_crit(i,j) = p->S30*s->reduce(i,j);
    
    
    
    s->MOB(i,j) = s->shields_eff(i,j)/(fabs(s->shields_crit(i,j))>1.0e-10?s->shields_crit(i,j):1.0e10);
    }
}