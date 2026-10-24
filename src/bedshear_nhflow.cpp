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
    
    
    if(p->S16==9)
    tau_wave_step(p);
    
    SEDSLICELOOP
    {
        k=0;
        tau_eff = 0.0;   // unsupported S 16 options: no shear instead of an uninitialised value
        
        if(p->S16==1 && p->B22==1)
        {
        U = d->U[IJK];
        V = d->V[IJK];
        W = d->W[IJK];
        
        dist = 0.5*p->DZN[KP]*d->WL(i,j);
        
        uabs = sqrt(U*U + V*V + W*W);
        
        u_plus = uplus(dist,s->ks_eff(i,j));

        tau_eff=density*(uabs*uabs)/(u_plus*u_plus);
        
        us = sqrt(tau_eff/density);
        
            for(int qn=0; qn<10; ++qn)
            {
            usn = us;
            
            us = uabs/((1.0/kappa) * log(MAX(dist*(us>1.0e-10?us:1.0e-10)/visc,1.0)) + 5.5);
            
            us += 0.5*(usn - us);
            }
        tau_eff = density * us * us;
        }
        
        if((p->S16==1 && p->B22==2) || p->S16==9)
        {
        U = d->U[IJK];
        V = d->V[IJK];
        W = d->W[IJK];
        
        dist = 0.5*p->DZN[KP]*d->WL(i,j);
        
        uabs = sqrt(U*U + V*V + W*W);
        
        u_plus = uplus(dist,s->ks_eff(i,j));

        tau_eff = density*(uabs*uabs)/(u_plus*u_plus);
        
        // tau_wave updates the running mean: called once (MAX is a macro)
        if(p->S16==9 && p->S59<1)
        {
        double tau_w = tau_wave(p,U,V,s->ks_eff(i,j),density);
        tau_eff = MAX(tau_eff, tau_w);
        }
        
        // S 59 1: max(tau_log, Soulsby DATA13), S 59 2: Soulsby DATA13
        if(p->S16==9 && p->S59>=1)
        {
        double tau_wc = tau_wave_current(p,U,V,s->ks_eff(i,j),density,dist);
        tau_eff = p->S59==1 ? MAX(tau_eff, tau_wc) : tau_wc;
        }
        }
        
        
        if(p->S16==2)
        {
        U = d->U[IJK];
        V = d->V[IJK];
        W = d->W[IJK];
        
        dist = 0.5*p->DZN[KP]*d->WL(i,j);
        
        uabs = sqrt(U*U + V*V + W*W);
        
        u_plus = uplus(dist,s->ks_eff(i,j));

        tau_eff = MAX(density*(uabs*uabs)/(u_plus*u_plus), density*d->KIN[IJK]*0.3);
        }
        
        // S 16 10: the bed shear stress of the CFD S 16 10 and the CPM bedload layer (bedshear::column_tau):
        // log law at the first cell centre, or S 63 cells above the bed, horizontal velocity, and
        // max(log law, 0.3 rho k) with k of the first cell; S 16 11: sqrt(log law * 0.3 rho k)
        if(p->S16==10 || p->S16==11)
        {
        double tx=0.0, ty=0.0;
        
        bedshear::column_tau(p,d,i,j,s->ks_eff(i,j),p->S63,p->S16==11?2:1,density,tx,ty);
        
        tau_eff = sqrt(tx*tx + ty*ty);
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
        
        k=0;
        
        U=wh>1.0e-20?U/wh:0.0;
        V=wh>1.0e-20?V/wh:0.0;
        
        u_abs = sqrt(U*U + V*V);

        //Cval=18.0*log10((12.0*wh)/s->ks_eff(i,j));	
        //tau_eff = density*pow(sqrt(9.81)*(u_abs/Cval),2.0);
        
        manning = pow(s->ks(i,j),1.0/6.0)/20.0;
        cf = pow(manning,2.0)/pow(MAX(d->WL(i,j),1.0e-6),1.0/3.0);
        
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
        
        u_plus = uplus(dist,s->ks(i,j));

        tau_i = density*(uabs*uabs)/(u_plus*u_plus);
        
        
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
        
            k=0;
            
            U=wh>1.0e-20?U/wh:0.0;
            V=wh>1.0e-20?V/wh:0.0;

            
            u_abs = sqrt(U*U + V*V);

            manning = pow(s->ks(i,j),1.0/6.0)/20.0;
            cf = pow(manning,2.0)/pow(MAX(d->WL(i,j),1.0e-6),1.0/3.0);

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
            
            // the inertia term is signed: no negative shear stress (sqrt below)
            tau_eff = MAX(tau_eff + tau_acc, 0.0);
            tau_i = MAX(tau_i + tau_acc, 0.0);
            
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
    }
    
    SEDSLICELOOP
    s->MOB(i,j) = s->shields_eff(i,j)/(fabs(s->shields_crit(i,j))>1.0e-10?s->shields_crit(i,j):1.0e10);
}
/*--------------------------------------------------------------------
bed shear stress of a bed column on the NHFLOW sigma grid (S 16 10), the same as bedshear::column_tau of CFD:
the first cell centre is half a cell above the bed; href > 0 takes the velocity href cells (of the first cell)
above the bed, at most half the water depth, linear between the cell centres (below the second cell centre
the first cell, as in CFD); horizontal velocity; k of the first cell
thin water: z_r at least e k_s/30 (u+ >= 1/kappa, as uplus and the NHFLOW wall drag), no singular shear
stress when the depth approaches zero; the CFD column has z_r >= half a cell
--------------------------------------------------------------------*/
bool bedshear::column_tau(lexer *p, fdm_nhf *d, int i, int j, double ks, double href, int kmode,
                          double rho, double &tx, double &ty)
{
    double wl = d->WL(i,j);
    int k=0;
    
    tx = ty = 0.0;
    ks = MAX(ks, 1.0e-6);
    
    if(wl<=1.0e-20)
    return false;
    
    double ur = d->U[IJK];
    double vr = p->j_dir==1 ? d->V[IJK] : 0.0;
    double zr = p->ZP[KP]*wl;
    double kin = d->KIN[IJK];
    
    if(href>0.0)
    {
        double zref = href*p->DZN[KP]*wl;
        
        zref = MIN(zref, 0.5*wl);
        zref = MAX(zref, zr);
        
        int kr = 0;
        while(kr+1<p->knoz && p->ZP[kr+1+marge]*wl<=zref)
        ++kr;
        
        if(kr>0 && kr+1<p->knoz)
        {
            double z0 = p->ZP[kr+marge]*wl, z1 = p->ZP[kr+1+marge]*wl;
            double f = (zref - z0)/(z1 - z0);
            
            k=kr;
            double u0 = d->U[IJK], v0 = p->j_dir==1 ? d->V[IJK] : 0.0;
            k=kr+1;
            double u1 = d->U[IJK], v1 = p->j_dir==1 ? d->V[IJK] : 0.0;
            
            ur = (1.0-f)*u0 + f*u1;
            vr = (1.0-f)*v0 + f*v1;
            zr = zref;
        }
        else if(kr>0)
        {
            k=kr;
            ur = d->U[IJK];
            vr = p->j_dir==1 ? d->V[IJK] : 0.0;
            zr = p->ZP[KP]*wl;
        }
    }
    
    zr = MAX(zr, 2.718281828459045*ks/30.0);
    
    double um = sqrt(ur*ur + vr*vr);
    double tau = tau_logk(um, zr, ks, kin, kmode, rho);
    
    tx = um>1.0e-12 ? tau*ur/um : 0.0;
    ty = um>1.0e-12 ? tau*vr/um : 0.0;
    
    return true;
}
