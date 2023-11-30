/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

bedshear::bedshear(lexer *p, turbulence *ppturb) : norm_vec(p), ks(p->S20*p->S21), kappa(0.4), taueff_loc(p), taucrit_loc(p)
{
    tau=0.0;
    tauc=0.0;
    pturb=ppturb;
}

bedshear::~bedshear()
{
}

void bedshear::taubed(lexer *p, fdm * a, ghostcell *pgc, sediment_fdm *s)
{
	int count;
	double zval,fac,topoval,taukin,tauvel,density;
    
    SLICELOOP4
    {
    
    k = s->bedk(i,j);
    
    xip= p->XP[IP];
	yip= p->YP[JP];
    dist = p->DZN[KP];
		
    density = a->ro(i,j,k);
	
    if(p->S16==1)
    {
    zval = s->bedzh(i,j) + p->S116*p->DZN[KP];
    
        if(p->S33==1)
        {
        uvel=p->ccipol1(a->u,xip,yip,zval);
        vvel=p->ccipol2(a->v,xip,yip,zval);
        wvel=p->ccipol3(a->w,xip,yip,zval);
        }
        
        if(p->S33==2)
        {
        uvel=p->ccipol1_a(a->u,xip,yip,zval);
        vvel=p->ccipol2_a(a->v,xip,yip,zval);
        wvel=p->ccipol3_a(a->w,xip,yip,zval);
        }
        
        
    u_abs = sqrt(uvel*uvel + vvel*vvel  + wvel*wvel);

    u_plus = (1.0/kappa)*log(30.0*(dist/ks));

    tau=density*(u_abs*u_abs)/pow((u_plus>0.0?u_plus:1.0e20),2.0);
    }
	
	if(p->S16==2)
    {
    double E,B,dB,ustar,y_plus,visc,ks_plus,tau0;
    
    visc=p->W2;
    
    zval = s->bedzh(i,j) + p->S116*p->DZN[KP];
    
        if(p->S33==1)
        {
        uvel=p->ccipol1(a->u,xip,yip,zval);
        vvel=p->ccipol2(a->v,xip,yip,zval);
        wvel=p->ccipol3(a->w,xip,yip,zval);
        }
        
        if(p->S33==2)
        {
        uvel=p->ccipol1_a(a->u,xip,yip,zval);
        vvel=p->ccipol2_a(a->v,xip,yip,zval);
        wvel=p->ccipol3_a(a->w,xip,yip,zval);
        }
        
    // predictor    
    u_abs = sqrt(uvel*uvel + vvel*vvel);
    u_plus = (1.0/kappa)*log(30.0*(dist/ks));
    tau0=tau=density*(u_abs*u_abs)/pow((u_plus>0.0?u_plus:1.0e20),2.0);
    ustar=sqrt(tau/density);
    
    
        // corrector
        for(int qn=0; qn<5;++qn)
        {
        y_plus = dist*ustar/visc;
        
        ks_plus = ks*ustar/visc;
        
        B = 5.2;
        
        if(ks_plus<2.25)
        dB = 0.0;
        
        if(ks_plus>=2.25 && ks_plus<90.0)
        dB = B - 8.5 + (1.0/kappa)*log(ks_plus)*sin(0.4258*(log(ks_plus)-0.811));
        
        if(ks_plus>=90.0)
        dB = B - 8.5 + (1.0/kappa)*log(ks_plus);

        E = exp(kappa*(B-dB));

        if(y_plus>11.53)
        u_plus = (1.0/kappa)*log(MAX(E*y_plus,1.0));
        
        if(y_plus<=11.53)
        {
        u_plus = u_abs/(y_plus>1.0e-10?y_plus:1.0e20);
        //cout<<"y_plus: "<<y_plus<<" ustar: "<<ustar<<" u_plus: "<<u_plus<<" u_abs: "<<u_abs<<" tau: "<<density*(u_abs*u_abs)/pow((u_plus>1.0e-10?u_plus:1.0e20),2.0)<<endl;
        }
        
        tau=MIN(density*(u_abs*u_abs)/pow((u_plus>1.0e-4?u_plus:1.0e20),2.0), tau0*3.5);

        ustar=sqrt(tau/density);
        }
    }
    
    if(p->S16==3)
    {
	double v_t,v_d;
		
	xip= p->XP[IP];
	yip= p->YP[JP];

    zval = s->bedzh(i,j) + p->S116*p->DZN[KP];
	
        if(p->S33==1)
        {
        uvel=p->ccipol1(a->u,xip,yip,zval);
        vvel=p->ccipol2(a->v,xip,yip,zval);
        wvel=p->ccipol3(a->w,xip,yip,zval);
        
        v_d=p->ccipol4(a->visc,xip,yip,zval);
        v_t=p->ccipol4(a->eddyv,xip,yip,zval);
        }
        
        if(p->S33==2)
        {
        uvel=p->ccipol1_a(a->u,xip,yip,zval);
        vvel=p->ccipol2_a(a->v,xip,yip,zval);
        wvel=p->ccipol3_a(a->w,xip,yip,zval);
        
        v_d=p->ccipol4_a(a->visc,xip,yip,zval);
        v_t=p->ccipol4_a(a->eddyv,xip,yip,zval);
        }

    u_abs = sqrt(uvel*uvel + vvel*vvel);
    
    tau=density*(v_d + v_t)*(u_abs/dist);
    }
    
	if(p->S16==4)
    {
    zval = s->bedzh(i,j) + 0.5*p->DZN[KP];
    
    if(p->S33==1)
    tau=density*pturb->kinval(i,j,k)*0.3;
    
    if(p->S33==2)
    tau=density*pturb->ccipol_a_kinval(p,pgc,xip,yip,zval)*0.3;
    
    //tau=density*pturb->kinval(i,j,k)*0.3;
    }
    
	
	if(p->S16==5)
    {
    zval = s->bedzh(i,j) + 0.5*p->DZN[KP];
    
        uvel=p->ccipol1(a->u,xip,yip,zval);
        vvel=p->ccipol2(a->v,xip,yip,zval);
        wvel=p->ccipol3(a->w,xip,yip,zval);

        u_abs = sqrt(uvel*uvel + vvel*vvel);
			   
    u_plus = (1.0/kappa)*log(30.0*(dist/ks));


    tauvel=density*(u_abs*u_abs)/pow((u_plus)>(0.0)?(u_plus):(1.0e20),2.0);
	
	taukin=density*pturb->ccipol_a_kinval(p,pgc,xip,yip,zval)*0.3;
	
	tau = sqrt(fabs(tauvel))*sqrt(fabs(taukin));
    }
    
    if(p->S16==6)
    {
        double Cval,wh;
        double bedlevel,waterlevel;
        count=0;
        uvel=vvel=wh=0.0;
        KLOOP
        PCHECK
        {
            if(a->phi(i,j,k)>=0.0)
            {
            uvel+=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
            vvel+=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
            ++count;
            wh+=p->DZN[KP];
            }
        }

        uvel=uvel/double(count);
        vvel=vvel/double(count);

        Cval=18.0*log10((12.0*wh)/ks);

        u_abs = sqrt(uvel*uvel + vvel*vvel);
	
    tau = density*pow(sqrt(9.81)*(u_abs/Cval),2.0);
    }
    
    if(p->S16==7)
    {
        if(p->S33==1)
        {
        uvel=p->ccipol1(a->u,xip,yip,zval);
        vvel=p->ccipol2(a->v,xip,yip,zval);
        }
        
        if(p->S33==2)
        {
        uvel=p->ccipol1_a(a->u,xip,yip,zval);
        vvel=p->ccipol2_a(a->v,xip,yip,zval);
        }

    u_abs = sqrt(uvel*uvel + vvel*vvel);

    u_plus = (1.0/kappa)*log(30.0*(dist/ks));

    tau=density*(u_abs*u_abs)/pow((u_plus>0.0?u_plus:1.0e20),2.0);
    }
    
    s->tau_eff(i,j) = taueff_loc(i,j) = tau;
    s->shearvel_eff(i,j) = sqrt(tau/p->W1);
    s->shields_eff(i,j) = tau/(p->W1*((p->S22-p->W1)/p->W1)*fabs(p->W22)*p->S20);
    }
}

void bedshear::taucritbed(lexer *p, fdm * a, ghostcell *pgc, sediment_fdm *s)
{
	double r,density;
    
    SLICELOOP4
    {
    k = s->bedk(i,j);
    
    density = a->ro(i,j,k);
    
    tauc = (p->S30*fabs(p->W22)*(p->S22-density))*p->S20*s->reduce(i,j);
  
    s->tau_crit(i,j) = taucrit_loc(i,j) = tauc;
    s->shearvel_crit(i,j) = sqrt(tauc/density);
    s->shields_crit(i,j) = tauc/(density*((p->S22-density)/density)*fabs(p->W22)*p->S20);
    }
}

void bedshear::taubed(lexer*, fdm*, ghostcell*, double &tau_eff)
{
    tau_eff = taueff_loc(i,j);
}

void bedshear::taucritbed(lexer*, fdm*, ghostcell*, double &tau_crit)
{
    tau_crit = taucrit_loc(i,j);
}

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
    
    //cout<<"MANNING: "<<manning<<" uabs: "<<uabs<<endl;
    
    if(p->S16==1)
    {
    cf = pow(manning,2.0)/pow(HP,1.0/3.0);
    
    tau = p->W1*9.81*cf*uabs*uabs; 
    }
    
    if(p->S16==2)
    {    
    cf = 2.5*log(12.0*b->hp(i,j)/ p->S21*s->ks(i,j));
    
    tau = p->W1*9.81*uabs*uabs/(fabs(cf*cf)>1.0e-20?(cf*cf):1.0e20); 
    }
    
    if(p->S16==3)
    {    
    cf = 2.5*log(12.0*b->hp(i,j)/ p->S21*s->ks(i,j));
    
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
