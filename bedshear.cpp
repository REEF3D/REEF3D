/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"bedshear.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"reduction_void.h"
#include"reduction_parker.h"
#include"reduction_deyemp.h"
#include"reduction_deyana.h"
#include"reduction_FD.h"
#include"sliceint.h"

bedshear::bedshear(lexer *p, turbulence *ppturb) : norm_vec(p), ks(p->S20*p->S21), kappa(0.4)
{
    tau=0.0;
    tauc=0.0;
    pturb=ppturb;

    if(p->S80==0)
    preduce=new reduction_void(p);

    if(p->S80==1)
    preduce=new reduction_parker(p);

    if(p->S80==2)
    preduce=new reduction_deyemp(p);

    if(p->S80==3)
    preduce=new reduction_deyana(p);
	
	if(p->S80==4)
    preduce=new reduction_FD(p);
}

bedshear::~bedshear()
{
}

void bedshear::taubed(lexer *p, fdm * a, ghostcell *pgc, double &tau_eff, double &shearvel_eff, double &shields_eff)
{
	int count;
	double zval,fac,topoval,taukin,tauvel,density;
    
    k = a->bedk(i,j)+1;
    
    xip= p->XP[IP];
	yip= p->YP[JP];
    zval = a->bedzh(i,j) + p->S116*p->DZN[k];
    dist = p->S117*p->DZN[k];
		
    density = p->W1;
	
    if(p->S16==1)
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
	
	if(p->S16==2)
    {    
        double wh;
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
        
        u_abs = sqrt(uvel*uvel + vvel*vvel);
	
	fac = 0.125/pow(log(12.0*(wh/ks)),2.0);

	tau = density*fac*pow(u_abs,2.0);
    }

	
	if(p->S16==3)
    {
	double v_t,v_d;
		
	xip= p->XP[IP];
	yip= p->YP[JP];
	
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
        
	v_d=p->ccipol4_a(a->visc,xip,yip,zval);
	v_t=p->ccipol4_a(a->eddyv,xip,yip,zval);

    u_abs = sqrt(uvel*uvel + vvel*vvel);
    

    tau=density*(v_d + v_t)*(u_abs/dist);
    }
	
	if(p->S16==4)
    tau=density*pturb->ccipol_kinval(p,pgc,xip,yip,zval)*0.3;
    
    if(p->S16==5)
    tau=density*pturb->ccipol_kinval(p,pgc,xip,yip,zval)*0.3;
	
	if(p->S16==6)
    {
		pip=1;
        uvel=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
        pip=0;

        pip=2;
        vvel=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
        pip=0;

        pip=3;
        wvel=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
        pip=0;

        u_abs = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);
			   
    u_plus = (1.0/kappa)*log(30.0*(dist/ks));


    tauvel=density*(u_abs*u_abs)/pow((u_plus)>(0.0)?(u_plus):(1.0e20),2.0);
	
	taukin=density*pturb->kinval(i,j,k)*0.3;
	
	tau = sqrt(fabs(tauvel))*sqrt(fabs(taukin));
    }
    
    if(p->S16==7)
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
    
    if(p->S16==8)
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
    
    tau_eff = tau;
    shearvel_eff = sqrt(tau/p->W1);
    shields_eff = tau/(p->W1*((p->S22-p->W1)/p->W1)*fabs(p->W22)*p->S20);
}

void bedshear::taucritbed(lexer *p, fdm * a, ghostcell *pgc, double &tau_crit, double &shearvel_crit, double &shields_crit)
{
	double r;

	k = a->bedk(i,j);
	
	r = preduce->start(p,a,pgc);
    
    tauc = (p->S30*fabs(p->W22)*(p->S22-p->W1))*p->S20*r;
  
    tau_crit = tauc;
    shearvel_crit = sqrt(tauc/p->W1);
    shields_crit = tauc/(p->W1*((p->S22-p->W1)/p->W1)*fabs(p->W22)*p->S20);

}

double bedshear::shear_reduction(lexer *p, fdm *a, ghostcell *pgc)
{
    double r=1.0;

    r=preduce->start(p,a,pgc);

    return r;
}

