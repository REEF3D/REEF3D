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
#include"sliceint.h"

void bedshear::taubedx(lexer *p, fdm * a, ghostcell *pgc, double &tau_eff, double &shearvel_eff, double &shields_eff)
{
	int count;
	double zval,fac,topoval,taukin,tauvel,density;
    
    k=a->bedk(i,j)+1;
    
		zval = a->bedzh(i,j) + p->S116*p->DZN[k];
		
		dist = p->S117*p->DZN[k];
		
	
    
    density = p->W1;
	
    if(p->S16==1)
    {
	xip= p->XN[IP1];
	yip= p->YP[JP];
	zip= p->ZP[KP];
    
	uvel=p->ccipol1(a->u,xip,yip,zval);
    vvel=p->ccipol2(a->v,xip,yip,zval);
    
    u_abs = sqrt(uvel*uvel + vvel*vvel);

    u_plus = (1.0/kappa)*log(30.0*(dist/ks));

    tau=density*(u_abs*u_abs)/pow((u_plus>0.0?u_plus:1.0e20),2.0);
    }
	
	
	
	if(p->S16==3)
    {
	double v_t,v_d;
		
	xip= p->XN[IP1];
	yip= p->YP[JP];
	zip= p->ZP[KP];
	
	uvel=p->ccipol1_a(a->u,xip,yip,zval);
    
	v_d=p->ccipol4_a(a->visc,xip,yip,zval);
	v_t=p->ccipol4_a(a->eddyv,xip,yip,zval);

    
    tau=density*(v_d + v_t)*(uvel/dist);
    }
    
    
    tau_eff = tau;
    shearvel_eff = sqrt(tau/p->W1);
    shields_eff = tau/(p->W1*((p->S22-p->W1)/p->W1)*fabs(p->W22)*p->S20);
}

void bedshear::taubedy(lexer *p, fdm * a, ghostcell *pgc, double &tau_eff, double &shearvel_eff, double &shields_eff)
{
	int count;
	double zval,fac,topoval,taukin,tauvel,density;

		zval = a->bedzh(i,j) + p->S116*p->dx-p->originz;
		
		dist = p->S117*p->dx;
		
	k=a->bedk(i,j)+1;
    
    density = p->W1;
	
    if(p->S16==1)
    {
	xip= p->XP[IP];
	yip= p->YN[JP1];
	zip= p->ZP[KP];
    
    uvel=p->ccipol1(a->u,xip,yip,zval);
	vvel=p->ccipol2(a->v,xip,yip,zval);
	
     u_abs = sqrt(uvel*uvel + vvel*vvel);
 
    u_plus = (1.0/kappa)*log(30.0*(dist/ks));

    tau=density*(u_abs*u_abs)/pow((u_plus>0.0?u_plus:1.0e20),2.0);
    }
	
	
	if(p->S16==3)
    {
	double v_t,v_d;
		
	xip= p->XP[IP];
	yip= p->YN[JP1];
	zip= p->ZP[KP];
	
	vvel=p->ccipol2_a(a->v,xip,yip,zval);

	v_d=p->ccipol4_a(a->visc,xip,yip,zval);
	v_t=p->ccipol4_a(a->eddyv,xip,yip,zval);

    tau=density*(v_d + v_t)*(vvel/dist);
    }
    
    
    tau_eff = tau;
    shearvel_eff = sqrt(tau/p->W1);
    shields_eff = tau/(p->W1*((p->S22-p->W1)/p->W1)*fabs(p->W22)*p->S20);
}