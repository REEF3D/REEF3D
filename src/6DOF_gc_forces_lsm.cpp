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

#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_gc::forces_lsm(lexer *p,fdm* a, ghostcell *pgc)
{
	Xe=Ye=Ze=Ke=Me=Ne=0.0;
	
	double ddnorm,dirac;
	double area=0.0;
	double volume;
	double areasum=0.0;
//	double epsi = p->X41*p->DXM;
	double x_norm,y_norm,z_norm;
	double xloc,yloc,zloc;
	double Fx,Fy,Fz,V;
	double Fzl,Fzr;
	double veff, txx,tyy,tzz,txy,txz,tyz;
	double Fxvisc,Fyvisc,Fzvisc; 
	double pval,fbval;
    	
	Fzl=Fzr=0.0;
	Fxvisc=Fyvisc=Fzvisc=0.0;
	
	pgc->start4(p,a->press,40);
	pgc->start4(p,a->press,401);
    pgc->start4(p,a->press,401);

    pgc->gcfb_update_extra_gcb(p,a,a->press);
    
    pgc->dgcpol(p,a->press,p->dgc4,p->dgc4_count,14);
    a->press.ggcpol(p);
    
	
    ALOOP
	{	
		double epsi = p->X41*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);	
		
	dstx=0.0;
	dsty=0.0;
	dstz=0.0;
	
	fbval = a->fb(i,j,k);
	lsSig=fbval/(fabs(fbval)>1.0e-10?sqrt(fbval*fbval):1.0e20);

	
	//dstx=(a->fb(i+1,j,k)-a->fb(i-1,j,k))/(2.0*p->DXM);
	//dsty=(a->fb(i,j+1,k)-a->fb(i,j-1,k))/(2.0*p->DXM);
	//dstz=(a->fb(i,j,k+1)-a->fb(i,j,k-1))/(2.0*p->DXM);
	
	dstx = (a->fb(i+1,j,k)-a->fb(i-1,j,k))/(p->DXP[IM1] + p->DXP[IP]);
	dsty = (a->fb(i,j+1,k)-a->fb(i,j-1,k))/(p->DYP[JM1] + p->DYP[JP]);
	dstz = (a->fb(i,j,k+1)-a->fb(i,j,k-1))/(p->DZP[KM1] + p->DZP[KP]);
			
	dnorm=sqrt(dstx*dstx + dsty*dsty + dstz*dstz);
	
	
	dirac=0.0;
	area=0.0;
	volume=0.0;
	Fx=Fy=Fz=0.0;
	
	
	
		if(fabs(fbval)<epsi)
		{
		dirac = (0.5/epsi)*(1.0 + cos((PI*fbval)/epsi));
		
		if(fbval<epsifb)
		H=0.0;

		if(fbval>epsifb)
		H=1.0;

		if(fabs(fbval)<=epsifb)
		H=0.5*(1.0 + fbval/epsifb + (1.0/PI)*sin((PI*(fbval))/epsifb));


		area =  p->DXN[IP]*p->DYN[JP]*p->DZN[KP] * dirac *dnorm;
		
		volume = p->DXN[IP]*p->DYN[JP]*p->DZN[KP] * H;
				
		x_norm = dstx/(dnorm>1.0e-20?dnorm:1.0e20);
		y_norm = dsty/(dnorm>1.0e-20?dnorm:1.0e20);
		z_norm = dstz/(dnorm>1.0e-20?dnorm:1.0e20);
		
		
		xloc = p->pos_x() - p->originx;
		yloc = p->pos_y() - p->originy;
		zloc = p->pos_z() - p->originz;
				
		pval=a->press(i,j,k);
		Fx = -area*pval*x_norm;
		Fy = -area*pval*y_norm;
		Fz = -area*pval*z_norm;
		

		veff = a->visc(i,j,k) + a->eddyv(i,j,k);
		txx = veff*pudx(p,a);
		tyy = veff*pvdy(p,a);
		tzz = veff*pwdz(p,a);
		txy = veff*pudy(p,a);
		txz = veff*pudz(p,a);
		tyz = veff*pvdz(p,a);
		
		
		Fx +=area* (txx*x_norm + txy*y_norm + txz*z_norm);
		Fy +=area* (txy*x_norm + tyy*y_norm + tyz*z_norm);
		Fz +=area* (txz*x_norm + tyz*y_norm + tzz*z_norm);
		
		
		Fxvisc +=area* (txx*x_norm + txy*y_norm + txz*z_norm);
		Fyvisc +=area* (txy*x_norm + tyy*y_norm + tyz*z_norm);
		Fzvisc +=area* (txz*x_norm + tyz*y_norm + tzz*z_norm);
		
		Xe += Fx;
		Ye += Fy;
		Ze += Fz;
		
		xloc+=p->originx;
		yloc+=p->originy;
		zloc+=p->originz;
		
		Ke += (yloc-yg)*Fz - (zloc-zg)*Fy;
		Me += (zloc-zg)*Fx - (xloc-xg)*Fz;
		Ne += (xloc-xg)*Fy - (yloc-yg)*Fx;
		
		}
		
	areasum += area;
	}
	
	
	areasum = pgc->globalsum(areasum);
	
	Xe = pgc->globalsum(Xe);
	Ye = pgc->globalsum(Ye);
	Ze = pgc->globalsum(Ze);
	Ke = pgc->globalsum(Ke);
	Me = pgc->globalsum(Me);
	Ne = pgc->globalsum(Ne);

	Xe += a->gi*Mfb;
	Ye += a->gj*Mfb;
	Ze += a->gk*Mfb;

	if(p->mpirank==0)
	cout<<"area: "<<areasum<<endl;
	
	if(p->mpirank==0)
	cout<<"Xe: "<<Xe<<" Ye: "<<Ye<<" Ze: "<<Ze<<" Ke: "<<Ke<<" Me: "<<Me<<" Ne: "<<Ne<<endl;
}
