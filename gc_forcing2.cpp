/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"density_f.h"

void ghostcell::forcing2(lexer *p, fdm *a, field& f, field &uvel, field &vvel, field &wvel, double alpha)
{
    // loop over gcb2
    // if sixdof gcb
    // calculate simple convection C
    // calculate simple diffusion D
    // interpolate u_gamma, from solid and fluid velcoity
    // fi = (u_gamma - uvel)/dt - C - D - g
    // add fi to F
	
	double dx,dy,dz,L;
	double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;
	double solidvel, ugamma;
    int bc;
	
	if(p->X38==1 || p->B19==1)
	GC2LOOP
	if(((p->gcb2[n][4]==41 || p->gcb2[n][4]==42 || p->gcb2[n][4]==43)&&p->X38==1)
        || ((p->gcb2[n][4]==21 || p->gcb2[n][4]==22) && p->B19==1))
	{
		i = p->gcb2[n][0];
		j = p->gcb2[n][1];
		k = p->gcb2[n][2];
        
        bc = p->gcb2[n][4];
		
		//convection
		dx=dy=dz=0.0;	

		pip=1;
		ivel1= 0.5*(uvel(i-1,j,k)+uvel(i-1,j+1,k))*(1.0/(0.25*(a->porosity(i-1,j,k)+a->porosity(i-1,j+1,k)+a->porosity(i,j,k)+a->porosity(i,j+1,k))));
		ivel2= 0.5*(uvel(i,j,k)+uvel(i,j+1,k))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i,j+1,k)+a->porosity(i+1,j,k)+a->porosity(i+1,j+1,k))));
		pip=0;
	
		pip=2;
		jvel1= 0.5*(vvel(i,j,k)+vvel(i,j-1,k))*(1.0/a->porosity(i,j-1,k));
		jvel2= 0.5*(vvel(i,j,k)+vvel(i,j+1,k))*(1.0/a->porosity(i,j,k));
		pip=0;

		pip=3;
		kvel1= 0.5*(wvel(i,j,k-1)+wvel(i,j+1,k-1))*(1.0/(0.25*(a->porosity(i,j,k-1)+a->porosity(i,j+1,k-1)+a->porosity(i,j,k)+a->porosity(i,j+1,k))));
		kvel2= 0.5*(wvel(i,j,k)+wvel(i,j+1,k))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i,j+1,k)+a->porosity(i,j,k+1)+a->porosity(i,j+1,k+1))));
		pip=0;
			
		dx = (ivel2*0.5*(f(i,j,k) + f(i+1,j,k))  -  ivel1*0.5*(f(i-1,j,k) +  f(i,j,k)))/(p->dx);
		dy = (jvel2*0.5*(f(i,j,k) + f(i,j+1,k))  -  jvel1*0.5*(f(i,j-1,k) +  f(i,j,k)))/(p->dx);		
		dz = (kvel2*0.5*(f(i,j,k) + f(i,j,k+1))  -  kvel1*0.5*(f(i,j,k-1) +  f(i,j,k)))/(p->dx);

		L = -dx-dy-dz;
		
		
		// Diffusion
		
		double visc_ddx_p,visc_ddx_m,visc_ddz_p,visc_ddz_m;
		double b_ijk,ev_ijk,visc_ijk;
		double vfm, vft,sqd;
		
		sqd = (1.0/(p->dx*p->dx));
	
		vfm=vft=0.0;
		
		if(p->D22==1)
		vfm=1.0;
		
		if(p->D23==1)
		vft=1.0;

		b_ijk=f(i,j,k);
		ev_ijk=a->eddyv(i,j,k);
		visc_ijk=a->visc(i,j,k);
		visc_ddx_p = (vfm*visc_ijk+ev_ijk + vfm*a->visc(i,j+1,k)+a->eddyv(i,j+1,k) + vfm*a->visc(i+1,j,k)+a->eddyv(i+1,j,k) + vfm*a->visc(i+1,j+1,k)+a->eddyv(i+1,j+1,k))*0.25;
		visc_ddx_m = (vfm*a->visc(i-1,j,k)+a->eddyv(i-1,j,k) + vfm*a->visc(i-1,j+1,k)+a->eddyv(i-1,j+1,k) + vfm*visc_ijk+ev_ijk + vfm*a->visc(i,j+1,k)+a->eddyv(i,j+1,k))*0.25;
		visc_ddz_p = (vfm*visc_ijk+ev_ijk + vfm*a->visc(i,j+1,k)+a->eddyv(i,j+1,k) + vfm*a->visc(i,j,k+1)+a->eddyv(i,j,k+1) + vfm*a->visc(i,j+1,k+1)+a->eddyv(i,j+1,k+1))*0.25;
		visc_ddz_m = (vfm*a->visc(i,j,k-1)+a->eddyv(i,j,k-1) + vfm*a->visc(i,j+1,k-1)+a->eddyv(i,j+1,k-1) + vfm*visc_ijk+ev_ijk + vfm*a->visc(i,j+1,k)+a->eddyv(i,j+1,k))*0.25;
		
		L	 +=	  sqd*((f(i+1,j,k)-b_ijk)*visc_ddx_p
					- (b_ijk-f(i-1,j,k))*visc_ddx_m)

			+ 2.0*sqd*((f(i,j+1,k)-b_ijk)*(vfm*a->visc(i,j+1,k)+a->eddyv(i,j+1,k))
				 -(b_ijk-f(i,j-1,k))*(vfm*visc_ijk+ev_ijk))

			+ sqd*((f(i,j,k+1)-b_ijk)*visc_ddz_p
			  -(b_ijk-f(i,j,k-1))*visc_ddz_m)

			+ sqd*((a->u(i,j+1,k)-a->u(i,j,k))*visc_ddx_p - (a->u(i-1,j+1,k)-a->u(i-1,j,k))*visc_ddx_m)
			+ sqd*((a->w(i,j+1,k)-a->w(i,j,k))*visc_ddz_p - (a->w(i,j+1,k-1)-a->w(i,j,k-1))*visc_ddz_m);



		// Gravity
		L += a->gj*PORVAL2;
        
        // Pressure
        L -= PORVAL2*(a->press(i,j+1,k)-a->press(i,j,k))/(p->dx*pdens->roface(p,a,0,1,0));
		
		
		// Interpolate u_gamma
		double dist = p->gcd2[n];

		
		double x,x0,x1,x2;
		double y0,y1,y2;
		int aa,bb,cc;
		
		aa=bb=cc=0;
		if(p->gcb2[n][3]==1)
		aa=1;
		
		if(p->gcb2[n][3]==2)
		bb=-1;
		
		if(p->gcb2[n][3]==3)
		bb=1;
		
		if(p->gcb2[n][3]==4)
		aa=-1;
		
		if(p->gcb2[n][3]==5)
		cc=1;
		
		if(p->gcb2[n][3]==6)
		cc=-1;
		
		x = 0.0;
		/*
		x0 = -2.0*p->dx;
		x1 = - p->dx;
		x2 = dist;
		
		y0 = f(i+2*aa,j+2*bb,k+2*cc);
		y1 = f(i+aa,j+bb,k+cc);

        if(bc==41 || bc==42 || bc==43)
		y2 = p->vfbi + (p->pos_x()-p->xg)*p->rfbi - (p->pos_z()-p->zg)*p->pfbi;
        
        if(bc==21)
        y2=0.0;*/

		/*ugamma =  y0*(x-x1)*(x-x2)/((x0-x1)*(x0-x2))
				+ y1*(x-x0)*(x-x2)/((x1-x0)*(x1-x2))
				+ y2*(x-x0)*(x-x1)/((x2-x0)*(x2-x1));*/

		x0 = - p->dx;
		x1 = dist;
		

		y0 = f(i+aa,j+bb,k+cc);

        if(bc==41 || bc==42 || bc==43)
		y1 = p->vfbi + (p->pos_x()-p->xg)*p->rfbi - (p->pos_z()-p->zg)*p->pfbi;
        
        if(bc==21)
        y1=0.0;
        
        ugamma =  y0*(x-x1)/((x0-x1))
				+ y1*(x-x0)/((x1-x0));
		
        //if(dist<p->dx)
		a->G(i,j,k) += (ugamma-f(i,j,k))/(alpha*p->dt) - L;
	}
}