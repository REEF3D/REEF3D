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

#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"density_f.h"

void ghostcell::forcing1(lexer *p, fdm *a, field& f, field &uvel, field &vvel, field &wvel, double alpha)
{
    // loop over gcb1
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
	GC1LOOP
	if(((p->gcb1[n][4]==41 || p->gcb1[n][4]==42 || p->gcb1[n][4]==43)&&p->X38==1)
        || ((p->gcb1[n][4]==21 || p->gcb1[n][4]==22) && p->B19==1))
	{
		i = p->gcb1[n][0];
		j = p->gcb1[n][1];
		k = p->gcb1[n][2];
        
        bc = p->gcb1[n][4];      
		
		//convection
		dx=dy=dz=0.0;	

		pip=1;
		ivel1= 0.5*(uvel(i,j,k)+uvel(i-1,j,k))*(1.0/a->porosity(i-1,j,k));
		ivel2= 0.5*(uvel(i,j,k)+uvel(i+1,j,k))*(1.0/a->porosity(i,j,k));	
		pip=0;
		
		pip=2;
		jvel1= 0.5*(vvel(i,j-1,k)+vvel(i+1,j-1,k))*(1.0/(0.25*(a->porosity(i,j-1,k)+a->porosity(i+1,j-1,k)+a->porosity(i,j,k)+a->porosity(i+1,j,k))));
		jvel2= 0.5*(vvel(i,j,k)+vvel(i+1,j,k))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i+1,j,k)+a->porosity(i,j+1,k)+a->porosity(i+1,j+1,k))));
		pip=0;
		
		pip=3;
		kvel1= 0.5*(wvel(i,j,k-1)+wvel(i+1,j,k-1))*(1.0/(0.25*(a->porosity(i,j,k-1)+a->porosity(i+1,j,k-1)+a->porosity(i,j,k)+a->porosity(i+1,j,k))));
		kvel2= 0.5*(wvel(i,j,k)+wvel(i+1,j,k))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i+1,j,k)+a->porosity(i,j,k+1)+a->porosity(i+1,j,k+1))));
		pip=0;
		
		dx = (ivel2*0.5*(f(i,j,k) + f(i+1,j,k))  -  ivel1*0.5*(f(i-1,j,k) +  f(i,j,k)))/(p->DXM);
		dy = (jvel2*0.5*(f(i,j,k) + f(i,j+1,k))  -  jvel1*0.5*(f(i,j-1,k) +  f(i,j,k)))/(p->DXM);		
		dz = (kvel2*0.5*(f(i,j,k) + f(i,j,k+1))  -  kvel1*0.5*(f(i,j,k-1) +  f(i,j,k)))/(p->DXM);

		L = -dx-dy-dz;
		
		
		// Diffusion
		
		double visc_ddy_p,visc_ddy_m,visc_ddz_p,visc_ddz_m;
		double b_ijk,ev_ijk,visc_ijk;
		double sqd;
		
		sqd = (1.0/(p->DXM*p->DXM));
	


		b_ijk=f(i,j,k);
		ev_ijk=a->eddyv(i,j,k);
		visc_ijk=a->visc(i,j,k);
		visc_ddy_p = (visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k) + a->visc(i,j+1,k)+a->eddyv(i,j+1,k) + a->visc(i+1,j+1,k)+a->eddyv(i+1,j+1,k))*0.25;
		visc_ddy_m = (a->visc(i,j-1,k)+a->eddyv(i,j-1,k)  +a->visc(i+1,j-1,k)+a->eddyv(i+1,j-1,k) + visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k))*0.25;
		visc_ddz_p = (visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k) + a->visc(i,j,k+1)+a->eddyv(i,j,k+1) + a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1))*0.25;
		visc_ddz_m = (a->visc(i,j,k-1)+a->eddyv(i,j,k-1) + a->visc(i+1,j,k-1)+a->eddyv(i+1,j,k-1) + visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k))*0.25;
		
		L += 2.0*sqd*((f(i+1,j,k)-b_ijk)*(a->visc(i+1,j,k)+a->eddyv(i+1,j,k))
					   -(b_ijk-f(i-1,j,k))*(visc_ijk+ev_ijk))

			+   sqd*((f(i,j+1,k)-b_ijk)*visc_ddy_p
				-(b_ijk-f(i,j-1,k))*visc_ddy_m)

			+   sqd*((f(i,j,k+1)-b_ijk)*visc_ddz_p
				-(b_ijk-f(i,j,k-1))*visc_ddz_m)

			+ sqd*((a->v(i+1,j,k)-a->v(i,j,k))*visc_ddy_p - (a->v(i+1,j-1,k)-a->v(i,j-1,k))*visc_ddy_m)
			+ sqd*((a->w(i+1,j,k)-a->w(i,j,k))*visc_ddz_p - (a->w(i+1,j,k-1)-a->w(i,j,k-1))*visc_ddz_m);



		// Gravity
		L += a->gi*PORVAL1;
        
        // Pressure
        L -= PORVAL1*(a->press(i+1,j,k)-a->press(i,j,k))/(p->DXM*pdens->roface(p,a,1,0,0));
		
		
		// Interpolate u_gamma
		double dist = p->gcd1[n];

		
		double x,x0,x1,x2;
		double y0,y1,y2;
		int aa,bb,cc;
		
		aa=bb=cc=0;
		if(p->gcb1[n][3]==1)
		aa=1;
		
		if(p->gcb1[n][3]==2)
		bb=-1;
		
		if(p->gcb1[n][3]==3)
		bb=1;
		
		if(p->gcb1[n][3]==4)
		aa=-1;
		
		if(p->gcb1[n][3]==5)
		cc=1;
		
		if(p->gcb1[n][3]==6)
		cc=-1;
		
		x = 0.0;
		
		x0 = -p->DXM;
		x1 = dist;
		
		y0 = f(i+aa,j+bb,k+cc);
        
        if(bc==41 || bc==42 || bc==43)
		y1 = p->ufbi + (p->pos_z()-p->zg)*p->qfbi - (p->pos_y()-p->yg)*p->rfbi;
        
        if(bc==21)
        y1=0.0;
                
        ugamma =  y0*(x-x1)/((x0-x1))
               +  y1*(x-x0)/((x1-x0));

        		
        //if(dist<p->DXM)
		a->F(i,j,k) += (ugamma-f(i,j,k))/(alpha*p->dt) - L;
	}
}