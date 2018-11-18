/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"6DOF_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include <math.h>

void sixdof_f::forces_stl(lexer* p, fdm *a, ghostcell *pgc)
{
	double x0,x1,x2,y0,y1,y2,z0,z1,z2;
	double xc,yc,zc;
	double at,bt,ct,st;
	double nx,ny,nz,norm;
	double A_triang,A;
	double p_int,rho_int,nu_int,u_int,v_int,w_int;
	double du,dv,dw;
	double xlocvel,ylocvel,zlocvel;
	double Fx,Fy,Fz;
	

    pgc->start4(p,a->press,40);
	pgc->start4(p,a->press,401);
    pgc->start4(p,a->press,401);

    pgc->gcfb_update_extra_gcb(p,a,a->press);
    
    pgc->dgcpol(p,a->press,p->dgc4,p->dgc4_count,14);
    a->press.ggcpol(p);
    
	
    A=0.0;
    Xe=Ye=Ze=Ke=Me=Ne=0.0;

    for (int n = 0; n < tricount; ++n)
    {     
		// Vertices of triangle
	
		x0 = tri_x[n][0];
		x1 = tri_x[n][1];
		x2 = tri_x[n][2];
		
		y0 = tri_y[n][0];
		y1 = tri_y[n][1];
		y2 = tri_y[n][2];
		
		z0 = tri_z[n][0];
		z1 = tri_z[n][1];
		z2 = tri_z[n][2];  
           
		   
		// Center of triangle
		
		xc = (x0 + x1 + x2)/3.0;
		yc = (y0 + y1 + y2)/3.0;
		zc = (z0 + z1 + z2)/3.0;
            
		if 
		(
			xc >= p->originx && xc < p->endx &&
			yc >= p->originy && yc < p->endy &&
			zc >= p->originz && zc < p->endz
		)
		{
			// Area of triangle using Heron's formula

			at = sqrt(pow(x1-x0,2.0) + pow(y1-y0,2.0) + pow(z1-z0,2.0));
			bt = sqrt(pow(x1-x2,2.0) + pow(y1-y2,2.0) + pow(z1-z2,2.0));
			ct = sqrt(pow(x2-x0,2.0) + pow(y2-y0,2.0) + pow(z2-z0,2.0));
				
			st = 0.5*(at+bt+ct);
				
			A_triang = sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
				

			// Normal vectors (always pointing outwards)      
				
			nx = (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0);
			ny = (x2 - x0) * (z1 - z0) - (x1 - x0) * (z2 - z0); 
			nz = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);

			norm = sqrt(nx*nx + ny*ny + nz*nz);
			
			nx /= norm > 1.0e-20 ? norm : 1.0e20;
			ny /= norm > 1.0e-20 ? norm : 1.0e20;
			nz /= norm > 1.0e-20 ? norm : 1.0e20;
			
			   
			// Interpolate values at centre of triangle

			p_int = p->ccipol4_a(a->press,xc,yc,zc);
			nu_int = p->ccipol4_a(a->visc,xc,yc,zc);
			rho_int = p->ccipol4_a(a->ro,xc,yc,zc);
			
			
			// Interpolate velocities in fluid near centre of triangle
			
			i = p->posc_i(xc);
			j = p->posc_j(yc);
			k = p->posc_k(zc);
	 
			xlocvel = xc + nx*p->DXP[IP];
			ylocvel = yc + ny*p->DYP[JP];
			zlocvel = zc + nz*p->DZP[KP];
			
			u_int = p->ccipol1_a(a->u,xlocvel,ylocvel,zlocvel);
			v_int = p->ccipol2_a(a->v,xlocvel,ylocvel,zlocvel);
			w_int = p->ccipol3_a(a->w,xlocvel,ylocvel,zlocvel);
				
			du = u_int/p->DXP[IP];
			dv = v_int/p->DYP[JP];
			dw = w_int/p->DZP[KP];
				

			// Calculate forces on triangle
				
			Fx = -nx*p_int*A_triang + rho_int*nu_int*A_triang*(du*ny+du*nz);
			Fy = -ny*p_int*A_triang + rho_int*nu_int*A_triang*(dv*nx+dv*nz);           
			Fz = -nz*p_int*A_triang + rho_int*nu_int*A_triang*(dw*nx+dw*ny);  
			

			// Add forces to global forces
			
			Xe += Fx;
			Ye += Fy;
			Ze += Fz;

			Ke += (yc - yg)*Fz - (zc - zg)*Fy;
			Me += (zc - zg)*Fx - (xc - xg)*Fz;
			Ne += (xc - xg)*Fy - (yc - yg)*Fx;
							
			A += A_triang;
		}
	}		
 
	// Communication with other processors
	
    A = pgc->globalsum(A);
	
	Xe = pgc->globalsum(Xe);
	Ye = pgc->globalsum(Ye);
	Ze = pgc->globalsum(Ze);
	Ke = pgc->globalsum(Ke);
	Me = pgc->globalsum(Me);
	Ne = pgc->globalsum(Ne);


	// Add gravity forces
	
	Xe += a->gi*Mfb;
	Ye += a->gj*Mfb;
	Ze += a->gk*Mfb;

/*
if (p->mpirank == 0) cout<<"Changed force calc"<<endl;
	double r = 0.0762;
	double d = 0.01;

	double hprime = p->zg - 0.8;

	double h = r - fabs(hprime);

	double s = 2*sqrt(2*r*h-h*h);

	double b = 2*r*asin(s/(2*r));

	if (hprime > 0.0)
	{
		double A = r*b/2 - s*(r-h)/2;
	}
	else
	{
		double A = PI*r*r - (r*b/2 - s*(r-h)/2);
	}

	double Fg = PI*r*r*d*9.81*500;
	double Fb = 1000*9.81*A*d;

	Xe = 0.0;
	Ze = Fb - Fg;
	Ye = 0.0;
	Ke = 0.0;
	Me = 0.0;
	Ne = 0.0;
*/


	// Print results
	
	if (p->mpirank == 0)
	{
		cout<<"area: "<<A<<endl;
		cout<<"Xe: "<<Xe<<" Ye: "<<Ye<<" Ze: "<<Ze<<" Ke: "<<Ke<<" Me: "<<Me<<" Ne: "<<Ne<<endl;
	}


	// Transform to body coordinate system
	
	transform_vec_ES(Xe,Ye,Ze,Xs,Ys,Zs);
	transform_vec_ES(Ke,Me,Ne,Ks,Ms,Ns);	
}