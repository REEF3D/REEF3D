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
#include <math.h>

void sixdof_gc::forces_stl(lexer* p, fdm *a, ghostcell *pgc)
{
    double x0,x1,x2,y0,y1,y2,z0,z1,z2;
	double xc,yc,zc;
	double at,bt,ct,st;
	double nx,ny,nz,norm;
	double A_triang,A;
	double p_int,rho_int,nu_int,enu_int,u_int,v_int,w_int;
	double du,dv,dw, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
	double dudxf, dudyf, dudzf, dvdxf, dvdyf, dvdzf, dwdxf, dwdyf, dwdzf;
	double dudxb, dudyb, dudzb, dvdxb, dvdyb, dvdzb, dwdxb, dwdyb, dwdzb;
	double xlocvel,ylocvel,zlocvel;
	double Fx,Fy,Fz;


    pgc->start4(p,a->press,40);
	//pgc->start4(p,a->press,401);
    //pgc->start4(p,a->press,401);

    //pgc->gcfb_update_extra_gcb(p,a,a->press);
    
    
	
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
			
            xlocvel = xc + nx*p->DXP[IP];
			ylocvel = yc + ny*p->DYP[JP];
			zlocvel = zc + nz*p->DZP[KP];
			   
			// Add normal stress contributions

			i = p->posc_i(xlocvel);
			j = p->posc_j(ylocvel);
			k = p->posc_k(zlocvel);
           
			p_int = p->ccipol4(a->press,xc,yc,zc);
			
			Fx = -nx*p_int*A_triang;
			Fy = -ny*p_int*A_triang;
            Fz = -nz*p_int*A_triang;
			

            // Add tangential stress contributions
			
			nu_int = p->ccipol4(a->visc,xlocvel,ylocvel,zlocvel);
			enu_int = p->ccipol4(a->eddyv,xlocvel,ylocvel,zlocvel);
			rho_int = p->ccipol4(a->ro,xlocvel,ylocvel,zlocvel);
			
            dudx = (a->u(i+1,j,k) - a->u(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
            dudy = (a->u(i,j+1,k) - a->u(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
            dudz = (a->u(i,j,k+1) - a->u(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);
                                                                           
            dvdx = (a->v(i+1,j,k) - a->v(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
            dvdy = (a->v(i,j+1,k) - a->v(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
            dvdz = (a->v(i,j,k+1) - a->v(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);
                                                                            
            dwdx = (a->w(i+1,j,k) - a->w(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
            dwdy = (a->w(i,j+1,k) - a->w(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
            dwdz = (a->w(i,j,k+1) - a->w(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);
            
            Fx += rho_int*(nu_int + enu_int)*A_triang*(2.0*dudx*nx + (dudy + dvdx)*ny + (dudz + dwdx)*nz);
            Fy += rho_int*(nu_int + enu_int)*A_triang*((dudy + dvdx)*nx + 2.0*dvdy*ny + (dvdz + dwdy)*nz);
            Fz += rho_int*(nu_int + enu_int)*A_triang*((dudz + dwdx)*nx + (dvdz + dwdy)*ny + 2.0*dwdz*nz);
			

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

    //if(p->mpirank==0)
    //cout<<"Fx: "<<Ze<<"  A: "<<A<<endl;

	// Add gravity force
	
	Xe += a->gi*Mfb;
	Ye += a->gj*Mfb;
	Ze += a->gk*Mfb;

	// Print results
	
	if (p->mpirank == 0)
	{
		cout<<"area: "<<A<<endl;
		cout<<"Xe: "<<Xe<<" Ye: "<<Ye<<" Ze: "<<Ze<<" Ke: "<<Ke<<" Me: "<<Me<<" Ne: "<<Ne<<endl;
	}
}
