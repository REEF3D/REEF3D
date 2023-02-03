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

#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include <math.h>

void sixdof_df_object::forces_stl
(
    lexer* p, fdm *a, ghostcell *pgc, double alpha,
    field& uvel, field& vvel, field& wvel
)
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
	double xlocvel,ylocvel,zlocvel,xlocp,ylocp,zlocp;
	double Fx,Fy,Fz,Fp_x,Fp_y,Fp_z,Fv_x,Fv_y,Fv_z;
    double Xe_p,Ye_p,Ze_p,Xe_v,Ye_v,Ze_v;

    A=0.0;
    Xe=Ye=Ze=Ke=Me=Ne=0.0;
    Xe_p=Ye_p=Ze_p=Xe_v=Ye_v=Ze_v=0.0;
    
    // Set new time
    curr_time += alpha*p->dt; 


    for (int n = 0; n < tricount; ++n)
    {     
		// Vertices of triangle
        x0 = tri_x[n][0];
        y0 = tri_y[n][0];
        z0 = tri_z[n][0];
        
        x1 = tri_x[n][1];
        y1 = tri_y[n][1];
        z1 = tri_z[n][1];
        
        x2 = tri_x[n][2];
        y2 = tri_y[n][2];
        z2 = tri_z[n][2];  
		   
		// Center of triangle
		xc = (x0 + x1 + x2)/3.0;
		yc = (y0 + y1 + y2)/3.0;
		zc = (z0 + z1 + z2)/3.0;
        
             at = sqrt(pow(x1-x0,2.0) + pow(y1-y0,2.0) + pow(z1-z0,2.0));
			bt = sqrt(pow(x1-x2,2.0) + pow(y1-y2,2.0) + pow(z1-z2,2.0));
			ct = sqrt(pow(x2-x0,2.0) + pow(y2-y0,2.0) + pow(z2-z0,2.0));
				
			st = 0.5*(at+bt+ct);
				
			A_triang = sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
				

			// Normal vectors (always pointing outwards)      
            nx = (y1 - y0)*(z2 - z0) - (y2 - y0)*(z1 - z0);
            ny = (x2 - x0)*(z1 - z0) - (x1 - x0)*(z2 - z0); 
            nz = (x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0);

			norm = sqrt(nx*nx + ny*ny + nz*nz);
			
			nx /= norm > 1.0e-20 ? norm : 1.0e20;
			ny /= norm > 1.0e-20 ? norm : 1.0e20;
			nz /= norm > 1.0e-20 ? norm : 1.0e20;
        
 
		if 
		(
			xc >= p->originx && xc < p->endx &&
			yc >= p->originy && yc < p->endy &&
			zc >= p->originz && zc < p->endz
		)
		{
            // Position of triangle
            i = p->posc_i(xc);
            j = p->posc_j(yc);
            k = p->posc_k(zc);
			
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
            
			// Add normal stress contributions
             xlocp = xc + p->X42*nx*p->DXP[IP];
			ylocp = yc + p->X42*ny*p->DYP[JP];
			zlocp = zc + p->X42*nz*p->DZP[KP];

			p_int = p->ccipol4_a(a->press,xlocp,ylocp,zlocp);
            
            //p_int =1.0;
            
			
			Fp_x = -nx*p_int*A_triang;
			Fp_y = -ny*p_int*A_triang;
             Fp_z = -nz*p_int*A_triang;
             
		   
        // Add viscous stress contributions
        
            double ustar, uplus, dist, value;
            double vstar,vplus,wstar,wplus;
            double uval, vval, wval;
            double kappa = 0.4;
            double xip, yip, zip, zval;
            double u_abs, tau, density;
            double ks;


            if(p->X39==0)
            {
            xlocvel = xc + p->X43*nx*p->DXP[IP];
            ylocvel = yc + p->X43*ny*p->DYP[JP];
            zlocvel = zc + p->X43*nz*p->DZP[KP];
            
             nu_int = p->ccipol4_a(a->visc,xlocvel,ylocvel,zlocvel);
			enu_int = 0.0; //p->ccipol4_a(a->eddyv,xlocvel,ylocvel,zlocvel);
			rho_int = p->ccipol4_a(a->ro,xlocvel,ylocvel,zlocvel);
	        
             i = p->posc_i(xlocvel);
			j = p->posc_j(ylocvel);
			k = p->posc_k(zlocvel);
            
            dudx = (uvel(i+1,j,k) - uvel(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
            dudy = (uvel(i,j+1,k) - uvel(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
            dudz = (uvel(i,j,k+1) - uvel(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);
                                                                           
            dvdx = (vvel(i+1,j,k) - vvel(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
            dvdy = (vvel(i,j+1,k) - vvel(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
            dvdz = (vvel(i,j,k+1) - vvel(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);
                                                                            
            dwdx = (wvel(i+1,j,k) - wvel(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
            dwdy = (wvel(i,j+1,k) - wvel(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
            dwdz = (wvel(i,j,k+1) - wvel(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);

            Fv_x = rho_int*(nu_int + enu_int)*A_triang*(2.0*dudx*nx + (dudy + dvdx)*ny + (dudz + dwdx)*nz);
            Fv_y = rho_int*(nu_int + enu_int)*A_triang*((dudy + dvdx)*nx + 2.0*dvdy*ny + (dvdz + dwdy)*nz);
            Fv_z = rho_int*(nu_int + enu_int)*A_triang*((dudz + dwdx)*nx + (dvdz + dwdy)*ny + 2.0*dwdz*nz);
            }
            
            
            if(p->X39==1)
            {
            xlocvel = xc + p->X43*nx*p->DXP[IP];
            ylocvel = yc + p->X43*ny*p->DYP[JP];
            zlocvel = zc + p->X43*nz*p->DZP[KP];
            
             nu_int  = p->ccipol4_a(a->visc,xlocvel,ylocvel,zlocvel);
			enu_int = p->ccipol4_a(a->eddyv,xlocvel,ylocvel,zlocvel);
			rho_int = p->ccipol4_a(a->ro,xlocvel,ylocvel,zlocvel);
	        
             i = p->posc_i(xlocvel);
			j = p->posc_j(ylocvel);
			k = p->posc_k(zlocvel);
            
            dudx = (uvel(i+1,j,k) - uvel(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
            dudy = (uvel(i,j+1,k) - uvel(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
            dudz = (uvel(i,j,k+1) - uvel(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);
                                                                           
            dvdx = (vvel(i+1,j,k) - vvel(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
            dvdy = (vvel(i,j+1,k) - vvel(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
            dvdz = (vvel(i,j,k+1) - vvel(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);
                                                                            
            dwdx = (wvel(i+1,j,k) - wvel(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
            dwdy = (wvel(i,j+1,k) - wvel(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
            dwdz = (wvel(i,j,k+1) - wvel(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);

            Fv_x = rho_int*(nu_int + enu_int)*A_triang*(2.0*dudx*nx + (dudy + dvdx)*ny + (dudz + dwdx)*nz);
            Fv_y = rho_int*(nu_int + enu_int)*A_triang*((dudy + dvdx)*nx + 2.0*dvdy*ny + (dvdz + dwdy)*nz);
            Fv_z = rho_int*(nu_int + enu_int)*A_triang*((dudz + dwdx)*nx + (dvdz + dwdy)*ny + 2.0*dwdz*nz);
            }
            
            
            if(p->X39==2)
            {
            xlocvel = xc + p->X43*nx*p->DXP[IP];
            ylocvel = yc + p->X43*ny*p->DYP[JP];
            zlocvel = zc + p->X43*nz*p->DZP[KP];
            
            nu_int  = p->ccipol4_a(a->visc,xlocvel,ylocvel,zlocvel);
            enu_int = p->ccipol4_a(a->eddyv,xlocvel,ylocvel,zlocvel);
            rho_int = p->ccipol4_a(a->ro,xlocvel,ylocvel,zlocvel);
	        
            i = p->posc_i(xlocvel);
            j = p->posc_j(ylocvel);
            k = p->posc_k(zlocvel);
            
            uval=p->ccipol1(uvel,xlocvel,ylocvel,zlocvel);
            vval=p->ccipol2(vvel,xlocvel,ylocvel,zlocvel);
            wval=p->ccipol3(wvel,xlocvel,ylocvel,zlocvel);
            
            dist = fabs(p->ccipol4(a->fb,xlocvel,ylocvel,zlocvel));
            
            ks = 0.01*dist;
            
            /*
            if(p->B5==2)
            {
            ustar=sqrt(fabs((0.5*a->visc(i,j,k)+0.5*a->visc(i,j-1,k))*a->u(i,j,k)/dist));
            value=((9.0*dist*ustar)/(0.5*a->visc(i,j,k)+0.5*a->visc(i,j-1,k)));
            
            if(value>1.0)
            uplus=(1.0/kappa)*log(value);
            
            if(value<=1.0)
            uplus=1.0e20;
            }*/
        
            // x-dir
            //ustar = sqrt(fabs(nu_int*uval)/(dist>0.0?dist:1.0e20));
            
            ustar = uval/((1.0/kappa)*log(30.0*(dist/ks)));
            
            //for(
            
            
            
            
            value = ((9.0*dist*ustar)/(nu_int>0.0?nu_int:1.0e20));
            
            if(value>1.0)
            uplus=(1.0/kappa)*log(value);
            
            if(value<=1.0)
            uplus=1.0e20;
  
            Fv_x =  (ny + nz)*rho_int*(uval*uval)/pow((uplus>0.0?uplus:1.0e20),2.0);
            
            cout<<"Fv_x: "<<Fv_x<<" ustar: "<<ustar<<" nx: "<<nx<<" uplus: "<<uplus<<" dist: "<<dist<<" value: "<<value<<" log(value): "<<log(value)<<endl;
            
            // y-dir
            vstar = sqrt(fabs(nu_int*vval)/(dist>0.0?dist:1.0e20));
            value=((9.0*dist*vstar)/(nu_int>0.0?nu_int:1.0e20));
            
            if(value>1.0)
            vplus=(1.0/kappa)*log(value);
            
            if(value<=1.0)
            vplus=1.0e20;
  
            Fv_y = 0.0;// (1.0-fabs(ny))*rho_int*(vval*vval)/pow((vplus>0.0?vplus:1.0e20),2.0);
            
            // z-dir
            //wstar = sqrt(fabs(nu_int*wval)/(dist>0.0?dist:1.0e20));
            
            wstar = wval/((1.0/kappa)*log(30.0*(dist/ks)));
            
            value = ((9.0*dist*wstar)/(nu_int>0.0?nu_int:1.0e20));
            
            if(value>1.0)
            wplus=(1.0/kappa)*log(value);
            
            if(value<=1.0)
            wplus=1.0e20;
  
            Fv_z = (nx + ny)*rho_int*(wval*wval)/pow((wplus>0.0?wplus:1.0e20),2.0);
            }
            
            if(p->X39==3)
            {
            double v_t,v_d;
                
            xip= p->XP[IP];
            yip= p->YP[JP];

            //zval = a->bed(i,j) + p->S116*p->DZN[KP];
            
                if(p->S33==1)
                {
                uval=p->ccipol1(a->u,xip,yip,zval);
                vval=p->ccipol2(a->v,xip,yip,zval);
                wval=p->ccipol3(a->w,xip,yip,zval);
                
                v_d=p->ccipol4(a->visc,xip,yip,zval);
                v_t=p->ccipol4(a->eddyv,xip,yip,zval);
                }
                
                if(p->S33==2)
                {
                uval=p->ccipol1_a(a->u,xip,yip,zval);
                vval=p->ccipol2_a(a->v,xip,yip,zval);
                wval=p->ccipol3_a(a->w,xip,yip,zval);
                
                v_d=p->ccipol4_a(a->visc,xip,yip,zval);
                v_t=p->ccipol4_a(a->eddyv,xip,yip,zval);
                }

            u_abs = sqrt(uval*uval + vval*vval);
            
            tau=density*(v_d + v_t)*(u_abs/dist);
            }
            
            
            if(p->X39==4)
            {
            //zval = s->bedzh(i,j) + 0.5*p->DZN[KP];
            
           /* if(p->S33==1)
            tau=density*pturb->ccipol_a_kinval(p,pgc,xip,yip,zval)*0.3;
            
            if(p->S33==2)
            tau=density*pturb->ccipol_kinval(p,pgc,xip,yip,zval)*0.3;*/
            }


            // Total forces
            
            Fx = Fp_x + Fv_x;
            Fy = Fp_y + Fv_y;
            Fz = Fp_z + Fv_z;
            

			// Add forces to global forces
			
			Xe += Fx;
			Ye += Fy;
			Ze += Fz;

			Ke += (yc - c_(1))*Fz - (zc - c_(2))*Fy;
			Me += (zc - c_(2))*Fx - (xc - c_(0))*Fz;
			Ne += (xc - c_(0))*Fy - (yc - c_(1))*Fx;
            
            Xe_p += Fp_x;
            Ye_p += Fp_y;
            Ze_p += Fp_z;
            
            Xe_v += Fv_x;
            Ye_v += Fv_y;
            Ze_v += Fv_z;
							
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

	Xe_p = pgc->globalsum(Xe_p);
	Ye_p = pgc->globalsum(Ye_p);
	Ze_p = pgc->globalsum(Ze_p);
	Xe_v = pgc->globalsum(Xe_v);
	Ye_v = pgc->globalsum(Ye_v);
	Ze_v = pgc->globalsum(Ze_v);

	// Add gravity force
	
    if(p->mpirank==0)
    cout<<"Viscous Force Fx_v: "<<Xe_v<<" Fz_v: "<<Ze_v<<endl<<endl;
    
	Xe += a->gi*Mass_fb;
	Ye += a->gj*Mass_fb;
	Ze += a->gk*Mass_fb;


    // Print results
	
    if (p->mpirank == 0) 
    {
        ofstream print;
        char str[1000];
       
        if(p->P14==0)
        sprintf(str,"REEF3D_6DOF_forces_%i.dat",n6DOF);
        if(p->P14==1)
        sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_forces_%i.dat",n6DOF);

        print.open(str, std::ofstream::out | std::ofstream::app);
        print<<curr_time<<" \t "<<Xe<<" \t "<<Ye<<" \t "<<Ze<<" \t "<<Ke
        <<" \t "<<Me<<" \t "<<Ne<<" \t "<<Xe_p<<" \t "<<Ye_p<<" \t "<<Ze_p<<" \t "<<Xe_v<<" \t "<<Ye_v<<" \t "<<Ze_v<<endl;   
        print.close();
    }

	if (p->mpirank == 0)
	cout<<"Xe: "<<Xe<<" Ye: "<<Ye<<" Ze: "<<Ze<<" Ke: "<<Ke<<" Me: "<<Me<<" Ne: "<<Ne<<endl;
}
