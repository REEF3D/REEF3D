/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include <math.h>

void sixdof_obj::forces_stl(lexer* p, fdm *a, ghostcell *pgc,field& uvel, field& vvel, field& wvel, int iter)
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
    curr_time += alpha[iter]*p->dt; 


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
        
             /*at = sqrt(pow(x1-x0,2.0) + pow(y1-y0,2.0) + pow(z1-z0,2.0));
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
			nz /= norm > 1.0e-20 ? norm : 1.0e20;*/
        
 
		if (xc >= p->originx && xc < p->endx &&
			yc >= p->originy && yc < p->endy &&
			zc >= p->originz && zc < p->endz)
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

            p_int = p->ccipol4_a(a->press,xlocp,ylocp,zlocp) - p->pressgage;
            
            Fp_x = -nx*p_int*A_triang;
            Fp_y = -ny*p_int*A_triang;
            Fp_z = -nz*p_int*A_triang;
             
            if(p->j_dir==0)
            Fp_y = 0.0;
             
		   
    // Add viscous stress contributions
        
            double ustar, uplus, dist, value;
            double vstar,vplus,wstar,wplus;
            double uval, vval, wval;
            double kappa = 0.4;
            double xip, yip, zip, zval;
            double u_abs, tau, density;
            double ks;
            double dir;


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
            
            uval=p->ccipol1(uvel,xlocvel,ylocvel,zlocvel);
            vval=p->ccipol2(vvel,xlocvel,ylocvel,zlocvel);
            wval=p->ccipol3(wvel,xlocvel,ylocvel,zlocvel);
	        
            i = p->posc_i(xlocvel);
            j = p->posc_j(ylocvel);
            k = p->posc_k(zlocvel);
            
            double delta = sqrt(pow(xc-xlocvel,2.0) + pow(yc-ylocvel,2.0) + pow(zc-zlocvel,2.0));
            
            /*
            dudx = (uval)/(p->DXP[IP]);
            dudy = (uval)/(p->DYP[JP]);
            dudz = (uval)/(p->DZP[KP]);
                                                                           
            dvdx = (vval)/(p->DXP[IP]);
            dvdy = (vval)/(p->DYP[JP]);
            dvdz = (vval)/(p->DZP[KP]);
                                                                            
            dwdx = (wval)/(p->DXP[IP]);
            dwdy = (wval)/(p->DYP[JP]);
            dwdz = (wval)/(p->DZP[KP]);
            */
            dudx = (uval)/delta;
            dudy = (uval)/delta;
            dudz = (uval)/delta;
                                                                           
            dvdx = (vval)/delta;
            dvdy = (vval)/delta;
            dvdz = (vval)/delta;
                                                                            
            dwdx = (wval)/delta;
            dwdy = (wval)/delta;
            dwdz = (wval)/delta;
            
            u_abs = sqrt(uval*uval + vval*vval + wval*wval);

            Fv_x = -rho_int*(nu_int + enu_int)*A_triang*(2.0*dudx*nx + (dudy + dvdx)*ny + (dudz + dwdx)*nz);
            Fv_y = -rho_int*(nu_int + enu_int)*A_triang*((dudy + dvdx)*nx + 2.0*dvdy*ny + (dvdz + dwdy)*nz);
            Fv_z = -rho_int*(nu_int + enu_int)*A_triang*((dudz + dwdx)*nx + (dvdz + dwdy)*ny + 2.0*dwdz*nz);
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
            
            ks = 0.00001;
            
        
        // x-dir            
            uplus = (1.0/kappa)*log(30.0*(dist/ks));
            
            dir = uval/(fabs(uval)>0.0?uval:1.0e20);
  
            Fv_x =  dir*fabs(sqrt(ny*ny + nz*nz)*A_triang*rho_int*(uval*uval)/pow((uplus>0.0?uplus:1.0e20),2.0));
            
            //cout<<"Fv_x: "<<Fv_x<<" ustar: "<<ustar<<" nx: "<<nx<<" uplus: "<<uplus<<" dist: "<<dist<<" value: "<<value<<" log(value): "<<log(value)<<endl;
            
        // y-dir
            vplus = (1.0/kappa)*log(30.0*(dist/ks));
            
            dir = vval/(fabs(vval)>0.0?uval:1.0e20);
  
            Fv_y =  dir*fabs(sqrt(nx*nx + nz*nz)*A_triang*rho_int*(vval*vval)/pow((vplus>0.0?vplus:1.0e20),2.0));
            
        // z-dir
            wplus = (1.0/kappa)*log(30.0*(dist/ks));
            
            dir = wval/(fabs(wval)>0.0?wval:1.0e20);
  
            Fv_z = dir*fabs(sqrt(nx*nx + ny*ny)*A_triang*rho_int*(wval*wval)/pow((wplus>0.0?wplus:1.0e20),2.0));
            }
            
            
            
            if(p->X39==3)
            {
            xlocvel = xc + p->X43*nx*p->DXP[IP];
            ylocvel = yc + p->X43*ny*p->DYP[JP];
            zlocvel = zc + p->X43*nz*p->DZP[KP];
            
            nu_int  = p->ccipol4_a(a->visc,xlocvel,ylocvel,zlocvel);
            enu_int = p->ccipol4_a(a->eddyv,xlocvel,ylocvel,zlocvel);
            rho_int = p->ccipol4_a(a->ro,xlocvel,ylocvel,zlocvel);
            
            uval=p->ccipol1(uvel,xlocvel,ylocvel,zlocvel);
            vval=p->ccipol2(vvel,xlocvel,ylocvel,zlocvel);
            wval=p->ccipol3(wvel,xlocvel,ylocvel,zlocvel);
	        
            i = p->posc_i(xlocvel);
            j = p->posc_j(ylocvel);
            k = p->posc_k(zlocvel);
            
            double delta = sqrt(pow(xc-xlocvel,2.0) + pow(yc-ylocvel,2.0) + pow(zc-zlocvel,2.0));
            
                
            nu_int  = p->ccipol4_a(a->visc,xlocvel,ylocvel,zlocvel);
            enu_int = p->ccipol4_a(a->eddyv,xlocvel,ylocvel,zlocvel);
            rho_int = p->ccipol4_a(a->ro,xlocvel,ylocvel,zlocvel);


            u_abs = sqrt(uval*uval + vval*vval + wval*wval);
            
            tau=density*(nu_int + enu_int)*(u_abs/delta);
            
            Fv_x = -rho_int*(nu_int + enu_int)*A_triang*(2.0*dudx*nx + (dudy + dvdx)*ny + (dudz + dwdx)*nz);
            Fv_y = -rho_int*(nu_int + enu_int)*A_triang*((dudy + dvdx)*nx + 2.0*dvdy*ny + (dvdz + dwdy)*nz);
            Fv_z = -rho_int*(nu_int + enu_int)*A_triang*((dudz + dwdx)*nx + (dvdz + dwdy)*ny + 2.0*dwdz*nz);
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
    cout<<"Hydrodynamic Forces:  Fx_p: "<<Xe_p<<" Fy_p: "<<Ye_p<<" Fz_p: "<<Ze_p<<"  |  Fx_v: "<<Xe_v<<" Fy_v: "<<Ye_v<<" Fz_v: "<<Ze_v<<endl;
    
	Xe += a->gi*Mass_fb;
	Ye += a->gj*Mass_fb;
	Ze += a->gk*Mass_fb;


    // Print results
	
    if (p->mpirank == 0) 
    {
        ofstream print;
        char str[1000];
       
        sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_forces_%i.dat",n6DOF);

        print.open(str, std::ofstream::out | std::ofstream::app);
        print<<curr_time<<" \t "<<Xe<<" \t "<<Ye<<" \t "<<Ze<<" \t "<<Ke
        <<" \t "<<Me<<" \t "<<Ne<<" \t "<<Xe_p<<" \t "<<Ye_p<<" \t "<<Ze_p<<" \t "<<Xe_v<<" \t "<<Ye_v<<" \t "<<Ze_v<<endl;   
        print.close();
    }

	if (p->mpirank == 0)
	cout<<"Xe: "<<Xe<<" Ye: "<<Ye<<" Ze: "<<Ze<<" Ke: "<<Ke<<" Me: "<<Me<<" Ne: "<<Ne<<endl;
}
