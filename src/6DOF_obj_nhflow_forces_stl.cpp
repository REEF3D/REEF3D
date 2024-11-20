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
#include"gradient.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

void sixdof_obj::force_calc_stl(lexer* p, fdm_nhf *d, ghostcell *pgc, bool finalize)
{
    double x0,x1,x2,y0,y1,y2,z0,z1,z2;
    double xs0,xs1,xs2,ys0,ys1,ys2,zs0,zs1,zs2;
	double xc,yc,zc;
    double xp,yp,zp;
	double at,bt,ct,st;
	double nx,ny,nz,norm;
	double A_triang,A,A_red;
    double f;
	double pval,rho_int,nu_int,enu_int,u_int,v_int,w_int;
	double du,dv,dw, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
	double dudxf, dudyf, dudzf, dvdxf, dvdyf, dvdzf, dwdxf, dwdyf, dwdzf;
	double dudxb, dudyb, dudzb, dvdxb, dvdyb, dvdzb, dwdxb, dwdyb, dwdzb;
	double xlocvel,ylocvel,zlocvel,xlocp,ylocp,zlocp;
	double Fx,Fy,Fz,Fp_x,Fp_y,Fp_z,Fv_x,Fv_y,Fv_z;
    double Xe_p,Ye_p,Ze_p,Xe_v,Ye_v,Ze_v;
    double fsf_z;
    double f_jdir;

    A=0.0;
    Xe=Ye=Ze=Ke=Me=Ne=0.0;
    Xe_p=Ye_p=Ze_p=Xe_v=Ye_v=Ze_v=0.0;
    
    // Set new time
    curr_time = p->simtime;

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
        
        xp = xc;
        yp = yc;
        zp = zc;
    
 
		if (xc >= p->originx && xc < p->endx &&
			yc >= p->originy && yc < p->endy &&
			zc >= p->originz && zc < p->endz)
		{
            
            
            // Normal vectors (always pointing outwards)     
			nx = (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0);
			ny = (x2 - x0) * (z1 - z0) - (x1 - x0) * (z2 - z0); 
			nz = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);

			norm = sqrt(nx*nx + ny*ny + nz*nz);
			
			nx /= norm > 1.0e-20 ? norm : 1.0e20;
			ny /= norm > 1.0e-20 ? norm : 1.0e20;
			nz /= norm > 1.0e-20 ? norm : 1.0e20;
            
            f_jdir=1.0;
            
            if(fabs(ny)>0.9 && p->j_dir==0)
            f_jdir=0.0;
            
            
            // Position of triangle
            i = p->posc_i(xc);
            j = p->posc_j(yc);
            k = p->posc_k(zc);
            
            etaval = p->ccslipol4(d->eta,xc,yc);  
            
            fsf_z = p->wd + etaval;
            

            //if(zc<fsf_z)
            if(z0<fsf_z || z1<fsf_z || z2<fsf_z)
            {
            // Area of triangle using Heron's formula
			A_triang = triangle_area(p,x0,y0,z0,x1,y1,z1,x2,y2,z2);
            

        // peak up 0
            if(z0>fsf_z && z1<fsf_z  && z2<fsf_z)
            {
            // Ps1
             f = fabs(fsf_z - z1)/(fabs(z0-z1)>1.0e-10?fabs(z0-z1):1.0e10);

             xs1 = x1 + f*(x0-x1);
             ys1 = y1 + f*(y0-y1);
             zs1 = z1 + f*(z0-z1);
             
             // Ps2
             f = fabs(fsf_z - z2)/(fabs(z0-z2)>1.0e-10?fabs(z0-z2):1.0e10);

             xs2 = x2 + f*(x0-x2);
             ys2 = y2 + f*(y0-y2);
             zs2 = z2 + f*(z0-z2);
             
             xp = (x1 + x2 + xs1 + xs2)/4.0;
             yp = (y1 + y2 + ys1 + ys2)/4.0;
             zp = (z1 + z2 + zs1 + zs2)/4.0;
             
             //cout<<"A_triang: "<<A_triang<<" A: "<<triangle_area(p,x0,y0,z0,xs1,ys1,zs1,xs2,ys2,zs2)<<endl;
                
            A_triang = A_triang - triangle_area(p,x0,y0,z0,xs1,ys1,zs1,xs2,ys2,zs2);
            }
            
        // peak up 1
            if(z1>fsf_z && z0<fsf_z  && z2<fsf_z)
            {
            // Ps0
             f = fabs(fsf_z - z0)/(fabs(z1-z0)>1.0e-10?fabs(z1-z0):1.0e10);
            
             xs0 = x0 + f*(x1-x0);
             ys0 = y0 + f*(y1-y0);
             zs0 = z0 + f*(z1-z0);
             
             // Ps2
             f = fabs(fsf_z - z2)/(fabs(z1-z2)>1.0e-10?fabs(z1-z2):1.0e10);
             
             xs2 = x2 + f*(x1-x2);
             ys2 = y2 + f*(y1-y2);
             zs2 = z2 + f*(z1-z2);
             
             xp = (x0 + x2 + xs0 + xs2)/4.0;
             yp = (y0 + y2 + ys0 + ys2)/4.0;
             zp = (z0 + z2 + zs0 + zs2)/4.0;
                
            A_triang = A_triang - triangle_area(p,xs0,ys0,zs0,x1,y1,z1,xs2,ys2,zs2);
            }
            
        // peak up 2
            if(z2>fsf_z && z0<fsf_z  && z1<fsf_z)
            {
            // Ps0
             f = fabs(fsf_z - z0)/(fabs(z2-z0)>1.0e-10?fabs(z2-z0):1.0e10);

             xs0 = x0 + f*(x2-x0);
             ys0 = y0 + f*(y2-y0);
             zs0 = z0 + f*(z2-z0);
             
             // Ps1
             f = fabs(fsf_z - z1)/(fabs(z1-z2)>1.0e-10?fabs(z1-z2):1.0e10);
             
             xs1 = x1 + f*(x2-x1);
             ys1 = y1 + f*(y2-y1);
             zs1 = z1 + f*(z2-z1);
             
             xp = (x0 + x1 + xs0 + xs1)/4.0;
             yp = (y0 + y1 + ys0 + ys1)/4.0;
             zp = (z0 + z1 + zs0 + zs1)/4.0;
                
            A_triang = A_triang - triangle_area(p,xs0,ys0,zs0,xs1,ys1,zs1,x2,y2,z2);
            }
        
        // -----
        
         // peak down 0
            if(z0<fsf_z && z1>fsf_z  && z2>fsf_z)
            {
            // Ps1
             f = fabs(fsf_z - z0)/(fabs(z0-z1)>1.0e-10?fabs(z0-z1):1.0e10);
             
             xs1 = x0 + f*(x1-x0);
             ys1 = y0 + f*(y1-y0);
             zs1 = z0 + f*(z1-z0);
             
             // Ps2
             f = fabs(fsf_z - z0)/(fabs(z0-z2)>1.0e-10?fabs(z0-z2):1.0e10);
             
             xs2 = x0 + f*(x2-x0);
             ys2 = y0 + f*(y2-y0);
             zs2 = z0 + f*(z2-z0);
             
             xp = (x0 + xs1 + xs2)/3.0;
             yp = (y0 + ys1 + ys2)/3.0;
             zp = (z0 + zs1 + zs2)/3.0;
                
            A_triang = triangle_area(p,x0,y0,z0,xs1,ys1,zs1,xs2,ys2,zs2);
            }
            
            
            // peak down 1
            if(z1<fsf_z && z0>fsf_z  && z2>fsf_z)
            {
            // Ps0
             f = fabs(fsf_z - z1)/(fabs(z1-z0)>1.0e-10?fabs(z1-z0):1.0e10);
             
             xs0 = x1 + f*(x0-x1);
             ys0 = y1 + f*(y0-y1);
             zs0 = z1 + f*(z0-z1);
             
             // Ps2
             f = fabs(fsf_z - z1)/(fabs(z1-z2)>1.0e-10?fabs(z1-z2):1.0e10);
             
             xs2 = x1 + f*(x2-x1);
             ys2 = y1 + f*(y2-y1);
             zs2 = z1 + f*(z2-z1);
             
             xp = (xs0 + x1 + xs2)/3.0;
             yp = (ys0 + y1 + ys2)/3.0;
             zp = (zs0 + z1 + zs2)/3.0;
                
            //cout<<"A_0: "<<A_triang<<" A_triang: "<<triangle_area(p,xs0,ys0,zs0,x1,y1,z1,xs2,ys2,zs2)<<" f: "<<f<<endl;
            
            A_triang = triangle_area(p,xs0,ys0,zs0,x1,y1,z1,xs2,ys2,zs2);
            }
            
        // peak down 2
            if(z2<fsf_z && z0>fsf_z  && z1>fsf_z)
            {
            // Ps0
             f = fabs(fsf_z - z2)/(fabs(z2-z0)>1.0e-10?fabs(z2-z0):1.0e10);
               
             xs0 = x2 + f*(x0-x2);
             ys0 = y2 + f*(y0-y2);
             zs0 = z2 + f*(z0-z2);
             
             // Ps1
             f = fabs(fsf_z - z2)/(fabs(z1-z2)>1.0e-10?fabs(z1-z2):1.0e10);

             xs1 = x2 + f*(x1-x2);
             ys1 = y2 + f*(y1-y2);
             zs1 = z2 + f*(z1-z2);
             
             xp = (xs0 + xs1 + x2)/3.0;
             yp = (ys0 + ys1 + y2)/3.0;
             zp = (zs0 + zs1 + z2)/3.0;
                
            A_triang = triangle_area(p,xs0,ys0,zs0,xs1,ys1,zs1,x2,y2,z2);
            }
   
            if(p->j_dir==0)
            ny=0.0;
            
            // Add normal stress contributions
            xlocp = xc + p->X42*nx*p->DXP[IP];
            ylocp = yc + p->X42*ny*p->DYP[JP];
            zlocp = zc + p->X42*nz*p->DZP[KP];
            
            /*
            double p0,p1,p2,pc;
            
            p0   = p->ccipol7P(d->P, d->WL, d->bed, x0, y0, z0);
            p1   = p->ccipol7P(d->P, d->WL, d->bed, x1, y1, z1);
            p2   = p->ccipol7P(d->P, d->WL, d->bed, x2, y2, z2);
            
            pc   = p->ccipol7P(d->P, d->WL, d->bed, xc, yc, zc);
            
            pval = (1.0/4.0)*(p0 + p1 + p2 + pc);*/
            
          

            // pressure
            pval   = p->ccipol7P(d->P, d->WL, d->bed, xp, yp, zp);// - p->pressgage;
            etaval = p->ccslipol4(d->eta,xp,yp);  
            hspval = (p->wd + etaval - zp)*p->W1*fabs(p->W22);

            Fp_x = -(pval + hspval)*A_triang*nx*f_jdir;
            Fp_y = -(pval + hspval)*A_triang*ny*f_jdir;
            Fp_z = -(pval + hspval)*A_triang*nz*f_jdir;
             
            if(p->j_dir==0)
            Fp_y = 0.0;
             
            // Total forces
            Fx = Fp_x;// + Fv_x;
            Fy = Fp_y;// + Fv_y;
            Fz = Fp_z;// + Fv_z;

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
							
            A += A_triang*nz;
            }
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
    
    Fx = Xe;
    Fy = Ye;
    Fz = Ze;

	Xe_p = pgc->globalsum(Xe_p);
	Ye_p = pgc->globalsum(Ye_p);
	Ze_p = pgc->globalsum(Ze_p);
	Xe_v = pgc->globalsum(Xe_v);
	Ye_v = pgc->globalsum(Ye_v);
	Ze_v = pgc->globalsum(Ze_v);

	// Add gravity force
	Xe += p->W20*Mass_fb;
	Ye += p->W21*Mass_fb;
	Ze += p->W22*Mass_fb;
    
    if(p->mpirank==0)
    {
    cout<<"Mass_fb: "<<Mass_fb<<" G_fb: "<<p->W22*Mass_fb<<endl;
    cout<<"Fx: "<<Fx<<" Fy: "<<Fy<<" Fz: "<<Fz<<endl;
    cout<<"A_tot: "<<A<<endl;
    cout<<"Xe: "<<Xe<<" Ye: "<<Ye<<" Ze: "<<Ze<<" Ke: "<<Ke<<" Me: "<<Me<<" Ne: "<<Ne<<endl;
    }

    // Print results	
    if (p->mpirank==0 && finalize==1) 
    {
        printforce<<curr_time<<" \t "<<Xe<<" \t "<<Ye<<" \t "<<Ze<<" \t "<<Ke
        <<" \t "<<Me<<" \t "<<Ne<<" \t "<<Xe_p<<" \t "<<Ye_p<<" \t "<<Ze_p<<" \t "<<Xe_v<<" \t "<<Ye_v<<" \t "<<Ze_v<<endl;   
    }
}


double sixdof_obj::triangle_area(lexer *p, double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2)
{
    double at,bt,ct,st,A;
    
    at = sqrt(pow(x1-x0,2.0) + pow(y1-y0,2.0) + pow(z1-z0,2.0));
    bt = sqrt(pow(x1-x2,2.0) + pow(y1-y2,2.0) + pow(z1-z2,2.0));
    ct = sqrt(pow(x2-x0,2.0) + pow(y2-y0,2.0) + pow(z2-z0,2.0));
				
    st = 0.5*(at+bt+ct);
				
    A = sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
    
    return A;
            
}
  