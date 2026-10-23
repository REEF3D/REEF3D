/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"6DOF_obj_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"amr_grid_lookup.h"
#include<vector>
#include <math.h>

void sixdof_obj_cfd::forces_stl(lexer* p, fdm *a, ghostcell *pgc,field& uvel, field& vvel, field& wvel, int iter, bool finalize)
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
    
    // distributed loads (REEF3D::CORAL fluid_member): the hull forces binned along the member, this rank's part
    external_strips_zero();
    
    // Set new time
    curr_time = p->simtime;

    // the triangles of this rank: in its level-0 box; with the placement of the CFD patches (G 40)
    // the rank of the finest grid at the triangle
    std::vector<char> mine;
    if(amr_loads!=nullptr && amr_loads->placed())
    {
        std::vector<double> cx(tricount), cy(tricount), cz(tricount);
        std::vector<int> fl(tricount);
        for(int n=0; n<tricount; ++n)
        {
            cx[n] = (tri_x[n][0]+tri_x[n][1]+tri_x[n][2])/3.0;
            cy[n] = (tri_y[n][0]+tri_y[n][1]+tri_y[n][2])/3.0;
            cz[n] = (tri_z[n][0]+tri_z[n][1]+tri_z[n][2])/3.0;
            fl[n] = (cx[n]>=p->originx && cx[n]<p->endx && cy[n]>=p->originy && cy[n]<p->endy && cz[n]>=p->originz && cz[n]<p->endz) ? 1 : 0;
        }
        amr_cfd_mine(amr_loads,tricount,cx.data(),cy.data(),cz.data(),fl.data(),mine);
    }


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
        
        
		if (mine.empty() ? (xc >= p->originx && xc < p->endx &&
			yc >= p->originy && yc < p->endy &&
			zc >= p->originz && zc < p->endz) : mine[n]!=0)
		{
            // the grid of the loads: with the mesh refinement the finest grid at the triangle
            // (its pressure, density, viscosity, velocities and spacing), else level 0
            lexer *q = p;
            fdm *fa = a;
            field *U = &uvel, *V = &vvel, *W = &wvel;
            if(amr_loads!=nullptr && amr_loads->finest_grid(xc,yc,zc,q,fa))
            {
                U = &fa->u;
                V = &fa->v;
                W = &fa->w;
            }

            // Position of triangle
            i = q->posc_i(xc);
            j = q->posc_j(yc);
            k = q->posc_k(zc);
			
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
            
            if(q->j_dir==0)
            ny=0.0;
            
            // Add normal stress contributions
            xlocp = xc + q->X42*nx*q->DXP[IP];
            ylocp = yc + q->X42*ny*q->DYP[JP];
            zlocp = zc + q->X42*nz*q->DZP[KP];

            p_int = q->ccipol4a(fa->press,xlocp,ylocp,zlocp) - p->pressgage;
            
            Fp_x = -nx*p_int*A_triang;
            Fp_y = -ny*p_int*A_triang;
            Fp_z = -nz*p_int*A_triang;
             
            if(q->j_dir==0)
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


            if(q->X38==0)
            {
            xlocvel = xc + q->X43*nx*q->DXP[IP];
            ylocvel = yc + q->X43*ny*q->DYP[JP];
            zlocvel = zc + q->X43*nz*q->DZP[KP];
            
             nu_int = q->ccipol4a(fa->visc,xlocvel,ylocvel,zlocvel);
			enu_int = 0.0; //q->ccipol4a(fa->eddyv,xlocvel,ylocvel,zlocvel);
			rho_int = q->ccipol4a(fa->ro,xlocvel,ylocvel,zlocvel);
	        
             i = q->posc_i(xlocvel);
			j = q->posc_j(ylocvel);
			k = q->posc_k(zlocvel);
            
            dudx = ((*U)(i+1,j,k) - (*U)(i-1,j,k))/(q->DXP[IP] + q->DXP[IM1]);
            dudy = ((*U)(i,j+1,k) - (*U)(i,j-1,k))/(q->DYP[JP] + q->DYP[JM1]);
            dudz = ((*U)(i,j,k+1) - (*U)(i,j,k-1))/(q->DZP[KP] + q->DZP[KM1]);
                                                                           
            dvdx = ((*V)(i+1,j,k) - (*V)(i-1,j,k))/(q->DXP[IP] + q->DXP[IM1]);
            dvdy = ((*V)(i,j+1,k) - (*V)(i,j-1,k))/(q->DYP[JP] + q->DYP[JM1]);
            dvdz = ((*V)(i,j,k+1) - (*V)(i,j,k-1))/(q->DZP[KP] + q->DZP[KM1]);
                                                                            
            dwdx = ((*W)(i+1,j,k) - (*W)(i-1,j,k))/(q->DXP[IP] + q->DXP[IM1]);
            dwdy = ((*W)(i,j+1,k) - (*W)(i,j-1,k))/(q->DYP[JP] + q->DYP[JM1]);
            dwdz = ((*W)(i,j,k+1) - (*W)(i,j,k-1))/(q->DZP[KP] + q->DZP[KM1]);

            Fv_x = rho_int*(nu_int + enu_int)*A_triang*(2.0*dudx*nx + (dudy + dvdx)*ny + (dudz + dwdx)*nz);
            Fv_y = rho_int*(nu_int + enu_int)*A_triang*((dudy + dvdx)*nx + 2.0*dvdy*ny + (dvdz + dwdy)*nz);
            Fv_z = rho_int*(nu_int + enu_int)*A_triang*((dudz + dwdx)*nx + (dvdz + dwdy)*ny + 2.0*dwdz*nz);
            }
            
            
            if(q->X38==1)
            {
            xlocvel = xc + q->X43*nx*q->DXP[IP];
            ylocvel = yc + q->X43*ny*q->DYP[JP];
            zlocvel = zc + q->X43*nz*q->DZP[KP];
            
            nu_int  = q->ccipol4a(fa->visc,xlocvel,ylocvel,zlocvel);
            enu_int = q->ccipol4a(fa->eddyv,xlocvel,ylocvel,zlocvel);
            rho_int = q->ccipol4a(fa->ro,xlocvel,ylocvel,zlocvel);
            
            uval=q->ccipol1(*U,xlocvel,ylocvel,zlocvel);
            vval=q->ccipol2(*V,xlocvel,ylocvel,zlocvel);
            wval=q->ccipol3(*W,xlocvel,ylocvel,zlocvel);
	        
            i = q->posc_i(xlocvel);
            j = q->posc_j(ylocvel);
            k = q->posc_k(zlocvel);
            
            double delta = sqrt(pow(xc-xlocvel,2.0) + pow(yc-ylocvel,2.0) + pow(zc-zlocvel,2.0));
            
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
            
            
            if(q->X38==2)
            {
            xlocvel = xc + q->X43*nx*q->DXP[IP];
            ylocvel = yc + q->X43*ny*q->DYP[JP];
            zlocvel = zc + q->X43*nz*q->DZP[KP];
            
            nu_int  = q->ccipol4a(fa->visc,xlocvel,ylocvel,zlocvel);
            enu_int = q->ccipol4a(fa->eddyv,xlocvel,ylocvel,zlocvel);
            rho_int = q->ccipol4a(fa->ro,xlocvel,ylocvel,zlocvel);
	        
            i = q->posc_i(xlocvel);
            j = q->posc_j(ylocvel);
            k = q->posc_k(zlocvel);
            
            uval=q->ccipol1(*U,xlocvel,ylocvel,zlocvel);
            vval=q->ccipol2(*V,xlocvel,ylocvel,zlocvel);
            wval=q->ccipol3(*W,xlocvel,ylocvel,zlocvel);
            
            dist = fabs(q->ccipol4(fa->fb,xlocvel,ylocvel,zlocvel));
            
            ks = 0.00001;
            
        
        // x-dir            
            uplus = (1.0/kappa)*log(30.0*(dist/ks));
            
            dir = uval/(fabs(uval)>0.0?uval:1.0e20);
  
            Fv_x =  dir*fabs(sqrt(ny*ny + nz*nz)*A_triang*rho_int*(uval*uval)/pow((uplus>0.0?uplus:1.0e20),2.0));

        // y-dir
            vplus = (1.0/kappa)*log(30.0*(dist/ks));
            
            dir = vval/(fabs(vval)>0.0?vval:1.0e20);
  
            Fv_y =  dir*fabs(sqrt(nx*nx + nz*nz)*A_triang*rho_int*(vval*vval)/pow((vplus>0.0?vplus:1.0e20),2.0));
            
        // z-dir
            wplus = (1.0/kappa)*log(30.0*(dist/ks));
            
            dir = wval/(fabs(wval)>0.0?wval:1.0e20);
  
            Fv_z = dir*fabs(sqrt(nx*nx + ny*ny)*A_triang*rho_int*(wval*wval)/pow((wplus>0.0?wplus:1.0e20),2.0));
            }
            
            
            if(q->X38==3)
            {
            xlocvel = xc + q->X43*nx*q->DXP[IP];
            ylocvel = yc + q->X43*ny*q->DYP[JP];
            zlocvel = zc + q->X43*nz*q->DZP[KP];
            
            nu_int  = q->ccipol4a(fa->visc,xlocvel,ylocvel,zlocvel);
            enu_int = q->ccipol4a(fa->eddyv,xlocvel,ylocvel,zlocvel);
            rho_int = q->ccipol4a(fa->ro,xlocvel,ylocvel,zlocvel);
            
            uval=q->ccipol1(*U,xlocvel,ylocvel,zlocvel);
            vval=q->ccipol2(*V,xlocvel,ylocvel,zlocvel);
            wval=q->ccipol3(*W,xlocvel,ylocvel,zlocvel);
	        
            i = q->posc_i(xlocvel);
            j = q->posc_j(ylocvel);
            k = q->posc_k(zlocvel);
            
            double delta = sqrt(pow(xc-xlocvel,2.0) + pow(yc-ylocvel,2.0) + pow(zc-zlocvel,2.0));
            
                
            nu_int  = q->ccipol4a(fa->visc,xlocvel,ylocvel,zlocvel);
            enu_int = q->ccipol4a(fa->eddyv,xlocvel,ylocvel,zlocvel);
            rho_int = q->ccipol4a(fa->ro,xlocvel,ylocvel,zlocvel);


            u_abs = sqrt(uval*uval + vval*vval + wval*wval);
            
            tau=density*(nu_int + enu_int)*(u_abs/delta);
            
            Fv_x = -rho_int*(nu_int + enu_int)*A_triang*(2.0*dudx*nx + (dudy + dvdx)*ny + (dudz + dwdx)*nz);
            Fv_y = -rho_int*(nu_int + enu_int)*A_triang*((dudy + dvdx)*nx + 2.0*dvdy*ny + (dvdz + dwdy)*nz);
            Fv_z = -rho_int*(nu_int + enu_int)*A_triang*((dudz + dwdx)*nx + (dvdz + dwdy)*ny + 2.0*dwdz*nz);
            }
            
            
            if(q->X38==4)
            {
            //zval = s->bedzh(i,j) + 0.5*q->DZN[KP];
            
            //tau=density*pturb->ccipol_a_kinval(p,pgc,xip,yip,zval)*0.3;
            }


            // Total forces
            Fx = Fp_x + Fv_x;
            Fy = Fp_y + Fv_y;
            Fz = Fp_z + Fv_z;
            

			// Add forces to global forces
			Xe += Fx;
			Ye += Fy;
			Ze += Fz;
            
            // strips: the pressure of the triangle (constant) split at the strip boundaries, the viscous
            // force at its centre
            if(external_strips_on())
            {
                const double ta[3] = {x0,y0,z0}, tb[3] = {x1,y1,z1}, tc[3] = {x2,y2,z2};
                const double pm[3] = {p_int,p_int,p_int};
                external_strips_tri(ta,tb,tc,pm,nx,ny,nz);
                external_strips_add(xc,yc,zc,Fv_x,Fv_y,Fv_z);
            }

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
	
	Xe += a->gi*Mass_fb;
	Ye += a->gj*Mass_fb;
	Ze += a->gk*Mass_fb;


    // Print results
	
    if (p->mpirank==0 && finalize==1) 
    {

        printforce<<curr_time<<" \t "<<Xe<<" \t "<<Ye<<" \t "<<Ze<<" \t "<<Ke
        <<" \t "<<Me<<" \t "<<Ne<<" \t "<<Xe_p<<" \t "<<Ye_p<<" \t "<<Ze_p<<" \t "<<Xe_v<<" \t "<<Ye_v<<" \t "<<Ze_v<<endl;   

    }

}
