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

void sixdof_obj::forces_lsm(lexer* p, fdm *a, ghostcell *pgc,field& uvel, field& vvel, field& wvel, int iter)
{
    triangulation(p,a,pgc,a->fb);
	reconstruct(p,a,a->fb);
    
    //print_vtp(p,a,pgc);
    
    forces_lsm_calc(p,a,pgc);
    
    
    p->del_Iarray(tri,numtri,4);
    p->del_Darray(pt,numvert,3);
    p->del_Darray(ls,numvert);
    p->del_Iarray(facet,numtri,4);
    p->del_Iarray(confac,numtri);
    p->del_Iarray(numfac,numtri);
	p->del_Iarray(numpt,numtri);
    p->del_Darray(ccpt,numtri*4,3);
    
}

void sixdof_obj::forces_lsm_calc(lexer* p, fdm *a, ghostcell *pgc)
{
    double ux,vy,wz,vel,pressure,density,viscosity;
    double du,dv,dw;
    double xloc,yloc,zloc;
	double xlocvel,ylocvel,zlocvel;
    double sgnx,sgny,sgnz;
    double Ax=0.0;
    double Ay=0.0;
    double A_tot,A;
    double Px=0.0;
    double xp1,xp2,yp1,yp2,zp1,zp2;
    double Fx,Fy,Fz,Fp_x,Fp_y,Fp_z,Fv_x,Fv_y,Fv_z;
    double Xe_p,Ye_p,Ze_p,Xe_v,Ye_v,Ze_v;
    

    Fx=Fy=Fz=0.0;
    A_tot=0.0;
    
    polygon_num=facount;
	
	polygon_sum=0;
	for(n=0;n<polygon_num;++n)
	polygon_sum+=numpt[n];
	

	
	vertice_num = ccptcount;
    
    if(numtri>0)
    cout<<"FORCES"<<endl;
    

    pgc->dgcpol(p,a->press,p->dgc4,p->dgc4_count,14);
    a->press.ggcpol(p);
    
    for(n=0;n<polygon_num;++n)
    {       
            // triangle
            if(numpt[n]==3)
            {
            // tri1
            x1 = ccpt[facet[n][0]][0];
            y1 = ccpt[facet[n][0]][1];
            z1 = ccpt[facet[n][0]][2];
            
            x2 = ccpt[facet[n][1]][0];
            y2 = ccpt[facet[n][1]][1];
            z2 = ccpt[facet[n][1]][2];
            
            x3 = ccpt[facet[n][2]][0];
            y3 = ccpt[facet[n][2]][1];
            z3 = ccpt[facet[n][2]][2];
            
            xc = (1.0/3.0)*(x1 + x2 + x3);
            yc = (1.0/3.0)*(y1 + y2 + y3);
            zc = (1.0/3.0)*(z1 + z2 + z3);
            
            at = sqrt(pow(x2-x1,2.0) + pow(y2-y1,2.0) + pow(z2-z1,2.0));
            bt = sqrt(pow(x2-x3,2.0) + pow(y2-y3,2.0) + pow(z2-z3,2.0));
            ct = sqrt(pow(x3-x1,2.0) + pow(y3-y1,2.0) + pow(z3-z1,2.0));
            
            st = 0.5*(at+bt+ct);
            
            A = sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
            }
            
            //quadrilidral
            if(numpt[n]==4)
            {
            x1 = ccpt[facet[n][0]][0];
            y1 = ccpt[facet[n][0]][1];
            z1 = ccpt[facet[n][0]][2];
            
            x2 = ccpt[facet[n][1]][0];
            y2 = ccpt[facet[n][1]][1];
            z2 = ccpt[facet[n][1]][2];
            
            x3 = ccpt[facet[n][3]][0];
            y3 = ccpt[facet[n][3]][1];
            z3 = ccpt[facet[n][3]][2];
            
            x4 = ccpt[facet[n][2]][0];
            y4 = ccpt[facet[n][2]][1];
            z4 = ccpt[facet[n][2]][2];
            
            xc = (1.0/4.0)*(x1 + x2 + x3 + x4);
            yc = (1.0/4.0)*(y1 + y2 + y3 + y4);
            zc = (1.0/4.0)*(z1 + z2 + z3 + z4);
            
            
            //tri1
            at = sqrt(pow(x2-x1,2.0) + pow(y2-y1,2.0) + pow(z2-z1,2.0));
            bt = sqrt(pow(x2-x3,2.0) + pow(y2-y3,2.0) + pow(z2-z3,2.0));
            ct = sqrt(pow(x3-x1,2.0) + pow(y3-y1,2.0) + pow(z3-z1,2.0));
            
            st = 0.5*(at+bt+ct);

            A = sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
            
            //tri2
            at = sqrt(pow(x3-x1,2.0) + pow(y3-y1,2.0) + pow(z3-z1,2.0));
            bt = sqrt(pow(x4-x3,2.0) + pow(y4-y3,2.0) + pow(z4-z3,2.0));
            ct = sqrt(pow(x4-x1,2.0) + pow(y4-y1,2.0) + pow(z4-z1,2.0));
            
            st = 0.5*(at+bt+ct);

            A += sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
            }
            
            xp1 = x2-x1;
            yp1 = y2-y1;
            zp1 = z2-z1;
            
            xp2 = x3-x1;
            yp2 = y3-y1;
            zp2 = z3-z1;
        
            
            nx=fabs(yp1*zp2 - zp1*yp2);
            ny=fabs(zp1*xp2 - xp1*zp2);
            nz=fabs(xp1*yp2 - yp1*xp2);
            
            norm = sqrt(nx*nx + ny*ny + nz*nz);
            
            nx/=norm>1.0e-20?norm:1.0e20;
            ny/=norm>1.0e-20?norm:1.0e20;
            nz/=norm>1.0e-20?norm:1.0e20;
            //----

            
            i = p->posc_i(xc);
            j = p->posc_j(yc);
            k = p->posc_k(zc);
            
            sgnx = (a->solid(i+1,j,k)-a->solid(i-1,j,k))/(p->DXP[IM1] + p->DXP[IP]);
            sgny = (a->solid(i,j+1,k)-a->solid(i,j-1,k))/(p->DYP[JM1] + p->DYP[JP]);
            sgnz = (a->solid(i,j,k+1)-a->solid(i,j,k-1))/(p->DZP[KM1] + p->DZP[KP]);
            
            nx = nx*sgnx/fabs(fabs(sgnx)>1.0e-20?sgnx:1.0e20);
            ny = ny*sgny/fabs(fabs(sgny)>1.0e-20?sgny:1.0e20);
            nz = nz*sgnz/fabs(fabs(sgnz)>1.0e-20?sgnz:1.0e20);
            
            
            xloc = xc + nx*p->DXP[IP]*p->P91;
            yloc = yc + ny*p->DYP[JP]*p->P91;
            zloc = zc + nz*p->DZP[KP]*p->P91;
            
            xlocvel = xc + nx*p->DXP[IP];
            ylocvel = yc + ny*p->DYP[JP];
            zlocvel = zc + nz*p->DZP[KP];
            
            uval = p->ccipol1_a(a->u,xlocvel,ylocvel,zlocvel);
            vval = p->ccipol2_a(a->v,xlocvel,ylocvel,zlocvel);
            wval = p->ccipol3_a(a->w,xlocvel,ylocvel,zlocvel);
            
            du = uval/p->DXN[IP];
            dv = vval/p->DYN[JP];
            dw = wval/p->DZN[KP];
            
            pval =      p->ccipol4_a(a->press,xloc,yloc,zloc) - p->pressgage;
            density =   p->ccipol4_a(a->ro,xloc,yloc,zloc);
            viscosity = p->ccipol4_a(a->visc,xloc,yloc,zloc);
            phival =    p->ccipol4_a(a->phi,xloc,yloc,zloc);
            
            if(p->P82==1)
            viscosity += p->ccipol4_a(a->eddyv,xloc,yloc,zloc);   
            
            i = p->posc_i(xloc);
            j = p->posc_j(yloc);
            k = p->posc_k(zloc);
            
            // Force
            if(phival>-1.6*p->DXM || p->P92==1)
            {
            Fx = -(pval)*A*nx
                       + density*viscosity*A*(du*ny+du*nz);
                       
            Fy = -(pval)*A*ny
                       + density*viscosity*A*(dv*nx+dv*nz);
                    
            Fz = -(pval)*A*nz
                       + density*viscosity*A*(dw*nx+dw*ny);   
                       
            }
    Ax+=A*nx;
    Ay+=A*ny;
    
    Px += pval*nx/fabs(nx);
    
    
    
    
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
							
			A_tot+=A;
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
    
    Ax = pgc->globalsum(Ax);
    Ay = pgc->globalsum(Ay);
    A_tot = pgc->globalsum(A_tot);
    Px = pgc->globalsum(Px);
    
    //if(p->mpirank==0)
    //cout<<"Ax : "<<Ax<<" Ay: "<<Ay<<" A_tot: "<<A_tot<<endl;
    
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