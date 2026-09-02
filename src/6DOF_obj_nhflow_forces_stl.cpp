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

#include"6DOF_obj.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

void sixdof_obj::force_calc_stl(lexer* p, fdm_nhf *d, ghostcell *pgc, slice &WL, bool finalize)
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
    double A0;
    int cut;

    A0=A=0.0;
    Xe=Ye=Ze=Ke=Me=Ne=0.0;
    Xe_p=Ye_p=Ze_p=Xe_v=Ye_v=Ze_v=0.0;
    
    Fv_x = 0.0;
    Fv_y = 0.0;
    Fv_z = 0.0;
    
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

		// Ownership on the raw centroid: pure geometry, so every rank agrees.
		// Half-open intervals partition the domain exactly. No z test - NHFLOW
		// holds the full sigma column on every rank.
		if (xc >= p->originx && xc < p->endx &&
		   (p->j_dir==0 || (yc >= p->originy && yc < p->endy)))
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
            
            if(p->j_dir==0)
            ny=0.0;

            // ---- predictor clip: free surface at the raw centroid
            fsf_z = p->wd + p->ccslipol4(d->eta,xc,yc);

            if(clip_facet(p,x0,y0,z0,x1,y1,z1,x2,y2,z2,fsf_z,A_triang,xp,yp,zp)==false)
            continue;

            // ---- corrector clip: resample eta at the wetted centroid and
            //      re-clip the ORIGINAL facet against the updated plane
            //if(p->X61==1)
            {
            /*fsf_z = p->wd + p->ccslipol4(d->eta,xp,yp);

                if(clip_facet(p,x0,y0,z0,x1,y1,z1,x2,y2,z2,fsf_z,A_triang,xp,yp,zp)==false)
                continue;*/
            }

            // pressure: one free surface for both the clip and the head
            double off = 1.1;
            double xprobe = xp + off*nx*p->DXP[IP];
            double yprobe = yp + off*ny*p->DYP[JP]*p->y_dir;
            double zprobe = zp + off*nz*p->DZP[KP]*WL(i,j);

            //pval  = p->ccipol7V(d->P, WL, d->bed, xprobe, yprobe, zprobe);

            pval   = p->ccipol7V(d->P, WL, d->bed, xp, yp, zp);
            hspval = (fsf_z - zp)*p->W1*fabs(p->W22);

            Fp_x = -(pval + hspval)*A_triang*nx*f_jdir;
            Fp_y = -(pval + hspval)*A_triang*ny*f_jdir;
            Fp_z = -(pval + hspval)*A_triang*nz*f_jdir;
             
            // Total forces
            Fx = Fp_x;// + Fv_x;
            Fy = Fp_y;// + Fv_y;
            Fz = Fp_z;// + Fv_z;

            // Add forces to global forces
            Xe += Fx;
            Ye += Fy;
            Ze += Fz;

            Ke += (yp - c_(1))*Fz - (zp - c_(2))*Fy;
            Me += (zp - c_(2))*Fx - (xp - c_(0))*Fz;
            Ne += (xp - c_(0))*Fy - (yp - c_(1))*Fx;
            
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
    double ax,ay,az;
    double bx,by,bz;
    double cx,cy,cz;
    
    ax = x1-x0;
    ay = y1-y0;
    az = z1-z0;
    
    bx = x2-x0;
    by = y2-y0;
    bz = z2-z0;
    
    // twice the area is the norm of the cross product
    cx = ay*bz - az*by;
    cy = az*bx - ax*bz;
    cz = ax*by - ay*bx;
    
    return 0.5*sqrt(cx*cx + cy*cy + cz*cz);
}

// Parameter along the segment a->b where z crosses fsf_z.
double sixdof_obj::clip_edge(double za, double zb, double fsf_z)
{
    double dz,f;
    
    dz = zb - za;
    
    if(fabs(dz) < 1.0e-12)
    return 0.0;
    
    f = (fsf_z - za)/dz;
    
    return MIN(MAX(f,0.0),1.0);
}


// Clip a facet against the horizontal plane z = fsf_z and return the wetted
// area together with its area-weighted centroid.
// Returns false when the facet is entirely dry.
bool sixdof_obj::clip_facet(lexer *p,
                            double x0,double y0,double z0,
                            double x1,double y1,double z1,
                            double x2,double y2,double z2,
                            double fsf_z,
                            double &A_triang, double &xp, double &yp, double &zp)
{
    double xA,yA,zA,xB,yB,zB,xC,yC,zC;
    double xSb,ySb,zSb,xSc,ySc,zSc;
    double A1,A2,f;
    int nabove;

    A_triang = 0.0;
    xp = yp = zp = 0.0;

    nabove = 0;

    if(z0 > fsf_z)
    ++nabove;

    if(z1 > fsf_z)
    ++nabove;

    if(z2 > fsf_z)
    ++nabove;

    // fully dry
    if(nabove==3)
    return false;

    // fully wet
    if(nabove==0)
    {
    A_triang = triangle_area(p,x0,y0,z0,x1,y1,z1,x2,y2,z2);

    xp = (x0 + x1 + x2)/3.0;
    yp = (y0 + y1 + y2)/3.0;
    zp = (z0 + z1 + z2)/3.0;

    return true;
    }

    // Relabel so that A is the odd vertex out: the single dry one when
    // nabove==1, the single wet one when nabove==2. The cut points then
    // always lie on edges A->B and A->C. The cyclic order 0,1,2 -> 1,2,0
    // -> 2,0,1 preserves the facet winding.
    if(nabove==1)
    {
        if(z0 > fsf_z)
        {xA=x0;yA=y0;zA=z0; xB=x1;yB=y1;zB=z1; xC=x2;yC=y2;zC=z2;}

        else if(z1 > fsf_z)
        {xA=x1;yA=y1;zA=z1; xB=x2;yB=y2;zB=z2; xC=x0;yC=y0;zC=z0;}

        else
        {xA=x2;yA=y2;zA=z2; xB=x0;yB=y0;zB=z0; xC=x1;yC=y1;zC=z1;}
    }

    else
    {
        if(z0 <= fsf_z)
        {xA=x0;yA=y0;zA=z0; xB=x1;yB=y1;zB=z1; xC=x2;yC=y2;zC=z2;}

        else if(z1 <= fsf_z)
        {xA=x1;yA=y1;zA=z1; xB=x2;yB=y2;zB=z2; xC=x0;yC=y0;zC=z0;}

        else
        {xA=x2;yA=y2;zA=z2; xB=x0;yB=y0;zB=z0; xC=x1;yC=y1;zC=z1;}
    }

    // cut points on A->B and A->C
    f = clip_edge(zA,zB,fsf_z);

    xSb = xA + f*(xB-xA);
    ySb = yA + f*(yB-yA);
    zSb = zA + f*(zB-zA);

    f = clip_edge(zA,zC,fsf_z);

    xSc = xA + f*(xC-xA);
    ySc = yA + f*(yC-yA);
    zSc = zA + f*(zC-zA);


    // two vertices dry: the wet part is the triangle A-Sb-Sc
    if(nabove==2)
    {
    A_triang = triangle_area(p,xA,yA,zA,xSb,ySb,zSb,xSc,ySc,zSc);

        if(A_triang < 1.0e-20)
        return false;

    xp = (xA + xSb + xSc)/3.0;
    yp = (yA + ySb + ySc)/3.0;
    zp = (zA + zSb + zSc)/3.0;

    return true;
    }


    // one vertex dry: the wet part is the quad B-C-Sc-Sb, fanned from B
    A1 = triangle_area(p,xB,yB,zB,xC,yC,zC,xSc,ySc,zSc);
    A2 = triangle_area(p,xB,yB,zB,xSc,ySc,zSc,xSb,ySb,zSb);

    A_triang = A1 + A2;

    if(A_triang < 1.0e-20)
    return false;

    xp = (A1*(xB + xC + xSc) + A2*(xB + xSc + xSb))/(3.0*A_triang);
    yp = (A1*(yB + yC + ySc) + A2*(yB + ySc + ySb))/(3.0*A_triang);
    zp = (A1*(zB + zC + zSc) + A2*(zB + zSc + zSb))/(3.0*A_triang);

    return true;
}