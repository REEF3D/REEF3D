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

void sixdof_gc::forces_triang(lexer* p, fdm *a, ghostcell *pgc)
{
    double ux,vy,wz,vel,pressure,density,viscosity;
    double du,dv,dw;
	double xlocvel,ylocvel,zlocvel;
    double xloc,yloc,zloc;
	double Fx,Fy,Fz,V;
	double Fzl,Fzr;
	double veff, txx,tyy,tzz,txy,txz,tyz;
	double Fxvisc,Fyvisc,Fzvisc; 
	double pval,fbval;
    double A_tot,A;
    
    double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    double xc,yc,zc;
    double at,bt,ct;
    double st;
    double nx,ny,nz,norm;
    double uval,vval,wval;
    double phival;

    double Ax=0.0;
    double Ay=0.0;
    
    pgc->start4(p,a->press,40);
	pgc->start4(p,a->press,401);
    pgc->start4(p,a->press,401);

    pgc->gcfb_update_extra_gcb(p,a,a->press);
    
    
    pgc->dgcpol(p,a->press,p->dgc4,p->dgc4_count,14);
    a->press.ggcpol(p);
    

    Fx=Fy=Fz=0.0;
    A_tot=0.0;
    
    Xe=Ye=Ze=Ke=Me=Ne=0.0;

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
            
            //----
            /*
			 xc-=p->originx;
            yc-=p->originy;
            zc-=p->originz;
            */
            i = p->posc_i(xc);
            j = p->posc_j(yc);
            k = p->posc_k(zc);
            
            nx = (a->fb(i+1,j,k)-a->fb(i-1,j,k))/(p->DXP[IM1] + p->DXP[IP]);
            ny = (a->fb(i,j+1,k)-a->fb(i,j-1,k))/(p->DYP[JM1] + p->DYP[JP]);
            nz = (a->fb(i,j,k+1)-a->fb(i,j,k-1))/(p->DZP[KM1] + p->DZP[KP]);
            
            norm = sqrt(nx*nx + ny*ny + nz*nz);
            
            nx/=norm>1.0e-20?norm:1.0e20;
            ny/=norm>1.0e-20?norm:1.0e20;
            nz/=norm>1.0e-20?norm:1.0e20;
            
            
            xloc = xc;
            yloc = yc;
            zloc = zc;
            
            xlocvel = xc + nx*p->DXP[IP];
            ylocvel = yc + ny*p->DYP[JP];
            zlocvel = zc + nz*p->DZP[KP];
            
            uval = p->ccipol1_a(a->u,xlocvel,ylocvel,zlocvel);
            vval = p->ccipol2_a(a->v,xlocvel,ylocvel,zlocvel);
            wval = p->ccipol3_a(a->w,xlocvel,ylocvel,zlocvel);
            
            du = uval/p->DXP[IP];
            dv = vval/p->DYP[JP];
            dw = wval/p->DZP[KP];
            
            pval =      p->ccipol4_a(a->press,xloc,yloc,zloc);
            density =   p->ccipol4_a(a->ro,xloc,yloc,zloc);
            viscosity = p->ccipol4_a(a->visc,xloc,yloc,zloc) + p->ccipol4_a(a->eddyv,xloc,yloc,zloc);
            phival =    p->ccipol4_a(a->phi,xloc,yloc,zloc);
     
   
            //cout<<"A: "<<A<<" nx: "<<nx<<" ny: "<<ny<<" nz: "<<nz<<" pval: "<<pval<<" Fx: "<<-(pval)*A*nx<<" Fy: "<<-(pval)*A*ny<<" Fz: "<<-(pval)*A*nz<<endl;
        
            // Force
            Fx = -(pval)*A*nx
                       + density*viscosity*A*(du*ny+du*nz);
                       
            Fy = -(pval)*A*ny
                       + density*viscosity*A*(dv*nx+dv*nz);
                    
            Fz = -(pval)*A*nz
                       + density*viscosity*A*(dw*nx+dw*ny);  
                       
            Xe += Fx;
            Ye += Fy;
            Ze += Fz;
	

            //xloc+=p->originx;
            //yloc+=p->originy;
            //zloc+=p->originz;
            
            Ke += (yloc-yg)*Fz - (zloc-zg)*Fy;
            Me += (zloc-zg)*Fx - (xloc-xg)*Fz;
            Ne += (xloc-xg)*Fy - (yloc-yg)*Fx;
       
    Ax+=A*nx;
    Ay+=A*ny;
                       
    A_tot+=A;
        
    }
        
        
    A_tot = pgc->globalsum(A_tot);
	
	Xe = pgc->globalsum(Xe);
	Ye = pgc->globalsum(Ye);
	Ze = pgc->globalsum(Ze);
	Ke = pgc->globalsum(Ke);
	Me = pgc->globalsum(Me);
	Ne = pgc->globalsum(Ne);

	Xe += a->gi*Mfb;
	Ye += a->gj*Mfb;
	Ze += a->gk*Mfb;

	if(p->mpirank==0)
	cout<<"area: "<<A_tot<<endl;
	
	if(p->mpirank==0)
	cout<<"Xe: "<<Xe<<" Ye: "<<Ye<<" Ze: "<<Ze<<" Ke: "<<Ke<<" Me: "<<Me<<" Ne: "<<Ne<<endl;
}


