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
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
/*
void sixdof_obj::ray_cast_io_zcorr(lexer *p, fdm_nhf *d, ghostcell *pgc, int ts, int te)
{
	double ys,ye,zs,ze;
	double Px,Py,Pz;
	double Qx,Qy,Qz;
	double Rx,Ry,Rz;
	double Ax,Ay,Az;
	double Bx,By,Bz;
	double Cx,Cy,Cz;
	double PQx,PQy,PQz;
	double PAx,PAy,PAz;
	double PBx,PBy,PBz;
	double PCx,PCy,PCz;
	double Mx,My,Mz;
	int is,ie,js,je,ks,ke;
	int ir;
    int checkin;
	double u,v,w;
	double denom;
	double psi = 1.0e-8*p->DXM;
    double margin = 2.0*p->DXM;
    
    LOOP
	{
	CL[IJK]=0;
	CR[IJK]=0;
	}


	for(n=ts; n<te; ++n)
	{ 
	Ax = tri_x[n][0];
	Ay = tri_y[n][0];
	Az = tri_z[n][0];
		
	Bx = tri_x[n][1];
	By = tri_y[n][1];
	Bz = tri_z[n][1];
		
	Cx = tri_x[n][2];
	Cy = tri_y[n][2];
	Cz = tri_z[n][2];
    
    checkin = 0;
    
	if(Ax>=p->global_xmin && Ax<=p->global_xmax 
    && ((Ay>=p->global_ymin && Ay<=p->global_ymax) || p->j_dir==0)
    && Az>=p->global_zmin && Az<=p->global_zmax)
    checkin=1;
    
    if(Bx>=p->global_xmin && Bx<=p->global_xmax 
    && ((By>=p->global_ymin && By<=p->global_ymax) || p->j_dir==0)
    && Bz>=p->global_zmin && Bz<=p->global_zmax)
    checkin=1;
    
    if(Cx>=p->global_xmin && Cx<=p->global_xmax 
    && ((Cy>=p->global_ymin && Cy<=p->global_ymax) || p->j_dir==0)
    && Cz>=p->global_zmin && Cz<=p->global_zmax)
    checkin=1;
    
    checkin=1;
        
    if(checkin==1)
    {
	xs = MIN3(Ax,Bx,Cx); 
	xe = MAX3(Ax,Bx,Cx);
	
	ys = MIN3(Ay,By,Cy);
	ye = MAX3(Ay,By,Cy);
	
	is = p->posc_i(xs);
	ie = p->posc_i(xe);
	
	js = p->posc_j(ys);
	je = p->posc_j(ye);
		
	
    xs = MIN3(Ax,Bx,Cx) - epsi*p->DXP[is + marge];
	xe = MAX3(Ax,Bx,Cx) + epsi*p->DXP[ie + marge];
	
	ys = MIN3(Ay,By,Cy) - epsi*p->DYP[js + marge];
	ye = MAX3(Ay,By,Cy) + epsi*p->DYP[je + marge];

	
	is = p->posc_i(xs);
	ie = p->posc_i(xe);
	
	js = p->posc_j(ys);
	je = p->posc_j(ye);
	
	is = MAX(is,0);
	ie = MIN(ie,p->knox);
	
	js = MAX(js,0);
	je = MIN(je,p->knoy);
	
		for(i=is;i<ie;i++)
		for(j=js;j<je;j++)
		{
		Px = p->XP[IP]-psi;
		Py = p->YP[JP]+psi;
		Pz = zmin-10.0*p->DXM ;
		
		Qx = p->XP[IP]+psi;
		Qy = p->YP[JP]-psi;
		Qz = zmax+10.0*p->DXM ;
		
		PQx = Qx-Px;
		PQy = Qy-Py;
		PQz = Qz-Pz;
		
		PAx = Ax-Px;
		PAy = Ay-Py;
		PAz = Az-Pz;
		
		PBx = Bx-Px;
		PBy = By-Py;
		PBz = Bz-Pz;
		
		PCx = Cx-Px;
		PCy = Cy-Py;
		PCz = Cz-Pz;
		
		// uvw
		Mx = PQy*Pz - PQz*Py;
		My = PQz*Px - PQx*Pz;
		Mz = PQx*Py - PQy*Px;

		
		u = PQx*(Cy*Bz - Cz*By) + PQy*(Cz*Bx - Cx*Bz) + PQz*(Cx*By - Cy*Bx)
		  + Mx*(Cx-Bx) + My*(Cy-By) + Mz*(Cz-Bz);
		  
		v = PQx*(Ay*Cz - Az*Cy) + PQy*(Az*Cx - Ax*Cz) + PQz*(Ax*Cy - Ay*Cx)
		  + Mx*(Ax-Cx) + My*(Ay-Cy) + Mz*(Az-Cz);
		  
		w = PQx*(By*Az - Bz*Ay) + PQy*(Bz*Ax - Bx*Az) + PQz*(Bx*Ay - By*Ax)
		  + Mx*(Bx-Ax) + My*(By-Ay) + Mz*(Bz-Az);
    
         int check=1;
		if(u==0.0 && v==0.0 && w==0.0)
		check = 0;

			if(((u>0.0 && v>0.0 && w>0.0) || (u<0.0 && v<0.0 && w<0.0)) && check==1)
			{
			denom = 1.0/(u+v+w);
			u *= denom;
			v *= denom;
			w *= denom;
			
			Rz = u*Az + v*Bz + w*Cz;
			
            for(k=0;k<p->knoz;++k)
            {
				if(p->ZSP[IJK]<Rz)
				CL[IJK] += 1;
				
				if(p->ZSP[IJK]>=Rz)
				CR[IJK] += 1;
            }
            }
		}
	}
    }
    
    LOOP
    WETDRY
	if((CL[IJK]+1)%2==0  && (CR[IJK]+1)%2==0)
	IO[IJK]=-1;
}
*/


void sixdof_obj::ray_cast_io_zcorr(lexer *p, fdm_nhf *d, ghostcell *pgc, int ts, int te)
{
	double xs,xe,ys,ye;
	double Ax,Ay,Az;
	double Bx,By,Bz;
	double Cx,Cy,Cz;
	double abx,aby,bcx,bcy,cax,cay;
	double e0,e1,e2;
	double area2,sgn,denom;
	double Px,Py,Rz;
	int is,ie,js,je;

    const double eps_area = 1.0e-20*DSM*DSM;
    const double margin   = 2.0*DSM;
    const double psi      = 1.0e-8*DSM;
    

    LOOP
	{
	CL[IJK]=0;
	CR[IJK]=0;
	}


	for(n=ts; n<te; ++n)
	{
	Ax = tri_x[n][0];
	Ay = tri_y[n][0];
	Az = tri_z[n][0];

	Bx = tri_x[n][1];
	By = tri_y[n][1];
	Bz = tri_z[n][1];

	Cx = tri_x[n][2];
	Cy = tri_y[n][2];
	Cz = tri_z[n][2];

    // triangle footprint in x-y
	xs = MIN3(Ax,Bx,Cx);
	xe = MAX3(Ax,Bx,Cx);

	ys = MIN3(Ay,By,Cy);
	ye = MAX3(Ay,By,Cy);

    // AABB reject: footprint does not overlap this subdomain
    if(xe < p->originx-margin || xs > p->endx+margin)
    continue;

    if(p->j_dir==1 && (ye < p->originy-margin || ys > p->endy+margin))
    continue;

    // projected edge vectors, loop invariant
    abx = Bx-Ax;
    aby = By-Ay;

    bcx = Cx-Bx;
    bcy = Cy-By;

    cax = Ax-Cx;
    cay = Ay-Cy;

    // twice the signed area of the projected triangle: (B-A) x (C-A)
    area2 = -abx*cay + aby*cax;

    // vertical triangle: cannot be crossed by a vertical ray
    if(fabs(area2) < eps_area)
    continue;

    // normalise the winding to counter-clockwise
    sgn   = (area2 > 0.0 ? 1.0 : -1.0);
    denom = 1.0/(sgn*area2);

	is = MAX(p->posc_i(xs)-1, 0);
	ie = MIN(p->posc_i(xe)+2, p->knox);

    js = 0;
    je = 1;

    if(p->j_dir==1)
    {
	js = MAX(p->posc_j(ys)-1, 0);
	je = MIN(p->posc_j(ye)+2, p->knoy);
    }


		for(i=is;i<ie;i++)
		for(j=js;j<je;j++)
		{

		Px = p->XP[IP] - psi;
		Py = p->YP[JP] + psi;

        // projected edge functions: e0 opposite C, e1 opposite A, e2 opposite B
        e0 = sgn*(abx*(Py-Ay) - aby*(Px-Ax));
        e1 = sgn*(bcx*(Py-By) - bcy*(Px-Bx));
        e2 = sgn*(cax*(Py-Cy) - cay*(Px-Cx));

        if(e0 < 0.0 || e1 < 0.0 || e2 < 0.0)
        continue;

        // barycentric interpolation of the crossing height
        Rz = (e1*Az + e2*Bz + e0*Cz)*denom;

            for(k=0;k<p->knoz;++k)
            {
				if(p->ZSP[IJK] < Rz)
				CL[IJK] += 1;

                else
				CR[IJK] += 1;
            }
		}
	}


    LOOP
	if((CL[IJK]+1)%2==0  && (CR[IJK]+1)%2==0)
	IO[IJK]=-1;
    

}
