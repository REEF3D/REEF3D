/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
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
#include"fieldint.h"

void sixdof_obj::ray_cast_io_zcorr(lexer *p, fdm *a, ghostcell *pgc, int ts, int te)
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
    
    ALOOP
	{
	cutl(i,j,k)=0;
	cutr(i,j,k)=0;
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
    && Ay>=p->global_ymin && Ay<=p->global_ymax
    && Az>=p->global_zmin && Az<=p->global_zmax)
    checkin=1;
    
    if(Bx>=p->global_xmin && Bx<=p->global_xmax 
    && By>=p->global_ymin && By<=p->global_ymax
    && Bz>=p->global_zmin && Bz<=p->global_zmax)
    checkin=1;
    
    if(Cx>=p->global_xmin && Cx<=p->global_xmax 
    && Cy>=p->global_ymin && Cy<=p->global_ymax
    && Cz>=p->global_zmin && Cz<=p->global_zmax)
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
		Pz = p->global_zmin-10.0*p->DXM ;
		
		Qx = p->XP[IP]+psi;
		Qy = p->YP[JP]-psi;
		Qz = p->global_zmax+10.0*p->DXM ;
		
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
				if(p->ZP[KP]<Rz)
				cutr(i,j,k) += 1;
				
				if(p->ZP[KP]>=Rz)
				cutl(i,j,k) += 1;
            }
            }
		
		}
	}
    }
    
    ALOOP
	if((cutl(i,j,k)+1)%2==0  && (cutr(i,j,k)+1)%2==0)
	fbio(i,j,k)=-1;

    /*
    count=0;
	ALOOP
	if(a->fb(i,j,k)>0)
	++count;
    
    count=pgc->globalisum(count);
    
    if(p->mpirank==0)
    cout<<"Number of active cells after fb_ray_io_z_corr: "<<count<<endl;*/
}
