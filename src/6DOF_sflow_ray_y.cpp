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

#include"6DOF_sflow.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"

void sixdof_sflow::ray_cast_y(lexer *p, ghostcell *pgc, int ts, int te)
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
	int ir,insidecheck;
	double u,v,w;
	double denom;
	double psi = 1.0e-8*p->DXM;	

	for(int n=ts; n<te; ++n)
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
	
	
	xs = MIN3(Ax,Bx,Cx);
	xe = MAX3(Ax,Bx,Cx);
	
	zs = MIN3(Az,Bz,Cz);
	ze = MAX3(Az,Bz,Cz);	
	
	is = p->posf_i(xs);
	ie = p->posf_i(xe);
	
	ks = p->posf_k(zs);
	ke = p->posf_k(ze);
    
	xs = MIN3(Ax,Bx,Cx) - epsi*p->DXP[is +marge];
	xe = MAX3(Ax,Bx,Cx) + epsi*p->DXP[ie +marge];
	
	zs = MIN3(Az,Bz,Cz) - epsi*p->DZP[ks +marge];
	ze = MAX3(Az,Bz,Cz) + epsi*p->DZP[ke +marge];
	
	
	is = p->posf_i(xs);
	ie = p->posf_i(xe);
	
	ks = p->posf_k(zs);
	ke = p->posf_k(ze);

	
	is = MAX(is,0);
	ie = MIN(ie,p->knox);
	
	ks = MAX(ks,0);
	ke = MIN(ke,p->knoz);
	
	
		for(i=is;i<ie;i++)
		{
		Px = p->XP[IP]+psi;
		Py = p->global_ymin-10.0*p->DXM;
		Pz = p->wd+psi;
		
		Qx = p->XP[IP]+psi;
		Qy = p->global_ymax+10.0*p->DXM;
		Qz = p->wd+psi;
		
		
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
		
			if((u>=0.0 && v>=0.0 && w>=0.0) || (u<0.0 && v<0.0 && w<0.0) && check==1)
			{
			denom = 1.0/(u+v+w);
			u *= denom;
			v *= denom;
			w *= denom;
			
			Ry = u*Ay + v*By + w*Cy;
			
			
            j = p->posf_j(Ry);
            
            int distcheck=1;
  
            
            if(Ry<p->YP[JP])
            if(j>=0 && j<p->knoy)
            if(fbio(i,j)<0 && fbio(i,j-1)<0)
            distcheck=0;
            
            if(Ry>=p->YP[JP])
            if(j>=0 && j<p->knoy)
            if(fbio(i,j)<0 && fbio(i,j+1)<0)
            distcheck=0;

            if(distcheck==1)
			for(j=0;j<p->knoy;++j)
			fb(i,j)=MIN(fabs(Ry-p->YP[JP]),fabs(fb(i,j)));
            
			for(j=0;j<p->knoy;++j)
            {
            Rymin(i,j) = MIN(Ry,Rymin(i,j));
            Rymax(i,j) = MAX(Ry,Rymax(i,j));
            }
			}
		}
	}

}

