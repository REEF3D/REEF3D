/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"6DOF_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fieldint.h"

void sixdof_f::ray_cast_z(lexer *p, fdm *a, ghostcell *pgc, int ts, int te)
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
	double u,v,w;
	double denom;
	double psi = 1.0e-6*p->DXM;


	for(n=ts; n<te; ++n)
	{ 
		
	Ax = tri_x[n][0]-p->originx;
	Ay = tri_y[n][0]-p->originy;
	Az = tri_z[n][0]-p->originz;
		
	Bx = tri_x[n][1]-p->originx;
	By = tri_y[n][1]-p->originy;
	Bz = tri_z[n][1]-p->originz;
		
	Cx = tri_x[n][2]-p->originx;
	Cy = tri_y[n][2]-p->originy;
	Cz = tri_z[n][2]-p->originz;
	
	
	xs = MIN3(Ax,Bx,Cx); 
	xe = MAX3(Ax,Bx,Cx);
	
	ys = MIN3(Ay,By,Cy);
	ye = MAX3(Ay,By,Cy);
	
	is = p->posf_i(xs+p->originx);
	ie = p->posf_i(xe+p->originx);
	
	js = p->posf_j(ys+p->originy);
	je = p->posf_j(ye+p->originy);
		
	
    xs = MIN3(Ax,Bx,Cx) - epsi*p->DXP[is + marge];
	xe = MAX3(Ax,Bx,Cx) + epsi*p->DXP[ie + marge];
	
	ys = MIN3(Ay,By,Cy) - epsi*p->DYP[js + marge];
	ye = MAX3(Ay,By,Cy) + epsi*p->DYP[je + marge];

	
	is = p->posf_i(xs+p->originx);
	ie = p->posf_i(xe+p->originx);
	
	js = p->posf_j(ys+p->originy);
	je = p->posf_j(ye+p->originy);
	
	is = MAX(is,0);
	ie = MIN(ie,p->knox);
	
	js = MAX(js,0);
	je = MIN(je,p->knoy);
	
		for(i=is;i<ie;i++)
		for(j=js;j<je;j++)
		{
		Px = p->XP[IP]+psi-p->originx;
		Py = p->YP[JP]+psi-p->originy;
		Pz = p->global_zmin-10.0*p->DXM ;
		
		Qx = p->XP[IP]+psi-p->originx;
		Qy = p->YP[JP]+psi-p->originy;
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


			if((u>0.0 && v>0.0 && w>0.0) || (u<0.0 && v<0.0 && w<0.0))
			{
			denom = 1.0/(u+v+w);
			u *= denom;
			v *= denom;
			w *= denom;
			
			Rz = u*Az + v*Bz + w*Cz;
			
			Rz+=p->originz;
            
            k = p->posf_k(Rz);
			
            int distcheck=1;
  
            
            if(Rz<p->ZP[KP])
            if(k>=0 && k<p->knoz)
            if(a->fb(i,j,k)<0 && a->fb(i,j,k-1)<0)
            distcheck=0;
            
            if(Rz>=p->ZP[KP])
            if(k>=0 && k<p->knoz)
            if(a->fb(i,j,k)<0 && a->fb(i,j,k+1)<0)
            distcheck=0;

            if(distcheck==1)
			for(k=0;k<p->knoz;++k)
            {
            if(a->fb(i,j,k)<0.0)
			a->fb(i,j,k)=-MIN(fabs(Rz-p->ZP[KP]),fabs(a->fb(i,j,k)));
            
            if(a->fb(i,j,k)>=0.0)
			a->fb(i,j,k)=MIN(fabs(Rz-p->ZP[KP]),fabs(a->fb(i,j,k)));
            
            }
			}
		
		}
	}



}