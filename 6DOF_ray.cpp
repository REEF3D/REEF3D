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

void sixdof_f::ray_cast(lexer *p, fdm *a, ghostcell *pgc)
{
	ALOOP
	{
	a->fb(i,j,k)=1.0;
	}
	
    for(rayiter=0; rayiter<2; ++rayiter)
    {

        for(int qn=0;qn<entity_sum;++qn)
        {
            if(rayiter==0)
            {
            ray_cast_io_x(p,a,pgc,tstart[qn],tend[qn]);
            ray_cast_io_ycorr(p,a,pgc,tstart[qn],tend[qn]);
            ray_cast_io_zcorr(p,a,pgc,tstart[qn],tend[qn]);
            }
        
            if(rayiter==1)
            {
            ray_cast_x(p,a,pgc,tstart[qn],tend[qn]);
            ray_cast_y(p,a,pgc,tstart[qn],tend[qn]);
            ray_cast_z(p,a,pgc,tstart[qn],tend[qn]);
            }
        

        }
    }
    
    
	
	ALOOP
	{
		if(a->fb(i,j,k)>10.0*p->DXM)
		a->fb(i,j,k)=10.0*p->DXM;
		
		if(a->fb(i,j,k)<-10.0*p->DXM)
		a->fb(i,j,k)=-10.0*p->DXM;
	}
	
	pgc->start4a(p,a->fb,50);
}





void sixdof_f::ray_cast_io(lexer *p, fdm *a, ghostcell *pgc, int ts, int te)
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
	int js,je,ks,ke;
	int ir;
	double u,v,w;
	double denom;
    double psi = 1.0e-6*p->DXM;
    
    MALOOP
	{
	cutl(i,j,k)=0;
	cutr(i,j,k)=0;
	}
	
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
	
	ys = MIN3(Ay,By,Cy);
	ye = MAX3(Ay,By,Cy);
	
	zs = MIN3(Az,Bz,Cz);
	ze = MAX3(Az,Bz,Cz);
    
	js = p->posf_j(ys+p->originy);
	je = p->posf_j(ye+p->originy);
	
	ks = p->posf_k(zs+p->originz);
	ke = p->posf_k(ze+p->originz);	
	
    ys = MIN3(Ay,By,Cy) - epsi*p->YP[js + marge-1];
	ye = MAX3(Ay,By,Cy) + epsi*p->YP[je + marge+1];
	
	zs = MIN3(Az,Bz,Cz) - epsi*p->ZP[ks + marge-1];
	ze = MAX3(Az,Bz,Cz) + epsi*p->ZP[ke + marge+1];
    
    
	js = p->posf_j(ys+p->originy);
	je = p->posf_j(ye+p->originy);
	
	ks = p->posf_k(zs+p->originz);
	ke = p->posf_k(ze+p->originz);	
	
	js = MAX(js,0);
	je = MIN(je,p->knoy);
	
	ks = MAX(ks,0);
	ke = MIN(ke,p->knoz);
	

		for(j=js;j<je;j++)
		for(k=ks;k<ke;k++)
		{
		Px = -10.0*p->DXM;
		Py = p->YP[JP]+psi-p->originy;
		Pz = p->ZP[KP]+psi-p->originz;
		
		Qx = 10.0*p->DXM;
		Qy = p->YP[JP]+psi-p->originy;
		Qz = p->ZP[KP]+psi-p->originz;
		
		
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
			
			Rx = u*Ax + v*Bx + w*Cx;
			Ry = u*Ay + v*By + w*Cy;
			Rz = u*Az + v*Bz + w*Cz;

			
				for(i=0;i<=p->knox;++i)
				{
				if(p->XP[IP]-p->originx<=Rx)
				cutr(i,j,k) += 1;
				
				if(p->XP[IP]-p->originx>Rx)
				cutl(i,j,k) += 1;
				}
			}
		}
	}
	
    
	ALOOP
	if((cutl(i,j,k)+1)%2==0  && (cutr(i,j,k)+1)%2==0)
	a->fb(i,j,k)=-fabs(a->fb(i,j,k));
	
	
	count=0;
	ALOOP
	if(a->fb(i,j,k)>0)
	++count;
}

