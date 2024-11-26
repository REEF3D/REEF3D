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
#include"fdm_nhf.h"
#include"ghostcell.h"

void sixdof_obj::ray_cast_x(lexer *p, fdm_nhf *d, ghostcell *pgc, int ts, int te)
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
    double zs_x,zs_y,ze_x,ze_y;
    int checkin,ii;
    int i_zs,i_ze,j_zs,j_ze;
	double u,v,w;
	double denom;
	double psi = 1.0e-8*p->DXM;
    int margin=3;
    
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
	
    // zs
    if(Az<=Bz)
    {
    zs = Az;
    zs_x = Ax;
    zs_y = Ay;
    }
    
    if(Bz<Cz)
    {
    zs = Bz;
    zs_x = Bx;
    zs_y = By;
    }
    
    if(Cz<Bz)
    {
    zs = Cz;
    zs_x = Cx;
    zs_y = Cy;
    }
    
    // ze
    if(Az>=Bz)
    {
    ze = Az;
    ze_x = Ax;
    ze_y = Ay;
    }
    
    if(Bz>Cz)
    {
    ze = Bz;
    ze_x = Bx;
    ze_y = By;
    }
    
    if(Cz>Bz)
    {
    ze = Cz;
    ze_x = Cx;
    ze_y = Cy;
    }
    
    // i,j
    is = p->posc_i(xs);
	ie = p->posc_i(xe);
	
	js = p->posc_j(ys);
	je = p->posc_j(ye);
	
    
    // k
    i_zs = p->posc_i(zs_x);
	i_ze = p->posc_i(ze_x);
	
	j_zs = p->posc_j(zs_y);
	j_ze = p->posc_j(ze_y);
    
    i_zs = MAX(is,0);
	i_ze = MIN(ie,p->knox-1);
    
	j_zs = MAX(js,0);
	j_ze = MIN(je,p->knoy-1);
    
    
    ks = p->posc_sig(i_zs,j_zs,zs);
	ke = p->posc_sig(i_ze,j_ze,ze);
    
    i=i_zs;
    j=j_zs;
    k=ks;
    
    
    is-=margin;
    ie+=margin;

    js-=margin;
    je+=margin;

    ks-=margin;
    ke+=margin;
   

	is = MAX(is,0);
	ie = MIN(ie,p->knox);
    
	js = MAX(js,0);
	je = MIN(je,p->knoy);
	
	ks = MAX(ks,0);
	ke = MIN(ke,p->knoz);
    

    for(i=is;i<ie;i++)
    for(j=js;j<je;j++)
    for(k=ks;k<ke;k++)
    if((IO[IJK] != IO[Im1JK]) || (IO[IJK] != IO[Ip1JK]))
    {
        
        if(IO[IJK] != IO[Ip1JK])
        {
		Px = p->global_xmin-10.0*p->DXM;
		Py = p->YP[JP]-psi;
		Pz = p->ZSP[IJK]+psi;
		
		Qx = p->global_xmax+10.0*p->DXM;
		Qy = p->YP[JP]+psi;
		Qz = p->ZSP[Ip1JK]-psi;
        }
        
        if(IO[IJK] != IO[Im1JK])
        {
		Px = p->global_xmax+10.0*p->DXM;
		Py = p->YP[JP]-psi;
		Pz = p->ZSP[IJK]+psi;
		
		Qx = p->global_xmin-10.0*p->DXM;
		Qy = p->YP[JP]+psi;
		Qz = p->ZSP[Im1JK]-psi;
        }
		
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
		if(fabs(u)<=1.0e-20 && fabs(v)<=1.0e-20 && fabs(w)<=1.0e-20)
		check = 0;

			if(((u>1.0e-20 && v>1.0e-20 && w>1.0e-20) || (u<-1.0e-20 && v<-1.0e-20 && w<-1.0e-20)) && check==1)
			{
			denom = 1.0/(u+v+w);
            
			u *= denom;
			v *= denom;
			w *= denom;
			
			Rx = u*Ax + v*Bx + w*Cx;
            

            ii=i;

            for(i=0;i<p->knox;++i)
            d->FB[IJK]=MIN(fabs(Rx-p->XP[IP]),d->FB[IJK]);
            
            i=ii;
            }
		
		}
	}
    }

}
