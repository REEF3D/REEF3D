/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_forcing::ray_cast_direct(lexer *p, fdm_nhf *d, ghostcell *pgc, int ts, int te)
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
    double xc,yc,zc;
	int is,ie,js,je,ks,ke;
	int ir;
	double u,v,w;
	double denom;	
	int checkin;
	double psi = 1.0e-8*p->DXM;
    int margin = 5;
    double dist;

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
    && Az>=zmin && Az<=zmax)
    checkin=1;
    
    if(Bx>=p->global_xmin && Bx<=p->global_xmax 
    && By>=p->global_ymin && By<=p->global_ymax
    && Bz>=zmin && Bz<=zmax)
    checkin=1;
    
    if(Cx>=p->global_xmin && Cx<=p->global_xmax 
    && Cy>=p->global_ymin && Cy<=p->global_ymax
    && Cz>=zmin && Cz<=zmax)
    checkin=1;
    
    checkin=1;
        
    if(checkin==1)
    {
    xs = MIN3(Ax,Bx,Cx);
	xe = MAX3(Ax,Bx,Cx);
    
	ys = MIN3(Ay,By,Cy);
	ye = MAX3(Ay,By,Cy);
	
	zs = MIN3(Az,Bz,Cz);
	ze = MAX3(Az,Bz,Cz);
    
    

	is = p->posc_i(xs)-margin;
	ie = p->posc_i(xe)+margin;
    
    js = p->posc_j(ys)-margin;
	je = p->posc_j(ye)+margin;
	
	ks = p->posc_sig(is+margin,js+margin,zs)-margin;
	ke = p->posc_sig(is+margin,js+margin,ze)+margin;	

   // ks=0;
   // ke=p->knoz;
    
	is = MAX(is,0);
	ie = MIN(ie,p->knox);
    
	js = MAX(js,0);
	je = MIN(je,p->knoy);
	
	ks = MAX(ks,0);
	ke = MIN(ke,p->knoz);		
        
        //cout<<"xs: "<<xs<<" xe: "<<xe<<" ys: "<<ys<<" ye: "<<ye<<" zs: "<<zs<<" ze: "<<ze<<endl;
        //cout<<"is: "<<is<<" ie: "<<ie<<" js: "<<js<<" je: "<<je<<" ks: "<<ks<<" ke: "<<ke<<endl;    
        
         for(i=is;i<ie;++i)
		for(j=js;j<je;++j) 
		for(k=ks;k<ke;++k)
		{
        xc = p->XP[IP];
        yc = p->YP[JP];
        zc = p->ZSP[IJK];
        
        dist = sqrt(pow(xc-Ax,2.0) + pow(yc-Ay,2.0) + pow(zc-Az,2.0));

        d->SOLID[IJK]=MIN(dist,fabs(d->SOLID[IJK]));
        
        dist = sqrt(pow(xc-Bx,2.0) + pow(yc-By,2.0) + pow(zc-Bz,2.0));
        
        d->SOLID[IJK]=MIN(dist,fabs(d->SOLID[IJK]));
        
        dist = sqrt(pow(xc-Cx,2.0) + pow(yc-Cy,2.0) + pow(zc-Cz,2.0));
        
        d->SOLID[IJK]=MIN(dist,fabs(d->SOLID[IJK]));
		}
	}
    }
	
}
