/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"ghostcell.h"

void nhflow_forcing::jacketmember(lexer *p, ghostcell *pgc, int id)
{
    double U,ds,eta;
	double rmax;
	int snum;
	int vertice_mem, center1_num,center2_num;
    double alpha,beta,gamma;
    double dX,dY,dZ,ddX,ddY,ddZ;
    double length;
    double xvec1,yvec1,zvec1;
    double alpha1,beta1,gamma1;
    double x1,y1,z1;
    double x2,y2,z2;
    double a1,b1,c1;
    double a2,b2,c2;
    double xm1,ym1,zm1,r1,xm2,ym2,zm2,r2;
    
    xrot=0.0;
	yrot=0.0;
	zrot=0.0;
    
    
    xm1 = p->A585_xm1[id];
    ym1 = p->A585_ym1[id];
    zm1 = p->A585_zm1[id];
    r1  = p->A585_r1[id];
    
    xm2 = p->A585_xm2[id];
    ym2 = p->A585_ym2[id];
    zm2 = p->A585_zm2[id];
    r2  = p->A585_r2[id];

    dX = xm2-xm1;
    dY = ym2-ym1;
    dZ = zm2-zm1;

    length = sqrt(dX*dX + dY*dY + dZ*dZ);
    
    alpha=beta=gamma=0.0;
    phi=theta=psi=0.0;

    // alpha
    angle_calc(dX,dY,dZ,alpha,beta,gamma);

    a1=0.0;
    b1=0.0;
    c1=0.0;
    a2=0.0;
    b2=0.0;
    c2=0.0;
 
    double ee=1.0e-4;
    int count=0;
    do
    {
        x1=length;
        y1=0.0;
        z1=0.0;

      rotation(x1,y1,z1,a1,b1,c1);
      angle_calc(x1,y1,z1,a2,b2,c2);

      if(a2>alpha+ee || a2<alpha-ee)
      a1 = a1 - 0.1*(a2-alpha);

      if(b2>beta+ee || b2<beta-ee)
      b1 = b1 - 0.1*(b2-beta);

      if(c2>gamma+ee || c2<gamma-ee)
      c1 = c1 - 0.1*(c2-gamma);
    

      if(a2<=alpha+ee && a2>=alpha-ee)
      if(b2<=beta+ee  && b2>=beta-ee)
      if(c2<=gamma+ee && c2>=gamma-ee)
      break;

     ++count;
    }while(count<2500);
    
    if(p->mpirank==0)
    {
    cout<<"ID: "<<id<<endl;
    cout<<"iteration: "<<count<<endl;
    cout<<"alpha: "<<alpha*(180.0/PI)<<" beta: "<<beta*(180.0/PI)<<" gamma: "<<gamma*(180.0/PI)<<endl;
    cout<<"a1: "<<a1*(180.0/PI)<<" b1: "<<b1*(180.0/PI)<<" c1: "<<c1*(180.0/PI)<<endl;

    cout<<"dX: "<<dX<<" dY: "<<dY<<" dZ: "<<dZ<<endl;
    cout<<"x1: "<<x1<<" y1: "<<y1<<" z1: "<<z1<<endl;
    cout<<"length: "<<length<<endl<<endl;
    }

	rmax = MAX(p->A585_r1[id],p->A585_r2[id]);

	U = 2.0 * PI * rmax;
	ds = 0.75*(U*p->DXM);
	snum = int(U/ds);

    dX = xm2-xm1;
    dY = ym2-ym1;
    dZ = zm2-zm1;
    
    xm2=xm1+length;
    ym2=ym1;
    zm2=zm1;

// Vertices
	ds = (2.0*PI)/double(snum);
	eta=0.0;
    
    tstart[entity_count]=tricount;

	for(n=0;n<snum;++n)
	{
	//bottom circle
	tri_x[tricount][0] = xm1;
	tri_y[tricount][0] = ym1;
	tri_z[tricount][0] = zm1;

	tri_x[tricount][1] = xm1;
	tri_y[tricount][1] = ym1 + r1*sin(eta);
	tri_z[tricount][1] = zm1 + r1*cos(eta);

	tri_x[tricount][2] = xm1;
	tri_y[tricount][2] = ym1 + r1*sin(eta+ds);
	tri_z[tricount][2] = zm1 + r1*cos(eta+ds);
	++tricount;

	//top circle
	tri_x[tricount][0] = xm2;
	tri_y[tricount][0] = ym2;
	tri_z[tricount][0] = zm2;

	tri_x[tricount][1] = xm2;
	tri_y[tricount][1] = ym2 + r2*sin(eta);
	tri_z[tricount][1] = zm2 + r2*cos(eta);

	tri_x[tricount][2] = xm2;
	tri_y[tricount][2] = ym2 + r2*sin(eta+ds);
	tri_z[tricount][2] = zm2 + r2*cos(eta+ds);
	++tricount;

	//side
	// 1st triangle
	tri_x[tricount][0] = xm1;
	tri_y[tricount][0] = ym1 + r1*sin(eta);
	tri_z[tricount][0] = zm1 + r1*cos(eta);

	tri_x[tricount][1] = xm2;
	tri_y[tricount][1] = ym2 + r2*sin(eta+ds);
	tri_z[tricount][1] = zm2 + r2*cos(eta+ds);

	tri_x[tricount][2] = xm1;
	tri_y[tricount][2] = ym1 + r1*sin(eta+ds);
	tri_z[tricount][2] = zm1 + r1*cos(eta+ds);
	++tricount;

	// 2nd triangle
	tri_x[tricount][0] = xm1;
	tri_y[tricount][0] = ym1 + r1*sin(eta);
	tri_z[tricount][0] = zm1 + r1*cos(eta);

	tri_x[tricount][1] = xm2;
	tri_y[tricount][1] = ym2 + r2*sin(eta+ds);
	tri_z[tricount][1] = zm2 + r2*cos(eta+ds);

	tri_x[tricount][2] = xm2;
	tri_y[tricount][2] = ym2 + r2*sin(eta);
	tri_z[tricount][2] = zm2 + r2*cos(eta);
	++tricount;

	eta+=ds;
	}
    tend[entity_count]=tricount;


    xrot=xm1;
	yrot=ym1;
	zrot=zm1;

    psi=c1;
    theta=b1;
    phi=a1;


    rotate_triangle(p,tstart[entity_count],tend[entity_count]);


}
