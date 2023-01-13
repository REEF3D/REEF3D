/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void sixdof_gc::read_stl(lexer *p, fdm *a, ghostcell *pgc)
{
	string word;
	int count, vert_count;
    double trivec_x,trivec_y,trivec_z;
	
	// read and count number of triangles
	
	ifstream stl("floating.stl", ios_base::in);
    
    tstart[entity_count]=tricount;
	
	count=tricount;
    
    int chk=0;
	while(!stl.eof())
	{
		stl>>word;
		
		if(word=="facet")
		++count;

        if(word=="solid")
        chk=1;

	}
	
	stl.close();
	stl.clear();
	
    if(chk==0)
	{
	cout<<"Please convert STL file to ASCII format!"<<endl<<endl;
	cout<<"See User's Guide for more information!"<<endl<<endl<<endl;
    pgc->final();
	exit(0);
	}
	
	// create vecs
	p->Dresize(tri_x,tricount,count,3,3);
	p->Dresize(tri_y,tricount,count,3,3);
	p->Dresize(tri_z,tricount,count,3,3);
	p->Dresize(tri_xn,tricount,count,3,3);
	p->Dresize(tri_yn,tricount,count,3,3);
	p->Dresize(tri_zn,tricount,count,3,3);		
	
	tricount=count;
	
	// reopen and read triangles
	stl.open("floating.stl", ios_base::in);
	
	count=-1;
	while(!stl.eof())
	{
	
		stl>>word;
		
		if(word=="facet")
		{
		++count;
		vert_count=0;
		}
		
		if(word=="normal")
		stl>>trivec_x>>trivec_y>>trivec_z;
		
		if(word=="vertex")
		{
		stl>>tri_x[count][vert_count]>>tri_y[count][vert_count]>>tri_z[count][vert_count];
		++vert_count;
		}
	}
	stl.close();
	
	tricount = count + 1;
    tend[entity_count] = tricount;
	
	// scale STL model
	if (p->X181 == 1)
	for(n=0; n<tricount; ++n)
	for(int q=0; q<3; ++q)
	{
        tri_x[n][q] *= p->X181_x;
		tri_y[n][q] *= p->X181_y;
		tri_z[n][q] *= p->X181_z;
	}
	
    // change orgin
	if(p->X182==1)
	for(n=0; n<tricount; ++n)
	for(q=0; q<3; ++q)
	{
		tri_x[n][q]+=p->X182_x;
		tri_y[n][q]+=p->X182_y;
		tri_z[n][q]+=p->X182_z;
	}


    // rotate STL model
    
    p->X183_phi *= -(PI/180.0);
    p->X183_theta *= -(PI/180.0);
    p->X183_psi *= -(PI/180.0);
    
	if (p->X13 == 0)
	{
		for(int qr=0;qr<tricount;++qr)
		{
			rotation_stl(p,tri_x[qr][0],tri_y[qr][0],tri_z[qr][0]);
			rotation_stl(p,tri_x[qr][1],tri_y[qr][1],tri_z[qr][1]);
			rotation_stl(p,tri_x[qr][2],tri_y[qr][2],tri_z[qr][2]);
		}
	}
	else if (p->X13 == 1)
	{
		for(int qr=0;qr<tricount;++qr)
		{
            rotation_stl_quaternion(p,p->X183_phi,p->X183_theta,p->X183_psi,tri_x[qr][0],tri_y[qr][0],tri_z[qr][0],p->X183_x,p->X183_y,p->X183_z);
            rotation_stl_quaternion(p,p->X183_phi,p->X183_theta,p->X183_psi,tri_x[qr][1],tri_y[qr][1],tri_z[qr][1],p->X183_x,p->X183_y,p->X183_z);
            rotation_stl_quaternion(p,p->X183_phi,p->X183_theta,p->X183_psi,tri_x[qr][2],tri_y[qr][2],tri_z[qr][2],p->X183_x,p->X183_y,p->X183_z);
		}
	}
    
}

void sixdof_gc::rotation_stl_quaternion
(
    lexer *p,
    double phi_,double theta_,double psi_, 
    double &xvec,double &yvec,double &zvec, 
    const double& x0, const double& y0, const double& z0
)
{
	// Distance to origin
    double dx = xvec - x0;
    double dy = yvec - y0;
    double dz = zvec - z0;

	// Rotation using Goldstein page 603 (but there is wrong result)
    xvec = dx*(cos(psi_)*cos(theta_)) + dy*(cos(theta_)*sin(psi_)) - dz*sin(theta_);
    yvec = dx*(cos(psi_)*sin(phi_)*sin(theta_)-cos(phi_)*sin(psi_)) + dy*(cos(phi_)*cos(psi_)+sin(phi_)*sin(psi_)*sin(theta_)) + dz*(cos(theta_)*sin(phi_));
    zvec = dx*(sin(phi_)*sin(psi_)+cos(phi_)*cos(psi_)*sin(theta_)) + dy*(cos(phi_)*sin(psi_)*sin(theta_)-cos(psi_)*sin(phi_)) + dz*(cos(phi_)*cos(theta_));
    
	// Moving back
    xvec += x0;
    yvec += y0;
    zvec += z0;
}	

void sixdof_gc::rotation_stl(lexer *p,double &xvec,double &yvec,double &zvec)
{
	double a,b,c;

	// p->X183_phi
	a = xvec-p->X183_x;
	
	b = (yvec-p->X183_y)*cos(p->X183_phi) - (zvec-p->X183_z)*sin(p->X183_phi); 
	
	c = (yvec-p->X183_y)*sin(p->X183_phi) + (zvec-p->X183_z)*cos(p->X183_phi); 
	
	xvec=a+p->X183_x;
	yvec=b+p->X183_y;
	zvec=c+p->X183_z;	
	
	// p->X183_theta
	a = (xvec-p->X183_x)*cos(p->X183_theta) + (zvec-p->X183_z)*sin(p->X183_theta); 
	
	b = yvec-p->X183_y;
	
	c = -(xvec-p->X183_x)*sin(p->X183_theta) + (zvec-p->X183_z)*cos(p->X183_theta); 
	
	xvec=a+p->X183_x;
	yvec=b+p->X183_y;
	zvec=c+p->X183_z;	
	
	// p->X183_psi
	a = (xvec-p->X183_x)*cos(p->X183_psi) - (yvec-p->X183_y)*sin(p->X183_psi); 
	
	b = (xvec-p->X183_x)*sin(p->X183_psi) + (yvec-p->X183_y)*cos(p->X183_psi);
	
	c = zvec-p->X183_z;
	
	xvec=a+p->X183_x;
	yvec=b+p->X183_y;
	zvec=c+p->X183_z;	

}
