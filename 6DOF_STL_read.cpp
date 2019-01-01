/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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


void sixdof_f::read_stl(lexer *p, fdm *a, ghostcell *pgc)
{
	string word;
	int count, vert_count;
    double trivec_x,trivec_y,trivec_z;
	
	// read and count number of triangles
	
	ifstream stl("floating.stl", ios_base::in);
    
    tstart[entity_count]=tricount;
	
	count=tricount;
	while(!stl.eof())
	{
	
		stl>>word;
		
		if(word=="facet")
		++count;
	}
	
	stl.close();
	stl.clear();
	
	
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
	
	int chk=0;
	while(!stl.eof())
	{
		stl>>word;
		
		if(word=="ascii")
		chk=1;

	}
	
	if(chk==0)
	{
	cout<<"Please convert STL file to ASCII format!"<<endl<<endl;
	cout<<"See User's Guide for more information!"<<endl<<endl<<endl;
    pgc->final();
	exit(0);
	}
	
	stl.close();
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
	for(n=0; n<tricount; ++n)
	for(q=0; q<3; ++q)
	{
		tri_x[n][q]*=p->X181;
		tri_y[n][q]*=p->X181;
		tri_z[n][q]*=p->X181;
	}
	
	
	
    // rotate STL model
    
    p->X183_phi *= (PI/180.0);
    p->X183_theta *= (PI/180.0);
    p->X183_psi *= (PI/180.0);
    
	if (p->X13 == 0)
	{
		for(int qr=0;qr<tricount;++qr)
		{
			rotation_stl(p,tri_x[qr][0],tri_y[qr][0],tri_z[qr][0]);
			rotation_stl(p,tri_x[qr][1],tri_y[qr][1],tri_z[qr][1]);
			rotation_stl(p,tri_x[qr][2],tri_y[qr][2],tri_z[qr][2]);
		}
	}
	else if (p->X13 > 0)
	{
		for(int qr=0;qr<tricount;++qr)
		{
			rotation_stl_quaternion(p,tri_x[qr][0],tri_y[qr][0],tri_z[qr][0]);
			rotation_stl_quaternion(p,tri_x[qr][1],tri_y[qr][1],tri_z[qr][1]);
			rotation_stl_quaternion(p,tri_x[qr][2],tri_y[qr][2],tri_z[qr][2]);
		}
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
}


void sixdof_f::rotation_stl_quaternion(lexer *p,double &xvec,double &yvec,double &zvec)
{
	double a,b,c;

	// Distance to origin
    a = xvec - p->X183_x;
    b = yvec - p->X183_y;
    c = zvec - p->X183_z;

	// Rotation using matrices from Fossen
    xvec = a*(cos(psi)*cos(theta)) + b*(cos(theta)*sin(psi)) - c*sin(theta);
    yvec = a*(cos(psi)*sin(phi)*sin(theta)-cos(phi)*sin(psi)) + b*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta)) + c*(cos(theta)*sin(phi));
    zvec = a*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta)) + b*(cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi)) + c*(cos(phi)*cos(theta));
    
	// Moving back
    xvec += p->X183_x;
    yvec += p->X183_y;
    zvec += p->X183_z;
}	


void sixdof_f::rotation_stl(lexer *p,double &xvec,double &yvec,double &zvec)
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