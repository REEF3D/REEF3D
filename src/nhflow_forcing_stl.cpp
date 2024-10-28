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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"ghostcell.h"

void nhflow_forcing::read_stl(lexer *p, ghostcell *pgc)
{
	string word;
	int count, vert_count;
    double trivec_x,trivec_y,trivec_z;
	
	// read and count number of triangles
    ifstream stl;

    stl.open("solid.stl", ios_base::in);

    
    tstart[entity_count]=tricount;
	
	count=tricount;
	
    int chk=0;
	while(!stl.eof())
	{
		stl>>word;
		
		if(word=="facet")
		++count;

        if(word=="ascii" || word=="solid")
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
	p->Dresize(tri_x0,tricount,count,3,3);
	p->Dresize(tri_y0,tricount,count,3,3);
	p->Dresize(tri_z0,tricount,count,3,3);		
	
	tricount=count;
	
	// reopen and read triangles
    stl.open("solid.stl", ios_base::in);

	
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
	if (p->A591==1)
	for(n=0; n<tricount; ++n)
	for(int q=0; q<3; ++q)
	{
         tri_x[n][q] *= p->A591_x;
		tri_y[n][q] *= p->A591_y;
		tri_z[n][q] *= p->A591_z;
	}
    
    // change orgin
	if(p->A592==1)
	for(n=0; n<tricount; ++n)
	for(int q=0; q<3; ++q)
	{
		tri_x[n][q] += p->A592_x;
		tri_y[n][q] += p->A592_y;
		tri_z[n][q] += p->A592_z;
	}
    
    // rotate STL model
    p->A593_phi *= -(PI/180.0);
    p->A593_theta *= -(PI/180.0);
    p->A593_psi *= -(PI/180.0);

    for(int qr=0;qr<tricount;++qr)
    {
        rotation_tri(p,p->A593_phi,p->A593_theta,p->A593_psi,tri_x[qr][0],tri_y[qr][0],tri_z[qr][0],p->A593_x,p->A593_y,p->A593_z);
        rotation_tri(p,p->A593_phi,p->A593_theta,p->A593_psi,tri_x[qr][1],tri_y[qr][1],tri_z[qr][1],p->A593_x,p->A593_y,p->A593_z);
        rotation_tri(p,p->A593_phi,p->A593_theta,p->A593_psi,tri_x[qr][2],tri_y[qr][2],tri_z[qr][2],p->A593_x,p->A593_y,p->A593_z);
    }
}

void nhflow_forcing::rotation_tri(lexer *p,double phi_,double theta_,double psi_, 
                                    double &xvec,double &yvec,double &zvec, 
                                    const double& x0, const double& y0, const double& z0)
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

