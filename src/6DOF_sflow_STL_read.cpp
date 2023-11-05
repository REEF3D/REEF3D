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
--------------------------------------------------------------------*/

#include"6DOF_sflow.h"
#include"lexer.h"
#include"ghostcell.h"

void sixdof_sflow::read_stl(lexer *p, ghostcell *pgc)
{
	string word;
	int count, vert_count;
    double trivec_x,trivec_y,trivec_z;
    
    trisum=1;
    p->Darray(tri_x,trisum,3);
	p->Darray(tri_y,trisum,3);
	p->Darray(tri_z,trisum,3);
    p->Darray(tri_x0,trisum,3);
	p->Darray(tri_y0,trisum,3);
	p->Darray(tri_z0,trisum,3); 
	
	// read and count number of triangles
    if(p->mpirank==0)
	cout<<"reading 6DOF STL "<<endl;
    
	ifstream stl("floating.stl", ios_base::in);
    
    //tstart[entity_count]=tricount;
	
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
    
    if(p->mpirank==0)
	cout<<"6DOF STL trisum: "<<count<<endl;
	++count;
	// create vecs
	p->Dresize(tri_x,trisum,count,3,3);
	p->Dresize(tri_y,trisum,count,3,3);
	p->Dresize(tri_z,trisum,count,3,3);
	p->Dresize(tri_x0,trisum,count,3,3);
	p->Dresize(tri_y0,trisum,count,3,3);
	p->Dresize(tri_z0,trisum,count,3,3);		
    
	
	trisum=count;
    
    if(p->mpirank==0)
	cout<<"re-reading 6DOF STL "<<endl;
	
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
    //if(p->mpirank==0)
	//cout<<"6DOF STL count "<<count<<endl;
        
	}
	stl.close();
    
    if(p->mpirank==0)
	cout<<"6DOF STL closed "<<endl;
	
	tricount = count + 1;
    //tend[entity_count] = tricount;
	
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
    

    for(int qr=0;qr<tricount;++qr)
    {
        rotation_stl_quaternion(p,p->X183_phi,p->X183_theta,p->X183_psi,tri_x[qr][0],tri_y[qr][0],tri_z[qr][0],p->X183_x,p->X183_y,p->X183_z);
        rotation_stl_quaternion(p,p->X183_phi,p->X183_theta,p->X183_psi,tri_x[qr][1],tri_y[qr][1],tri_z[qr][1],p->X183_x,p->X183_y,p->X183_z);
        rotation_stl_quaternion(p,p->X183_phi,p->X183_theta,p->X183_psi,tri_x[qr][2],tri_y[qr][2],tri_z[qr][2],p->X183_x,p->X183_y,p->X183_z);
    }

    // xg analysis
    STL_xmin=1.0e10;
    STL_xmax=-1.0e10;
    STL_ymin=1.0e10;
    STL_ymax=-1.0e10;
    
	for(n=0; n<tricount; ++n)
	for(q=0; q<3; ++q)
	{
		STL_xmin = MIN(STL_xmin,tri_x[n][q]);
        STL_xmax = MAX(STL_xmax,tri_x[n][q]);
        
        STL_ymin = MIN(STL_ymin,tri_y[n][q]);
        STL_ymax = MAX(STL_ymax,tri_y[n][q]);
	}
    
    p->xg = STL_xmin + 0.5*(STL_xmax-STL_xmin);
    p->yg = STL_ymin + 0.5*(STL_ymax-STL_ymin);
    
    //cout<<p->mpirank<<" STL_ymin: "<<STL_ymin<<" STL_ymax: "<<STL_ymax<<" yg: "<<p->yg<<endl;
}

void sixdof_sflow::rotation_stl_quaternion
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

void sixdof_sflow::rotation_stl(lexer *p,double &xvec,double &yvec,double &zvec)
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
