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
#include"fdm.h"
#include"ghostcell.h"

void sixdof_obj::read_stl(lexer *p, fdm *a, ghostcell *pgc)
{
	string word;
	int count, vert_count;
    double trivec_x,trivec_y,trivec_z;
	
	// read and count number of triangles
    ifstream stl;
    if (n6DOF == 0)
    {
	    stl.open("floating.stl", ios_base::in);
    }
    else
    {
        char str[1000];
        sprintf(str,"floating-%i.stl",n6DOF);
	    stl.open(str, ios_base::in);
    }
    
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
	p->Dresize(tri_x0,tricount,count,3,3);
	p->Dresize(tri_y0,tricount,count,3,3);
	p->Dresize(tri_z0,tricount,count,3,3);		
	
	tricount=count;
	
	// reopen and read triangles
    if (n6DOF == 0)
    {
	    stl.open("floating.stl", ios_base::in);
    }
    else
    {
        char str[1000];
        sprintf(str,"floating-%i.stl",n6DOF);
	    stl.open(str, ios_base::in);
    }
	
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
	for(int q=0; q<3; ++q)
	{
		tri_x[n][q] += p->X182_x;
		tri_y[n][q] += p->X182_y;
		tri_z[n][q] += p->X182_z;
	}
    
    // rotate STL model
    p->X183_phi *= -(PI/180.0);
    p->X183_theta *= -(PI/180.0);
    p->X183_psi *= -(PI/180.0);

    for(int qr=0;qr<tricount;++qr)
    {
        rotation_tri(p,p->X183_phi,p->X183_theta,p->X183_psi,tri_x[qr][0],tri_y[qr][0],tri_z[qr][0],p->X183_x,p->X183_y,p->X183_z);
        rotation_tri(p,p->X183_phi,p->X183_theta,p->X183_psi,tri_x[qr][1],tri_y[qr][1],tri_z[qr][1],p->X183_x,p->X183_y,p->X183_z);
        rotation_tri(p,p->X183_phi,p->X183_theta,p->X183_psi,tri_x[qr][2],tri_y[qr][2],tri_z[qr][2],p->X183_x,p->X183_y,p->X183_z);
    }
}


