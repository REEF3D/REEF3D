/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"vrans.h"
   

void sixdof_sflow::initialize(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
}
    
void sixdof_sflow::ini(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // Initialise parameters
	ini_parameter(p,b,pgc);
    
    // Initialise folder structure
	print_ini(p,b,pgc);

    // Initialise object and distance field
    cylinder(p,b,pgc);

    ray_cast(p,pgc);
    time_preproc(p); 
	reini(p,pgc,fb);

SLICELOOP1
{
    b->test(i,j) = fb(i,j);
}
pgc->gcsl_start4(p,b->test,50);

    // Initialise position of bodies
    iniPosition_RBM(p,b,pgc);
    
	// Recalculate distances
	ray_cast(p,pgc);
	reini(p,pgc,fb);
    
    // Print initial body 
    print_stl(p,pgc);



/*
    // Ini pressure field ship
	SLICELOOP1
    {
        xpos = p->pos1_x();
        ypos = p->pos1_y();
       
        if (-Ls/2.0 <= xpos && xpos <= Ls/2.0 && -Bs/2.0 <= ypos && ypos <= Bs/2.0)
        {
            press_x(i,j) = -4.0*press0*cl/Ls*pow(xpos/Ls,3)*(1.0 - cb*pow(ypos/Bs,2))*exp(-as*ypos*ypos/(Bs*Bs));
        }
        else
        {
             press_x(i,j) = 0.0;
        }
    }
    
	SLICELOOP2
    {
        xpos = p->pos2_x();
        ypos = p->pos2_y();
       
        if (-Ls/2.0 <= xpos && xpos <= Ls/2.0 && -Bs/2.0 <= ypos && ypos <= Bs/2.0)
        {
            press_y(i,j) = -2.0*press0*as/Bs*(1.0 - cl*pow(xpos/Ls,4))*(cb/as + 1.0 - cb*pow(ypos/Bs,2))*ypos/Bs*exp(-as*ypos*ypos/(Bs*Bs));
        }
        else
        {
             press_y(i,j) = 0.0;
        }
    }
    */
}



void sixdof_sflow::ini_parameter(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // Prescribed motions
    Uext = Vext = Wext = Pext = Qext = Rext = 0.0; 
    
    if (p->X210 == 1)
    {
        Uext = p->X210_u;
        Vext = p->X210_v;
        Wext = p->X210_w;
    }
    if (p->X211 == 1)
    {
        Pext = p->X211_p;
        Qext = p->X211_q;
        Rext = p->X211_r;
    }
    if (p->X221==1)
    {
        //motion_vec(p,a,pgc);
        cout<<"not implemented yet"<<endl;
    }

    // Position
    p->ufbi=p->vfbi=p->wfbi=0.0;
	p->pfbi=p->qfbi=p->rfbi=0.0;
    p->xg=p->yg=p->zg=0.0;
    phi = theta = psi = 0.0;
    quatRotMat << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
	
    // Printing
    printtime = 0.0;
    p->printcount_sixdof = 0;

    n6DOF = 0;
}

void sixdof_sflow::iniPosition_RBM(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // Initial position
    if(p->X23==1)
    {
        p->xg = p->X23_x; 
        p->yg = p->X23_y; 
        p->zg = p->X23_z; 
    }
    else
    {
         cout<<"Please provide centre of floating body using X 23!"<<endl;
    }

    // Store initial position of triangles
	for(n=0; n<tricount; ++n)
	{
        for(int q=0; q<3; q++)
        {        
            tri_x0[n][q] = tri_x[n][q] - p->xg;
            tri_y0[n][q] = tri_y[n][q] - p->yg;
            tri_z0[n][q] = tri_z[n][q] - p->zg;
        }
    }
	
	// Initial rotation
	if (p->X101==1)
	{	
        phi = p->X101_phi*(PI/180.0);
        theta = p->X101_theta*(PI/180.0);
        psi = p->X101_psi*(PI/180.0);	
	
		for (n=0; n<tricount; ++n)
		{
			rotation_tri
				(p,-phi,-theta,-psi,tri_x[n][0],tri_y[n][0],tri_z[n][0],p->xg,p->yg,p->zg);
			rotation_tri
				(p,-phi,-theta,-psi,tri_x[n][1],tri_y[n][1],tri_z[n][1],p->xg,p->yg,p->zg);
			rotation_tri
				(p,-phi,-theta,-psi,tri_x[n][2],tri_y[n][2],tri_z[n][2],p->xg,p->yg,p->zg);
		}
	}
	
	// Initialise quaternions (Goldstein p. 604)
	e_(0) = 
		 cos(0.5*phi)*cos(0.5*theta)*cos(0.5*psi) 
		+ sin(0.5*phi)*sin(0.5*theta)*sin(0.5*psi);
	e_(1) = 
		 sin(0.5*phi)*cos(0.5*theta)*cos(0.5*psi) 
		- cos(0.5*phi)*sin(0.5*theta)*sin(0.5*psi);
	e_(2) = 
		 cos(0.5*phi)*sin(0.5*theta)*cos(0.5*psi) 
		+ sin(0.5*phi)*cos(0.5*theta)*sin(0.5*psi);
	e_(3) = 
		 cos(0.5*phi)*cos(0.5*theta)*sin(0.5*psi) 
		- sin(0.5*phi)*sin(0.5*theta)*cos(0.5*psi);   

    // Initialise rotation matrices
    quat_matrices(e_);
}


void sixdof_sflow::cylinder(lexer *p, fdm2D *b, ghostcell *pgc)
{
	double U,ds,angle;
	double xm,ym,zm,z1,z2,r;
	int snum, trisum;
	
	xm = p->xg;
	ym = p->yg;
    zm = p->phimean;
	z1 = zm - 1.0;
	z2 = zm + 1.0;
    r = 40.0;

    // Prepare fields
    U = 2.0 * PI * r;
	ds = p->DXM;
	snum = int(U/ds);
    trisum=5*(snum+1);
    p->Darray(tri_x,trisum,3);
	p->Darray(tri_y,trisum,3);
	p->Darray(tri_z,trisum,3);
    p->Darray(tri_x0,trisum,3);
	p->Darray(tri_y0,trisum,3);
	p->Darray(tri_z0,trisum,3);    	

    // Vertices	
	ds = (2.0*PI)/double(snum);
	angle=0.0;
    tricount = 0;

	for(int n=0; n<snum; ++n)
	{
        //bottom circle	
        tri_x[tricount][0] = xm;
        tri_y[tricount][0] = ym;
        tri_z[tricount][0] = z1;
        
        tri_x[tricount][1] = xm + r*cos(angle);
        tri_y[tricount][1] = ym + r*sin(angle);
        tri_z[tricount][1] = z1;
        
        tri_x[tricount][2] = xm + r*cos(angle+ds);
        tri_y[tricount][2] = ym + r*sin(angle+ds);
        tri_z[tricount][2] = z1;
        ++tricount;
            
        //top circle
        tri_x[tricount][0] = xm;
        tri_y[tricount][0] = ym;
        tri_z[tricount][0] = z2;
        
        tri_x[tricount][1] = xm + r*cos(angle);
        tri_y[tricount][1] = ym + r*sin(angle);
        tri_z[tricount][1] = z2;
        
        tri_x[tricount][2] = xm + r*cos(angle+ds);
        tri_y[tricount][2] = ym + r*sin(angle+ds);
        tri_z[tricount][2] = z2;
        ++tricount;
        
        //side		
        // 1st triangle
        tri_x[tricount][0] = xm + r*cos(angle);
        tri_y[tricount][0] = ym + r*sin(angle);
        tri_z[tricount][0] = z1;
        
        tri_x[tricount][1] = xm + r*cos(angle+ds);
        tri_y[tricount][1] = ym + r*sin(angle+ds);
        tri_z[tricount][1] = z2;
        
        tri_x[tricount][2] = xm + r*cos(angle+ds);
        tri_y[tricount][2] = ym + r*sin(angle+ds);
        tri_z[tricount][2] = z1;
        ++tricount;
        
        // 2nd triangle
        tri_x[tricount][0] = xm + r*cos(angle);
        tri_y[tricount][0] = ym + r*sin(angle);
        tri_z[tricount][0] = z1;
        
        tri_x[tricount][1] = xm + r*cos(angle+ds);
        tri_y[tricount][1] = ym + r*sin(angle+ds);
        tri_z[tricount][1] = z2;
        
        tri_x[tricount][2] = xm + r*cos(angle);
        tri_y[tricount][2] = ym + r*sin(angle);
        tri_z[tricount][2] = z2;
        ++tricount;
            
        angle+=ds;
	}
}
