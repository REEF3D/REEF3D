/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2020 Tobias Martin

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

#include<sys/stat.h>
#include"net_sheet.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void net_sheet::ini(lexer *p, fdm *a, ghostcell *pgc)
{
    l_c = p->X321_lambda[nNet];  // Length of twine
    d_c = p->X321_d[nNet];       // Diameter of twine
    rho_c = p->X321_rho[nNet];   // Density of material
    knot_d = p->X321_dk[nNet];   // Knot diameter

/*---------------------------------------*/

// Havfarm 2
    if (p->mpirank == 0) cout<<"Havfarm 2 sheets"<<endl;
    // Cylinder
    double xm = p->X322_x0[nNet];
    double ym = p->X322_y0[nNet];
    double zb = p->X322_z0[nNet];
    double r = p->X322_D[nNet]/2.0; 
	double dalpha = 2.0*PI/p->X321_nd[nNet];
    double dz = p->X322_L[nNet]/p->X321_nl[nNet];
    
    int triD = int(2.0*PI/dalpha)*2;
    int triL = int(p->X322_L[nNet]/dz); 
    tend = triD*triL;    
    
    vector<double> vec3(3,0.0);
    tri_x.resize(tend, vec3);
    tri_y.resize(tend, vec3);
    tri_z.resize(tend, vec3);
    tri_x0.resize(tend, vec3);
    tri_y0.resize(tend, vec3);
    tri_z0.resize(tend, vec3);

    int tricount = 0;
    for(int n = 0; n < triD/2; ++n)
	{
        for (int q = 0; q < triL; ++q)
        {
            // 1st triangle
            tri_x[tricount][0] = xm + r*cos(n*dalpha);
            tri_y[tricount][0] = ym + r*sin(n*dalpha);
            tri_z[tricount][0] = zb + q*dz;
            
            tri_x[tricount][1] = xm + r*cos((n+1)*dalpha);
            tri_y[tricount][1] = ym + r*sin((n+1)*dalpha);
            tri_z[tricount][1] = zb + (q+1)*dz;
            
            tri_x[tricount][2] = xm + r*cos((n+1)*dalpha);
            tri_y[tricount][2] = ym + r*sin((n+1)*dalpha);
            tri_z[tricount][2] = zb + q*dz;

            ++tricount;
            
            // 2nd triangle
            tri_x[tricount][0] = xm + r*cos(n*dalpha);
            tri_y[tricount][0] = ym + r*sin(n*dalpha);
            tri_z[tricount][0] = zb + q*dz;
            
            tri_x[tricount][1] = xm + r*cos((n+1)*dalpha);
            tri_y[tricount][1] = ym + r*sin((n+1)*dalpha);
            tri_z[tricount][1] = zb + (q+1)*dz;
            
            tri_x[tricount][2] = xm + r*cos(n*dalpha);
            tri_y[tricount][2] = ym + r*sin(n*dalpha);
            tri_z[tricount][2] = zb + (q+1)*dz;
            
            ++tricount;
        }
    }


// OceanFarm 1
/*
    if (p->mpirank == 0) cout<<"OceanFarm 1 sheets"<<endl;
    // Top cylinder
    double xm = p->X322_x0[nNet];
    double ym = p->X322_y0[nNet];
    double zb = p->X322_z0[nNet];
    double r = p->X322_D[nNet]/2.0; 
	double dalpha = 2.0*PI/p->X321_nd[nNet];
    double dz = p->X322_L[nNet]/p->X321_nl[nNet];
    
    int triD = int(2.0*PI/dalpha)*2;
    int triL = int(p->X322_L[nNet]/dz); 
    tend = triD*triL;    
    
    // Bottom cone
    double h = 0.064;
    double zb2 = p->X322_z0[nNet] - h;
    double r2 = 0.04; 
    double dz2 = h/1;
    
    int triL2 = int(h/dz2);
    tend += triL2*triD;

    vector<double> vec3(3,0.0);
    tri_x.resize(tend, vec3);
    tri_y.resize(tend, vec3);
    tri_z.resize(tend, vec3);
    tri_x0.resize(tend, vec3);
    tri_y0.resize(tend, vec3);
    tri_z0.resize(tend, vec3);

    // Top cylinder
    int tricount = 0;
    for(int n = 0; n < triD/2; ++n)
	{
        for (int q = 0; q < triL; ++q)
        {
            // 1st triangle
            tri_x[tricount][0] = xm + r*cos(n*dalpha);
            tri_y[tricount][0] = ym + r*sin(n*dalpha);
            tri_z[tricount][0] = zb + q*dz;
            
            tri_x[tricount][1] = xm + r*cos((n+1)*dalpha);
            tri_y[tricount][1] = ym + r*sin((n+1)*dalpha);
            tri_z[tricount][1] = zb + (q+1)*dz;
            
            tri_x[tricount][2] = xm + r*cos((n+1)*dalpha);
            tri_y[tricount][2] = ym + r*sin((n+1)*dalpha);
            tri_z[tricount][2] = zb + q*dz;

            ++tricount;
            
            // 2nd triangle
            tri_x[tricount][0] = xm + r*cos(n*dalpha);
            tri_y[tricount][0] = ym + r*sin(n*dalpha);
            tri_z[tricount][0] = zb + q*dz;
            
            tri_x[tricount][1] = xm + r*cos((n+1)*dalpha);
            tri_y[tricount][1] = ym + r*sin((n+1)*dalpha);
            tri_z[tricount][1] = zb + (q+1)*dz;
            
            tri_x[tricount][2] = xm + r*cos(n*dalpha);
            tri_y[tricount][2] = ym + r*sin(n*dalpha);
            tri_z[tricount][2] = zb + (q+1)*dz;
            
            ++tricount;
        }
    }

    // Bottom cone
    for(int n = 0; n < triD/2; ++n)
	{
        for (int q = 0; q < triL2; ++q)
        {
            // 1st triangle
            tri_x[tricount][0] = xm + r2*cos(n*dalpha);
            tri_y[tricount][0] = ym + r2*sin(n*dalpha);
            tri_z[tricount][0] = zb2 + q*dz2;
            
            tri_x[tricount][1] = xm + r*cos((n+1)*dalpha);
            tri_y[tricount][1] = ym + r*sin((n+1)*dalpha);
            tri_z[tricount][1] = zb2 + (q+1)*dz2;
            
            tri_x[tricount][2] = xm + r2*cos((n+1)*dalpha);
            tri_y[tricount][2] = ym + r2*sin((n+1)*dalpha);
            tri_z[tricount][2] = zb2 + q*dz2;

            ++tricount;
          
            // 2nd triangle
            tri_x[tricount][0] = xm + r2*cos(n*dalpha);
            tri_y[tricount][0] = ym + r2*sin(n*dalpha);
            tri_z[tricount][0] = zb2 + q*dz2;
            
            tri_x[tricount][1] = xm + r*cos((n+1)*dalpha);
            tri_y[tricount][1] = ym + r*sin((n+1)*dalpha);
            tri_z[tricount][1] = zb2 + (q+1)*dz2;
            
            tri_x[tricount][2] = xm + r*cos(n*dalpha);
            tri_y[tricount][2] = ym + r*sin(n*dalpha);
            tri_z[tricount][2] = zb2 + (q+1)*dz2;
            
            ++tricount;
        }
    }

*/



/*-----------------------------------*/

    //- Triangulate net sheet
    triangulation(p,a,pgc);
    nK = tend;             
    
    //- Initialise fields    
    p->Darray(coupledField, nK, 4);		// fluid coupling matrix (velocity 1,2,3 + phi 4)
    p->Darray(coupledFieldn, nK, 4);
    x0_ = MatrixXd::Zero(nK,3); 
    x_ = MatrixXd::Zero(nK,3);
    xdot_ = MatrixXd::Zero(nK,3);
    forces_knot = MatrixXd::Zero(nK, 3);
    mass_knot = VectorXd::Zero(nK);
    weight_knot = VectorXd::Zero(nK);
    added_mass = VectorXd::Zero(nK);
    lagrangePoints.resize(nK);
    lagrangeForces.resize(nK);

    //- Define knots as centres of triangles
    double xc,yc,zc,x0,x1,x2,y0,y1,y2,z0,z1,z2;
    double A_panel, A_solid, l_solid;
    Vector3d side1, side2, normalVec;
    
    for (int i = 0; i < nK; i++)
    {
        // Coordinates
		x0 = tri_x[i][0];
		x1 = tri_x[i][1];
		x2 = tri_x[i][2];
		
		y0 = tri_y[i][0];
		y1 = tri_y[i][1];
		y2 = tri_y[i][2];
		
		z0 = tri_z[i][0];
		z1 = tri_z[i][1];
		z2 = tri_z[i][2];  
        
        xc = (x0 + x1 + x2)/3.0;
        yc = (y0 + y1 + y2)/3.0;
        zc = (z0 + z1 + z2)/3.0;

        x0_.row(i) << xc - p->xg, yc - p->yg, zc - p->zg;
        x_.row(i) << xc, yc, zc;
    
        // Mass lumping at each knot
        side1 << x1-x0, y1-y0, z1-z0;
        side2 << x2-x0, y2-y0, z2-z0;
        normalVec = side1.cross(side2);
        
        A_panel = 0.5*normalVec.norm();    
        A_solid = p->X321_Sn[nNet]*A_panel;
        l_solid = A_solid/d_c;    
           

        // Mass (in air)
        mass_knot(i) = rho_c*PI/4.0*d_c*d_c*l_solid; 
        
        // Weight (in water)
        weight_knot(i) = p->W1*PI/4.0*d_c*d_c*l_solid; 
        
        // Added mass assuming ca = 1.0
        added_mass(i) = p->W1*PI/4.0*d_c*d_c*1.0*l_solid; 
    }

    // Initialise probe points
    Eigen::Vector3d ppI;
    double dist, dist_new;
    
    probeKnot.resize(p->X324);

    for (int pp = 0; pp < p->X324; pp++)
    {
        ppI << p->X324_x[pp], p->X324_y[pp], p->X324_z[pp];
        
        dist = 1e10;

        for (int kI = 0; kI < nK; kI++)
        {
            dist_new = (x_.row(kI).transpose() - ppI).norm();

            if (dist_new < dist) 
            {
                dist = dist_new;
                probeKnot(pp) = kI;
            }
        }
        
        // Initialise print
        if(p->mpirank==0 && p->P14==1)
        {
            char str[1000];
            sprintf(str,"./REEF3D_6DOF_Net/REEF3D_6DOF_Net_%i_Point_Probe_%i.dat",nNet,pp+1);
            ofstream header_out;
            header_out.open(str);
            header_out<<"Knot point probe located near ("<<ppI.transpose()<<")"<<endl;
            header_out<<"time [s] \t x [m] \t y [m] \t z [m]"<<endl;
            header_out.close();
        }		
    }
    
    // Initialise force print
    if(p->mpirank==0 && p->P14==1)
    {
        char str[1000];
        sprintf(str,"./REEF3D_6DOF_Net/REEF3D_6DOF_Net_Forces_%i.dat",nNet);
        ofstream header_out;
        header_out.open(str);
        header_out<<"time [s] \t Fx [N] \t Fy [N] \t Fz [N]"<<endl;
        header_out.close();
    }		
    
    printtime = 0.0;


    // Initialise communication 
    p->Darray(xstart, p->mpi_size);
    p->Darray(xend, p->mpi_size);
    p->Darray(ystart, p->mpi_size);
    p->Darray(yend, p->mpi_size);
    p->Darray(zstart, p->mpi_size);
    p->Darray(zend, p->mpi_size);
    
    xstart[p->mpirank] = p->originx;
    ystart[p->mpirank] = p->originy;
    zstart[p->mpirank] = p->originz;
    xend[p->mpirank] = p->endx;
    yend[p->mpirank] = p->endy;
    zend[p->mpirank] = p->endz;
    
    for (int i = 0; i < p->mpi_size; i++)
    {
        MPI_Bcast(&xstart[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
        MPI_Bcast(&xend[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
        MPI_Bcast(&ystart[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
        MPI_Bcast(&yend[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
        MPI_Bcast(&zstart[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
        MPI_Bcast(&zend[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
    }
}



void net_sheet::triangulation(lexer *p, fdm *a, ghostcell *pgc)
{
    // Refine according to cell size DXM

	double x01,x02,x12,y01,y02,y12,z01,z02,z12,mag;
	double at,bt,ct,st;
	double nx,ny,nz;	
  
    vector<double> vec3(3,0.0);
	
    for (int n = 0; n < tend; n++)
	{
		double x0 = tri_x[n][0];
		double x1 = tri_x[n][1];
		double x2 = tri_x[n][2];
		
		double y0 = tri_y[n][0];
		double y1 = tri_y[n][1];
		double y2 = tri_y[n][2];
		
		double z0 = tri_z[n][0];
		double z1 = tri_z[n][1];
		double z2 = tri_z[n][2];  
           
		at = sqrt(pow(x1 - x0, 2.0) + pow(y1 - y0, 2.0) + pow(z1 - z0, 2.0));
		bt = sqrt(pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0) + pow(z1 - z2, 2.0));
		ct = sqrt(pow(x2 - x0, 2.0) + pow(y2 - y0, 2.0) + pow(z2 - z0, 2.0));   
		   

		// Check size of triangle and split into 4 triangles if too big
		
		if ((at + bt + ct)/3.0 > 1.5*p->DXM)
		{
			// Half points
			x01 = x0 + (x1 - x0)/2.0;
			y01 = y0 + (y1 - y0)/2.0;
			z01 = z0 + (z1 - z0)/2.0;

			x02 = x0 + (x2 - x0)/2.0;
			y02 = y0 + (y2 - y0)/2.0;
			z02 = z0 + (z2 - z0)/2.0;			

			x12 = x1 + (x2 - x1)/2.0;
			y12 = y1 + (y2 - y1)/2.0;
			z12 = z1 + (z2 - z1)/2.0;
			
			
			// Old normal vector    
			nx = (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0);
			ny = (x2 - x0) * (z1 - z0) - (x1 - x0) * (z2 - z0); 
			nz = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
           

			// Delete old triangles
			tri_x.erase(tri_x.begin() + n); 
			tri_y.erase(tri_y.begin() + n); 
			tri_z.erase(tri_z.begin() + n); 
			n--;
            

            // Create new triangles
			create_triangle(tri_x,tri_y,tri_z,x0,y0,z0,x01,y01,z01,x02,y02,z02,nx,ny,nz);
			create_triangle(tri_x,tri_y,tri_z,x01,y01,z01,x12,y12,z12,x02,y02,z02,nx,ny,nz);
			create_triangle(tri_x,tri_y,tri_z,x01,y01,z01,x1,y1,z1,x12,y12,z12,nx,ny,nz);
			create_triangle(tri_x,tri_y,tri_z,x02,y02,z02,x12,y12,z12,x2,y2,z2,nx,ny,nz);
		}

		if (tri_x.size() > 100000) break;
        
        tend = tri_x.size(); 
	}


    // Store initial position of triangles
    tri_x0.resize(tend, vec3);
    tri_y0.resize(tend, vec3);
    tri_z0.resize(tend, vec3);
	
    for(int n = 0;  n < tend; ++n)
	{
        for(int q = 0; q < 3; q++)
        {      
            tri_x0[n][q] = tri_x[n][q] - p->xg;
            tri_y0[n][q] = tri_y[n][q] - p->yg;
            tri_z0[n][q] = tri_z[n][q] - p->zg;
        }
    }

    if (p->mpirank == 0)
    {
        cout<<"Number of net triangles: "<<tend<<endl;
    }
}


void net_sheet::create_triangle
(
    MatrixVd& tri_x_, MatrixVd& tri_y_, MatrixVd& tri_z_, 
	const double& x0, const double& y0, const double& z0,
	const double& x1, const double& y1, const double& z1,
	const double& x2, const double& y2, const double& z2,
	const double& nx_old, const double& ny_old, const double& nz_old
)
{
	double nx,ny,nz,temp;

	vector<double> tri_x_new(3,0.0);
	vector<double> tri_y_new(3,0.0);
	vector<double> tri_z_new(3,0.0); 

	// Calculate new normal vector
	
	nx = (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0);
	ny = (x2 - x0) * (z1 - z0) - (x1 - x0) * (z2 - z0); 
	nz = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);		

	nx = nx > 1.0e-5 ? nx : nx_old;
	ny = ny > 1.0e-5 ? ny : ny_old;
	nz = nz > 1.0e-5 ? nz : nz_old;	
	
	
	// Arrange triangle such that normal vector points outward
	
	if 
	(
		   SIGN(nx) != SIGN(nx_old) 
		|| SIGN(ny) != SIGN(ny_old) 
		|| SIGN(nz) != SIGN(nz_old)
	)
	{
		tri_x_new[0] = x2;
		tri_x_new[1] = x1;
		tri_x_new[2] = x0;

		tri_y_new[0] = y2;
		tri_y_new[1] = y1;
		tri_y_new[2] = y0;

		tri_z_new[0] = z2;
		tri_z_new[1] = z1;
		tri_z_new[2] = z0;				
	}
	else
	{	
		tri_x_new[0] = x0;
		tri_x_new[1] = x1;
		tri_x_new[2] = x2;

		tri_y_new[0] = y0;
		tri_y_new[1] = y1;
		tri_y_new[2] = y2;

		tri_z_new[0] = z0;
		tri_z_new[1] = z1;
		tri_z_new[2] = z2;	
	}
	
	
	// Add triangle to list
	
	tri_x_.push_back(tri_x_new);
	tri_y_.push_back(tri_y_new);
	tri_z_.push_back(tri_z_new);
}
