/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2024 Tobias Martin

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

#include"net_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_fsf.h"
#include"vrans.h"	


void net_barDyn::vransCoupling(lexer *p, fdm *a, ghostcell *pgc)
{
    
    //- Coupling to collar  
    collarPoints.resize(0);
    collarVel.resize(0);

    for (int knotI = 0; knotI < nbK; knotI++)
    {
        const Eigen::Vector3d& coordI = x_.row(knotI);
        
        if 
        (
            coordI(0) >= xstart[p->mpirank] && coordI(0) < xend[p->mpirank] &&
            coordI(1) >= ystart[p->mpirank] && coordI(1) < yend[p->mpirank] &&
            coordI(2) >= zstart[p->mpirank] && coordI(2) < zend[p->mpirank]
        )
        {
            collarPoints.push_back(coordI);
            collarVel.push_back(xdot_.row(knotI));
        }
    }


    //- Triangulate net
    triangulation(p,a,pgc);
    
    
    //- Save Lagrangian coordinates and forces
     
    double xc,yc,zc,x0,x1,x2,y0,y1,y2,z0,z1,z2,nx,ny,nz,mag,area; 
   
    lagrangePoints.resize(tend);
    lagrangeForces.resize(tend);
    
    int index = 0;
    int indexSide = 0;
    
    Vector3d side1, side2, normalVec;
    Vector3d v_rel, n_d, n_s, n_l;
 
    for (int i = 0; i < tend; i++)
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
        
        lagrangePoints[i] << xc, yc, zc;
        
        // Forces
        side1 << x1-x0, y1-y0, z1-z0;
        side2 << x2-x0, y2-y0, z2-z0;
        
        normalVec = side1.cross(side2);
        mag = normalVec.norm();

        area = 0.5*mag;
        
        normalVec /= mag;

        const Eigen::Vector3d& coordI = lagrangePoints[i];
        
        if 
        (
            coordI(0) >= xstart[p->mpirank] && coordI(0) < xend[p->mpirank] &&
            coordI(1) >= ystart[p->mpirank] && coordI(1) < yend[p->mpirank] &&
            coordI(2) >= zstart[p->mpirank] && coordI(2) < zend[p->mpirank]
        )
        {
            //- Calculate relative velocity at knot
            
            v_rel << 
                p->ccipol1_a(a->u,coordI(0),coordI(1),coordI(2)) - tri_vel[i][0],
                p->ccipol2_a(a->v,coordI(0),coordI(1),coordI(2)) - tri_vel[i][1],
                p->ccipol3_a(a->w,coordI(0),coordI(1),coordI(2)) - tri_vel[i][2];
            
            
            // Calculate normal velocity vector
            
            double v_mag = v_rel.norm();
            
            n_d = v_rel/(v_mag + 1e-10);
           

            // Correct direciton of normal vector of triangle
    
            n_s = SIGN(n_d.dot(normalVec))*normalVec;
            
            
            // Angle between velocity and normal vector
            
            double thetan = acos(n_d.dot(n_s));     


            // Normal vector of lift force
            
            n_l = (n_d.cross(n_s)).cross(n_d).normalized();
            
            
            //- Get drag and lift force coefficients
            
            double cd, cl;

            double v_mag_corr = v_mag;
            
            double error = 1.0;
            int nIt = 0;

            while (error > 1e-3 && nIt < 10)
            {
                error = v_mag_corr;    
                
                screenForceCoeff(p,cd,cl,v_mag_corr,thetan,p->X321_Sn[nNet]);
                
                // Froude momentum theory
                v_mag_corr = v_mag*cd/(2.0*(sqrt(1.0 + cd) - 1.0)); 

                error = fabs(v_mag_corr - error);
                
                nIt++;
            }
            
            if (std::isnan(v_mag_corr))
            {
                v_mag_corr = v_mag;
                screenForceCoeff(p,cd,cl,v_mag_corr,thetan,p->X321_Sn[nNet]);
            }            
            

            // Save directional forces at lagrangian points (w/o density since multiplied later)

            lagrangeForces[i] = 0.5*area*pow(v_mag_corr,2.0)*(cd*n_d + cl*n_l);

            //lagrangeForces[i] << tri_forces[i][0], tri_forces[i][1], tri_forces[i][2];
        }
        else
        {
            lagrangeForces[i] << 0.0, 0.0, 0.0;   
        }
    }    
    
    for (int pI = 0; pI < tend; pI++)
    {
        Eigen::Vector3d& forceI = lagrangeForces[pI];
        forceI << pgc->globalsum(forceI(0)), pgc->globalsum(forceI(1)), pgc->globalsum(forceI(2));
    }

}


void net_barDyn::triangulation(lexer *p, fdm *a, ghostcell *pgc)
{
    // Triangulate rectangular meshes into 2 triangles
    
    int index = 0.0;
    vector<double> vec3(3,0.0);
    vector<double> vec3f(3,0.0);

    tri_x.resize(meshID.size()*2, vec3);
    tri_y.resize(meshID.size()*2, vec3);
    tri_z.resize(meshID.size()*2, vec3);
    
    tri_vel.resize(meshID.size()*2, vec3);
    tri_forces.resize(meshID.size()*2, vec3);

    for (int i = 0; i < meshID.size(); i++)
    {
        // Tri 1
        tri_x[index][0] = x_(meshID[i][0], 0);
        tri_x[index][1] = x_(meshID[i][1], 0);
        tri_x[index][2] = x_(meshID[i][2], 0);

        tri_y[index][0] = x_(meshID[i][0], 1);
        tri_y[index][1] = x_(meshID[i][1], 1);
        tri_y[index][2] = x_(meshID[i][2], 1);

        tri_z[index][0] = x_(meshID[i][0], 2);
        tri_z[index][1] = x_(meshID[i][1], 2);
        tri_z[index][2] = x_(meshID[i][2], 2);
    
        tri_vel[index][0] = 1.0/3.0*(xdot_(meshID[i][0],0) + xdot_(meshID[i][1],0) + xdot_(meshID[i][2],0));
        tri_vel[index][1] = 1.0/3.0*(xdot_(meshID[i][0],1) + xdot_(meshID[i][1],1) + xdot_(meshID[i][2],1));
        tri_vel[index][2] = 1.0/3.0*(xdot_(meshID[i][0],2) + xdot_(meshID[i][1],2) + xdot_(meshID[i][2],2));
       
        tri_forces[index][0] = forces_knot(meshID[i][0], 0)/4.0 + forces_knot(meshID[i][1], 0)/8.0 + forces_knot(meshID[i][2], 0)/8.0;
        tri_forces[index][1] = forces_knot(meshID[i][0], 1)/4.0 + forces_knot(meshID[i][1], 1)/8.0 + forces_knot(meshID[i][2], 1)/8.0;
        tri_forces[index][2] = forces_knot(meshID[i][0], 2)/4.0 + forces_knot(meshID[i][1], 2)/8.0 + forces_knot(meshID[i][2], 2)/8.0;
        
        index++;

        // Tri 2
        tri_x[index][0] = x_(meshID[i][1], 0);
        tri_x[index][1] = x_(meshID[i][3], 0);
        tri_x[index][2] = x_(meshID[i][2], 0);      
        
        tri_y[index][0] = x_(meshID[i][1], 1);
        tri_y[index][1] = x_(meshID[i][3], 1);
        tri_y[index][2] = x_(meshID[i][2], 1); 

        tri_z[index][0] = x_(meshID[i][1], 2);
        tri_z[index][1] = x_(meshID[i][3], 2);
        tri_z[index][2] = x_(meshID[i][2], 2); 

        tri_vel[index][0] = 1.0/3.0*(xdot_(meshID[i][0],0) + xdot_(meshID[i][1],0) + xdot_(meshID[i][2],0));
        tri_vel[index][1] = 1.0/3.0*(xdot_(meshID[i][0],1) + xdot_(meshID[i][1],1) + xdot_(meshID[i][2],1));
        tri_vel[index][2] = 1.0/3.0*(xdot_(meshID[i][0],2) + xdot_(meshID[i][1],2) + xdot_(meshID[i][2],2));

        tri_forces[index][0] = forces_knot(meshID[i][3], 0)/4.0 + forces_knot(meshID[i][1], 0)/8.0 + forces_knot(meshID[i][2], 0)/8.0;
        tri_forces[index][1] = forces_knot(meshID[i][3], 1)/4.0 + forces_knot(meshID[i][1], 1)/8.0 + forces_knot(meshID[i][2], 1)/8.0;
        tri_forces[index][2] = forces_knot(meshID[i][3], 2)/4.0 + forces_knot(meshID[i][1], 2)/8.0 + forces_knot(meshID[i][2], 2)/8.0;
        
        index++;
    }

    tend = meshID.size()*2;
   

    // Refine according to cell size DXM

	double x01,x02,x12,y01,y02,y12,z01,z02,z12,mag;
	double at,bt,ct,st;
	double nx,ny,nz;	
  
    tri_x.reserve(tend*4);
    tri_y.reserve(tend*4);
    tri_z.reserve(tend*4); 

    tri_vel.reserve(tend*4);    


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
		
		if ((at + bt + ct)/3.0 > p->DXM)
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
            
            
            // Store old velocity
           
			vec3 = tri_vel[n];


            // Store forces, split in four new triangles

            vec3f[0] = tri_forces[n][0]/4.0;
            vec3f[1] = tri_forces[n][1]/4.0;
            vec3f[2] = tri_forces[n][2]/4.0;

			
			// Delete old triangles
		
			tri_x.erase(tri_x.begin() + n); 
			tri_y.erase(tri_y.begin() + n); 
			tri_z.erase(tri_z.begin() + n); 
            tri_vel.erase(tri_vel.begin() + n); 
            tri_forces.erase(tri_forces.begin() + n); 
			n--;
            

            // Create new triangles
            
			create_triangle(tri_x,tri_y,tri_z,x0,y0,z0,x01,y01,z01,x02,y02,z02,nx,ny,nz);
            tri_vel.push_back(vec3);
            tri_forces.push_back(vec3f);
			create_triangle(tri_x,tri_y,tri_z,x01,y01,z01,x12,y12,z12,x02,y02,z02,nx,ny,nz);
            tri_vel.push_back(vec3);
            tri_forces.push_back(vec3f);
			create_triangle(tri_x,tri_y,tri_z,x01,y01,z01,x1,y1,z1,x12,y12,z12,nx,ny,nz);
            tri_vel.push_back(vec3);
            tri_forces.push_back(vec3f);
			create_triangle(tri_x,tri_y,tri_z,x02,y02,z02,x12,y12,z12,x2,y2,z2,nx,ny,nz);
            tri_vel.push_back(vec3); 
            tri_forces.push_back(vec3f);
		}

		if (tri_x.size() > 5000) break;
        
        tend = tri_x.size(); 
	}

    // Save net as .stl
/*
    ofstream result;
    result.open("REEF3D_CFD_6DOF_Net/REEF3D_net.stl", ios::binary);

	result<<"solid"<<" "<<"ascii"<<endl;

    double x0, x1, x2, y0, y1, y2, z0, z1, z2;

	for(int n = 0; n < tend; ++n)
	{
		x0 = tri_x[n][0];
		x1 = tri_x[n][1];
		x2 = tri_x[n][2];
		
		y0 = tri_y[n][0];
		y1 = tri_y[n][1];
		y2 = tri_y[n][2];
		
		z0 = tri_z[n][0];
		z1 = tri_z[n][1];
		z2 = tri_z[n][2];          
        
        nx = (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0);
        ny = (x2 - x0) * (z1 - z0) - (x1 - x0) * (z2 - z0); 
        nz = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
            
        mag = sqrt(nx*nx + ny*ny + nz*nz);
        
        nx /= mag;
        ny /= mag;
        nz /= mag;
        
        result<<" facet normal "<<nx<<" "<<ny<<" "<<nz<<endl;
        result<<"  outer loop"<<endl;
        result<<"   vertex "<<tri_x[n][0]<<" "<<tri_y[n][0]<<" "<<tri_z[n][0]<<endl;
        result<<"   vertex "<<tri_x[n][1]<<" "<<tri_y[n][1]<<" "<<tri_z[n][1]<<endl;
        result<<"   vertex "<<tri_x[n][2]<<" "<<tri_y[n][2]<<" "<<tri_z[n][2]<<endl;
        result<<"  endloop"<<endl;
        result<<" endfacet"<<endl;
	}

	result<<"endsolid"<<endl;

	result.close();  
*/
}


void net_barDyn::create_triangle
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
		   SIGN(nx)!=SIGN(nx_old) 
		|| SIGN(ny)!=SIGN(ny_old) 
		|| SIGN(nz)!=SIGN(nz_old)
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
