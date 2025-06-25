/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2025 Tobias Martin

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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"net_sheet.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_nhf.h"
#include"ghostcell.h"	

net_sheet::net_sheet(int number, lexer *p):nNet(number){}

net_sheet::~net_sheet(){}


void net_sheet::initialize_cfd(lexer *p, fdm *a, ghostcell *pgc)
{    
    //- Initialise net model
    ini(p,pgc);
    
    //- Initialise printing
    printtime = 0.0;
    print(p);
}

void net_sheet::initialize_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{    
    //- Initialise net model
    ini(p,pgc);
    
    //- Initialise printing
    printtime = 0.0;
    print(p);
}

void net_sheet::start_cfd(lexer *p, fdm *a, ghostcell *pgc, double alpha, Eigen::Matrix3d quatRotMat)
{
    double starttime1 = pgc->timer();    
    dt_ = alpha*p->dt;

    //- Store old velocities
    for (int i = 0; i < nK; i++)
    {
        coupledFieldn[i][0] = coupledField[i][0];
        coupledFieldn[i][1] = coupledField[i][1];
        coupledFieldn[i][2] = coupledField[i][2];
    }

    //- Get velocities at knots
    updateField_cfd(p,a,pgc,0);
    updateField_cfd(p,a,pgc,1);	
    updateField_cfd(p,a,pgc,2);
    
    //- Get density at knots
    updateField_cfd(p,a,pgc,3);       
    
    //- Calculate velocities from rigid body motion
    for (int knotI = 0; knotI < nK; knotI++)
    {
        xdot_(knotI,0) = p->ufbi + (x_(knotI,2) - p->zg)*p->qfbi - (x_(knotI,1) - p->yg)*p->rfbi;
        xdot_(knotI,1) = p->vfbi + (x_(knotI,0) - p->xg)*p->rfbi - (x_(knotI,2) - p->zg)*p->pfbi;
        xdot_(knotI,2) = p->wfbi + (x_(knotI,1) - p->yg)*p->pfbi - (x_(knotI,0) - p->xg)*p->qfbi;
    }
    
    //- Calculate force vector
    forces_knot *= 0.0; 
    gravityForce(p);
    inertiaForce(p);
    dragForce(p);
    
    Fx = 0.0;
    Fy = 0.0;
    Fz = 0.0;
    for (int knotI = 0; knotI < nK; knotI++)
    {
        Fx += forces_knot(knotI, 0); 
        Fy += forces_knot(knotI, 1); 
        Fz += forces_knot(knotI, 2); 
    }

	// Update position of triangles
    Eigen::Vector3d point;
	for(int n = 0; n < nK; ++n)
	{
        for(int q = 0; q < 3; q++)
        {
            // (tri_x0 is initial vector between tri_x and xg)
            point << tri_x0[n][q], tri_y0[n][q], tri_z0[n][q];
					
            point = quatRotMat*point;
        
            tri_x[n][q] = point(0) + p->xg;
            tri_y[n][q] = point(1) + p->yg;
            tri_z[n][q] = point(2) + p->zg;
        }
        
        // (x0_ is initial vector between x_ and xg)
        point = quatRotMat*(x0_.row(n)).transpose();
    
        x_(n,0) = point(0) + p->xg;
        x_(n,1) = point(1) + p->yg;
        x_(n,2) = point(2) + p->zg;

	}

    //- Coupling to vrans model
    coupling_dlm_cfd(p,a,pgc);
	
    //- Build and save net
	print(p);	

    //- Print output
    double endtime1 = pgc->timer() - starttime1; 
    if (p->mpirank==0)
    {
        cout<<"Net time: "<<endtime1<<endl;    
    }
}

void net_sheet::start_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, Eigen::Matrix3d quatRotMat)
{
    double starttime1 = pgc->timer();    
    dt_ = alpha*p->dt;

    //- Store old velocities
    for (int i = 0; i < nK; i++)
    {
        coupledFieldn[i][0] = coupledField[i][0];
        coupledFieldn[i][1] = coupledField[i][1];
        coupledFieldn[i][2] = coupledField[i][2];
    }

    //- Get velocities at knots
    updateField_nhflow(p,d,pgc,0);
    updateField_nhflow(p,d,pgc,1);
    updateField_nhflow(p,d,pgc,2);
    
    //- Get density at knots
    updateField_nhflow(p,d,pgc,3);       
    
    //- Calculate velocities from rigid body motion
    for (int knotI = 0; knotI < nK; knotI++)
    {
        xdot_(knotI,0) = p->ufbi + (x_(knotI,2) - p->zg)*p->qfbi - (x_(knotI,1) - p->yg)*p->rfbi;
        xdot_(knotI,1) = p->vfbi + (x_(knotI,0) - p->xg)*p->rfbi - (x_(knotI,2) - p->zg)*p->pfbi;
        xdot_(knotI,2) = p->wfbi + (x_(knotI,1) - p->yg)*p->pfbi - (x_(knotI,0) - p->xg)*p->qfbi;
    }
    
    //- Calculate force vector
    forces_knot *= 0.0; 
    gravityForce(p);
    inertiaForce(p);
    dragForce(p);
    
    Fx = 0.0;
    Fy = 0.0;
    Fz = 0.0;
    for (int knotI = 0; knotI < nK; knotI++)
    {
        Fx += forces_knot(knotI, 0); 
        Fy += forces_knot(knotI, 1); 
        Fz += forces_knot(knotI, 2); 
    }

	// Update position of triangles
    Eigen::Vector3d point;
	for(int n = 0; n < nK; ++n)
	{
        for(int q = 0; q < 3; q++)
        {
            // (tri_x0 is initial vector between tri_x and xg)
            point << tri_x0[n][q], tri_y0[n][q], tri_z0[n][q];
					
            point = quatRotMat*point;
        
            tri_x[n][q] = point(0) + p->xg;
            tri_y[n][q] = point(1) + p->yg;
            tri_z[n][q] = point(2) + p->zg;
        }
        
        // (x0_ is initial vector between x_ and xg)
        point = quatRotMat*(x0_.row(n)).transpose();
    
        x_(n,0) = point(0) + p->xg;
        x_(n,1) = point(1) + p->yg;
        x_(n,2) = point(2) + p->zg;

	}

    //- Coupling to vrans model
    coupling_dlm_nhflow(p,d,pgc);
	
    //- Build and save net
	print(p);	

    //- Print output
    double endtime1 = pgc->timer() - starttime1; 
    if (p->mpirank==0)
    {
        cout<<"Net time: "<<endtime1<<endl;    
    }
}



