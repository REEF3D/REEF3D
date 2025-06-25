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

#include"net_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void net_barDyn::coupling_dlm_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
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
    triangulation(p,pgc);
    
    
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
                p->ccipol4V(d->U,d->WL,d->bed,coordI(0),coordI(1),coordI(2)) - tri_vel[i][0],
                p->ccipol4V(d->V,d->WL,d->bed,coordI(0),coordI(1),coordI(2)) - tri_vel[i][1],
                p->ccipol4V(d->W,d->WL,d->bed,coordI(0),coordI(1),coordI(2)) - tri_vel[i][2];

            
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

