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
Author: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_obj::update_position_3D(lexer *p, fdm *a, ghostcell *pgc, bool finalise)
{
    // Calculate new position
    update_Euler_angles(p,pgc,finalise);
    
    // Update STL mesh
    update_trimesh_3D(p,a,pgc,finalise);

    // Update angular velocities 
    omega_B = I_.inverse()*h_;
    omega_I = R_*omega_B;
      
    // Global body variables
    update_fbvel(p);    
    maxvel(p,pgc);
}

void sixdof_obj::update_Euler_angles(lexer *p, ghostcell *pgc, bool finalise)
{
	// Calculate Euler angles from quaternion
	
	// around z-axis
	psi = atan2(2.0*(e_(1)*e_(2) + e_(3)*e_(0)), 1.0 - 2.0*(e_(2)*e_(2) + e_(3)*e_(3))); 
	
	// around new y-axis
	double arg = 2.0*(e_(0)*e_(2) - e_(1)*e_(3));
	
	if (fabs(arg) >= 1.0)
	theta = SIGN(arg)*PI/2.0;
    
	else
	theta = asin(arg);														
	
		
	// around new x-axis
	phi = atan2(2.0*(e_(2)*e_(3) + e_(1)*e_(0)), 1.0 - 2.0*(e_(1)*e_(1) + e_(2)*e_(2)));


	if(p->mpirank==0 && finalise == true)
    {
        cout<<"XG: "<<c_(0)<<" YG: "<<c_(1)<<" ZG: "<<c_(2)<<" phi: "<<phi*(180.0/PI)<<" theta: "<<theta*(180.0/PI)<<" psi: "<<psi*(180.0/PI)<<endl;
        cout<<"Ue: "<<u_fb(0)<<" Ve: "<< u_fb(1)<<" We: "<< u_fb(2)<<" Pe: "<<omega_I(0)<<" Qe: "<<omega_I(1)<<" Re: "<<omega_I(2)<<endl;
    }
}

void sixdof_obj::update_trimesh_3D(lexer *p, fdm *a, ghostcell *pgc, bool finalise)
{
	// Update position of triangles 
	for(n=0; n<tricount; ++n)
	{
        for(int q=0; q<3; q++)
        {
            // Update coordinates of triangles 
            // (tri_x0 is vector between tri_x and xg)
  
            Eigen::Vector3d point(tri_x0[n][q], tri_y0[n][q], tri_z0[n][q]);
					
            point = R_*point;
        
            tri_x[n][q] = point(0) + c_(0);
            tri_y[n][q] = point(1) + c_(1);
            tri_z[n][q] = point(2) + c_(2);

			// 2D
			if(p->X11_v != 1 && p->X11_p != 1 && p->X11_r != 1) 
			{
				tri_y[n][q] = tri_y0[n][q] + c_(1);	
			}
        }
	}
	
    // Update floating level set function
	ray_cast(p,a,pgc);
	reini_RK2(p,a,pgc,a->fb);
    pgc->start4a(p,a->fb,50);   
}


