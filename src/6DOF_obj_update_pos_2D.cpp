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

void sixdof_obj::update_position_2D(lexer *p, ghostcell *pgc, slice &fsglobal)
{
    // Calculate new position
    update_Euler_angles(p,pgc);
    
    if(p->mpirank==0)
    {
        cout<<"XG: "<<c_(0)<<" YG: "<<c_(1)<<" ZG: "<<c_(2)<<" phi: "<<phi*(180.0/PI)<<" theta: "<<theta*(180.0/PI)<<" psi: "<<psi*(180.0/PI)<<endl;
        cout<<"Ue: "<<u_fb(0)<<" Ve: "<< u_fb(1)<<" We: "<< u_fb(2)<<" Pe: "<<omega_I(0)<<" Qe: "<<omega_I(1)<<" Re: "<<omega_I(2)<<endl;
    }
    
    // Update STL mesh
    update_trimesh_2D(p,pgc);

    // Update angular velocities 
    omega_B = I_.inverse()*h_;
    omega_I = R_*omega_B;
    
    SLICELOOP4
    fsglobal(i,j) = fs(i,j);
    
    pgc->gcsl_start4(p,fsglobal,50);
}


void sixdof_obj::update_trimesh_2D(lexer *p, ghostcell *pgc)
{
    double starttime, endtime;
    starttime=pgc->timer();
    
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
            
            //cout<<tri_x0[n][q]<<" "<<tri_x[n][q]<<" "<<point(0)<<" "<<c_(0)<<endl;

			// 2D
			if(p->X11_v!=1 && p->X11_p!=1 && p->X11_r!=1) 
			{
				tri_y[n][q] = tri_y0[n][q] + c_(1);	
			}
        }
	}
    
    endtime=pgc->timer();
    
    //if(p->mpirank==0)
    //cout<<"6DOF update time 1: "<<endtime-starttime<<endl;
    
    starttime=pgc->timer();
	
    // Update floating level set function
	ray_cast_2D(p,pgc);
    
    endtime=pgc->timer();
    
    //if(p->mpirank==0)
    //cout<<"6DOF update time 2: "<<endtime-starttime<<endl;
    
    starttime=pgc->timer();
    
    
	reini_2D(p,pgc,fs);
    
    endtime=pgc->timer();
    
    //if(p->mpirank==0)
    //cout<<"6DOF update time 3: "<<endtime-starttime<<endl;
    
    pgc->gcsl_start4(p,fs,50);  
}

