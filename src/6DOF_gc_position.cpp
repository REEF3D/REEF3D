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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_gc::fb_position(lexer *p, fdm *a, ghostcell *pgc)
{
	double fbval;
	int check;
	int si,sj,sk;
	
	xg += dxg;
	yg += dyg;
	zg += dzg;
	phi += dphi;
	theta += dtheta;
	psi += dpsi;
	
	if(p->mpirank==0)
	cout<<"XG: "<<xg<<" YG: "<<yg<<" ZG: "<<zg<<" phi: "<<phi*(180.0/PI)<<" theta: "<<theta*(180.0/PI)<<" psi: "<<psi*(180.0/PI)<<endl;

	for(n=0; n<tricount; ++n)
	{
	rotation(tri_x[n][0],tri_y[n][0],tri_z[n][0],dphi,dtheta,dpsi);
	rotation(tri_x[n][1],tri_y[n][1],tri_z[n][1],dphi,dtheta,dpsi);
	rotation(tri_x[n][2],tri_y[n][2],tri_z[n][2],dphi,dtheta,dpsi);
	
	tri_x[n][0] += dxg;
	tri_y[n][0] += dyg;
	tri_z[n][0] += dzg;
	
	tri_x[n][1] += dxg;
	tri_y[n][1] += dyg;
	tri_z[n][1] += dzg;
	
	tri_x[n][2] += dxg;
	tri_y[n][2] += dyg;
	tri_z[n][2] += dzg;
	}
	
}

void sixdof_gc::fb_position_quaternion(lexer *p, fdm *a, ghostcell *pgc, const std::vector<double>& e)
{
	if(p->X210==1 || p->X211==1 || p->X221==1)
    {
		xg = p->xg + dxg;
		yg = p->yg + dyg;
		zg = p->zg + dzg;
	}

	// Calculate Euler angles from quaternion
	
	// around z-axis
	psi = atan2(2.0*(e[1]*e[2] + e[3]*e[0]), 1.0 - 2.0*(e[2]*e[2] + e[3]*e[3])); 
	
	// around new y-axis
	double arg = 2.0*(e[0]*e[2] - e[1]*e[3]);
	
	if (fabs(arg) >= 1.0)
	{
		theta = SIGN(arg)*PI/2.0;
	}
	else
	{
		theta = asin(arg);															
	}	
		
	// around new x-axis
	phi = atan2(2.0*(e[2]*e[3] + e[1]*e[0]), 1.0 - 2.0*(e[1]*e[1] + e[2]*e[2]));

	if(p->mpirank==0)
	cout<<"XG: "<<xg<<" YG: "<<yg<<" ZG: "<<zg<<" phi: "<<phi*(180.0/PI)<<" theta: "<<theta*(180.0/PI)<<" psi: "<<psi*(180.0/PI)<<endl;


	// Update position of triangles 
	
    std::vector<double> point(3,0.0);
	
	for(n=0; n<tricount; ++n)
	{
        for(int q=0; q<3; q++)
        {
			// tri_xn is vector between tri_x and xg
            point[0] = tri_xn[n][q];
            point[1] = tri_yn[n][q];
            point[2] = tri_zn[n][q];

            point = rotation_R(point);
        
            tri_x[n][q] = point[0] + xg;
            tri_y[n][q] = point[1] + yg;
            tri_z[n][q] = point[2] + zg;
			
			// 2D
			if(p->X11_v != 1 && p->X11_p != 1 && p->X11_r != 1) 
			{
				tri_y[n][q] = tri_yn[n][q] + yg;	
			}
        }
	}
 }
