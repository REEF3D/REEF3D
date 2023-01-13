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

void sixdof_gc::position_ini(lexer *p, fdm *a, ghostcell *pgc)
{	
	dxg = p->X100_x;
	dyg = p->X100_y;
	dzg = p->X100_z;
	dphi = p->X101_phi*(PI/180.0);
	dtheta = p->X101_theta*(PI/180.0);
	dpsi = p->X101_psi*(PI/180.0);
	
	if(p->X100==1 || p->X101==1)
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
	
	xg += dxg;
	yg += dyg;
	zg += dzg;
	phi += dphi;
	theta += dtheta;
	psi += dpsi;
	
	dxg=dyg=dzg=0.0;
	dphi=dtheta=dpsi=0.0;
}


void sixdof_gc::position_ini_quaternion(lexer *p, fdm *a, ghostcell *pgc)
{
	// Translation
	
	dxg = p->X100_x;
	dyg = p->X100_y;
	dzg = p->X100_z;
	
	if (p->X100==1)
	{
		for(n=0; n<tricount; ++n)
		{
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
	
	xg += dxg;
	yg += dyg;
	zg += dzg;
	
	for(n=0; n<tricount; ++n)
	{
        for(int q=0; q<3; q++)
        {        
            tri_xn[n][q] = tri_x[n][q] - xg;
            tri_yn[n][q] = tri_y[n][q] - yg;
            tri_zn[n][q] = tri_z[n][q] - zg;
        }
    }
	
    e_[7] = xg;
    e_[8] = yg;
    e_[9] = zg;	
	
	dxg=dyg=dzg=0.0;
	
	
	// Rotation
	
	dphi = p->X101_phi*(PI/180.0);
	dtheta = p->X101_theta*(PI/180.0);
	dpsi = p->X101_psi*(PI/180.0);	
	
	if (p->X101==1)
	{
		for (n=0; n<tricount; ++n)
		{
			rotation_quaternion
				(tri_x[n][0],tri_y[n][0],tri_z[n][0],dphi,dtheta,dpsi);
			rotation_quaternion
				(tri_x[n][1],tri_y[n][1],tri_z[n][1],dphi,dtheta,dpsi);
			rotation_quaternion
				(tri_x[n][2],tri_y[n][2],tri_z[n][2],dphi,dtheta,dpsi);
		}
        
        // Rotate mooring end point
        if (p->X313==1)
        {
            for (int line=0; line < p->mooring_count; line++)
            {
			    rotation_quaternion(p->X311_xe[line],p->X311_ye[line],p->X311_ze[line],dphi,dtheta,dpsi);
            }
        }
	}
	
	phi += dphi;
	theta += dtheta;
	psi += dpsi;
	
	dphi=dtheta=dpsi=0.0;

	// Goldstein page 604
	e_[0] = 
		 cos(0.5*phi)*cos(0.5*theta)*cos(0.5*psi) 
		+ sin(0.5*phi)*sin(0.5*theta)*sin(0.5*psi);
	e_[1] = 
		 sin(0.5*phi)*cos(0.5*theta)*cos(0.5*psi) 
		- cos(0.5*phi)*sin(0.5*theta)*sin(0.5*psi);
	e_[2] = 
		 cos(0.5*phi)*sin(0.5*theta)*cos(0.5*psi) 
		+ sin(0.5*phi)*cos(0.5*theta)*sin(0.5*psi);
	e_[3] = 
		 cos(0.5*phi)*cos(0.5*theta)*sin(0.5*psi) 
		- sin(0.5*phi)*sin(0.5*theta)*cos(0.5*psi);  

    // Initialise rotation matrix
	R_[0][0] = e_[0]*e_[0] + e_[1]*e_[1] - e_[2]*e_[2] - e_[3]*e_[3]; 
	R_[0][1] = 2.0*e_[1]*e_[2] - 2.0*e_[0]*e_[3]; 
	R_[0][2] = 2.0*e_[0]*e_[2] + 2.0*e_[1]*e_[3];
	R_[1][0] = 2.0*e_[0]*e_[3] + 2.0*e_[1]*e_[2]; 
	R_[1][1] = e_[0]*e_[0] - e_[1]*e_[1] + e_[2]*e_[2] - e_[3]*e_[3];
	R_[1][2] = 2.0*e_[2]*e_[3] - 2.0*e_[0]*e_[1];
	R_[2][0] = 2.0*e_[1]*e_[3] - 2.0*e_[0]*e_[2]; 
	R_[2][1] = 2.0*e_[0]*e_[1] + 2.0*e_[2]*e_[3]; 
	R_[2][2] = e_[0]*e_[0] - e_[1]*e_[1] - e_[2]*e_[2] + e_[3]*e_[3];    
}
