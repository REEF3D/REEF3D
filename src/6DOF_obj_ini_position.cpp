/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"momentum.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>

void sixdof_obj::iniPosition_RBM(lexer *p, ghostcell *pgc)
{
    // Store initial position of triangles (body frame, relative to the CoG)
    geom.store_body_frame(c_);
	
	// Initial rotation

	if (p->X101==1)
	{	
        phi = p->X101_phi*(PI/180.0);
        theta = p->X101_theta*(PI/180.0);
        psi = p->X101_psi*(PI/180.0);	
	
        geom.rotate(-phi,-theta,-psi,c_);

        // Rotate mooring end point
        if (p->X313==1) 
        {
            for (int line=0; line < p->mooring_count; line++)
            {
			    sixdof_geometry::rotation_tri(-phi,-theta,-psi,p->X311_xe[line],p->X311_ye[line],p->X311_ze[line],c_(0),c_(1),c_(2));
            }
        }
	}
	

	// Initialise quaternions from the Euler angles
    rb.quaternion_from_euler();
    
    // Initial angular velocity (X 103, body-fixed p, q, r in rad/s): angular momentum h = I omega
    if(p->X103==1)
    h_ = I_*Eigen::Vector3d(p->X103_p, p->X103_q, p->X103_r);
    
    // stage and history copies
    rb.init_history();
    
    // Initialise rotation matrices
    quat_matrices(p);
}

