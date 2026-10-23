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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_obj_cfd::update_position_3D(lexer *p, fdm *a, ghostcell *pgc, bool finalize)
{
    // Calculate new position
    rb.euler_angles();

    // Update STL mesh
    update_trimesh_3D(p,a,pgc,finalize);
    
    // Update angular velocities 
    rb.update_omega();
    
    // the pose and velocity: in the 6DOF log files (REEF3D_*_6DOF), no longer on the terminal

}

void sixdof_obj_cfd::update_trimesh_3D(lexer *p, fdm *a, ghostcell *pgc, bool finalize)
{
	// Update position of triangles 
    geom.transform(R_,c_);

    // Update floating level set function
	ray_cast(p,a,pgc);
	reini_RK2(p,a,pgc,a->fb);
    
    pgc->start4a(p,a->fb,50);   
}

