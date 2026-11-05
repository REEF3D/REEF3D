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
#include"fdm.h"
#include"ghostcell.h"
#include"6DOF_motionext.h"

void sixdof_obj::get_trans(lexer *p, ghostcell *pgc)
{
    // dp = F, dc = p/m
    rb.derivatives_trans();
	
	// External motions
	pmotion->motionext_trans(p,pgc,rb.dp,rb.dc);
} 

void sixdof_obj::get_rot(lexer *p)
{
    // de, dh (updates the transformation matrices)
    rb.derivatives_rot();
	
	// External motions
    pmotion->motionext_rot(p,rb.dh,rb.h,rb.de,rb.G,rb.I);
} 

void sixdof_obj::quat_matrices(lexer *p)
{   
    rb.quat_matrices();
}
