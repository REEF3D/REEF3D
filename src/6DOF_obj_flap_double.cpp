/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"

void sixdof_obj::flap_double(lexer *p, ghostcell *pgc, int id)
{
    uwm1 = 0.0;
    
    // find correct time step    
    if((p->simtime>kinematics[timecount][0]))
    timecount_old=timecount;
    
	while(p->simtime>kinematics[timecount][0])
	++timecount;
    
    f0 = (p->simtime - kinematics[timecount_old][0])/(kinematics[timecount][0]-kinematics[timecount_old][0]);

    if(p->simtime>=ts && p->simtime<=te && timecount<ptnum-1 && timecount_old<ptnum)
    xwm1 = kinematics[timecount][1]*f0 + kinematics[timecount_old][1]*(1.0-f0); 
    
    if(p->simtime>=ts && p->simtime<=te && timecount<ptnum-1 && timecount_old<ptnum)
    uwm1 = (kinematics[timecount][1]-kinematics[timecount_old][1])/(kinematics[timecount][0]-kinematics[timecount_old][0]);
    
	xs = p->X172_xs + xwm1;
    xe = p->X172_xe + xwm1;
	
    ys = p->X172_ys;
    ye = p->X172_ye;

    zs = p->X172_zs;
    ze = p->X172_ze;    
	
	// Lower part 5
	// bottom
        // Tri 1
    count=0;
	
	tri_x[count][0] = xs;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe;

	tri_y[count][0] = ye;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = zs;
	tri_z[count][1] = zs;
	tri_z[count][2] = zs;
	++count;

        // Tri 2
	tri_x[count][0] = xs;
	tri_x[count][1] = xe;
	tri_x[count][2] = xs;

	tri_y[count][0] = ys;
	tri_y[count][1] = ye;
	tri_y[count][2] = ys;

	tri_z[count][0] = zs;
	tri_z[count][1] = zs;
	tri_z[count][2] = zs;
	++count;
    
    
    // Lower part 1
	// bottom
        // Tri 1
    count=0;
	
	tri_x[count][0] = xs;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ys;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = zs;
	tri_z[count][1] = zs;
	tri_z[count][2] = z1;
	++count;

        // Tri 2
	tri_x[count][0] = xs;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ys;
	tri_y[count][1] = ys;
	tri_y[count][2] = ye;

	tri_z[count][0] = zs;
	tri_z[count][1] = z1;
	tri_z[count][2] = z1;
	++count;
    
    
    // Lower part 4
	// bottom
        // Tri 1
    count=0;
	
	tri_x[count][0] = xe;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe;

	tri_y[count][0] = ys;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = zs;
	tri_z[count][1] = zs;
	tri_z[count][2] = z1;
	++count;

        // Tri 2
	tri_x[count][0] = xs;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe;

	tri_y[count][0] = ys;
	tri_y[count][1] = ys;
	tri_y[count][2] = ye;

	tri_z[count][0] = zs;
	tri_z[count][1] = z1;
	tri_z[count][2] = z1;
	++count;

	 
	
}

