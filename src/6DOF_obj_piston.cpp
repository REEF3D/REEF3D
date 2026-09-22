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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"

void sixdof_obj::piston(lexer *p, ghostcell *pgc, int id)
{
    KLOOP
    {
    uwm[k] = 0.0;
    wwm[k] = 0.0;
    }
    
    // find correct time step    
    if((p->simtime>kinematics[timecount][0]))
    timecount_old=timecount;
    
	while(p->simtime>kinematics[timecount][0])
	++timecount;
    
    //cout<<"timecount: "<<timecount <<" kin[timecount][0]: "<<kinematics[timecount][0]<<" kin[timecount+1][0]: "<<kinematics[timecount+1][0]<<endl;
    
    if(timecount>0)
    {
    f0 = (p->simtime - kinematics[timecount_old][0])/(kinematics[timecount][0]-kinematics[timecount_old][0]);

    if(p->simtime>=ts && p->simtime<=te && timecount<ptnum-1 && timecount_old<ptnum)
    xwm1 = p->B118*kinematics[timecount][1]*f0 + kinematics[timecount_old][1]*(1.0-f0); 
    
    if(p->simtime>=ts && p->simtime<=te && timecount<ptnum-1 && timecount_old<ptnum)
    KLOOP
    uwm[k] = p->B118*(kinematics[timecount][1]-kinematics[timecount_old][1])/(kinematics[timecount][0]-kinematics[timecount_old][0]);
    
    //cout <<"uwm[0]: "<<uwm[0]<<" f0: "<<f0<<endl;
    }
    
	xs = p->X170_xs + xwm1;
    xe = p->X170_xe + xwm1;
	
    ys = p->X170_ys;
    ye = p->X170_ye;

    zs = p->X170_zs;
    ze = p->X170_ze;    
	
	// Face 3
	// Tri 1
	
    count=0;
	
	tri_x[count][0] = xs;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe;

	tri_y[count][0] = ys;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = zs;
	tri_z[count][1] = zs;
	tri_z[count][2] = ze;
	++count;

	// Tri 2
	tri_x[count][0] = xe;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ys;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = ze;
	tri_z[count][1] = ze;
	tri_z[count][2] = zs;
	++count;

	// Face 4
	// Tri 3
	tri_x[count][0] = xe;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe;

	tri_y[count][0] = ye;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = ze;
	tri_z[count][1] = ze;
	tri_z[count][2] = zs;
	++count;

	// Tri 4
	tri_x[count][0] = xe;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe;

	tri_y[count][0] = ye;
	tri_y[count][1] = ye;
	tri_y[count][2] = ys;

	tri_z[count][0] = zs;
	tri_z[count][1] = ze;
	tri_z[count][2] = zs;
	++count;

	// Face 1
	// Tri 5
	tri_x[count][0] = xs;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ye;
	tri_y[count][1] = ye;
	tri_y[count][2] = ys;

	tri_z[count][0] = ze;
	tri_z[count][1] = zs;
	tri_z[count][2] = zs;
	++count;

	// Tri 6
	tri_x[count][0] = xs;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ys;
	tri_y[count][1] = ye;
	tri_y[count][2] = ys;

	tri_z[count][0] = ze;
	tri_z[count][1] = ze;
	tri_z[count][2] = zs;
	++count;
	
	// Face 2
	// Tri 7
	tri_x[count][0] = xe;
	tri_x[count][1] = xe;
	tri_x[count][2] = xs;

	tri_y[count][0] = ye;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = ze;
	tri_z[count][1] = zs;
	tri_z[count][2] = zs;
	++count;

	// Tri 8
	tri_x[count][0] = xs;
	tri_x[count][1] = xe;
	tri_x[count][2] = xs;

	tri_y[count][0] = ye;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = ze;
	tri_z[count][1] = ze;
	tri_z[count][2] = zs;
	++count;

	// Face 5
	// Tri 9
	tri_x[count][0] = xe;
	tri_x[count][1] = xe;
	tri_x[count][2] = xs;

	tri_y[count][0] = ye;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = zs;
	tri_z[count][1] = zs;
	tri_z[count][2] = zs;
	++count;

	// Tri 10
	tri_x[count][0] = xs;
	tri_x[count][1] = xe;
	tri_x[count][2] = xs;

	tri_y[count][0] = ye;
	tri_y[count][1] = ye;
	tri_y[count][2] = ys;

	tri_z[count][0] = zs;
	tri_z[count][1] = zs;
	tri_z[count][2] = zs;
	++count;

	// Face 6
	// Tri 11
	tri_x[count][0] = xs;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe;

	tri_y[count][0] = ys;
	tri_y[count][1] = ys;
	tri_y[count][2] = ye;

	tri_z[count][0] = ze;
	tri_z[count][1] = ze;
	tri_z[count][2] = ze;
	++count;

	// Tri 12
	tri_x[count][0] = xs;
	tri_x[count][1] = xe;
	tri_x[count][2] = xs;

	tri_y[count][0] = ys;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = ze;
	tri_z[count][1] = ze;
	tri_z[count][2] = ze;
	++count;
	 
	
}

