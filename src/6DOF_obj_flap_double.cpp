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

void sixdof_obj::flap_double(lexer *p, ghostcell *pgc, int id, double t)
{
    double z1,z2; 
    double xe2,xe3;
    
    xs = p->X172_xs;
    xe = p->X172_xe;
	
    ys = p->X172_ys;
    ye = p->X172_ye;

    zs = p->X172_zs;
    ze = p->X172_ze;   

    z1 = p->X172_z1;
    z2 = p->X172_z2;  
    
    // displacement of the upper hinge (column 1) and of the top relative to it (column 2), and
    // their rates, at the stage time t (B 118 on both)
    double d1,d2;
    wm_sample(p,t,1,d1,wm_v1);
    wm_sample(p,t,2,d2,wm_v2);
    
    // position
    xe2 = xe  + d1;
    xe3 = xe2 + d2;
    
    // velocity profile: evaluated at the cell height in the forcing (wm_u); uwm[k] is kept for
    // output only, from the lower and upper flap rates
    KLOOP
    {
    uwm[k] = 0.0;
    wwm[k] = 0.0;
    }

    //cout<<"xs: "<<xs<<" xe: "<<xe<<" ys: "<<ys<<" ye: "<<ye<<"zs: "<<zs<<" ze: "<<ze<<endl;
    //cout<<"xe2: "<<xe2<<" xe3: "<<xe3<<" z1: "<<z1<<" z2: "<<z2<<endl;
    
    //if(p->mpirank==0)
    //cout<<"xe2: "<<xe2<<" xe3: "<<xe3<<endl;
	
// Lower part 5
	// 3
        // Tri 1
    count=wm_tri0;
	
	tri_x[count][0] = xs;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe;

	tri_y[count][0] = ys;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = zs;
	tri_z[count][1] = zs;
	tri_z[count][2] = z1;
	++count;

        // Tri 2
	tri_x[count][0] = xe;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ys;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = z1;
	tri_z[count][1] = z1;
	tri_z[count][2] = zs;
	++count;
    
    
	// 2
        // Tri 1
	tri_x[count][0] = xs;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe;

	tri_y[count][0] = ye;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = zs;
	tri_z[count][1] = zs;
	tri_z[count][2] = z1;
	++count;

        // Tri 2
	tri_x[count][0] = xe;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ye;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = z1;
	tri_z[count][1] = z1;
	tri_z[count][2] = zs;
	++count;
    
    // 1 in
        // Tri 1
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

	tri_y[count][0] = ye;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = z1;
	tri_z[count][1] = z1;
	tri_z[count][2] = zs;
	++count;
    
    // 4 out
        // Tri 1
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
	tri_x[count][0] = xe;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe;

	tri_y[count][0] = ye;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = z1;
	tri_z[count][1] = z1;
	tri_z[count][2] = zs;
	++count;
    
    
    // 5 bed
        // Tri 1
	tri_x[count][0] = xs;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe;

	tri_y[count][0] = ys;
	tri_y[count][1] = ys;
	tri_y[count][2] = ye;

	tri_z[count][0] = zs;
	tri_z[count][1] = zs;
	tri_z[count][2] = zs;
	++count;

        // Tri 2
	tri_x[count][0] = xe;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ye;
	tri_y[count][1] = ye;
	tri_y[count][2] = ys;

	tri_z[count][0] = zs;
	tri_z[count][1] = zs;
	tri_z[count][2] = zs;
	++count;
    
    
// Middle part 4
	// 3
        // Tri 1
	tri_x[count][0] = xs;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe2;

	tri_y[count][0] = ys;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = z1;
	tri_z[count][1] = z1;
	tri_z[count][2] = z2;
	++count;

        // Tri 2
	tri_x[count][0] = xe2;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ys;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = z2;
	tri_z[count][1] = z2;
	tri_z[count][2] = z1;
	++count;
    
    
    
    // 2
        // Tri 1
	tri_x[count][0] = xs;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe2;

	tri_y[count][0] = ye;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = z1;
	tri_z[count][1] = z1;
	tri_z[count][2] = z2;
	++count;

        // Tri 2
	tri_x[count][0] = xe2;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ye;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = z2;
	tri_z[count][1] = z2;
	tri_z[count][2] = z1;
	++count;

	 // 1 in
        // Tri 1
	tri_x[count][0] = xs;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ys;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = z1;
	tri_z[count][1] = z1;
	tri_z[count][2] = z2;
	++count;

        // Tri 2
	tri_x[count][0] = xs;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ye;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = z2;
	tri_z[count][1] = z2;
	tri_z[count][2] = z1;
	++count;
    
    // 4 out
        // Tri 1
	tri_x[count][0] = xe;
	tri_x[count][1] = xe;
	tri_x[count][2] = xe2;

	tri_y[count][0] = ys;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = z1;
	tri_z[count][1] = z1;
	tri_z[count][2] = z2;
	++count;

        // Tri 2
	tri_x[count][0] = xe2;
	tri_x[count][1] = xe2;
	tri_x[count][2] = xe;

	tri_y[count][0] = ye;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = z2;
	tri_z[count][1] = z2;
	tri_z[count][2] = z1;
	++count;
    
    
// Top part 5
	// 3
        // Tri 1
	tri_x[count][0] = xs;
	tri_x[count][1] = xe2;
	tri_x[count][2] = xe3;

	tri_y[count][0] = ys;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = z2;
	tri_z[count][1] = z2;
	tri_z[count][2] = ze;
	++count;

        // Tri 2
	tri_x[count][0] = xe3;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ys;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = ze;
	tri_z[count][1] = ze;
	tri_z[count][2] = z2;
	++count;
    
    
    
    // 2
        // Tri 1
	tri_x[count][0] = xs;
	tri_x[count][1] = xe2;
	tri_x[count][2] = xe3;

	tri_y[count][0] = ye;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = z2;
	tri_z[count][1] = z2;
	tri_z[count][2] = ze;
	++count;

        // Tri 2
	tri_x[count][0] = xe3;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ye;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = ze;
	tri_z[count][1] = ze;
	tri_z[count][2] = z2;
	++count;

	 // 1 in
        // Tri 1
	tri_x[count][0] = xs;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ys;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = z2;
	tri_z[count][1] = z2;
	tri_z[count][2] = ze;
	++count;

        // Tri 2
	tri_x[count][0] = xs;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ye;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = ze;
	tri_z[count][1] = ze;
	tri_z[count][2] = z2;
	++count;
    
    // 4 out
        // Tri 1
	tri_x[count][0] = xe2;
	tri_x[count][1] = xe2;
	tri_x[count][2] = xe3;

	tri_y[count][0] = ys;
	tri_y[count][1] = ye;
	tri_y[count][2] = ye;

	tri_z[count][0] = z2;
	tri_z[count][1] = z2;
	tri_z[count][2] = ze;
	++count;

        // Tri 2
	tri_x[count][0] = xe3;
	tri_x[count][1] = xe3;
	tri_x[count][2] = xe2;

	tri_y[count][0] = ye;
	tri_y[count][1] = ys;
	tri_y[count][2] = ys;

	tri_z[count][0] = ze;
	tri_z[count][1] = ze;
	tri_z[count][2] = z2;
	++count;
    
    
    // 6 bed
        // Tri 1
	tri_x[count][0] = xs;
	tri_x[count][1] = xe3;
	tri_x[count][2] = xe3;

	tri_y[count][0] = ys;
	tri_y[count][1] = ys;
	tri_y[count][2] = ye;

	tri_z[count][0] = ze;
	tri_z[count][1] = ze;
	tri_z[count][2] = ze;
	++count;

        // Tri 2
	tri_x[count][0] = xe3;
	tri_x[count][1] = xs;
	tri_x[count][2] = xs;

	tri_y[count][0] = ye;
	tri_y[count][1] = ye;
	tri_y[count][2] = ys;

	tri_z[count][0] = ze;
	tri_z[count][1] = ze;
	tri_z[count][2] = ze;
	++count;
	
    
    
}

