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

void sixdof_obj::wm_sample(lexer *p, double t, int col, double &x, double &v)
{
    // linear interpolation of column col between the samples that bracket t; before the first and
    // after the last sample the paddle rests at the end value. timecount is moved in both
    // directions, so the stage times t^n + c dt (not monotone within a step) are fine, and the
    // interval is always [timecount-1, timecount], also when dt spans several samples.
    x = v = 0.0;
    
    if(ptnum<1)
    return;
    
    if(ptnum<2 || t<=kinematics[0][0])
    {
    x = p->B118*kinematics[0][col];
    return;
    }
    
    if(t>=kinematics[ptnum-1][0])
    {
    x = p->B118*kinematics[ptnum-1][col];
    return;
    }
    
    timecount = MIN(MAX(timecount,1),ptnum-1);
    
    while(timecount<ptnum-1 && kinematics[timecount][0]<t)
    ++timecount;
    
    while(timecount>1 && kinematics[timecount-1][0]>=t)
    --timecount;
    
    timecount_old = timecount-1;
    
    const double t0 = kinematics[timecount_old][0];
    const double t1 = kinematics[timecount][0];
    const double dt = t1-t0;
    
    if(dt<=1.0e-20)
    {
    x = p->B118*kinematics[timecount][col];
    return;
    }
    
    f0 = (t-t0)/dt;
    
    x = p->B118*(kinematics[timecount_old][col] + f0*(kinematics[timecount][col]-kinematics[timecount_old][col]));
    v = p->B118*(kinematics[timecount][col]-kinematics[timecount_old][col])/dt;
}

double sixdof_obj::wm_u(lexer *p, double x, double z) const
{
    // horizontal paddle velocity at (x,z)
    if(p->X170==1)
    return wm_v1;
    
    if(p->X172==1)
    {
        // double flap: lower flap hinged at z1, upper flap hinged at z2 on top of it, so above z2 the
        // velocity is the translation of the upper hinge (v1) plus the rotation of the upper flap
        const double z1 = p->X172_z1, z2 = p->X172_z2, ze = p->X172_ze;
        double u = 0.0;
        
        if(z>z1)
        u = wm_v1*MIN((z-z1)/(z2-z1),1.0);
        
        if(z>z2)
        u += wm_v2*MIN((z-z2)/(ze-z2),1.0);
        
        // ramp across the paddle thickness (back half at rest)
        const double x1 = p->X172_xs + 0.5*(p->X172_xe-p->X172_xs);
        const double x2 = p->X172_xs + 0.75*(p->X172_xe-p->X172_xs);
        double fac = (x-x1)/(x2-x1);
        fac = MIN(MAX(fac,0.0),1.0);
        
        return fac*u;
    }
    
    return 0.0;
}

void sixdof_obj::piston(lexer *p, ghostcell *pgc, int id, double t)
{
    // paddle displacement and velocity at the stage time t (B 118 on both)
    wm_sample(p,t,1,xwm1,wm_v1);
    wm_v2 = 0.0;
    
    KLOOP
    {
    uwm[k] = wm_v1;
    wwm[k] = 0.0;
    }
    
	xs = p->X170_xs + xwm1;
    xe = p->X170_xe + xwm1;
	
    ys = p->X170_ys;
    ye = p->X170_ye;

    zs = p->X170_zs;
    ze = p->X170_ze;    
	
	// Face 3
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

