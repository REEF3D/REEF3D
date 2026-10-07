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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_obj::update_fbvel(lexer *p, ghostcell *pgc)         
{
    // Determine floating body velocities (free: momentum, prescribed: motionext)
    rb.velocity(u_fb);
    
    // Velocities
	p->ufbi = p_(0)/Mass_fb;
	p->vfbi = p_(1)/Mass_fb;
	p->wfbi = p_(2)/Mass_fb;
    
    p->pfbi = omega_I(0);
    p->qfbi = omega_I(1);
    p->rfbi = omega_I(2);


    // Position
	p->xg = c_(0);
	p->yg = c_(1);
	p->zg = c_(2);

	p->phi_fb = phi;
	p->theta_fb = theta;
	p->psi_fb = psi;
    
    maxvel(p,pgc);
}

void sixdof_obj::saveTimeStep(lexer *p, int iter)
{
    rb.save_history(alpha[iter]*p->dt);
}

void sixdof_obj::maxvel(lexer *p, ghostcell *pgc)
{
    // Maximum rigid-body velocity |u + omega x r| over the body surface (STL vertices, which all
    // ranks hold). Previously evaluated over every cell of the domain, so a rotating body got a
    // lever arm up to the domain size, which inflated ufbmax and cut the time step; that value
    // was also rank-local.
	p->ufbmax = fabs(p->ufbi);
    p->vfbmax = fabs(p->vfbi); 
    p->wfbmax = fabs(p->wfbi);
    
    double uvel,vvel,wvel,rx,ry,rz;
    
	for(int n=0; n<tricount; ++n)
    for(int q=0; q<3; ++q)
	{
        rx = tri_x[n][q] - p->xg;
        ry = tri_y[n][q] - p->yg;
        rz = tri_z[n][q] - p->zg;
        
        uvel = p->ufbi + rz*p->qfbi - ry*p->rfbi;
        vvel = p->vfbi + rx*p->rfbi - rz*p->pfbi;
        wvel = p->wfbi + ry*p->pfbi - rx*p->qfbi;
        
        p->ufbmax = MAX(p->ufbmax, fabs(uvel));
        p->vfbmax = MAX(p->vfbmax, fabs(vvel));
        p->wfbmax = MAX(p->wfbmax, fabs(wvel));
	}
}
