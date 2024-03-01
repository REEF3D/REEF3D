/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

void sixdof_obj::update_fbvel(lexer *p)         
{
    // Determine floating body velocities
        // U
        if(p->X11_u==0)
        u_fb(0) = 0.0;
        
        if(p->X11_u>=1)
        u_fb(0) = p_(0)/Mass_fb;
        
        
        // V
        if(p->X11_v==0 || p->j_dir==0)
        u_fb(1) = 0.0;
        
        if(p->X11_u>=1 && p->j_dir==1)
        u_fb(1) = p_(1)/Mass_fb;
        
        
        // W
        if(p->X11_w==0)
        u_fb(2) = 0.0;
        
        if(p->X11_w>=1)
        u_fb(2) = p_(2)/Mass_fb;
        
        
        // rotation
        if(p->j_dir==0)
        {
        u_fb(3) = 0.0;
        u_fb(4) = omega_I(1);
        u_fb(5) = 0.0;
        }
        
        if(p->j_dir==1)
        {
        u_fb(3) = omega_I(0);
        u_fb(4) = omega_I(1);
        u_fb(5) = omega_I(2);
        }
        
        
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
}

void sixdof_obj::saveTimeStep(lexer *p, int iter)
{
    deltan3_ = deltan2_;
    deltan2_ = deltan1_;
    deltan1_ = delta_;

    dpn3_ = dpn2_;
    dpn2_ = dpn1_;
    dpn1_ = dp_;
    dcn3_ = dcn2_;
    dcn2_ = dcn1_;
    dcn1_ = dc_;
    dhn3_ = dhn2_;
    dhn2_ = dhn1_;
    dhn1_ = dh_;
    den3_ = den2_;
    den2_ = den1_;
    den1_ = de_;
    
    pn3_ = pn2_;
    pn2_ = pn1_;
    pn1_ = p_;
    cn3_ = cn2_;
    cn2_ = cn1_;
    cn1_ = c_;
    hn3_ = hn2_;
    hn2_ = hn1_;
    hn1_ = h_;
    en3_ = en2_;
    en2_ = en1_;   
    en1_ = e_;
    
    dtn3 = dtn2;
    dtn2 = dtn1;
    dtn1 = alpha[iter]*p->dt;   
}


void sixdof_obj::maxvel(lexer *p, ghostcell *pgc)
{
	p->ufbmax = p->ufbi;
    p->vfbmax = p->vfbi;
    p->wfbmax = p->wfbi;
    
    double uvel,vvel,wvel;

	ALOOP
	{
        uvel = p->ufbi + (p->pos_z() - p->zg)*p->qfbi - (p->pos_y() - p->yg)*p->rfbi;
        vvel = p->vfbi + (p->pos_x() - p->xg)*p->rfbi - (p->pos_z() - p->zg)*p->pfbi;
        wvel = p->wfbi + (p->pos_y() - p->yg)*p->pfbi - (p->pos_x() - p->xg)*p->qfbi;
        
        p->ufbmax = MAX(p->ufbmax, fabs(uvel));
        p->vfbmax = MAX(p->vfbmax, fabs(vvel));
        p->wfbmax = MAX(p->wfbmax, fabs(wvel));
	}
}

