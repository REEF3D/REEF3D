/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"momentum.h"
#include"ghostcell.h"
#include<sys/stat.h>
  
void sixdof_obj::ini_fbvel(lexer *p, ghostcell *pgc)
{
    // external velocity
      Uext = Vext = Wext = Pext = Qext = Rext = 0.0; 

    // Rigid body motion ini    
    R_ << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    e_ << 0.0, 0.0, 0.0, 0.0;
    p_ << 0.0, 0.0, 0.0;
    c_ << 0.0, 0.0, 0.0;
    h_ << 0.0, 0.0, 0.0;
    
    dp_   << 0.0, 0.0, 0.0;
    dpn1_ << 0.0, 0.0, 0.0;
    dpn2_ << 0.0, 0.0, 0.0;
    dpn3_ << 0.0, 0.0, 0.0;
    
    dc_   << 0.0, 0.0, 0.0;
    dcn1_ << 0.0, 0.0, 0.0;
    dcn2_ << 0.0, 0.0, 0.0;
    dcn3_ << 0.0, 0.0, 0.0;
    
    dh_   << 0.0, 0.0, 0.0;
    dhn1_ << 0.0, 0.0, 0.0;
    dhn2_ << 0.0, 0.0, 0.0;
    dhn3_ << 0.0, 0.0, 0.0;
    
    de_   << 0.0, 0.0, 0.0, 0.0;
    den1_ << 0.0, 0.0, 0.0, 0.0;
    den2_ << 0.0, 0.0, 0.0, 0.0;
    den3_ << 0.0, 0.0, 0.0, 0.0;
    
    omega_B << 0.0, 0.0, 0.0;
    omega_I << 0.0, 0.0, 0.0;
    
    
    
	if (p->X102==1)
	{
		p_(0) += p->X102_u*Mass_fb;
		p_(1) += p->X102_v*Mass_fb;
		p_(2) += p->X102_w*Mass_fb;
	} 
    
	if (p->X103==1)
	{
		h_(0) = p->X103_p;
		h_(1) = p->X103_q;
		h_(2) = p->X103_r;
	}  
	
    // Velocities
	p->ufb = p->vfb = p->wfb = 0.0;
	p->pfb = p->qfb = p->rfb = 0.0; 
	p->ufbi = p->vfbi = p->wfbi = 0.0;
	p->pfbi = p->qfbi = p->rfbi = 0.0; 
    
	if(p->X210==1)
	{
        p->ufbi = p->X210_u*ramp_vel(p);
        p->vfbi = p->X210_v*ramp_vel(p);
        p->wfbi = p->X210_w*ramp_vel(p);
	}
	
	if(p->X211==1)
	{
        p->pfbi = p->X211_p;
        p->qfbi = p->X211_q;
        p->rfbi = p->X211_r;
	}

    // Positions
    phi = theta = psi = 0.0;
    
    // Forces
    Xext = Yext = Zext = Kext = Mext = Next = 0.0;
    Ffb_ << 0.0, 0.0, 0.0;
    Mfb_ << 0.0, 0.0, 0.0;
    
    // Printing
	printtime = 0.0;
    printtimenormal = 0.0;
    p->printcount_sixdof = 0;
}

void sixdof_obj::ini_parameter_stl(lexer *p, fdm *a, ghostcell *pgc)
{
    
    
}

