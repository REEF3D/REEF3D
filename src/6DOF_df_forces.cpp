/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"mooring.h"
#include"net.h"
#include"vrans.h"

void sixdof_df_object::updateForces(fdm *a)
{
    // Forces in inertial system

    Ffb_ << 0.0, 0.0, 0.0;
    Mfb_ << 0.0, 0.0, 0.0;

    if(p->X11_u == 1)
    Ffb_(0) = Xext + Xe - p->X26_Cu*p_(0)/Mass_fb; 
    if(p->X11_v == 1)
    Ffb_(1) = Yext + Ye - p->X26_Cv*p_(1)/Mass_fb;
    if(p->X11_w == 1)
    Ffb_(2) = Zext + Ze - p->X26_Cw*p_(2)/Mass_fb; 
    
    if(p->X11_p == 1)
    Mfb_(0) = Kext + Ke - p->X25_Cp*omega_I(0); 
    if(p->X11_q == 1)
    Mfb_(1) = Mext + Me - p->X25_Cq*omega_I(1);
    if(p->X11_r == 1)
    Mfb_(2) = Next + Ne - p->X25_Cr*omega_I(2);
}
