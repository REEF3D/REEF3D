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
#include"momentum.h"
#include"ghostcell.h"
#include<sys/stat.h>
  
void sixdof_obj::ini_fbvel(lexer *p, ghostcell *pgc)
{
    // external velocity
      Uext = Vext = Wext = Pext = Qext = Rext = 0.0; 

    // Rigid body motion ini: state, derivatives, history, angles and loads
    rb.reset();
    
    for(int qn=0;qn<p->X102;++qn)
    {
            p_(0) += p->X102_u[qn]*Mass_fb;
            p_(1) += p->X102_v[qn]*Mass_fb;
            p_(2) += p->X102_w[qn]*Mass_fb;
    }
    
    // initial angular velocity (X 103): set in iniPosition_RBM, once the inertia tensor is known
    
	
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

    // Forces
    Xext = Yext = Zext = Kext = Mext = Next = 0.0;
    
    // Printing
	printtime = 0.0;
    printtimenormal = 0.0;
    p->printcount_sixdof = 0;
}

void sixdof_obj::ini_parameter_stl(lexer *p, fdm *a, ghostcell *pgc)
{
    
    
}

