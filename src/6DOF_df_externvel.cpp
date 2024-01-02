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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_df_object::prescribedMotion(lexer *p, fdm *a, ghostcell *pgc, Eigen::Vector3d& dp, Eigen::Vector3d& dc)
{
    
    if (p->X11_u == 2)
    {
        dp(0) = 0.0; 
        dc(0) = Uext*ramp_vel(p);
    }
    
    if (p->X11_v == 2)
    {
        dp(1) = 0.0; 
        dc(1) = Vext*ramp_vel(p);
    }

    if (p->X11_w == 2)
    {
        dp(2) = 0.0; 
        dc(2) = Wext*ramp_vel(p);
    }
    
    if(p->X11_p==2)
    {
        //Pext;
        cout<<"not implemented yet"<<endl;
    }
                
    if(p->X11_q==2)
    {
        //Qext;
        cout<<"not implemented yet"<<endl;
    }
    
    if(p->X11_r==2)
    {
        //Rext;
        cout<<"not implemented yet"<<endl;
    }
    
}

double sixdof_df_object::ramp_vel(lexer *p)
{
    double f=1.0;
    
    if(p->X205==1 && p->X206==1 && p->simtime>=p->X206_ts && p->simtime<p->X206_te)
    {
    f = (p->simtime-p->X206_ts)/(p->X206_te-p->X206_ts);
    }
    
    if(p->X205==2 && p->X206==1 && p->simtime>=p->X206_ts && p->simtime<p->X206_te)
    {
    f = (p->simtime-p->X206_ts)/(p->X206_te-p->X206_ts)-(1.0/PI)*sin(PI*(p->simtime-p->X206_ts)/(p->X206_te-p->X206_ts));
    }
    
    if(p->X206==1 && p->simtime<p->X206_ts)
    f=0.0;
    
    return f;
}




