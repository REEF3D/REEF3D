/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"6DOF_sflow.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans.h"

double sixdof_sflow::ramp_vel(lexer *p)
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

double sixdof_sflow::ramp_draft(lexer *p)
{
    double f=1.0;
    
    if(p->X205==1 && p->X207==1 && p->simtime>=p->X207_ts && p->simtime<p->X207_te)
    {
    f = p->simtime/(p->X207_te-p->X207_ts);
    }
    
    if(p->X205==2 && p->X207==1 && p->simtime>=p->X207_ts && p->simtime<p->X207_te)
    {
    f = p->simtime/(p->X207_te-p->X207_ts) - (1.0/PI)*sin(PI*(p->simtime/(p->X207_te-p->X207_ts)));
    }
    
    if(p->X207==1 && p->simtime<p->X207_ts)
    f=0.0;
    
    //cout<<"RAMP F: "<<f<<endl;

    return f;
}