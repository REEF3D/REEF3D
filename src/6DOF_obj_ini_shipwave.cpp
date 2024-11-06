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


void sixdof_obj::initialize_shipwave(lexer *p, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"6DOF_obj_ini "<<endl;
    
    // Initialise folder structure
    if(p->X50==1)
	print_ini_vtp(p,pgc);
    
    if(p->X50==2)
    print_ini_stl(p,pgc);
    
    // Initialise processor boundaries
    ini_parallel(p,pgc);
    
    // Initialise objects
	objects_create(p,pgc);
    
    // Initialise fbvel
	ini_fbvel(p,pgc);
    
    // Raycast
    ray_cast_2D(p,pgc);
	reini_2D(p,pgc,fs);
    pgc->gcsl_start4(p,fs,50);
    
    // Calculate geometrical properties
    geometry_parameters_2D(p,pgc);
    
    // Initialise position of bodies
    iniPosition_RBM(p,pgc);
    
    // Raycast
    ray_cast_2D(p,pgc);
	reini_2D(p,pgc,fs);
    pgc->gcsl_start4(p,fs,50);
    
    // Initialise global variables
	update_fbvel(p,pgc);
    
    
    // Print initial body 
    if(p->X50==1)
    print_vtp(p,pgc);
    
    if(p->X50==2)
    print_stl(p,pgc);
}

double sixdof_obj::ramp_vel(lexer *p)
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

double sixdof_obj::ramp_draft(lexer *p)
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
    
    return f;
}