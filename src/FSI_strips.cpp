/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
--------------------------------------------------------------------*/

#include"FSI_strips.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"FSI_strip.h"

fsi_strips::fsi_strips()
{
    numberStrips = 1;
	pstrip.reserve(numberStrips);

    for (int num = 0; num < numberStrips; num++)
	{
        pstrip.push_back(new fsi_strip(num));
    }
}
    
fsi_strips::~fsi_strips(){}

void fsi_strips::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
    if (p->mpirank == 0) cout<<"Initializing strips"<<endl;

    for (int num = 0; num < numberStrips; num++)
    {
        pstrip[num]->initialize(p, a, pgc);
    }
}

void fsi_strips::start(lexer*,fdm*,ghostcell*){}


void fsi_strips::forcing(lexer* p, fdm* a, ghostcell* pgc, double alpha, field& uvel, field& vvel, field& wvel, field1& fx, field2& fy, field3& fz, bool finalise)
{
    for (int num = 0; num < numberStrips; num++)
    {
        // Update Lagrangian points 
        pstrip[num]->update_points();
        
        // Get velocity at Lagrangian points
        pstrip[num]->interpolate_vel(p,a,pgc,uvel,vvel,wvel);

        // Advance strip in time
        pstrip[num]->start(p,a,pgc);

        // Get coupling velocities at Lagrangian points
        pstrip[num]->coupling_vel();
        
        // Get coupling forces at Lagrangian points
        pstrip[num]->coupling_force(p,alpha);
        
        // Distribute coupling forces on Eulerian grid 
        pstrip[num]->distribute_forces(p,a,pgc,fx,fy,fz);

        // Calculate body velocities
        //p_df_obj[nb]->calculate_fb_vel(p,a,pgc,alpha,uvel,vvel,wvel,uveln,vveln,wveln);
/*
// Calculate forces
//p_df_obj[nb]->forces_stl(p,a,pgc,alpha,uvel,vvel,wvel);
        // Calculate external forces
        p_df_obj[nb]->calculate_forces(p,a,pgc,alpha,uvel,vvel,wvel,uveln,vveln,wveln);
// Update position and fb level set
p_df_obj[nb]->updateFSI(p,a,pgc,finalise);
// Update forcing terms
p_df_obj[nb]->updateForcing(p,a,pgc,alpha,uvel,vvel,wvel,fx,fy,fz);

        //finalise=true;



        // Save and print
        p_df_obj[nb]->interface(p,true);

        if (finalise == true)
        {
            p_df_obj[nb]->saveTimeStep(p,alpha);
            p_df_obj[nb]->print_stl(p,a,pgc);
            p_df_obj[nb]->print_parameter(p, a, pgc);
        }*/
    }
};
