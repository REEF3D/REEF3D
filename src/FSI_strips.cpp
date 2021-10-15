/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2021 Tobias Martin

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

fsi_strips::fsi_strips(lexer *p, ghostcell *pgc)
{
	MPI_Bcast(&p->FSI_count,1,MPI_DOUBLE,0,pgc->mpi_comm);
    numberStrips = p->FSI_count;

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
        // Get velocity at Lagrangian points
        pstrip[num]->interpolate_vel(p,a,pgc,uvel,vvel,wvel);

        // Advance strip in time
        pstrip[num]->start(p,a,pgc,alpha);

        // Get coupling velocities at Lagrangian points
        pstrip[num]->coupling_vel();

        // Get coupling forces at Lagrangian points
        pstrip[num]->coupling_force(p,alpha);
        
        // Distribute coupling forces on Eulerian grid 
        pstrip[num]->distribute_forces(p,a,pgc,fx,fy,fz);
        
        // Update Lagrangian points 
        pstrip[num]->update_points();

        // Store variables
        pstrip[num]->store_variables(p);

        // Print
        if (finalise == true)
        {
            pstrip[num]->print_stl(p,a,pgc);
            pstrip[num]->print_parameter(p, a, pgc);
        }
    }
};
