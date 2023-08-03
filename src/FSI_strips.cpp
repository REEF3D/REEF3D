/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2023 Tobias Martin

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


void fsi_strips::forcing(lexer* p, fdm* a, ghostcell* pgc, double alpha, field &uvel, field &vvel, field &wvel, field &fx, field &fy, field &fz, bool finalise)
{
    starttime0=pgc->timer();
    
    for (int num = 0; num < numberStrips; num++)
    {
        // Get velocity at Lagrangian points
        pstrip[num]->interpolate_vel(p,a,pgc,uvel,vvel,wvel);
        
        endtime=pgc->timer();
        
        if(p->mpirank==0)
        cout<<"T0: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();

        // Advance strip in time
        pstrip[num]->start(p,a,pgc,alpha);
        
        endtime=pgc->timer();
        
        if(p->mpirank==0)
        cout<<"T1: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();

        // Get coupling velocities at Lagrangian points
        pstrip[num]->coupling_vel();
        
        endtime=pgc->timer();
        
        if(p->mpirank==0)
        cout<<"T2: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();

        // Get coupling forces at Lagrangian points
        pstrip[num]->coupling_force(p,alpha);
        
        endtime=pgc->timer();
        
        if(p->mpirank==0)
        cout<<"T3: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();
        
        // Distribute coupling forces on Eulerian grid 
        pstrip[num]->distribute_forces(p,a,pgc,fx,fy,fz);
        
        endtime=pgc->timer();
        
        if(p->mpirank==0)
        cout<<"T4: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();
        
        // Update Lagrangian points 
        pstrip[num]->update_points();
        
        endtime=pgc->timer();
        
        if(p->mpirank==0)
        cout<<"T5: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();

        // Store variables
        pstrip[num]->store_variables(p);
        
        endtime=pgc->timer();
        
        if(p->mpirank==0)
        cout<<"T6: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();

        // Print
        if (finalise == true)
        {
            pstrip[num]->print_stl(p,a,pgc);
            pstrip[num]->print_parameter(p, a, pgc);
        }
        
        endtime=pgc->timer();
        
        if(p->mpirank==0)
        cout<<"T2: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();
    }
    

    if(p->mpirank==0)
    cout<<"FSI time: "<<pgc->timer()-starttime0<<endl;
};
