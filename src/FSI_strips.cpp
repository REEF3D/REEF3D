/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2026 Tobias Martin

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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"FSI_strips.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"FSI_strip.h"
#include<mpi.h>

fsi_strips::fsi_strips(lexer *p, ghostcell *pgc) : eps0(p), pturb(nullptr)
{
	pgc->bcast_int(&p->FSI_count,1);
    numberStrips = p->FSI_count;

    pstrip.reserve(numberStrips);
    for (int num = 0; num < numberStrips; num++)
	{
        pstrip.push_back(new fsi_strip(p,num));
    }
}
    
fsi_strips::~fsi_strips(){}

void fsi_strips::initialize(lexer *p, fdm *a, ghostcell *pgc, turbulence *ppturb)
{
    if (p->mpirank==0) cout<<"Initializing strips"<<endl;
    
    pturb = ppturb;

    for (int num = 0; num < numberStrips; num++)
    {
        pstrip[num]->initialize(p, a, pgc, pturb);
        
        pstrip[num]->print_stl(p,a,pgc);
        pstrip[num]->print_parameter(p, a, pgc);
    }
    
    int ntot = 0;
    for (int num = 0; num < numberStrips; num++)
    ntot += 3*pstrip[num]->numLagrangePoints();
    
    velbuf.assign(ntot,0.0);
}

void fsi_strips::start(lexer*,fdm*,ghostcell*){}


void fsi_strips::forcing(lexer* p, fdm* a, ghostcell* pgc, double alpha, field &uvel, field &vvel, field &wvel, field &fx, field &fy, field &fz, bool finalize)
{
    starttime0=pgc->timer();
    
    // Get velocity at Lagrangian points: local interpolation for all strips,
    // then a single reduction (was 3 MPI_Allreduce per Lagrangian point and strip)
    int offset = 0;
    for (int num = 0; num < numberStrips; num++)
    {
        pstrip[num]->interpolate_vel(p,a,pgc,uvel,vvel,wvel);
        pstrip[num]->pack_vel(&velbuf[offset]);
        offset += 3*pstrip[num]->numLagrangePoints();
    }
    
    starttime=pgc->timer();
    if (offset>0)
    MPI_Allreduce(MPI_IN_PLACE, velbuf.data(), offset, MPI_DOUBLE, MPI_SUM, pgc->mpi_comm);
    double synctime = pgc->timer()-starttime;
    
    offset = 0;
    for (int num = 0; num < numberStrips; num++)
    {
        pstrip[num]->unpack_vel(&velbuf[offset]);
        offset += 3*pstrip[num]->numLagrangePoints();
    }
    
    // RANS forcing field, shared by all strips
    if(p->T10==2)
    {
        LOOP
        eps0(i,j,k) = 0.0;
        
        pgc->start4(p,eps0,30);
    }
    
    for (int num = 0; num < numberStrips; num++)
    {
        // Advance strip in time
        pstrip[num]->start(p,a,pgc,alpha);     // main time consumer
        
        // Get coupling velocities at Lagrangian points
        pstrip[num]->coupling_vel();

        // Get coupling forces at Lagrangian points
        pstrip[num]->coupling_force(p,alpha);
        
        // Distribute coupling forces on Eulerian grid 
        pstrip[num]->distribute_forces(p,a,pgc,fx,fy,fz,eps0);
        
        // Update Lagrangian points 
        pstrip[num]->update_points();

        // Store variables
        pstrip[num]->store_variables(p);

        // Print
        if (finalize==true)
        {
            pstrip[num]->print_stl(p,a,pgc);
            pstrip[num]->print_parameter(p, a, pgc);
        }
    }
    
    // RANS turbulence forcing and ghost-cell update, once for all strips
    if(p->T10==2)
    LOOP
    if(eps0(i,j,k)>1.0e-8)
    pturb->epsget(i,j,k,eps0(i,j,k));
    
    pgc->start1(p,fx,10);
    pgc->start2(p,fy,11);
    pgc->start3(p,fz,12);

    if(p->mpirank==0)
    cout<<"FSI time: "<<pgc->timer()-starttime0<<"  FSI_sync time: "<<synctime<<endl;
    
    // Beam solver statistics of the first strip (per Integrate call)
    if(p->mpirank==0 && finalize==true && numberStrips>0)
    {
        int steps, rejected, nfcn, njac;
        pstrip[0]->getSolverStats(steps, rejected, nfcn, njac);
        cout<<"FSI strip 0 RADAU5: steps "<<steps<<"  rejected "<<rejected<<"  RHS evals "<<nfcn<<"  Jacobians "<<njac<<endl;
    }
};
