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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "sedpart.h"

#include "lexer.h"
#include "ghostcell.h"
#include "fdm.h"
#include "vrans_f.h"
#include "reinitopo.h"
#include "ioflow.h"

/// @brief Initializes everything in the sediment for the CFD solver
/// Determines cell which should be filled with particles
/// Allocates memory for the particles
/// Seeds the particles
/// Prepares particles and particle related variables for simulation
/// Initializes VRANS
void sedpart::ini_cfd(lexer *p, fdm *a,ghostcell *pgc)
{
    // vrans
    pvrans->sed_update(p,a,pgc);
    if(p->I40!=1)
    {
        PLAINLOOP
        cellSumTopo[IJK]=maxParticlesPerCell(p,a,PP.d50);
        // seed
        seed_ini(p,a,pgc);
        PP.reserve(maxparticle);
        // if(p->mpirank==0)
        seed(p,a);
        make_stationary(p,a,&PP);
    }
    
    gparticle_active = pgc->globalisum(PP.size);

    if(gparticle_active>0)
    {
        particlesPerCell(p,a,pgc,&PP);
        particleStressTensor(p,a,pgc,&PP);
    }
    
    // print
    if((p->I40!=1)||(p->I40==1&&inicount>0))
    print_particles(p);
    
    if(p->mpirank==0)
        if(p->I40!=1)
            cout<<"Sediment particles: "<<gparticle_active<<endl;
        else if (inicount>0)
            cout<<"Loaded particles "<<gparticle_active<<" from state file."<<endl;
    
    
    ++inicount;
    if(p->mpirank==1)
    {
        i=8;j=12;k=2;
        std::cout<<maxParticlesPerCell(p,a,PP.d50,true)<<"|"<<maxParticlesPerCell(p,a,PP.d50,false,true)<<endl;
    }
    PLAINLOOP
    a->test(i,j,k)=theta_s(p,a,&PP,i,j,k);
    // particle_func::debug(p,a,pgc,&PP);
}

/// @brief CFD calculation function
/// @param a fdm object
/// @param pflow IO-flow object
/// @param preto topography reinitialization object
void sedpart::start_cfd(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow,
                                    reinitopo* preto, solver* psolv)
{
    double starttime=pgc->timer();
	int xchange=0;
	int removed=0;

	if (p->count>=p->Q43)
	{
        /// runtime seeding
		if(p->Q120==1&&p->count%p->Q121==0)
			posseed_suspended(p,a);
        point_source(p,a);
        if(p->Q101>0)
        {
            // topo_influx(p,a);
            // solid_influx(p,a);
            // seed_topo(p,a);
        }
        particlesPerCell(p,a,pgc,&PP);
        particleStressTensor(p,a,pgc,&PP);

        /// transport
        erode(p,a);
        transport(p,a,&PP);
		xchange=transfer(p,pgc,&PP,maxparticle);
		removed=remove(p,&PP);
        removed += deposit(p,a);

        /// topo update
        if(p->Q13==1)
            update_cfd(p,a,pgc,pflow,preto);

        /// cleanup
        if(p->count%p->Q20==0)
        {
            if(PP.size == 0)
                PP.erase_all();
            // PP.optimize();
            cleanup(p,a,&PP,0);
        }
	}

    /// print out
	print_particles(p);

	gparticle_active = pgc->globalisum(PP.size);
    gremoved = pgc->globalisum(removed);
    gxchange = pgc->globalisum(xchange);
	p->sedsimtime=pgc->timer()-starttime;

    PLAINLOOP
    a->test(i,j,k)=theta_s(p,a,&PP,i,j,k);
    // particle_func::debug(p,a,pgc,&PP);

    if(p->mpirank==0 && (p->count%p->P12==0))
    	cout<<"Sediment particles: "<<gparticle_active<<" | xch: "<<gxchange<<" rem: "<<gremoved<<" | sim. time: "<<p->sedsimtime<<" relative: "<<p->sedsimtime/double(gparticle_active)*(10^3)<<" ms\nTotal bed volume change: "<<std::setprecision(9)<<volumeChangeTotal<<endl;
}

/// @brief Updates the topography for the CFD solver
void sedpart::update_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reinitopo* preto)
{
    const double tolerance=5e-10;
    ILOOP
    JLOOP
    if(topoVolumeChange[IJ]>0)
    {
        double dh=topoVolumeChange[IJ]/p->DXN[IP]/p->DYN[JP];
        a->bed(i,j)+=dh;
        KLOOP
        {
            // Topo update
            a->topo(i,j,k) -= dh;

            // Seeding update
            // active_topo(i,j,k) = 0.0;
            // if( (a->topo(i,j,k)<0.5*p->DZN[KP]-tolerance) && (a->topo(i,j,k)>-p->DZN[KP]*ceil(p->Q102)-tolerance) && (a->solid(i,j,k)>=-p->DXM))
            // {
            // active_topo(i,j,k) = 1.0;
            // if(p->flag1[Im1JK]==SOLID_FLAG&&p->flag1[IJK]==WATER_FLAG)
            // active_topo(i,j,k) = 10.0;
            // }
        }
        volumeChangeTotal += topoVolumeChange[IJ];

        // Reset
        topoVolumeChange[IJ]=0;
    }

    pgc->start4a(p,a->topo,150);
    pgc->gcsl_start4(p,a->bed,50);
    preto->start(p,a,pgc,a->topo);
    if(p->mpirank==0)
        cout<<"Topo: update grid..."<<endl;
    pvrans->sed_update(p,a,pgc);
    pflow->gcio_update(p,a,pgc);
}