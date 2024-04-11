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
#include "looping.h"
#include "fdm.h"
#include "reinitopo.h"
#include "vrans_f.h"
#include "ioflow.h"
#include "turbulence.h"
#include "bedshear.h"

#include <sys/stat.h>
#include <string>

/// This class is enabled when using the options for Lagrangian particles and VRANS.\n
/// Initialization of the topography with particles, modification of topo values and print out.
/// @param p control object
/// @param pgc ghostcell object
/// @param pturb turbulence object
sedpart::sedpart(lexer* p, ghostcell* pgc, turbulence *pturb) : particle_func(p), PP(10,p->S20,p->S22,true), active_box(p), active_topo(p), irand(10000), drand(irand)
{
    pvrans = new vrans_f(p,pgc);
    if(p->I40!=1)
    {
        printcount = 0;
        p->partprinttime=0.0;
    }

    // Create Folder
	if(p->mpirank==0 && p->Q180>0 && (p->Q181>0||p->Q182>0))
	    mkdir("./REEF3D_CFD_SedPart",0777);

    // Output configuration to console
    if(p->mpirank==0)
    {
        string buff;
        buff.append("\nSedPart configuration\nParticles and VRANS active\n");
        buff.append("General configuration:\n\tTopo deformation: ");
        p->Q13>0?buff.append("True\n"):buff.append("False\n");
        buff.append("\tBox seeding: ");
        p->Q110>0?buff.append("True\n"):buff.append("False\n");
        buff.append("\tPoint seeding: ");
        p->Q61>0?buff.append("True\n"):buff.append("False\n");
        buff.append("\tTopo seeding: ");
        p->Q101>0?buff.append("True\n"):buff.append("False\n");
        buff.append("\tSuspension seeding: ");
        p->Q120>0?buff.append("True\n"):buff.append("False\n");
        buff.append("Particle properties:\n\td50: "+std::to_string(p->S20)+" m\n\tDensity: "+std::to_string(p->S22)+" kg/m/m/m\n\tPorosity: "+std::to_string(p->S24)+"\n");
        buff.append("Seeding properties:\n\tSeed: "+(p->Q29>0?std::to_string(p->Q29):"time dep.")+"\n\tParticles per cell: "+std::to_string(p->Q24)+"\n\tParticles represened by one: "+std::to_string(p->Q41)+"\n");
        cout<<buff<<endl;
    }
    inicount=0;
}

sedpart::~sedpart()
{
    delete pvrans;
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
        //     topo_influx(p,a);
        //     solid_influx(p,a);
        }
        particlesPerCell(p,pgc,&PP);
        particleStressTensor(p,a,pgc,&PP);

        /// transport
        erode(p,a,pgc);
        transport(p,a,&PP);
        // fixPos(p,a,&PP);
		xchange=transfer(p,pgc,&PP,maxparticle);
		removed=remove(p,&PP);
        make_stationary(p,a,&PP);

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

    if(p->mpirank==0 && (p->count%p->P12==0))
    	cout<<"Sediment particles: "<<gparticle_active<<" | xch: "<<gxchange<<" rem: "<<gremoved<<" | sim. time: "<<p->sedsimtime<<" relative: "<<p->sedsimtime/double(gparticle_active)*(10^3)<<" ms\nTotal bed volume change: "<<std::setprecision(9)<<volumeChangeTotal<<endl;
}

/// @brief Initializes everything in the sediment for the CFD solver
/// Determines cell which should be filled with particles
/// Allocates memory for the particles
/// Seeds the particles
/// Prepares particles and particle related variables for simulation
/// Initializes VRANS
void sedpart::ini_cfd(lexer *p, fdm *a,ghostcell *pgc)
{
    if(p->I40!=1)
    {
        // seed
        seed_ini(p,a,pgc);
        PP.reserve(maxparticle);
        seed(p,a);
        make_stationary(p,a,&PP);
    }
    
    gparticle_active = pgc->globalisum(PP.size);

    if(gparticle_active>0)
    {
        particlesPerCell(p,pgc,&PP);
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
    
    // vrans
    pvrans->sed_update(p,a,pgc);
    ++inicount;
    // if(p->mpirank==1)
    // {
    // i=8;j=12;
    // KLOOP
    // if(active_topo(i,j,k)>0)
    // cout<<k<<","<<p->ZN[KP]<<endl;
    // }
}

/// @brief SFLOW calculation function
void sedpart::start_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow* pflow, slice &P, slice &Q)
{

}

/// @brief SFLOW initialization function
void sedpart::ini_sflow(lexer *p, fdm2D *b, ghostcell *pgc)
{

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
            active_topo(i,j,k) = 0.0;
            if( (a->topo(i,j,k)<0.5*p->DZN[KP]-tolerance) && (a->topo(i,j,k)>-p->DZN[KP]*ceil(p->Q102)-tolerance) && (a->solid(i,j,k)>=-p->DXM))
            {
            active_topo(i,j,k) = 1.0;
            if(p->flag1[Im1JK]==SOLID_FLAG&&p->flag1[IJK]==WATER_FLAG)
            active_topo(i,j,k) = 10.0;
            }
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

/// @brief Updates the topography for the SFLOW solver
void sedpart::update_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow)
{
}

/// @brief Enables erosion of particles
void sedpart::erode(lexer* p, fdm* a, ghostcell* pgc)
{
    if(p->Q101>0)
        make_moving(p,a,&PP);
}

/// @brief Write out particle data to state file
/// @param result statefile
void sedpart::write_state_particles(ofstream& result)
{
    float ffn=num;
    result.write((char*)&ffn, sizeof (float));
    ffn=volumeChangeTotal;
    result.write((char*)&ffn, sizeof (float));
    size_t ffs=PP.capacity;
    result.write((char*)&ffs, sizeof (size_t));
    ffs=PP.size;
    result.write((char*)&ffs, sizeof (size_t));
    PARTICLELOOP
    {
        ffn=PP.X[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.Y[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.Z[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.Flag[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.U[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.V[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.W[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.PackingFactor[n];
        result.write((char*)&ffn, sizeof (float));
    }
}

/// @brief Read in particle data from state file
/// @param result statefile
void sedpart::read_state_particles(ifstream& result)
{
    float ffn;
    result.read((char*)&ffn, sizeof (float));
    printcount=ffn;
    result.read((char*)&ffn, sizeof (float));
    volumeChangeTotal=ffn;
    PP.erase_all();
    size_t ffs;
    result.read((char*)&ffs, sizeof (size_t));
    maxparticle=size_t(ffs);
    PP.reserve(maxparticle);
    result.read((char*)&ffs, sizeof (size_t));
    double x,y,z,flag,u,v,w,packing;
    for(size_t n=0; n<ffs;n++)
    {
        result.read((char*)&ffn, sizeof (float));
        x=double(ffn);
        result.read((char*)&ffn, sizeof (float));
        y=double(ffn);
        result.read((char*)&ffn, sizeof (float));
        z=double(ffn);
        result.read((char*)&ffn, sizeof (float));
        flag=int(ffn);
        result.read((char*)&ffn, sizeof (float));
        u=double(ffn);
        result.read((char*)&ffn, sizeof (float));
        v=double(ffn);
        result.read((char*)&ffn, sizeof (float));
        w=double(ffn);
        result.read((char*)&ffn, sizeof (float));
        packing=double(ffn);
        PP.add(x,y,z,flag,u,v,w,packing);
    } 
}