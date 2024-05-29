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
#include "sedpart_movement.h"

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
#include <vector>

/// This class is enabled when using the options for Lagrangian particles and VRANS.\n
/// Initialization of the topography with particles, modification of topo values and print out.
/// @param p control object
/// @param pgc ghostcell object
/// @param pturb turbulence object
sedpart::sedpart(lexer* p, ghostcell* pgc, turbulence *pturb) : particle_func(p), PP(10,p->S20,p->S22,true), active_box(p), active_topo(p), irand(10000), drand(irand)
{
    pvrans = new vrans_f(p,pgc);
    movement = new sediment_particle::movement::Tavouktsoglou(p);
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

/// @brief Enables erosion of particles
void sedpart::erode(lexer* p, fdm* a)
{
    if(p->Q101>0)
        make_moving(p,a,&PP);
}

/// @brief Deposits moving particles onto topo
int sedpart::deposit(lexer* p, fdm* a)
{
    if(p->Q101>0)
        make_stationary(p,a,&PP);
    return solid_clean(p,&PP);
}

void  sedpart::debug(lexer* p, fdm* a, ghostcell* pgc)
{
    std::vector<double> x,y,z;
    std::vector<double> u,v,w;
    x.push_back(0.2855);
    for(int n=0;n<100;n++)
    x.push_back(x[n]+0.0005);
    y.push_back(0.15);
    z.push_back(0);
    if(p->mpirank==2)
    {
        for(int n=0;n<x.size();n++)
        {
            u.push_back(p->ccipol1c(a->u,a->solid,x[n],y[0],z[0]));
        }
        // n=0;
        // double dist = p->ccipol4_b(a->solid,x[n],y[0],z[0]);
        // cout<<dist<<","<<x[n]+dist<<","<<p->ccipol1c(a->u,a->solid,x[n]+dist,y[0],z[0])<<endl;
        // u.push_back(p->ccipol1c(a->u,x[n]+dist,y[0],z[0]));
        string output;
        for(int n=0;n<u.size();n++)
        if(u[n]>0)
        output += "("+std::to_string(x[n])+")"+std::to_string(u[n])+",";
        else
        {
            output += "\n("+std::to_string(x[n])+")"+std::to_string(u[n]);
            break;
        }
        cout<<p->mpirank<<"|"<<output<<endl;
    }
}