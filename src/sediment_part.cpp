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

#include "sediment_part.h"
#include "partres.h"

#include "lexer.h"
#include "ghostcell.h"
#include "looping.h"
#include "fdm.h"
#include "reinitopo.h"
#include "vrans_f.h"
#include "vrans_v.h"
#include "ioflow.h"
#include "turbulence.h"
#include "bedshear.h"
#include "sediment_fdm.h"
#include "reduction_FD.h"

#include <sys/stat.h>
#include <string>
#include <vector>

/// This class is enabled when using the options for Lagrangian particles and VRANS.\n
/// Initialization of the topography with particles, modification of topo values and print out.
/// @param p control object
/// @param pgc ghostcell object
/// @param pturb turbulence object
sediment_part::sediment_part(lexer* p, ghostcell* pgc, turbulence *pturb) : particle_func(p), PP(10,p->S20,p->S22,true), active_box(p), active_topo(p), irand(10000), drand(irand), s(p), pbedshear(p,pturb), prelax(p), pslope(p)
{
    // pvrans = new vrans_f(p,pgc);
    pvrans =  new vrans_v(p,pgc);
    movement = new partres(p);
    preduce = new reduction_FD(p);
    sediment_part::pturb = pturb;

    prec = 6;

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
        buff.append("\nSedPart configuration\nParticles and !VRANS active\n");
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
        buff.append("Particle properties:\n\td50: "+std::to_string(p->S20)+" m\n\tDensity: "+std::to_string(p->S22)+" kg/m^3\n\tPorosity: "+std::to_string(p->S24)+"\n");
        buff.append("Seeding properties:\n\tSeed: "+(p->Q29>0?std::to_string(p->Q29):"time dep.")+"\n\tParticles per cell: "+std::to_string(p->Q24)+"\n\tParticles represened by one: "+std::to_string(p->Q41)+"\n");
        cout<<buff<<endl;
    }
    inicount=0;
}

sediment_part::~sediment_part()
{
    delete pvrans;
    delete movement;
}

/// @brief Enables erosion of particles
void sediment_part::erosion(lexer* p, fdm* a)
{
    if(p->Q101>0)
    {
        movement->erosion(p,*a,PP,s);
    }
}

/// @brief Deposits moving particles onto topo
void sediment_part::deposition(lexer* p, fdm* a)
{    
    if(p->Q101>0)
    {
        movement->deposition(p,*a,PP,s);
    }
}

void sediment_part::debug(lexer* p, fdm* a, ghostcell* pgc)
{
    movement->debug(p,*a,*pgc,PP,s);
}