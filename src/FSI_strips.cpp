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
    for (int num = 0; num < numberStrips; num++)
    {
        pstrip[num]->initialize(p, a, pgc);
    }
}

void fsi_strips::start(lexer*,fdm*,ghostcell*)
{
    for (int num = 0; num < numberStrips; num++)
    {
        // Calculate forces
       /* p_df_obj[nb]->forces_stl(p,a,pgc,alpha,uvel,vvel,wvel);

        // Advance body in time
        p_df_obj[nb]->start(p,a,pgc,alpha,pvrans,pnet);

        // Update position and fb level set
        p_df_obj[nb]->updateFSI(p,a,pgc,finalise);

        // Update forcing terms
        p_df_obj[nb]->updateForcing(p,a,pgc,alpha,uvel,vvel,wvel,fx,fy,fz);

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
