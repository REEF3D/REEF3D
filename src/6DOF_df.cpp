/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"6DOF_df.h"
#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

sixdof_df::sixdof_df
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc 
)
{
    number6DOF = 1;
    
    for (int nb = 0; nb < number6DOF; nb++)
    {
        p_df_obj.push_back(new sixdof_df_object(p,a,pgc,nb));
    }
}
    
sixdof_df::~sixdof_df(){}

void sixdof_df::initialize(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
    for (int nb = 0; nb < number6DOF; nb++)
    {
        p_df_obj[nb]->initialize(p, a, pgc, pnet);
    }
}
	
void sixdof_df::start(lexer*,fdm*,ghostcell*,double,vrans*,vector<net*>&){};

void sixdof_df::forcing(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, double alpha, field& uvel, field& vvel, field& wvel, field1& fx, field2& fy, field3& fz)
{
    for (int nb = 0; nb < number6DOF; nb++)
    {
        // Calculate forces
        p_df_obj[nb]->forces_stl(p,a,pgc,alpha,uvel,vvel,wvel);

        // Advance body in time
        p_df_obj[nb]->start(p,a,pgc,alpha,pvrans,pnet);

        // Update position and fb level set
        p_df_obj[nb]->updateFSI(p,a,pgc,alpha);

        // Update forcing terms
        p_df_obj[nb]->updateForcing(p,a,pgc,alpha,uvel,vvel,wvel,fx,fy,fz);

        // Save and print
        p_df_obj[nb]->saveTimeStep(p,alpha);
        p_df_obj[nb]->interface(p,true);
        p_df_obj[nb]->print_stl(p,a,pgc);
        p_df_obj[nb]->print_parameter(p, a, pgc);
    }
}
