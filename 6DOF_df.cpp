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

sixdof_df::sixdof_df
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc 
) 
{
    p_df_obj = new sixdof_df_object(p,a,pgc);
}
    
sixdof_df::~sixdof_df(){}

void sixdof_df::initialize(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
    p_df_obj->initialize(p, a, pgc, pnet);
}
	
void sixdof_df::sixdof_df::start(lexer *p, fdm *a, ghostcell *pgc, double alpha, vrans *pvrans, vector<net*>& pnet)
{
    p_df_obj->start(p, a, pgc, alpha, pvrans, pnet);
}
    
void sixdof_df::updateFSI(lexer *p, fdm *a, ghostcell *pgc, bool& conv)
{
    p_df_obj->updateFSI(p,a,pgc,conv);
};
    
void sixdof_df::updateForcing(lexer *p, fdm *a, ghostcell *pgc, double alpha, field& uvel, field& vvel, field& wvel, field1& fx, field2& fy, field3& fz)
{
    p_df_obj->updateForcing(p,a,pgc,alpha,uvel,vvel,wvel,fx,fy,fz);
};
	
void sixdof_df::forces_stl(lexer* p, fdm *a, ghostcell *pgc, double alpha)
{
    p_df_obj->forces_stl(p,a,pgc,alpha);
};
    
void sixdof_df::saveTimeStep(lexer *p, double alpha)
{
    p_df_obj->saveTimeStep(p,alpha);
};
	
void sixdof_df::interface(lexer* p, bool conv)
{
    p_df_obj->interface(p,true);
};
    
void sixdof_df::print_parameter(lexer *p,fdm *a,ghostcell *pgc)
{
    p_df_obj->print_parameter(p,a,pgc);
};
 
void sixdof_df::print_stl(lexer *p,fdm *a,ghostcell *pgc)
{
    p_df_obj->print_stl(p,a,pgc);
};
