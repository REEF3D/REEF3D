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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_df_void.h"
#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ddweno_f_nug.h"

sixdof_df_void::sixdof_df_void(lexer *p, fdm *a, ghostcell *pgc)
{
}
    
sixdof_df_void::~sixdof_df_void()
{
}

void sixdof_df_void::initialize(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
}
	
void sixdof_df_void::start(lexer*,fdm*,ghostcell*,double,vrans*,vector<net*>&)
{
}

void sixdof_df_void::start_forcing(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, int iter, field& uvel, field& vvel, field& wvel, field& fx, field& fy, field& fz, bool finalise)
{
}

void sixdof_df_void::isource(lexer*,fdm*,ghostcell*)
{
}

void sixdof_df_void::jsource(lexer*,fdm*,ghostcell*)
{
}

void sixdof_df_void::ksource(lexer*,fdm*,ghostcell*)
{
}
    
void sixdof_df_void::isource2D(lexer*,fdm2D*,ghostcell*)
{
}

void sixdof_df_void::jsource2D(lexer*,fdm2D*,ghostcell*)
{
}
