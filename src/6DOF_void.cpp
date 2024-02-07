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
--------------------------------------------------------------------*/

#include"6DOF_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans.h"
#include<sys/stat.h>

#include"mooring_void.h"
#include"mooring_barQuasiStatic.h"
#include"mooring_Catenary.h"
#include"mooring_Spring.h"
#include"mooring_dynamic.h"
#include"net.h"
#include"net_void.h"
#include"net_barDyn.h"
#include"net_barQuasiStatic.h"
#include"net_sheet.h"
    
sixdof_void::sixdof_void(lexer*,ghostcell*)
{
}

sixdof_void::~sixdof_void()
{
}

void sixdof_void::start_twoway(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, int iter, field& uvel, field& vvel, field& wvel, field& fx, field& fy, field& fz, bool finalise)
{
}

void sixdof_void::start_oneway(lexer *p, ghostcell *pgc)
{
}

void sixdof_void::initialize(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
}

void sixdof_void::ini(lexer *p, ghostcell *pgc)
{
}

void sixdof_void::isource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_void::jsource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_void::ksource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_void::isource(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
}

void sixdof_void::jsource(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
}

void sixdof_void::ksource(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
}


void sixdof_void::isource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
}

void sixdof_void::jsource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
}