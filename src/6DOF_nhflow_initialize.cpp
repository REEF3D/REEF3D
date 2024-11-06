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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_nhflow.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"vrans.h"
   
void sixdof_nhflow::ini(lexer *p, ghostcell *pgc)
{
}

void sixdof_nhflow::initialize(lexer *p, fdm_nhf *d, ghostcell *pgc, vector<net*>& pnet)
{
    if(p->X10==1 || p->X10==2)
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj[nb]->initialize_nhflow(p, d, pgc, pnet);
    
    if(p->X10==3)
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj[nb]->initialize_shipwave(p, pgc);
}

void sixdof_nhflow::initialize(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
}

void sixdof_nhflow::start_sflow(lexer *p, ghostcell *pgc, int iter, slice &fsglobal, slice &P, slice&Q, slice &fx, slice &fy, bool finalize)
{
}

void sixdof_nhflow::start_cfd(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, int iter, field &uvel, field &vvel, field &wvel, field &fx, field &fy, field &fz, bool finalize)
{
}

