
/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"6DOF_fsi.h"
#include"mooring_void.h"
#include"mooring_DGSEM.h"
#include"mooring_barQuasiStatic.h"
#include"mooring_Catenary.h"
#include"mooring_Spring.h"
#include"net.h"

#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans.h"
#include"reinidisc_fsf.h"

sixdof_fsi::sixdof_fsi
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc 
) : gradient(p), dt(p), L(p), f(p), frk1(p), cutl(p), cutr(p), fbio(p),epsifb(1.6*p->DXM), epsi(1.6)
{
    prdisc = new reinidisc_fsf(p);
}

sixdof_fsi::~sixdof_fsi(){}

void sixdof_fsi::initialize(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{}

void sixdof_fsi::start
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc, 
    vrans *pvrans,
    vector<net*>& pnet
)
{}
