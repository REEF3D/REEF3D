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

#include"6DOF_sflow.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"vrans.h"
   

sixdof_sflow::sixdof_sflow(lexer *p, fdm2D *b, ghostcell *pgc):press_x(p),press_y(p),ddweno_f_nug(p),frk1(p),frk2(p),L(p),dt(p),fb(p),fbio(p),cutr(p),cutl(p),epsi(1.6*p->DXM)
{
    // Initialise
    ini(p,b,pgc);
}

sixdof_sflow::~sixdof_sflow()
{
}

void sixdof_sflow::start
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc,
    double alpha,
    vrans *pvrans,
    vector<net*>& pnet
)
{}

void sixdof_sflow::start
(
	lexer *p, 
	fdm2D *b, 
	ghostcell *pgc
)
{
    // Move body
    p->xg += Uext*p->dt;
    p->yg += Vext*p->dt;

    updateFSI(p,pgc);
    updateForcing_hemisphere(p,pgc);

    print_parameter(p,pgc);
    print_stl(p,pgc);
    
    ++p->printcount_sixdof;
}
