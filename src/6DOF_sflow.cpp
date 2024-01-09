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

#include"6DOF_sflow.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"vrans.h"
   
sixdof_sflow::sixdof_sflow(lexer *p, ghostcell *pgc):press(p),ddweno_f_nug(p),frk1(p),frk2(p),L(p),dt(p),
                                                              fb(p),fbio(p),cutr(p),cutl(p),Ls(p),Bs(p),
                                                              Rxmin(p),Rxmax(p),Rymin(p),Rymax(p),draft(p),epsi(1.6*p->DXM)
{
    trisum=1;
    p->Darray(tri_xn,trisum,3);
	p->Darray(tri_yn,trisum,3);
	p->Darray(tri_zn,trisum,3);
}

sixdof_sflow::~sixdof_sflow()
{
}

void sixdof_sflow::start(lexer *p, ghostcell *pgc)
{

// FB/Ship location

    // Move body
    p->xg += ramp_vel(p)*Uext*p->dt;
    p->yg += ramp_vel(p)*Vext*p->dt;

    // Update position
    updateFSI(p,pgc);
    
// --------------------------

    // Update pressure field
    if (p->X400 == 1)
    {
        updateForcing_hemisphere(p,pgc);
    }
    
    else if (p->X400 == 2)
    {
        updateForcing_box(p,pgc);
    }
    
    else if (p->X400 == 3)
    {
        updateForcing_oned(p,pgc);
    }
    
    else if (p->X400 == 10)
    {
        updateForcing_ship(p,pgc);
    }

    // Print
    print_parameter(p,pgc);
    
    if(p->X50==1)
    print_vtp(p,pgc);
    
    if(p->X50==2)
    print_stl(p,pgc);
}
