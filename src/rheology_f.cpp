/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"rheology_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h" 

rheology_f::rheology_f(lexer *p, fdm *a) : strain(p,a), tau_x(p), tau_y(p), tau_z(p), epsi(p->F45*p->DXM)
{
    tanphi=0.0;
    if(p->W101>0)
    tanphi=tan(p->W102_phi*(PI/180.0));

}

rheology_f::~rheology_f()
{
}

double rheology_f::viscosity(lexer *p, fdm *a, ghostcell *pgc)
{
	val=0.0;
	
	if(p->W90==1)
	val = Herschel_Bulkley(p,a,pgc);
    
    return val;
}
