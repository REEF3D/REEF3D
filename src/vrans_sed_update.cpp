/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"vrans_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

// particle sediment (CPM, Q 10 >= 1): porosity from the solid volume fraction of the parcels.
// The grain drag on the fluid is in the particle coupling (CPM_coupling), so these cells get no
// Darcy-Forchheimer resistance; the porous structures B 270 - B 291 are kept and take precedence.
void vrans_f::sedpart_update(lexer *p, fdm *a, ghostcell *pgc, field &por, field &d50)
{
    structures(p,a);
    
    BASELOOP
    if(por(i,j,k)<1.0 && a->porosity(i,j,k)>=1.0)
    set_cell(p,a,por(i,j,k),d50(i,j,k),0.0,0.0);
    
    exchange(p,a,pgc);
}

// continuum sediment (S 10 2, Q 10 0): the bed below topo = 0 is a porous layer (S 24, S 20, S 26)
void vrans_f::sed_update(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->Q10>0)
    return;
    
    if(p->mpirank==0)
    cout<<"Update sediment for VRANS"<<endl;
	
    structures(p,a);
    sediment_bed(p,a);
    exchange(p,a,pgc);
}

void vrans_f::sediment_bed(lexer *p, fdm *a)
{
    LOOP
	if(a->topo(i,j,k)<0.0)
	set_cell(p,a,p->S24,p->S20,p->S26_a,p->S26_b);
}
