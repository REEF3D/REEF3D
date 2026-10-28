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

#include"vrans_nhflow_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void vrans_nhflow_f::initialize(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"VRANS ini"<<endl;
    
    // ************************
    geometry_ini(p, d, pgc);
    objects_create_vrans(p, pgc);
    
    layer_faces(p, d, pgc);
    
    // pressure projections per stage (B 205): the collocated correction (wide gradient, cell value
    // of CPORNH) and the compact Poisson matrix (face averages of CPORNH) differ where CPORNH jumps at
    // the structure surface; repeated passes converge the remaining divergence, as for membranes
    d->MPROJ = MAX(d->MPROJ, p->B205);
    
    update(p, d, pgc, 1.0, 0);
    // ************************
    
    
}

