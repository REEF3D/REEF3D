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

// B 295 2: suppress the eddy viscosity inside the static porous structures,
// blended with the structure Heaviside: nu_t -> (1-H) nu_t + H 0.001 nu
void vrans_nhflow_f::eddyv_func(lexer *p, fdm_nhf *d)
{
    if(p->B295==2)
    LOOP
    {
        H = Hporface(p,d,0,0,0);
        
        if(H>1.0e-12)
        {
        d->EV[IJK]  = (1.0-H)*d->EV[IJK]  + H*0.001*d->VISC[IJK];
        d->EV0[IJK] = (1.0-H)*d->EV0[IJK] + H*0.001*d->VISC[IJK];
        }
    }
}
