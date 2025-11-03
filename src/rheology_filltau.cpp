/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"rheology_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void rheology_f::filltau(lexer *p, fdm *a, ghostcell *pgc)
{
    LOOP
    {
        pressurePhi(p,a,1,0,0,true);
        
        // Yield Stress
        tau0 = yield_stress(p,a);
         
        if(p->W110==7)
        tau0 += ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);   
        
        tau_x(i,j,k) = tau0;
        tau_y(i,j,k) = tau0;
        tau_z(i,j,k) = tau0;
    }

    pgc->start4(p,tau_x,1);
    pgc->start4(p,tau_y,1);
    pgc->start4(p,tau_z,1);
}
