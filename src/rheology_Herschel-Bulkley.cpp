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
#include<algorithm>

double rheology_f::Herschel_Bulkley(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
	gamma = strainterm(p,u,v,w); 
    
    tau0=val=0.0;
    pressureval = a->press(i,j,k)-p->pressgage;
    
    if(p->W110==1)
    tau0 = yield_stress(p,a);
    
    if(p->W110!=3 && p->W110!=5)
    {
        val =  (tau0/(gamma>1.0e-20?gamma:1.0e-20) + (p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        
        val = MIN(val,p->W95);
    }
	
    return val;
}
