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

#include"vrans_veg.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

// S_eps = C_feps eps/k S_k, C_feps = 1.33 (Lopez & Garcia 2001, see vrans_veg_k.cpp)
void vrans_veg::eps_source(lexer *p, fdm *a, field &kin, field &eps)
{
	int count=0;
    const double C_feps=1.33;
    
	if(p->B295==1)
    LOOP
    {
    if(N(i,j,k)*D(i,j,k)>0.0)
    a->rhsvec.V[count] += C_feps*MAX(eps(i,j,k),0.0)/MAX(kin(i,j,k),1.0e-10)*veg_sk(p,a);
    
    ++count;
    }
}
