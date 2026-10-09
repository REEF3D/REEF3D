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

#include"interface_width.h"
#include"density_df.h"
#include"lexer.h"
#include"fdm.h"
#include"heaviside.h"

density_df::density_df(lexer* p)
{
    H=0.0;
}

density_df::~density_df()
{
}

double density_df::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{     
    phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));

    H = heaviside(phival,interface_width_face(p,a->phi,i,j,k,aa,bb,cc));
    
    
    roval = p->W1*H + p->W3*(1.0-H);

	return roval;	
}




