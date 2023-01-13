/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"density_vof.h"
#include"lexer.h"
#include"fdm.h"

density_vof::density_vof(lexer* p) : epsi(p->F45*p->DXM), eps(2.1*p->DXM)
{
    if(p->j_dir==0)        
    psi = p->F45*(1.0/2.0)*(p->DRM+p->DTM);
        
    if(p->j_dir==1)
    psi = p->F45*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);
}

density_vof::~density_vof()
{
}

double density_vof::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{
	H= 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));

		H=MAX(H,0.0);
		H=MIN(H,1.0);

    roval = p->W1*H +   p->W3*(1.0-H);
	
	return roval;		
}




