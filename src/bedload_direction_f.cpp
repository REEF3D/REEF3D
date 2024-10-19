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

#include"bedload_direction_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

bedload_direction_f::bedload_direction_f(lexer *p)
{
    beta = 1.3;
}

bedload_direction_f::~bedload_direction_f()
{
}

void bedload_direction_f::start(lexer* p, ghostcell* pgc, sediment_fdm *s)
{
	double cosa, sina;
    double uvel,vvel,u_abs;
    double bx0,by0;
    
	SLICELOOP4
    {
        uvel=0.5*(s->P(i,j)+s->P(i-1,j));
        vvel=0.5*(s->Q(i,j)+s->Q(i,j-1));
        
         u_abs = sqrt(uvel*uvel + vvel*vvel);
		cosa=fabs(u_abs)>1.0e-10?uvel/fabs(u_abs):0.0;
		sina=fabs(u_abs)>1.0e-10?vvel/fabs(u_abs):0.0;
        
        bx0 = (s->bedzh(i+1,j)-s->bedzh(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);
        by0 = (s->bedzh(i,j+1)-s->bedzh(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);

        s->qbe(i,j) = s->qbe(i,j)*(1.0 - beta*(cosa*bx0 + sina*by0));
	}
    
    pgc->gcsl_start4(p,s->qbe,1);
    
}
