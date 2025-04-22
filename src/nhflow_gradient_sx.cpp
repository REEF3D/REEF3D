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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"nhflow_gradient.h"
#include"fdm_nhf.h"
#include"lexer.h"
#include"slice.h"

double nhflow_gradient::sx(slice &f)
{
    grad = (f(i+1,j)-f(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);

	return grad;
}

double nhflow_gradient::sy(slice &f)
{
    grad = (f(i,j+1)-f(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);

	return grad;
}

double nhflow_gradient::limiter(double v1, double v2)
{
    denom = fabs(v1) + fabs(v2);
    
    denom = fabs(denom)>1.0e-10?denom:1.0e10;
    
    val =  (v1*fabs(v2) + fabs(v1)*v2)/denom;

    return val;	
}
