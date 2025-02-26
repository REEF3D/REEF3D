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

double nhflow_gradient::sxx(slice &f)
{
    return ((f(i+1,j)-f(i,j))/p->DXP[IP] - (f(i,j)-f(i-1,j))/p->DXP[IM1])/p->DXN[IP];
}

double nhflow_gradient::syy(slice &f)
{
    return ((f(i,j+1)-f(i,j))/p->DYP[JP] - (f(i,j)-f(i,j-1))/p->DYP[JM1])/p->DYN[JP];    
}