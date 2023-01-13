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

#include"sflow_fluxlim_minmod.h"
#include"lexer.h"
#include"slice.h"

sflow_fluxlim_minmod::sflow_fluxlim_minmod (lexer *p)
{
}

sflow_fluxlim_minmod::~sflow_fluxlim_minmod()
{
}

double sflow_fluxlim_minmod::iphi(slice& f,int n1, int n2, int q1, int q2)
{
    denom=(f(i+q1,j)-f(i+q2,j));
    r=(f(i+n1,j)-f(i+n2,j))/(fabs(denom)>1.0e-10?denom:1.0e20);

    phi = MAX(0.0, MIN(1.0,r));

    return phi;
}

double sflow_fluxlim_minmod::jphi(slice& f,int n1, int n2, int q1, int q2)
{
    denom=(f(i,j+q1)-f(i,j+q2));
    r=(f(i,j+n1)-f(i,j+n2))/(fabs(denom)>1.0e-10?denom:1.0e20);

    phi = MAX(0.0, MIN(1.0,r));

    return phi;
}
