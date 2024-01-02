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

#include"vanleer.h"
#include"lexer.h"
#include"fdm.h"

vanleer::vanleer (lexer *p)
{
}

vanleer::~vanleer()
{
}

double vanleer::iphi(field& b,int n1, int n2, int q1, int q2)
{
    denom=(b(i+q1,j,k)-b(i+q2,j,k));
    r=(b(i+n1,j,k)-b(i+n2,j,k))/(denom+1.0e-20);

    phi = (r+fabs(r))/(1.0+fabs(r));

    return phi;
}

double vanleer::jphi(field& b,int n1, int n2, int q1, int q2)
{
    denom=(b(i,j+q1,k)-b(i,j+q2,k));
    r=(b(i,j+n1,k)-b(i,j+n2,k))/(denom+1.0e-20);

    phi = (r+fabs(r))/(1.0+fabs(r));
	
    return phi;
}

double vanleer::kphi(field& b,int n1, int n2, int q1, int q2)
{
    denom=(b(i,j,k+q1)-b(i,j,k+q2));
    r=(b(i,j,k+n1)-b(i,j,k+n2))/(denom+1.0e-20);

    phi = (r+fabs(r))/(1.0+fabs(r));
	
    return phi;
}
