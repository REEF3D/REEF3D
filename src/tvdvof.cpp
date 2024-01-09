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

#include"tvdvof.h"
#include"lexer.h"
#include"fdm.h"

tvdvof::tvdvof(lexer *p)
{
}

tvdvof::~tvdvof()
{
}

double tvdvof::iphi(field& b,int n1, int n2, int q1, int q2)
{
    rp=b(i,j,k);
    rn=b(i+n1,j,k);

    phi = MIN( MAX(   1.0- pow(MAX(pow(1.0-4.0*rp*(1.0-rp),2.0), 1.0-(1.0-4.0*rn*(1.0-rn))),2.0), 0.0), 1.0);

    return phi;
}

double tvdvof::jphi(field& b,int n1, int n2, int q1, int q2)
{
    rp=b(i,j,k);
    rn=b(i,j+n1,k);

    phi = MIN( MAX(   1.0- pow(MAX(pow(1.0-4.0*rp*(1.0-rp),2.0), 1.0-(1.0-4.0*rn*(1.0-rn))),2.0), 0.0), 1.0);

    return phi;
}

double tvdvof::kphi(field& b,int n1, int n2, int q1, int q2)
{
    rp=b(i,j,k);
    rn=b(i,j,k+n1);

    phi = MIN( MAX(   1.0- pow(MAX(pow(1.0-4.0*rp*(1.0-rp),2.0), 1.0-(1.0-4.0*rn*(1.0-rn))),2.0), 0.0), 1.0);

    return phi;
}

