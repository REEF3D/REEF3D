/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"reduction_deyemp.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

reduction_deyemp::reduction_deyemp(lexer *p) : bedslope(p)
{
}

reduction_deyemp::~reduction_deyemp()
{
}

double reduction_deyemp::start(lexer *p, fdm * a, ghostcell *pgc)
{
    double r=1.0;

	slope(p,a,pgc,teta,alpha,gamma,phi);

	alpha = fabs(alpha);

	r = 0.954*pow(1.0-teta/phi, 0.745)*pow(1.0-alpha/phi,0.372);

	if( 1.0-teta/phi < 0.0 || 1.0-alpha/phi< 0.0)
	r = 0.1/(fabs(gamma) + 0.0000001)+0.1;


    if(r<0.0)
    r = 0.0001;

	if(p->pos_x()<p->S71)
	r=1.0;

	if(p->pos_x()>p->S72)
	r=10.0;
	
    return r;
}




