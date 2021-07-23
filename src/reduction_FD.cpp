/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"reduction_FD.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

reduction_FD::reduction_FD(lexer *p) : bedslope(p)
{
}

reduction_FD::~reduction_FD()
{
}

double reduction_FD::start(lexer *p, fdm * a, ghostcell *pgc)
{
    double r=1.0;

	slope(p,a,pgc,teta,alpha,gamma,phi);
    
    
    r = cos(teta)*(1.0 - tan(teta)/tan(phi));
    
    r*= cos(alpha)*(1.0 - pow(tan(alpha),2.0)/pow(tan(phi),2.0));
    
    r=MIN(r,1.25);
    r=MAX(r,0.01);
    

	if(p->pos_x()<p->S71)
	r=1.0;

	if(p->pos_x()>p->S72)
	r=10.0;
	
    return r;
}




