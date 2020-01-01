/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"reduction_deyana.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

reduction_deyana::reduction_deyana(lexer *p) : bedslope(p)
{
}

reduction_deyana::~reduction_deyana()
{
}

double reduction_deyana::start(lexer *p, fdm * a, ghostcell *pgc)
{
    double r=1.0;
    eta = 0.85;

	slope(p,a,pgc,teta,alpha,gamma,phi);


	alpha = fabs(alpha);

	r = (1.0/((1-eta*tan(phi))*tan(phi)))*( -sin(teta)  - eta*tan(phi)*tan(phi) * sqrt(cos(teta)*cos(teta)-sin(alpha)*sin(alpha))
		+ pow((pow(( sin(teta) + eta*tan(phi)*tan(phi)*sqrt(cos(teta)*cos(teta)-sin(alpha)*sin(alpha))),2.0) +(1.0 - eta*eta*tan(phi)*tan(phi))
		*(cos(teta)*cos(teta)*tan(phi)*tan(phi) - sin(alpha)*sin(alpha)*tan(phi)*tan(phi) - sin(teta)*sin(teta) - sin(alpha)*sin(alpha) ) ),0.5 ));

	if(  (pow(( sin(teta) + eta*tan(phi)*tan(phi)*sqrt(cos(teta)*cos(teta)-sin(alpha)*sin(alpha))),2.0) +(1.0 - eta*eta*tan(phi)*tan(phi))
		*(cos(teta)*cos(teta)*tan(phi)*tan(phi) - sin(alpha)*sin(alpha)*tan(phi)*tan(phi) - sin(teta)*sin(teta) - sin(alpha)*sin(alpha) ) )  < 0.0 || cos(teta)*cos(teta)-sin(alpha)*sin(alpha) < 0.0)
    {
        r = cos(teta)*(1.0 - tan(teta/tan(phi)));
        r*= cos(alpha)*(1.0 - pow(tan(alpha),2.0)/pow(tan(phi),2.0));
    }


    r = MAX(r,0.01);
    r = MIN(r,1.25);

	if(p->pos_x()<p->S71)
	r=1.0;

	if(p->pos_x()>p->S72)
	r=10.0;
    
    

    return r;
}




