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

#include"reduction_deyana.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

reduction_deyana::reduction_deyana(lexer *p) : bedslope(p)
{
}

reduction_deyana::~reduction_deyana()
{
}

void reduction_deyana::start(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double r=1.0;
    eta = 0.85;
    
    SLICELOOP4
    {
    
    alphaval = s->alpha(i,j);
    tetaval = s->teta(i,j);
    phival = s->phi(i,j);
    tanphi = tan(phival);
	alphaval = fabs(alphaval);

	r = (1.0/((1-eta*tanphi)*tanphi))*( -sin(tetaval)  - eta*tanphi*tanphi * sqrt(cos(tetaval)*cos(tetaval)-sin(alphaval)*sin(alphaval))
		+ pow((pow(( sin(tetaval) + eta*tanphi*tanphi*sqrt(cos(tetaval)*cos(tetaval)-sin(alphaval)*sin(alphaval))),2.0) +(1.0 - eta*eta*tanphi*tanphi)
		*(cos(tetaval)*cos(tetaval)*tanphi*tanphi - sin(alphaval)*sin(alphaval)*tanphi*tanphi - sin(tetaval)*sin(tetaval) - sin(alphaval)*sin(alphaval) ) ),0.5 ));

    // limiter
	if(  (pow(( sin(tetaval) + eta*tanphi*tanphi*sqrt(cos(tetaval)*cos(tetaval)-sin(alphaval)*sin(alphaval))),2.0) +(1.0 - eta*eta*tanphi*tanphi)
		*(cos(tetaval)*cos(tetaval)*tanphi*tanphi - sin(alphaval)*sin(alphaval)*tanphi*tanphi - sin(tetaval)*sin(tetaval) - sin(alphaval)*sin(alphaval) ) )  < 0.0 || cos(tetaval)*cos(tetaval)-sin(alphaval)*sin(alphaval) < 0.0)
    {
        if(p->S84==1)
        {
        r = cos(tetaval)*(1.0 - tan(tetaval/tanphi));
        r*= cos(alphaval)*(1.0 - pow(tan(alphaval),2.0)/pow(tanphi,2.0));
        }
        
        if(p->S84==2)
        r = 0.1/(fabs(s->gamma(i,j)) + 0.0000001)+0.1;
    }


    r = MAX(r,0.01);
    r = MIN(r,1.25);

	if(p->pos_x()<p->S71)
	r=1.0;

	if(p->pos_x()>p->S72)
	r=10.0;
    
    
    s->reduce(i,j)=r;
    }
}




