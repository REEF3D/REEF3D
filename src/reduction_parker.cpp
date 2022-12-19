/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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
#include"reduction_parker.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

reduction_parker::reduction_parker(lexer *p) : bedslope(p)
{
	eta = 0.85;
}

reduction_parker::~reduction_parker()
{
}


void reduction_parker::start(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double r=1.0;
	double r1,r2;
    
    SLICELOOP4
    {
	alphaval = s->alpha(i,j);
    tetaval = s->teta(i,j);
    phival = s->phi(i,j);
    tanphi = tan(phival);

	alphaval = fabs(alphaval);

	mu = atan(1.0/phival);
	d = (4.0/3.0)*mu*0.85*0.74;
	
	pval = (2.0/(1.0-d))*(d/sqrt(1.0 + tan(alphaval)*tan(alphaval) + tan(tetaval)*tan(tetaval)) + sin(tetaval)/mu);

	qval = ((1.0+d)/(1.0-d))*(1.0/(1.0 + tan(alphaval)*tan(alphaval) + tan(tetaval)*tan(tetaval)))*(-1.0 + ((tan(alphaval)*tan(alphaval) + tan(tetaval)*tan(tetaval))/mu));

	r1 = -0.5*pval - sqrt(pval*pval*0.25 - qval);
	
	r = -0.5*pval + sqrt(pval*pval*0.25 - qval);

	// limiter
	if(( (1.0 + tan(alphaval)*tan(alphaval) + tan(tetaval)*tan(tetaval))  < 0.0 || (pval*pval*0.25 - qval) < 0.0) || r<0.0)
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
	r=1.0;
    
    s->reduce(i,j)=r;
    }
}



