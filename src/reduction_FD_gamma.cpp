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

#include"reduction_FD_gamma.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

reduction_FD_gamma::reduction_FD_gamma(lexer *p) : bedslope(p)
{
}

reduction_FD_gamma::~reduction_FD_gamma()
{
}

void reduction_FD_gamma::start(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double r=1.0;
    double u0,v0,uvel,vvel,uabs,fx,fy;
    
    SEDSLICELOOP
    {
    u0=0.5*(s->P(i,j)+s->P(i-1,j));
    v0=0.5*(s->Q(i,j)+s->Q(i,j-1));
    
    uvel = (cos(s->beta(i,j))*u0-sin(s->beta(i,j))*v0);
	vvel = (sin(s->beta(i,j))*u0+cos(s->beta(i,j))*v0);
    
    uabs=sqrt(uvel*uvel + vvel*vvel);
    
    fx = fabs(uvel)/(fabs(uabs)>1.0e-10?uabs:1.0e10);
    fy = fabs(vvel)/(fabs(uabs)>1.0e-10?uabs:1.0e10);
    
    r = fx*cos(s->teta(i,j))*(1.0 - tan(s->teta(i,j))/tan(s->phi(i,j))) + (1.0-fx);
    
    r*= (fy*cos(s->alpha(i,j))*(1.0 - pow(tan(s->alpha(i,j)),2.0)/pow(tan(s->phi(i,j)),2.0))  + (1.0-fy));
        
        // limiter
        if( 1.0-s->gamma(i,j)/s->phi(i,j) < 0.0)
        if(p->S84==2)
        r = 0.1/(fabs(s->gamma(i,j)) + 0.0000001)+0.1;
    
    
    r=MAX(r,0.01);
    r=MIN(r,1.5);

	if(p->pos_x()<p->S71)
	r=1.0;

	if(p->pos_x()>p->S72)
	r=10.0;
	
    s->reduce(i,j)=r;
    }
}




