/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

double iowave::ra1(lexer *p, double x)
{
    double r=0.0;

    x=x/dist1;
    r = -2.0*x*x*x + 3.0*x*x;
	r*=ramp(p);    

    return r;
}

double iowave::ra2(lexer *p, double x)
{
    double r=0.0;

    x=(x-dist1)/(dist2-dist1);

    r = -2.0*(1.0-x)*(1.0-x)*(1.0-x) + 3.0*(1.0-x)*(1.0-x);

    if(dist2==dist1)
    r=0.0;

    r*=ramp(p);

    return r;
}

double iowave::ra3(lexer *p, double x)
{
    double r=0.0;

    x=(dist3-fabs(x))/dist3;

    r = 1.0-pow(x,p->B118);
	
	return r;
}

double iowave::rb1(lexer *p, double x)
{
    double r=0.0;

    x=1.0-x/dist1;
    x=MAX(x,0.0);
    
    r = 1.0 - (pow(1.0 + (pow(x,p->B119))/4000.0, 4000.0)-1.0)*expinverse;

    return r;
}

double iowave::rb3(lexer *p, double x)
{
    double r=0.0;

    x=(dist3-fabs(x))/dist3;
    x=MAX(x,0.0);
    
    r = 1.0 - (pow(1.0 + (pow(x,p->B119))/4000.0, 4000.0)-1.0)*expinverse;
	
	return r;
}

double iowave::ramp(lexer *p)
{
    double f=1.0;

    if(p->B101==1 && p->simtime<p->B102*p->wT)
    {
    f = p->simtime/(p->B102*p->wT) - (1.0/PI)*sin(PI*(p->simtime/(p->B102*p->wT)));
    }

    return f;
}
