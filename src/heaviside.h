/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef HEAVISIDE_H_
#define HEAVISIDE_H_

#include<cmath>
#include"definitions.h"

// smoothed Heaviside step of a signed distance x over the half width e:
// 0 for x < -e, 1 for x > e, 0.5 (1 + x/e + sin(pi x/e)/pi) in between
inline double heaviside(double x, double e)
{
    if(x>e)
    return 1.0;

    if(x<-e)
    return 0.0;

    return 0.5*(1.0 + x/e + (1.0/PI)*sin((PI*x)/e));
}

// derivative of heaviside(x,e) with respect to x (smoothed delta function)
inline double heaviside_delta(double x, double e)
{
    if(fabs(x)>e)
    return 0.0;

    return 0.5*(1.0 + cos((PI*x)/e))/e;
}

#endif
