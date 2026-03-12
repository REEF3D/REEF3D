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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

double iowave::cosh_func(double x)
{
    if(x<3.0)
    {
    double x2 = x * x;
    double num = 1.0 + x2 * (115.0/252.0 + x2 * (313.0/15120.0));
    double den = 1.0 + x2 * (-11.0/252.0 + x2 * (13.0/15120.0));
    return num / den;
    }
    
    else if(x>=3.0 && x<20.0)
    {
    double e = exp(x);
    return 0.5 * (e + 1.0/e);
    }
    
    else if(x>20.0)
    return 0.5 * exp(x);
}