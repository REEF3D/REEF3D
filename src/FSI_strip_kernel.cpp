/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2024 Tobias Martin

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"FSI_strip.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
    
double fsi_strip::kernel_roma(const double& dist)
{
    double D = 0.0;

    if (fabs(dist) <= 0.5)
    {
        D = 1.0/3.0*(1.0 + sqrt(-3*dist*dist + 1));
    }
    else if (fabs(dist) <= 1.5)
    {    
        D = 1.0/6.0*(5.0 - 3.0*fabs(dist) - sqrt(-3*(1 - fabs(dist))*(1 - fabs(dist)) + 1));
    }
    
    return D;
}