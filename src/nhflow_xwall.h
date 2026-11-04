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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#ifndef NHFLOW_XWALL_H_
#define NHFLOW_XWALL_H_

#include"lexer.h"
#include<cmath>

// NHFLOW: is the wall (IO 0) at the x side a slip wall with dp/dn = 0, as the walls at the y sides?
// side 1: x-, 4: x+.
// Not where waves travelling along x are generated (B 98 >= 2, x-) or absorbed (B 99 >= 1, x+): the relaxation
// zone or the wavemaker there imposes the wave up to the side, which is then open to it (zero-gradient U at x-,
// pressure correction 0), as before. Waves along x: generation direction B 105 within 60 deg of x.
inline bool nhflow_xwall_slip(lexer *p, int side)
{
    const bool along_x = fabs(cos(p->B105_1*3.14159265358979323846/180.0))>0.5;

    if(side==1 && p->B98>=2 && along_x)
    return false;

    if(side==4 && p->B99>=1 && along_x)
    return false;

    return true;
}

#endif
