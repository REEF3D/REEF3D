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

#include"nhflow_filter.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"
#include<cmath>
#include<cstdio>

nhflow_filter::nhflow_filter(lexer *p)
{
    int poly_order = 1; 
    int half_width = 1;
    
    order = poly_order;
    hw    = half_width;
    nw    = 2*hw + 1;

    build_coeffs(p);
}

nhflow_filter::~nhflow_filter()
{
}
