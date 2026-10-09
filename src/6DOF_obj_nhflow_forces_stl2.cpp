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

#include"6DOF_obj_nhflow.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

double sixdof_obj_nhflow::triangle_area(lexer *p, double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2)
{
    double ax,ay,az;
    double bx,by,bz;
    double cx,cy,cz;
    
    ax = x1-x0;
    ay = y1-y0;
    az = z1-z0;
    
    bx = x2-x0;
    by = y2-y0;
    bz = z2-z0;
    
    // twice the area is the norm of the cross product
    cx = ay*bz - az*by;
    cy = az*bx - ax*bz;
    cz = ax*by - ay*bx;
    
    return 0.5*sqrt(cx*cx + cy*cy + cz*cz);
}

// Parameter along the segment a->b where the height above the local free
// surface (already evaluated at each endpoint) crosses zero.


// Clip a facet against the free surface sampled at the facet's own three
// vertices (fsf0,fsf1,fsf2 = p->wd + eta at (x0,y0),(x1,y1),(x2,y2)) rather
// than a single per-facet constant. Facets sharing an edge in a conforming
// STL mesh share both endpoint vertices exactly, so they compute identical
// fsf values -- and identical cut points -- along that shared edge, which
// is what keeps the wetted patch closed across facet boundaries.
// Returns false when the facet is entirely dry.


