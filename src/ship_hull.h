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

#ifndef SHIP_HULL_H_
#define SHIP_HULL_H_

#include<vector>

//  Hull hydrostatics from the surface triangulation (solver independent, no lexer):
//  the parts of the triangles below a horizontal plane z = zw (the still water level in the
//  body frame) give the wetted surface, the waterline length and the draft distribution T(x)
//  along the body x-axis.
//
//  Triangles: tri_x/tri_y/tri_z[n][3], body frame relative to the centre of gravity.

class ship_hull
{
public:
    
    // area of the triangles below z = zw; with skip_y the faces normal to y are left out
    // (2D: the side walls of the one cell wide hull)
    static double wetted_surface(double**, double**, double**, int, double, bool skip_y=false);
    
    // x-range of the triangles below z = zw: aft and forward end of the waterline
    static void waterline_extent(double**, double**, double**, int, double, double&, double&);
    
    // draft distribution: nstrip strips between xa and xf, strip centres xs, widths dx and
    // drafts T (zw minus the lowest submerged point of the hull in the strip)
    static void draft_strips(double**, double**, double**, int, double, double, double, int,
                             std::vector<double>&, std::vector<double>&, std::vector<double>&);
    
    // part of the triangle (a,b,c) below z = zw: polygon with up to 4 vertices, returns their number
    static int clip_below(const double*, const double*, const double*, double, double (*)[3]);
    
    // area of a planar polygon
    static double polygon_area(const double (*)[3], int);
};

#endif
