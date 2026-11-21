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

#ifndef KS_MAP_H_
#define KS_MAP_H_

#include<string>
#include<vector>

// B 58 1: bed roughness ks (m) per cell from the raster ./ks-map.dat, in the coordinates of the grid (on a
// curvilinear grid G 50 the physical x, y of the nodes; the cell centroid is looked up), bilinear between the
// raster points. Outside the raster, and where a value is <= 0, ks = B 50.
//
//   # comment lines
//   x0 y0 dx dy nx ny          raster point (i,j) at x0 + i dx, y0 + j dy, i = 0..nx-1, j = 0..ny-1
//   ks(0,0) ks(1,0) ... ks(nx-1,0) ks(0,1) ...     nx*ny values, i fastest, free format
//
// Every rank reads the file once (ks_map::get). Used by NHFLOW (level 0 and the multi-block / AMR patches).

class lexer;

class ks_map
{
public:
    static const ks_map& get(lexer*);       // the map of ./ks-map.dat, read on the first call
    double at(double x, double y, double def) const;
    bool ok = false;

private:
    double x0=0.0, y0=0.0, dx=1.0, dy=1.0;
    int nx=0, ny=0;
    std::vector<double> v;
};

#endif
