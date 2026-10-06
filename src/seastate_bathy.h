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


#ifndef SEASTATE_BATHY_H_
#define SEASTATE_BATHY_H_

#include<string>
#include<vector>

/*--------------------------------------------------------------------
REEF3D::SEASTATE bathymetry raster (A 790 1, file seastate-bathy.dat)

The bed level of SEASTATE on every grid (level 0 and the REEFAMR
patches) from a regular raster, so that refined patches see the
bathymetry and coastline at their own resolution (the 2D grid of
DIVEMesh only carries the bed of the level-0 cells).

  line 1      any title
  line 2      nx ny                 raster nodes in x and y
  line 3      x0 y0 dx dy           first node and spacing [m], model frame
  then        ny lines of nx values: the bed level z [m] (positive up,
              the same vertical datum as F 60) of the nodes, row j at
              y0 + j dy from south to north, node i at x0 + i dx
  '$' starts a comment (also at the start of a line)

SWAN bottom files (depth positive down) are converted with
tools/seastate_forcing.py bathy.

  read       reads the file
  at(x,y)    bilinear in the raster, clamped at its edges
  cell(...)  bed of a cell [x0,x1]x[y0,y1] from the raster nodes inside
             the cell: with zdry (the bed level above which a node is dry,
             SEASTATE: still water level - A 705) the cell takes the
             majority of its nodes, wet or dry, and the mean bed of that
             majority (islands and narrow land stay land, the depth of a
             wet cell is not reduced by land nodes); without zdry the mean
             of all nodes. No node inside (cells finer than the raster):
             at() of the centre

No lexer or MPI dependency (unit test Regression/unit/seastate_test.cpp).
--------------------------------------------------------------------*/

class seastate_bathy
{
public:
    bool read(const std::string &file, std::string &err);

    double at(double x, double y) const;
    double cell(double x0, double x1, double y0, double y1, double zdry=1.0e30) const;

    int nx = 0, ny = 0;
    double x0 = 0.0, y0 = 0.0, dx = 1.0, dy = 1.0;
    std::vector<double> z;      // z[j*nx + i]

private:
    double node(int i, int j) const {return z[size_t(j)*nx + i];}
};

#endif
