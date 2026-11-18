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


#ifndef SFLOW_LEVEL_FILE_H_
#define SFLOW_LEVEL_FILE_H_

#include"increment.h"
#include<vector>

class lexer;
class fdm2D;
class ghostcell;

// F 73 1: the initial free-surface level of SFLOW from the file sflow-level.dat, a regular grid in the coordinates of
// the case (after the shift and rotation of the geodata of DIVEMesh, the coordinates REEF3D works in):
//   nx ny x0 y0 dx dy            (x0, y0: the first value; dx, dy: the spacing)
//   nx*ny levels [m], i (x) fastest, then j (y); a value < -9000 is no data
// Bilinear between the four surrounding values (with no data among them: the nearest valid one, within one spacing).
// Where the level lies above the bed the cell gets it: eta = level - wd. Applied after F 72 (an initial slide mass or
// lake from a terrain difference, e.g. a pre-event DEM minus the sliding surface).
class sflow_level_file : public increment
{
public:
    void start(lexer*, fdm2D*, ghostcell*);
};

// the regular grid of levels of such a file, shared by F 73 (sflow-level.dat) and W 161 (slide-level.dat, the top of the
// initial slide of the two-layer model); stops the run with a message when the file is missing or short
class level_grid
{
public:
    void read(lexer*, ghostcell*, const char *file, const char *key);
    double value(double x, double y) const;     // < -9000: no level there
    
    int nx=0, ny=0;
    double x0=0.0, y0=0.0, dx=1.0, dy=1.0;
    
private:
    std::vector<double> L;
};

#endif
