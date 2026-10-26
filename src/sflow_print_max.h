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

#ifndef SFLOW_PRINT_MAX_H_
#define SFLOW_PRINT_MAX_H_

#include"increment.h"
#include"vtp3D.h"
#include"lexer.h"
#include<functional>
#include<fstream>

class fdm2D;
class ghostcell;
class slice4;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::SFLOW maxima of a run (P 112 1, or points P 115): the peak water
level and what a storm surge or flood assessment needs from it, per cell,
updated every time step from P 113 on (after the spin-up):

  max_eta           highest water level while wet [m, relative to the
                    still water level F 60, as the gauges P 51]
  time_max_eta      its time [s]
  max_depth         largest water depth while wet [m]
  max_speed         largest depth-averaged speed while wet [m/s]
  time_max_speed    its time [s]
  flood_duration    time wet [s]
  first_flood_time  first time wet of a cell dry at P 113 [s]
  inundation        0 never wet, 1 wet at P 113, 2 dry at P 113 and
                    flooded later
  max_wind10        largest 10 m wind [m/s]          (W 200, wind)
  min_air_pressure  lowest air pressure [Pa]          (W 200, pressure)

A cell is wet when the solver has it wet (the wet/dry flag) and its water
depth exceeds P 114 (default 0.05 m): dry cells hold the depth A 244 and a
film of water does not flood a cell.  Values that do not exist (max_eta of
a cell never wet) are NaN.  Level 0 of the mesh refinement.

Output, overwritten at every print of the free surface (P 181/P 182), at an
emergency stop and at the end of the run:
  REEF3D_SFLOW_MAX/REEF3D-SFLOW-MAX.pvtp and the pieces -<rank>.vtp: the
     cell values (CellData of the two triangles of a cell), the points at the
     bed (bednode)
  REEF3D_SFLOW_MAX/REEF3D-SFLOW-MAX-points.dat (P 115 x y): the values of
     the cell holding each point, one line per point (high-water marks)

The fields are in the checkpoints (P 48): a continued run (I 45) carries
them on, on any decomposition.
--------------------------------------------------------------------*/

class sflow_print_max : public increment, private vtp3D
{
public:
    sflow_print_max(lexer*, fdm2D*, ghostcell*);
    virtual ~sflow_print_max();

    static bool active(lexer *p) {return p->P112>0 || p->P115>0;}

    // every time step, after it (simtime at its end)
    void update(lexer*, fdm2D*, ghostcell*);

    // the VTP files and the point table
    void write(lexer*, fdm2D*, ghostcell*);

    void ckpt_fields(const std::function<void(const char*, char, double*)>&, const std::function<void(const char*, char, int*)>&);

private:
    void start_window(lexer*, fdm2D*, ghostcell*);
    void cellvalues(lexer*, fdm2D*, int, double*);
    void write_vtp(lexer*, fdm2D*, ghostcell*);
    void write_pvtp(lexer*);
    void write_points(lexer*, fdm2D*, ghostcell*);
    bool wet(lexer*, fdm2D*);

    static const int nfield = 10;
    int nout;                   // fields written (8, + wind, + pressure)
    const char *fname[nfield];

    // the state: max_eta, time_max_eta, max_depth, max_speed, time_max_speed, flood_duration,
    // first_flood_time, the class (0 not started, 1 dry so far, 2 wet at the start, 3 flooded),
    // max_wind10, min_air_pressure
    slice4 *mx[nfield];

    bool started = false, checked = false;
    bool has_wind = false, has_pressure = false;
    int *ploc;                  // P 115: the point is in a fluid cell of this rank
    char name[250];
};

#endif
