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

#ifndef HYDRO_COUPLE_H_
#define HYDRO_COUPLE_H_

#include"increment.h"
#include<fstream>
#include<string>
#include<vector>

class lexer;
class fdm2D;
class ghostcell;
class hydro_core;

using namespace std;

/*
REEF3D::HYDRO M3: the catchment model (hydro_core) coupled with SFLOW in one run (R 100 2, A 10 2)

Overlap (built at the start from the cells of both grids, no file): every SFLOW cell belongs to the hydrology
cell its centre lies in (world coordinates, Xout/Yout); A_ov(h) is the SFLOW area of hydrology cell h.

Every hydrology step dt_h (R 101), at the SFLOW time t - t0 (t0: the time of the hot start, R 108):
 1. SFLOW -> HYDRO: the outflow of SFLOW over the open boundaries in the last interval enters the network as
    the inflow of cell ret, the cell below the most downstream cell whose flow path enters SFLOW
    (ret = -1: it leaves the catchment); the infiltration of SFLOW (R 30) on the cells of hydrology cell h
    goes into the upper zone of h.
 2. the catchment model steps from t to t + dt_h (hydro_core::step)
 3. HYDRO -> SFLOW for the interval, as source rates C of hydro_sflow (hydro_sflow::C):
      the share s = min(1, A_ov/A) of the runoff of every cell, spread evenly over its SFLOW cells,
      the discharge of every cell whose fine flow path (hydro_grid.dat version 2) enters SFLOW, at the
      first point of the path in a fluid SFLOW cell (the cell is not passed downstream)
 4. the SFLOW steps are cut to end on t + dt_h, so that SFLOW takes exactly the volumes HYDRO gives.

Budget of both models REEF3D_HYDRO/REEF3D-HYDRO-Coupled.dat at every exchange, cumulative [m3]:
precipitation, glacier melt, evaporation (HYDRO), outflow of the catchment, SFLOW rain, sources of R 40,
boundary inflow, evaporation, abstraction, infiltration outside the hydrology cells, outflow that leaves
(ret = -1), the storage of both models with the return volumes of the last interval, residual.

Not yet: mesh refinement and multi-block grids of SFLOW, checkpoints of the coupled run (R 107 / R 108 for
the catchment model alone), NHFLOW (M6).
*/

class hydro_couple : public increment
{
public:
    hydro_couple(lexer*, fdm2D*, ghostcell*);
    ~hydro_couple();

    void ini(lexer*, fdm2D*, ghostcell*);
    void step(lexer*, fdm2D*, ghostcell*);      // the start of an SFLOW step
    void final(lexer*, fdm2D*, ghostcell*);     // the end of the run

private:
    void stop(lexer*, ghostcell*, const string&) const;
    void exchange_in(lexer*, fdm2D*, ghostcell*);
    void sources(lexer*, fdm2D*, ghostcell*);
    void budget(lexer*, fdm2D*, ghostcell*);
    double xw(lexer*, int, int) const;
    double yw(lexer*, int, int) const;

    hydro_core *hc;
    double t0, tnext;
    vector<int> hcell;              // hydrology cell of every SFLOW cell (SLICELOOP4 order), -1 none
    vector<double> aov;             // SFLOW area of every hydrology cell [m2]
    struct entry {int n, ci, cj; double area;};
    vector<entry> ent;              // entry points in cells of this rank
    vector<double> Fp;              // (1 - sealed) F of every SFLOW cell at the last exchange [m]
    double Vqo0, Vret, Vleave, Vinfout, Vinf0, S0;
    int nexch;
    ofstream fout;
};

#endif
