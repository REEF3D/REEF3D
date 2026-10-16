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


#ifndef AMR_GRID_LOOKUP_H_
#define AMR_GRID_LOOKUP_H_

class lexer;
class fdm;
class turbulence;

// The finest grid of a mesh refinement at a point (rank local): the loads of a body (sixdof_obj_cfd
// forces_stl) take the pressure, density, viscosity and velocities there, the probes and wave
// gauges of the CFD printer (P 51, P 61, P 64) their values.  false: level 0.
class amr_grid_lookup
{
public:
    virtual ~amr_grid_lookup() {}
    virtual bool finest_grid(double x, double y, double z, lexer *&q, fdm *&a) const = 0;

    // the same with the turbulence model of that grid
    virtual bool finest_grid(double x, double y, double z, lexer *&q, fdm *&a, turbulence *&t) const = 0;
};

// The same for the 2D-column models of the refinement core reefamr (NHFLOW, FNPF, SFLOW): the
// finest patch of this rank for the output at (x,y), the point two cells inside the cells the
// patch takes from other grids (reefamr::finest_point); its level (0: none here, the output reads
// level 0 where this rank holds the point), its lexer and the model's fdm (fdm_nhf, fdm_fnpf,
// fdm2D).  The printers of the models (P 51 and their probes) read the finest grid with it.
class amr_column_lookup
{
public:
    virtual ~amr_column_lookup() {}
    virtual int finest_column(double x, double y, lexer *&q, void *&f) const = 0;
};

// The rank that writes a point of the output, over all ranks: per point the level of this rank's
// grid there (-1: none, 0: level 0); mine (1: this rank writes it).  A point on level 0 alone:
// every rank that holds it (the caller reduces the values as before); else the rank of the finest
// grid, the lowest of equal ones (with the placement of the patches, G 40, it need not be the rank
// of the level-0 cell).  Collective.
void amr_point_owners(const int *lev, int n, char *mine);

// the grids of the points of an output: find (collective) sets per point the grid to read on this
// rank (q, f; f nullptr: level 0, the printer's own fields) and whether this rank writes the point
// (mine); flag: the point on this rank's level-0 grid
#include<vector>
struct amr_points
{
    std::vector<int> lev;
    std::vector<lexer*> q;
    std::vector<void*> f;
    std::vector<char> mine;
    void find(const amr_column_lookup *amr, lexer *p, int n, const double *x, const double *y, const int *flag);
};

#endif
