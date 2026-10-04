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

#ifndef LAGOON_OUTPUT_H_
#define LAGOON_OUTPUT_H_

// LAGOON store output of a solver's VTU volume (P 18; FNPF and NHFLOW σ-grids, CFD
// Cartesian grids): every rank takes the arrays
// of its VTU piece, as the printer has just put them together, and writes them to
// its block of ./REEF3D_<SOLVER>.lagoon (see lagoon_store.h). The VTU points are
// columns of knoz+1 σ-levels, written x fastest, then y, then the level, which is
// exactly the store's block layout.

#include"lagoon_store.h"
#include"increment.h"
#include<string>
#include<vector>

class lexer;
class ghostcell;

class lagoon_output : public increment
{
public:
    lagoon_output(lexer*, ghostcell*, const char *solver);

    // one output of this rank's VTU piece: buffer holds the whole piece, its XML
    // header (point arrays, their offsets) before data_start, the appended data after
    void vtu_piece(lexer*, ghostcell*, const std::vector<char> &buffer, size_t data_start, int num);


    // P 18 2: the VTU files are left out
    static bool vtu_files(lexer*);

private:
    std::string solver;
    lagoon_store store;
    bool ready;
    bool usable;
    bool cartesian;  // CFD: levels at fixed heights
    int t;
    int nx, ny, nz;
    std::vector<double> sigma;

    bool start(lexer*, ghostcell*, const std::vector<lagoon_store::variable> &fields);
};

#endif
