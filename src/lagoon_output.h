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
#include<thread>
#include<utility>
#include<vector>

class lexer;
class ghostcell;

class lagoon_output : public increment
{
public:
    lagoon_output(lexer*, ghostcell*, const char *solver);

    // one output of this rank's VTU piece: buffer holds the whole piece, its XML
    // header (point arrays, their offsets) before data_start, the appended data after
    // The arrays are copied and compressed and written by a thread of this rank
    // while the solver goes on; an output is counted (its time committed) at the
    // next output, or at finish(), once every rank has written it.
    void vtu_piece(lexer*, ghostcell*, const std::vector<char> &buffer, size_t data_start, int num);

    // the end of the run: wait for the last output and count it
    void finish(lexer*, ghostcell*);

    // P 18 1: the VTU files are left out (P 18 2: written as well)
    static bool vtu_files(lexer*);

    // set when any LAGOON output of this run gave up (on every rank alike): from then
    // on the VTU/VTP files are written again, also with P 18 1
    static inline bool failed = false;

private:
    std::string solver;
    lagoon_store store;
    bool ready;
    bool usable;
    bool cartesian;  // CFD: levels at fixed heights
    int t;
    int nx, ny, nz;
    int rank;
    std::vector<double> sigma;

    bool start(lexer*, ghostcell*, const std::vector<lagoon_store::variable> &fields);

    // the output being written in the background
    struct job
    {
        int t = -1;
        int num = 0;
        double time = 0.0;
        std::vector<float> z;  // σ-grids: heights of the points, levels outermost
        std::vector<std::pair<std::string, std::vector<float> > > fields;
        bool ok = true;
    };
    job pending;
    std::thread worker;

    void write_job(job *j);
    bool settle(lexer*, ghostcell*);  // wait for the pending output, count it if all ranks wrote it
};

// LAGOON store output of a solver's free surface or bed VTP (P 18; FNPF, NHFLOW, CFD
// topography): each rank's piece, just written, goes into its block of the output
// "free_surface" or "bed" (grid "surface": the point heights z and the point arrays).
// The VTP points are the grid nodes, x outermost (TPSLICELOOP); the store's blocks
// have x fastest. Every rank writes its block, rank 0 then counts the output.
class lagoon_surface
{
public:
    lagoon_surface(lexer*, const char *solver, const char *output, const char *source);

    // after a rank's VTP piece of an output was written to file (all ranks call it):
    // into the store with P 18; with P 18 1 the file is removed again
    static void piece_written(lexer*, ghostcell*, lagoon_surface *&writer, const char *solver,
                              const char *output, const char *source, const char *file, int num);

    void vtp_piece(lexer*, ghostcell*, const std::string &buffer, int num);

private:
    std::string solver, output, source;
    lagoon_store store;
    bool ready, usable;
    int t, nx, ny, rank;

    bool start(lexer*, ghostcell*, const std::vector<lagoon_store::variable> &fields);
};

// LAGOON store output of an AMR solver's free surface (P 18; FNPF, NHFLOW, SFLOW with
// G 1): the grids of every rank (its level 0 and its patches, as its .vtr files have
// them) are gathered on rank 0, which writes them to the AMR set <solver>_amr of
// ./REEF3D_<SOLVER>.lagoon (lagoon_store.h, lagoon_amr): level 0 of all ranks first,
// then the patches, rank by rank, as the .vtm lists them.
class lagoon_amr_output
{
public:
    // solver: "FNPF", "NHFLOW", "SFLOW"; fields: the cell fields of every grid, in order
    lagoon_amr_output(const char *solver, const std::vector<lagoon_amr::field> &fields);

    // every rank, at every AMR output, with its grids. True on every rank when the output
    // is in the store; false (on every rank) when the .vtr and .vtm files are needed
    // instead: from then on nothing more goes to the store.
    bool write(lexer*, ghostcell*, const std::vector<lagoon_amr::grid> &grids, int printcount);

    // with P 18 1, an output in the store has no .vtr and .vtm files
    static bool files_needed(lexer*, bool stored);

private:
    std::string solver;
    std::vector<lagoon_amr::field> fields;
    lagoon_amr *writer;
    bool usable;
};

#endif
