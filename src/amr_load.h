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

#ifndef AMR_LOAD_H_
#define AMR_LOAD_H_

#include<string>
#include<fstream>
#include<vector>
#include<functional>

class lexer;

//  Load statistics of the mesh refinement on several ranks (G 43 N: a report every N steps, all
//  modules with refinement: SFLOW, FNPF, NHFLOW, CFD).
//
//  Per rank, over the steps since the last report:
//   wall    the wall-clock time
//   mpi     the time in blocking MPI calls (collectives, waits, sendrecv; mpi_clock.cpp wraps them
//           through the MPI profiling interface, so the waits of every exchange and of the
//           reductions of the composite solves are counted)
//   patch   the patch kernels: the time with the MPI exchange of the ghostcell class switched off
//           (reefamr_comms_off), without the MPI calls made within. This holds the rank-local
//           parts of the composite operations as well, so a rank without patches shows some
//           (FNPF: the level-0 rows of the composite Laplace solve). CFD: its own clock of the
//           loops over the patches and their V-cycles (wtp, tp[0])
//   regrid  the regrids (with their MPI time)
//   busy    wall - mpi: the work of the rank; the ranks with less work wait in the MPI calls
//  The load imbalance is max/mean of busy over the ranks, the parallel efficiency mean(busy)/wall.
//  The phases are reported as well (patch, and level 0: busy - patch - the work of the regrids):
//  the ranks synchronise between them, so a balanced sum with unbalanced phases still waits.
//  The predicted work of a rank is the model of the patch placement (G 40 1): its level-0 cells
//  plus the cells of its patches (EXT cells and layers included; CFD: times its weight wpatch).
//  The fit busy = a n0 + b np over the ranks gives the measured cost of a patch cell in level-0
//  cells (b/a).
//
//  Files (rank 0): <dir>/REEF3D_<M>_AMR_load.dat (one line per report) and
//  <dir>/REEF3D_<M>_AMR_load_ranks.dat (one line per rank and report).

// process-wide clocks
struct amr_clock
{
    static double mpi;          // blocking MPI calls (mpi_clock.cpp)
    static long   mpi_calls;
    static double patch;        // patch kernels (reefamr_comms_off, outside the regrids)
    static double regrid;       // regrids, their MPI time included
    static double regrid_mpi;   // the MPI time within the regrids
    static int    depth;        // nesting of reefamr_comms_off
    static int    in_regrid;    // nesting of the regrid clock
    static double t_patch;      // start of the outermost patch scope
    static double m_patch;      // MPI time at its start (the MPI calls within a patch scope are MPI time)

    static void patch_begin();
    static void patch_end();
};

// scope of a regrid
struct amr_regrid_clock
{
    double t, m;
    amr_regrid_clock();
    ~amr_regrid_clock();
};

// counts of this rank, given by the module at every step
struct amr_load_counts
{
    double n0 = 0.0;        // level-0 cells (multi-block: with its level-0 blocks)
    double n0a = 0.0;       // level-0 cells that are not solid
    double np = 0.0;        // cells of the patches of the rank (EXT cells and layers included)
    int npatch = 0;         // patches of the rank
    double wpred = 0.0;     // predicted work (the model of the placement)
    double tpatch = -1.0;   // the module's own clock of the patch kernels (cumulative), <0: amr_clock::patch
};

class amr_load
{
public:
    // module name as in the messages ("SFLOW AMR"); output into ./REEF3D_<name with _>
    void configure(const std::string &name);
    bool ready() const { return !name.empty(); }
    // once per step; a report every G 43 steps, with the counts of this rank (evaluated only then)
    void step(lexer*, const std::function<amr_load_counts()>&);

private:
    void report(lexer*, const amr_load_counts&, int);

    std::string name, dir, tag;
    bool started = false;
    int steps = 0;
    double w0 = 0.0, m0 = 0.0, p0 = 0.0, r0 = 0.0, rm0 = 0.0;
    std::ofstream sout, rout;
};

#endif
