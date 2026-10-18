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

#ifndef REEFRUN_CKPT_H_
#define REEFRUN_CKPT_H_

//  Checkpoints of a run on any decomposition (P 48: write every P 48 steps, I 45 1: restart from the
//  latest; the files: reefckpt.h) and the rebalance on the measured load (G 46, G 47), for every
//  module that hands its state over (hooks): FNPF and NHFLOW (driver_ckpt.cpp), SFLOW
//  (sflow_f_ckpt.cpp).  The run part is here: the counters and times of the run, the counters of the
//  outputs, the rebalance (load map, the weighted partition it gives, checkpoint, the marker
//  REEF3D_Checkpoint/rebalance, exit code 3; tools/reef3d_rebalance.sh starts again).

#include<functional>
#include<string>
#include<vector>

class lexer;
class ghostcell;
class reefckpt;

class reefrun_ckpt
{
public:
    struct hooks
    {
        std::function<std::string()> unsupported;                        // what the checkpoints do not cover, empty: all
        std::function<void(reefckpt&)> write;                            // the arrays (level 0 and the patches)
        std::function<bool(reefckpt&, std::string&)> read;               // the same, false with a message: no restart
        std::function<void(std::vector<double>&)> printer_state;         // the counters of the outputs
        std::function<void(const std::vector<double>&)> printer_restore;
        std::function<void()> flush;                                     // the outputs before the stop of a rebalance
        std::function<void(lexer*)> loadmap;                             // the load map of the mesh refinement (G 45)
    };

    reefrun_ckpt(lexer*, ghostcell*, hooks);

    void write();
    void read();                // I 45 1, exits on a failure
    void step();                // after the outputs of a step: P 48 and G 46
    void rebalance_check();

    double rebal_L = 0.0;       // the load max/mean that stopped the run before (restart), 0 none
    int rebal_n = 0;            // rebalance restarts so far
    bool rebal_off = false;     // no further rebalance (the last one did not help)

private:
    lexer *p;
    ghostcell *pgc;
    hooks h;
};

#endif
