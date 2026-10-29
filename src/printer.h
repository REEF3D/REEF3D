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

#ifndef PRINTER_H_
#define PRINTER_H_

#include<vector>
#include<cstddef>

class lexer;
class fdm;
class fdm_nhf;
class fdm_fnpf;
class ghostcell;
class turbulence;
class nhflow_turbulence;
class heat;
class ioflow;
class expdata;
class concentration;
class multiphase;
class sediment;
class amr_grid_lookup;
class amr_column_lookup;

class printer
{
public:
    // checkpoints (P 48, I 45): the counters of the outputs (file numbers, print windows)
    virtual void ckpt_state(std::vector<double>&) {};
    virtual void ckpt_restore(const std::vector<double>&) {};
    // the files written so far are complete (a rebalance stop, G 46)
    virtual void flush_outputs() {};
    // CFD
    virtual void start(lexer*,fdm*,ghostcell*,turbulence*,heat*,ioflow*,expdata*,concentration*,multiphase*,sediment*) {};
    virtual void print_stop(lexer*,fdm*,ghostcell*,turbulence*,heat*,ioflow*,expdata*,concentration*,multiphase*,sediment*) {};
    // the end of a run (N 45, N 41, S 19): the last state, unless the last step wrote it
    virtual void print_end(lexer*,fdm*,ghostcell*,turbulence*,heat*,ioflow*,expdata*,concentration*,multiphase*,sediment*) {};
    // CFD mesh refinement: the probes and wave gauges read the finest grid at their position
    virtual void amr_probes(const amr_grid_lookup*) {};
    // the 2D-column models (NHFLOW, FNPF): the gauges and probes read the finest grid
    virtual void amr_columns(lexer*, const amr_column_lookup*) {};
    // NHFLOW
    virtual void start(lexer*,fdm_nhf*,ghostcell*,ioflow*,nhflow_turbulence*,sediment*) {};
    virtual void print_stop(lexer*,fdm_nhf*,ghostcell*,ioflow*,nhflow_turbulence*,sediment*) {};
    virtual void print_end(lexer*,fdm_nhf*,ghostcell*,ioflow*,nhflow_turbulence*,sediment*) {};
    // FNPF
    virtual void start(lexer*,fdm_fnpf*,ghostcell*,ioflow*) {};
    virtual void print_stop(lexer*,fdm_fnpf*,ghostcell*) {};
    virtual void print_end(lexer*,fdm_fnpf*,ghostcell*) {};
protected:
    void writeFile(const char*, const size_t);
    
    // the volume output is rare (P 30): its buffer (the whole 3D grid with connectivity, about
    // 150 bytes per node) is given back after each output instead of being kept for the run
    void release_buffer() {std::vector<char>().swap(buffer);}

    std::vector<char> buffer;
};

#endif
