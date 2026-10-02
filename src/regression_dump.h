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

#ifndef REGRESSION_DUMP_H_
#define REGRESSION_DUMP_H_

#include"increment.h"
#include<fstream>
#include<string>
#include<vector>

class lexer;
class fdm;
class ghostcell;
class turbulence;
class concentration;
class field;

/*--------------------------------------------------------------------
regression_dump: exact (double precision) state output for the
regression test suite in tests/regression.

Inactive unless the environment variable REEF3D_REGRESSION_DIR is set,
so it never changes results or normal output.

  REEF3D_REGRESSION_DIR=<dir>   activate, write files into <dir>
  REEF3D_REGRESSION_EVERY=<n>   also dump the full state every n steps
                                (default 0: only initial + final state)

Files, one set per MPI rank r:
  steps_r<r>.txt          one line per time step: count, simtime, dt and
                          rank-local L2 sums of the main fields, in
                          C99 hexfloat (exact)
  state_<count>_r<r>.bin  full state at a step, see write_state()
--------------------------------------------------------------------*/

class regression_dump : public increment
{
public:
    regression_dump(lexer*);
    virtual ~regression_dump();

    bool active() const {return is_active;}

    // CFD
    void cfd_ini(lexer*, fdm*, ghostcell*, turbulence*, concentration*);     // after initialisation
    void cfd_step(lexer*, fdm*, ghostcell*, turbulence*, concentration*);    // end of each time step
    void cfd_final(lexer*, fdm*, ghostcell*, turbulence*, concentration*);   // after the main loop

private:
    void cfd_state(lexer*, fdm*, turbulence*, concentration*);
    void cfd_collect(lexer*, fdm*, turbulence*, concentration*);
    void add(const char*, int);
    void write_state(lexer*);

    bool is_active;
    int every;
    int last_written;
    std::string dir;
    std::ofstream steplog;

    std::vector<std::string> names;
    std::vector<std::vector<double>> data;
};

#endif
