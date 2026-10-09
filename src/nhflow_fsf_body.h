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

#ifndef NHFLOW_FSF_BODY_H_
#define NHFLOW_FSF_BODY_H_

#include"increment.h"
#include"slice4.h"
#include<vector>
#include<fstream>

class lexer;
class fdm_nhf;
class ghostcell;
class slice;

using namespace std;

// Water level in the columns pierced by a floating body (X 13 1, NHFLOW direct forcing X 10 1/2).
// See nhflow_fsf_body.cpp. Called by nhflow_momentum_RK2/RK3 once per RK stage:
//   level()     after the water level update of the stage (phase_F, before omega_update)
//   momentum()  before velcalc of the same stage (phase_P1): UH, VH, WH of the changed columns
//   print()     last stage, after the reforcing (phase_P2): log REEF3D_NHFLOW_6DOF/REEF3D_6DOF_pierced.dat (X 19)

class nhflow_fsf_body : public increment
{
public:
    nhflow_fsf_body(lexer*);
    virtual ~nhflow_fsf_body();

    void level(lexer*, fdm_nhf*, ghostcell*, slice&, double);
    void momentum(lexer*, double*, double*, double*);
    void print(lexer*, fdm_nhf*, ghostcell*, slice&);

    // mesh refinement patches (nhflow_amr): the level is only treated on level 0
    void off() {mode=0; patch=true;}

private:
    void solve(lexer*, fdm_nhf*, ghostcell*);

    slice4 Hb,Hbold,etaext,ratio;
    vector<int> fi,fj;              // interior columns with Hb > 0 on this rank

    int mode,gcval_eta;
    bool patch,told;
    int nfoot;                      // global number of columns with Hb > 0 (last stage)
    int itstep;                     // SOR sweeps since the last log line
    bool pending;                   // level() changed WL, momentum() not yet applied
    double omega;                   // SOR factor
    double dVcore,dVedge;           // level change x area since the last log line, Hb = 1 and 0 < Hb < 1
    double dVedge_sum;              // 0 < Hb < 1 summed over the run (upper bound of the change of real water)
    double tstep,ttot;              // wall time since the last log line, of the run
    double starttime;

    ofstream out;
};

#endif
