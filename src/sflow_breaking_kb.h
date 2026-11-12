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

#ifndef SFLOW_BREAKING_KB_H_
#define SFLOW_BREAKING_KB_H_

#include"increment.h"
#include"slice4.h"
#include"sliceint4.h"

class lexer;
class fdm2D;
class ghostcell;
class ioflow;
class slice;
class breaking_crest3d;

using namespace std;

// SFLOW A 246 4: the breaking model of NHFLOW A 550 4 for the depth-averaged equations.
//
// Onset from crest tracking (breaking_crest3d, Barthelemy et al. 2018: B = u_s/C > A 359, end below the front-face
// slope A 360) with u_s = (U, V), the horizontal velocity of the quadratic pressure is uniform over the depth. The
// tracked crests seed the eddy viscosity of Kennedy et al. (2000) (A 501-504, A 554), continuous in deta/dt, with the
// age of the event carried with the front, the length L_k = tanh(k h)/k of the breaking wave and the band of A 502 L_k
// behind the crest (nhflow_breaking::breaking_kennedy, A 550 4). The viscosity b->vb acts through the diffusion of
// the momentum (A 212 1 explicit, A 212 2 implicit); the non-hydrostatic pressure stays on in the breaking cells
// (switching it off there radiates spurious waves). b->breaking stays 0 (it switches the pressure off), the
// breaking columns are in b->breaking_print. No breaking in the relaxation zones. Not with mesh refinement.
class sflow_breaking_kb : public increment
{
public:
    sflow_breaking_kb(lexer*);
    virtual ~sflow_breaking_kb();

    // relaxation zones (once)
    void zones(lexer*, ioflow*);

    // eta: water level of the stage, eta_n: the level of the start of the step, (eta - eta_n)/(alpha dt) = deta/dt
    void start(lexer*, fdm2D*, ghostcell*, slice &eta, slice &eta_n, double alpha);

private:
    void tracking(lexer*, fdm2D*, ghostcell*, slice &eta);
    void viscosity(lexer*, fdm2D*, ghostcell*, slice &eta, slice &eta_n, double alpha);

    breaking_crest3d *bt = nullptr;
    int bt_count = -1, count_n = -1;
    bool zones_built = false;

    slice4 brkage, brkage_n, brkk, brkk_n, brkr, brkr_n, brkkseed, us, vs, vbs, vbk, vbp;
    sliceint4 brkseed, brkflag, brkzone, bd;
};

#endif
