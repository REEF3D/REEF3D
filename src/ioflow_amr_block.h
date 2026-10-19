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

#ifndef IOFLOW_AMR_BLOCK_H_
#define IOFLOW_AMR_BLOCK_H_

#include"ioflow_void.h"

// in- and outflow of a block 1.. of a multi-block grid (G 50 1, its level-0 grid of REEFAMR) through its open
// faces (G 54), as ioflow_f on the native grid with B 60 1: the discharge Q through the -i face (the velocity
// Q / wet area along the face normal in the ghost cells; NHFLOW: B 61 1 plain, else the log profile scaled
// to Q) and the water level at the +i face (start ramp F 64 from F 60).  A block lies on one rank: the sums
// over the face are local.  Everything else as ioflow_v.
class ioflow_amr_block : public ioflow_v
{
public:
    ioflow_amr_block(lexer*,ghostcell*,patchBC_interface*);
    virtual ~ioflow_amr_block();

    int in = 0;             // the -i face is open: discharge Q (m3/s)
    double Q = 0.0;
    int out = 0;            // the +i face is open: water level wsf (m)
    double wsf = 0.0;

    // SFLOW
    void discharge2D(lexer*,fdm2D*,ghostcell*) override;
    void waterlevel2D(lexer*,fdm2D*,ghostcell*,slice&) override;
    void inflow2D(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,slice&) override;

    // NHFLOW
    void discharge_nhflow(lexer*,fdm_nhf*,ghostcell*) override;
    void inflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*,slice&) override;
    void rkinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*,slice&) override;
    void rkinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*) override;
    void fsfinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,slice&) override;

private:
    double level(lexer*) const;
    void qin_nhflow(lexer*,fdm_nhf*,double&,double&);
    void curv_dir(lexer*,double*,double*,double*,double*);

    double Ui = 0.0;        // inflow velocity normal to the -i face (m/s)
    double Ai = 0.0;        // wet area of the -i face (m2)
};

#endif
