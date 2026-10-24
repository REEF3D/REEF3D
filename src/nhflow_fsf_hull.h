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

#ifndef NHFLOW_FSF_HULL_H_
#define NHFLOW_FSF_HULL_H_

#include"increment.h"
#include"slice4.h"
#include<vector>
#include<fstream>

class lexer;
class fdm_nhf;
class ghostcell;
class slice;
class sixdof;
class sixdof_obj_nhflow;
class ioflow;

using namespace std;

// Floating body as the top of the water columns (X 13 2, NHFLOW 6DOF X 10 1/2).
//
// In the columns under the hull (the footprint) the sigma grid runs from the bed to the hull bottom
// S(x,y,t): water only, no body cells, no direct forcing. The level of the hydrostatic pressure
// terms there is the piezometric head zeta, continuous with the surrounding free surface and
// determined by the depth-integrated continuity of the column (dS/dt + div Q = 0), solved
// implicitly after the momentum stage. The non-hydrostatic pressure keeps P = 0 at the top of
// every column, so the hull pressure is rho g (zeta - S). See nhflow_fsf_hull.cpp.
//
// Called by nhflow_momentum_RK3 once per RK stage:
//   level()     after the water level update of the stage (phase_F): footprint of the stage pose,
//               WL = S - bed, eta = S - bed - depth, detadt = dS/dt in the footprint; zeta
//   momentum()  before velcalc (phase_P1): UH, VH, WH of the columns whose depth level() changed
//   upgrad()    momentum stage (phase_M): hydrostatic source of the footprint cells
//   correct()   after velcalc (phase_P1), before the projection: implicit piezometric head
//   slope()     from nhflow_sigma::sigma_update: one-sided eta slopes at the footprint edge
//   forces()    from sixdof_obj_nhflow::hydrodynamic_forces_nhflow: the hull loads

class nhflow_fsf_hull : public increment
{
public:
    nhflow_fsf_hull(lexer*);
    virtual ~nhflow_fsf_hull();

    void level(lexer*, fdm_nhf*, ghostcell*, sixdof*, ioflow*, int, slice&, double);
    void momentum(lexer*, double*, double*, double*);
    void upgrad(lexer*, fdm_nhf*, slice&);
    void correct(lexer*, fdm_nhf*, ghostcell*, int, double, double*, double*, double*, slice&);
    void slope(lexer*, fdm_nhf*);
    bool edge_slope(lexer*, fdm_nhf*, int, int, double&, double&);
    void forces(lexer*, fdm_nhf*, ghostcell*, slice&, sixdof_obj_nhflow*, bool);
    void print(lexer*, fdm_nhf*, ghostcell*, slice&);

    bool active() const {return mode==1;}
    slice4& zeta_slice() {return zeta;}

private:
    void footprint(lexer*, fdm_nhf*, ghostcell*, slice&);
    void solve(lexer*, fdm_nhf*, ghostcell*, double, slice&);
    double gsum(ghostcell*, double);
    int gsum_i(ghostcell*, int);
    bool wall_x(lexer*, int) const;
    bool wall_y(lexer*, int) const;

    slice4 cover,S,St,zeta,dzeta,R,ratio,usk,vsk,wsk;
    vector<int> fi,fj;              // footprint columns on this rank
    vector<int> ei,ej;              // free columns next to the footprint
    vector<char> was;               // cover of the last stage

    int mode,gcval_eta;
    bool told,pending,fresh;
    int nfoot,nedge;                // global number of footprint columns, edge columns
    int itstep;                     // SOR sweeps since the last log line
    double omega;                   // SOR factor
    double rmax0,rmax1;             // continuity residual before and after the correction (max over the step)
    double dVstep,dVsum;            // volume replaced by the prescribed level: since the last line, over the run
    double Fz_p,Fx_p;               // last hull loads (pressure on the bottom, on the sides)
    double starttime,tstep,ttot;

    ofstream out;
};

#endif
