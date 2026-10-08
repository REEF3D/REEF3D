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
Author: Ahmet Soydan
--------------------------------------------------------------------*/

#ifndef ACOUSTICS_FWH_H_
#define ACOUSTICS_FWH_H_

#include"acoustics.h"
#include"increment.h"
#include<vector>
#include<fstream>

class fwh_permeable;

//  Permeable-surface FW-H for REEF3D::CFD (U 10 1), see acoustics_fwh_kernel.h for the formulation.
//
//  The permeable surface is the box U 20, snapped to the nearest grid nodes, so that the panels
//  are cell faces: the normal velocity is the staggered velocity on the face, the pressure and
//  the tangential velocities are averaged from the neighbour cells. Each face belongs to the rank
//  of the cell above it (half-open node range), so no panel is counted twice. The face U 21 can
//  be left open (end cap where the wake leaves the box).
//
//  p' is the pressure without the hydrostatic part rho0 g.(x - x_fs), rho0 = W 1. The observers
//  U 30 get image points across z = U 41 with weight -1 for U 40 1.
//
//  Output: REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-<n>.dat, time and p' of the observer
//  times that are complete (all panels and both neighbour samples have contributed), appended
//  every time step.
//  REEF3D-CFD-FWH-Box.dat: the box integrals per time step, pressure part Lp_i = int p' n_i dS,
//  momentum flux Lm_i = int rho u_i u_n dS, viscous traction V_i = int tau_ij n_j dS and
//  M_i = int x_i rho u_n dS. For incompressible flow and bodies at rest inside the box, the
//  momentum balance of the box gives the force of the fluid on the bodies
//  F = -(Lp + Lm - V + dM/dt); the FW-H integral (viscous stresses neglected) radiates the dipole
//  of -(Lp + Lm + dM/dt).

class acoustics_fwh final : public acoustics, public increment
{
public:
    acoustics_fwh(lexer*, fdm*, ghostcell*);
    virtual ~acoustics_fwh();
    
    void start(lexer*, fdm*, ghostcell*) override final;
    
private:
    struct face
    {
        int ii,jj,kk;   // cell on the low side, the face is its high side in direction dir
        int dir;
        double sgn;     // normal +-1 along dir, out of the box
        double x[3];
        double dS;
    };
    
    void ini(lexer*, fdm*, ghostcell*);
    double snap(ghostcell*, const double*, int, double);
    void box_faces(lexer*, fdm*, ghostcell*);
    void sample(lexer*, fdm*);
    void write(lexer*, ghostcell*);
    void box_integrals(lexer*, fdm*, ghostcell*);
    
    fwh_permeable *pfwh;
    std::vector<face> fc;
    std::vector<double> pv, uv;
    std::vector<double> dmin, dmax;     // delay bounds over all ranks, per observer
    std::vector<long> kwrite;           // next observer sample to write, per observer
    std::vector<double> buf;
    std::ofstream *pout;
    std::ofstream bout;
    
    double box[6];
    double zfs, rho0;
    double tprev;                       // first source time, for the default observer time step
    int nobs;
};

#endif
