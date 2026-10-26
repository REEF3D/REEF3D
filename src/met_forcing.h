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

#ifndef MET_FORCING_H_
#define MET_FORCING_H_

#include"increment.h"
#include"slice4.h"
#include"lexer.h"
#include<fstream>

class ghostcell;
class slice;
class met_provider;

using namespace std;

/*--------------------------------------------------------------------
REEF3D atmospheric forcing (W 200) and Coriolis force (W 211) of one grid
of SFLOW (A 10 2) or NHFLOW (A 10 5): level 0, every refined patch and
every block has its own, created by the momentum class of the grid (as
the wind forcing wind_f) and pointed to by fdm2D::met / fdm_nhf::met for
the output.

Fields of the cells (incl. ghost cells), evaluated from the provider
once per time step at the model time of the start of the step:
  U10, V10   wind 10 m above the sea along the model axes [m/s]
  Pa         air pressure at mean sea level [Pa]
  Tx, Ty     kinematic surface stress tau/rho_w = r (rho_a/rho_w) Cd |U10| U10
             [m^2/s^2], Cd = met_drag(W 207, |U10|), r the ramp W 204,
             rho_a = W 3, rho_w = W 1
  Gx, Gy     kinematic air pressure gradient r (dPa/dx, dPa/dy)/rho_w (W 208 1),
             central differences of the cell values
  Fc         Coriolis parameter f [1/s] (W 211; constant in time)

Momentum sources of the conserved variables (wet cells):
  SFLOW   d(UH)/dt += Tx - WL Gx + f VH,   d(VH)/dt += Ty - WL Gy - f UH
          (sflow_momentum_func::stage, the predictor of sflow_momentum_SI)
  NHFLOW  d(UH)/dt += Tx/dsigma (top layer) - WL Gx + f VH   (every layer)
          d(VH)/dt += Ty/dsigma (top layer) - WL Gy - f UH
          (nhflow_momentum_RK2, RK3, the predictor of nhflow_momentum_SI)
with the depth WL and the momentum UH, VH of the stage. The Coriolis
force is explicit in the Runge-Kutta stages (and in the Runge-Kutta
predictor of the semi-implicit schemes): SSP-RK3 is stable for
f dt < 1.7 and damps an inertial oscillation by (f dt)^4/24 per step,
SSP-RK2 amplifies it by (f dt)^4/8 per step; with f = 1e-4 1/s and
dt = 600 s that is 5e-7 and 2e-6 per step.

The inverse barometer of the open edges is in met_provider::ib (iowave).

Output in the free surface VTP files (REEF3D_NHFLOW_VTP_FSF,
REEF3D_SFLOW_VTP_FSF) with W 200 > 0: wind10 (vector), air_pressure,
wind_stress (vector, tau [N/m^2]).
--------------------------------------------------------------------*/

class met_forcing : public increment
{
public:
    met_forcing(lexer*, ghostcell*);
    virtual ~met_forcing() = default;

    // fields at the model time of the grid (once per time step)
    void update(lexer*, ghostcell*);

    // momentum sources, dir 0: x, 1: y; WL, UH, VH of the stage
    void source_sflow(lexer*, ghostcell*, slice &F, int dir, slice &WL, slice &UH, slice &VH);
    void source_nhflow(lexer*, ghostcell*, double *F, int dir, slice &WL, double *UH, double *VH);

    // free surface VTP output (W 200 > 0)
    void offset_vtp(lexer*, int*, int&);
    void name_vtp(ofstream&, int*, int&);
    void name_pvtp(ofstream&);
    void print_vtp(lexer*, ofstream&);

    static bool active(lexer *p) {return p->W200>0 || p->W211>0;}

    slice4 U10, V10, Pa, Tx, Ty, Gx, Gy, Fc;

    bool wind = false, pressure = false, coriolis = false;

private:
    met_provider *pm = nullptr;
    bool has_t = false;
    double tlast = 0.0;
};

#endif
