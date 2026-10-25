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
#include<fstream>

class lexer;
class ghostcell;
class slice;
class met_provider;

using namespace std;

/*--------------------------------------------------------------------
REEF3D atmospheric forcing (W 200) of one grid of SFLOW (A 10 2) or
NHFLOW (A 10 5): level 0, every refined patch and every block has its
own, created by the momentum class of the grid (as the wind forcing
wind_f) and pointed to by fdm2D::met / fdm_nhf::met for the output.

Fields of the cells (incl. ghost cells), evaluated from the provider
once per time step at the model time of the start of the step:
  U10, V10   wind 10 m above the sea along the model axes [m/s]
  Pa         air pressure at mean sea level [Pa]
  Tx, Ty     kinematic surface stress tau/rho_w = r (rho_a/rho_w) Cd |U10| U10
             [m^2/s^2], Cd = met_drag(W 207, |U10|), r the ramp W 204,
             rho_a = W 3, rho_w = W 1

Momentum sources (wet cells):
  SFLOW   d(UH)/dt += Tx                         (sflow_momentum_func, _SI)
  NHFLOW  d(UH)/dt += Tx/dsigma in the top layer (nhflow_momentum_RK2,
          RK3, SI; as A 573 1, without the forcing area)

The air pressure gradient, the inverse barometer at open boundaries and
the Coriolis force come with the next patch; Pa is evaluated and printed
already.

Output in the free surface VTP files (REEF3D_NHFLOW_VTP_FSF,
REEF3D_SFLOW_VTP_FSF): wind10 (vector), air_pressure, wind_stress
(vector, tau [N/m^2]).
--------------------------------------------------------------------*/

class met_forcing : public increment
{
public:
    met_forcing(lexer*, ghostcell*);
    virtual ~met_forcing() = default;

    // fields at the model time of the grid (once per time step)
    void update(lexer*, ghostcell*);

    // momentum sources, dir 0: x, 1: y
    void source_sflow(lexer*, ghostcell*, slice &F, int dir);
    void source_nhflow(lexer*, ghostcell*, double *F, int dir);

    // free surface VTP output
    void offset_vtp(lexer*, int*, int&);
    void name_vtp(ofstream&, int*, int&);
    void name_pvtp(ofstream&);
    void print_vtp(lexer*, ofstream&);

    slice4 U10, V10, Pa, Tx, Ty;

private:
    met_provider *pm;
    bool has_t = false;
    double tlast = 0.0;
};

#endif
