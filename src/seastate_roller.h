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


#ifndef SEASTATE_ROLLER_H_
#define SEASTATE_ROLLER_H_

#include"increment.h"
#include"slice4.h"
#include<vector>

class lexer;
class ghostcell;
class fdm_seastate;
class seastate_store;
class seastate_exchange;
class seastate_source;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::SEASTATE surfbeat - roller energy balance (A 748 1)

Per direction theta, in the variance units of the wave energy
(E(theta) = sig N(theta) [m^2/rad], roller energy per rho g):

  dR/dt + d((c cos(theta) + U) R)/dx + d((c sin(theta) + V) R)/dy
        = D_w(theta) - D_r(theta)

  D_w = (D/E) E(theta)      breaking dissipation of the waves (Roelvink,
                            A 740 2; seastate_source::brk_rate), the
                            source of the roller
  D_r = 2 g beta R(theta)/c roller dissipation (Reniers et al. 2004;
                            XBeach), beta = A 749, c = sig/k

First-order upwind, implicit (BDF1): an M-matrix as the wave action
solver, so R >= 0 for any time step. Gauss-Seidel sweeps in the
downstream order of the quadrant of each direction (exact upwinding
within a rank), halo exchange after each sweep. No roller enters
through the domain sides; land, dry cells and the sides absorb.
Refraction of the roller is neglected (the roller lives where the
waves break, over a short distance).

Integrals per cell: R = sum R dtheta [m^2], Dr = sum D_r dtheta and
Dw = sum D_w dtheta [m^2/s]. The coupling (seastate_coupling) adds the
roller to the radiation stress, g sum R (cos^2, cos sin, sin^2) dtheta,
and to the mass flux, 2 g sum R/c (cos, sin) dtheta, and with the
vortex force the roller dissipation replaces the breaking dissipation
in the dissipation force.
--------------------------------------------------------------------*/

class seastate_roller : public increment
{
public:
    seastate_roller(lexer*, fdm_seastate*, double beta);
    virtual ~seastate_roller();

    // one step of length dt after the wave step; src: the wave model's source terms (breaking)
    void step(lexer*, ghostcell*, fdm_seastate*, seastate_source *src, double dt, int iterations);

    // roller energy per direction of cell (i,j) [m^2/rad], nullptr: no storage
    const float *spec(int i, int j) const;

    // phase speed of cell (i,j)
    static double celerity(fdm_seastate*, int i, int j);

    slice4 R, Dr, Dw;
    const double beta;

private:
    void sweep(lexer*, fdm_seastate*, int q, double rdt);

    seastate_store *Rt, *R0, *Sw;
    seastate_exchange *pex;
    int ndir;
    vector<int> m0, m1;
    vector<double> P, D;
};

#endif
