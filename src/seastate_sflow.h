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

#ifndef SEASTATE_SFLOW_H_
#define SEASTATE_SFLOW_H_

#include"increment.h"
#include"slice4.h"

class lexer;
class fdm2D;
class ghostcell;
class slice;
class seastate_f;
class seastate_source;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::SEASTATE <-> SFLOW coupling (A 10 2, A 750 1)

SFLOW is the host: it owns the 2D grid, the time step and the time
loop. SEASTATE runs on the same grid and is advanced every A 706
seconds (the coupling interval, the wave time step); in between the
wave forcing is held constant.

  SFLOW -> SEASTATE (A 753)  water depth eta + still-water depth
                             (dries and re-wets cells of the initial
                             wave domain), depth change dd/dt, and the
                             depth-averaged Eulerian current:
                               A 751 1: U - M_St/h  (SFLOW velocity is
                                        the mass transport velocity)
                               A 751 2: U           (Eulerian)
  SEASTATE -> SFLOW          wave force per unit area and density
                             [m^2/s^2], added to the SFLOW momentum
                             right-hand side in every RK stage, ramped
                             up over A 752 seconds:

  A 751 1  radiation stress (Longuet-Higgins and Stewart 1964; Phillips
           1977), SFLOW velocity = Lagrangian (mass transport) velocity:
             F_i = -dS_ij/dx_j
             S_ij = g sum E [ n k_i k_j/k^2 + (n - 1/2) delta_ij ]

  A 751 2  vortex force (McWilliams et al. 2004, Smith 2006, depth-
           averaged), SFLOW velocity = Eulerian velocity u, the Stokes
           transport M_St = g sum E k/sig enters the continuity:
             dWL/dt + div(WL u) = -div(M_St)
             F_i = -g sum S_N k_i                 (dissipation and
                                                  nonlinear transfers,
                                                  without wind input)
                 + g sum N k sig/sinh(2kd) dd/dx_i (bottom slope)
                 - d/dx_i [g sum E (n - 1/2)]      (wave pressure)
                 + J_i - u_i div(M_St)             (vortex force)
             J = chi (M_St,y, -M_St,x),  chi = dv/dx - du/dy
           div(M_St) in flux form with no Stokes transport through the
           domain edges, walls and the shoreline (the mass source
           integrates to zero: no mass enters with the waves)
           For a steady wave field without currents this equals the
           radiation stress force exactly (wave action balance and
           ray equations); the rotational part comes only from the
           dissipation, so non-breaking waves over a sloping bed do
           not force spurious circulation.

  E = sig N dsig dtheta (m^2, per rho g), n = cg/c, S_N = P - D N of
  the source terms (seastate_source, wind input switched off: the
  momentum the waves gain from the wind is not taken from the current).
  Spatial derivatives: central differences between wave-active cells,
  one-sided next to dry cells and at the domain edge.
--------------------------------------------------------------------*/

class seastate_sflow : public increment
{
public:
    seastate_sflow(lexer*, fdm2D*, ghostcell*);
    virtual ~seastate_sflow();

    void ini(lexer*, fdm2D*, ghostcell*);

    // every SFLOW step, before the momentum step: wave step when the coupling interval is over
    void start(lexer*, fdm2D*, ghostcell*);

    // SFLOW momentum right-hand side (every RK stage) and continuity source (A 751 2)
    void u_source(lexer*, fdm2D*);
    void v_source(lexer*, fdm2D*);
    void mass_source(lexer*, fdm2D*, slice&);

private:
    void environment(lexer*, fdm2D*, ghostcell*);
    void forces(lexer*, fdm2D*, ghostcell*);
    double ramp(lexer*) const;
    double ddx(lexer*, slice&, int, int);
    double ddy(lexer*, slice&, int, int);

    seastate_f *pwave;
    seastate_source *pnet;      // net source terms without wind input (A 751 2)

    slice4 Fw,Gw,divM;          // wave force [m^2/s^2] and div(M_St) [m/s]
    slice4 Sxx,Sxy,Syy;         // radiation stress / (rho) [m^3/s^2]
    slice4 Mx,My;               // Stokes transport [m^2/s]
    slice4 Fdx,Fdy,Q,B;         // vortex force: dissipation force, wave pressure, bottom-slope coefficient
    slice4 uw,vw;               // current passed to SEASTATE

    double tlast, tnext;
    int nstep;
};

#endif
