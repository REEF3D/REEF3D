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


#ifndef SEASTATE_COUPLING_H_
#define SEASTATE_COUPLING_H_

#include"increment.h"
#include"slice4.h"

class lexer;
class ghostcell;
class slice;
class seastate_f;
class seastate_source;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::SEASTATE coupling with a depth-averaged or layered host model,
host-independent part (A 750 1). The hosts are SFLOW (seastate_sflow,
A 10 2) and NHFLOW (seastate_nhflow, A 10 5): they own the 2D grid, the
time step and the time loop, write the environment from their fields
and apply the wave force, the Stokes mass source and the long-wave
boundary.

Wave steps: every A 706 seconds (the coupling interval, the wave
forcing is held constant in between); with surfbeat (A 770 1) every
host step, because the wave groups and the long waves interact on the
host's time scale.

  host -> SEASTATE (A 753)  water depth eta + still-water depth (dries
                            and re-wets cells of the initial wave
                            domain), dd/dt, depth-averaged Eulerian
                            current:
                              A 751 1: U - M/h  (the host velocity is
                                       the mass transport velocity)
                              A 751 2: U        (Eulerian)
  SEASTATE -> host          wave force per unit area and density
                            [m^2/s^2], ramped up over A 752 seconds:

  A 751 1  radiation stress (Longuet-Higgins and Stewart 1964; Phillips
           1977), host velocity = Lagrangian (mass transport) velocity:
             F_i = -dS_ij/dx_j
             S_ij = g sum E [ n k_i k_j/k^2 + (n - 1/2) delta_ij ]
                  + g sum R k_i k_j/k^2                     (roller)

  A 751 2  vortex force (McWilliams et al. 2004, Smith 2006, depth-
           averaged), host velocity = Eulerian velocity u, the Stokes
           transport M enters the continuity:
             dWL/dt + div(WL u) = -div(M)
             F_i = -g sum S_N k_i                 (dissipation and
                                                  nonlinear transfers,
                                                  without wind input;
                                                  with the roller the
                                                  roller dissipation
                                                  instead of breaking)
                 + g sum N k sig/sinh(2kd) dd/dx_i (bottom slope)
                 - d/dx_i [g sum E (n - 1/2)]      (wave pressure)
                 + J_i - u_i div(M)                (vortex force)
             J = chi (M_y, -M_x),  chi = dv/dx - du/dy
           div(M) in flux form with no Stokes transport through the
           domain edges, walls and the shoreline.

  M = g sum E k/sig + 2 g sum R/c (cos, sin)   wave and roller mass flux
  E = sig N dsig dtheta (m^2, per rho g), n = cg/c, S_N = P - D N of the
  source terms (seastate_source, wind input switched off), R the roller
  energy (seastate_roller, A 748). Spatial derivatives: central
  differences between wave-active cells, one-sided next to dry cells
  and at the domain edge.

Long waves at the x- side (surfbeat, A 770 1): absorbing-generating
boundary (Van Dongeren and Svendsen 1997, as XBeach), the host's ghost
cells get
   eta_g = eta_1,  u_g = q_bx/h - sqrt(g/h) (eta_1 - zeta_b),  v_g = v_1
with zeta_b, q_b the incoming bound long wave of the row (A 774 1, zero
with A 774 0) and eta_1, v_1, h the first host cell: the incoming wave
enters, everything else leaves through the side.
--------------------------------------------------------------------*/

class seastate_coupling : public increment
{
public:
    seastate_coupling(lexer*);
    virtual ~seastate_coupling();

    seastate_f *wave() {return pwave;}

protected:
    // wave model set-up; the host then calls forces
    void ini_wave(lexer*, ghostcell*, const char *host);

    // a wave step is due at the host's simtime; dt: length of the wave step
    bool due(lexer*, double &dt) const;

    // environment from the host fields (2D, depth-averaged velocities)
    void environment(lexer*, slice &eta, slice &WL, slice &U, slice &V);

    // wave step of length dt, then the host computes the forces
    void wave_step(lexer*, ghostcell*, double dt);

    // wave force, div(M) and the integrals; U, V: host velocity (vortex force)
    void forces(lexer*, ghostcell*, slice &U, slice &V);

    // VTP output after a wave step
    void print(lexer*, ghostcell*);

    // long-wave ghost cell values of row j at the x- side; false: no long-wave boundary
    bool longwave_bc(lexer*, int j, double h, double eta1, double v1, double &etag, double &ug, double &vg);

    double ramp(lexer*) const;
    double ddx(lexer*, slice&, int, int);
    double ddy(lexer*, slice&, int, int);

    seastate_f *pwave;
    seastate_source *pnet;      // net source terms without wind input (A 751 2)

    slice4 Fw,Gw,divM;          // wave force [m^2/s^2] and div(M) [m/s]
    slice4 Sxx,Sxy,Syy;         // radiation stress / (rho) [m^3/s^2]
    slice4 Mx,My;               // wave (and roller) mass flux [m^2/s]
    slice4 Fdx,Fdy,Q,B;         // vortex force: dissipation force, wave pressure, bottom-slope coefficient

    double tlast, tnext, tprint;
    int nstep;
};

#endif
