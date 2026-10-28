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
#include<vector>

class lexer;
class ghostcell;
class slice;
class seastate_f;
class seastate_source;
class met_forcing;

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
                                                  nonlinear transfers;
                                                  with the roller the
                                                  roller dissipation
                                                  instead of breaking;
                                                  the wind input only
                                                  when the host gets
                                                  the full wind stress,
                                                  W 200 or A 570: minus
                                                  the wave-supported
                                                  stress)
                 + g sum N k sig/sinh(2kd) dd/dx_i (bottom slope)
                 - d/dx_i [g sum E (n - 1/2)]      (wave pressure)
                 + J_i - u_i div(M)                (vortex force)
             J = chi (M_y, -M_x),  chi = dv/dx - du/dy
           div(M) in flux form with no Stokes transport through the
           domain edges, walls and the shoreline.

  M = g sum E k/sig + 2 g sum R/c (cos, sin)   wave and roller mass flux
  E = sig N dsig dtheta (m^2, per rho g), n = cg/c, S_N = P - D N of the
  source terms (seastate_source, wind input as above), R the roller
  energy (seastate_roller, A 748). Spatial derivatives: central
  differences between wave-active cells, one-sided next to dry cells
  and at the domain edge.

Remote waves (A 750 2): SEASTATE runs as its own model (A 10 7, A 750 2)
on its own grid, in the same MPI job as the host (reef3d_mpmd, e.g.
mpirun -np 4 -wdir surge REEF3D : -np 2 -wdir waves REEF3D), the same
model frame and clock. Every A 706 s of host time the host sends eta,
WL, U, V (global fields, root to root, seastate_remote) and receives the
wave force, div(M), Hs, Tp and the wet cells of the wave grid:
  waves -> host   bilinear between the wet wave cell centres
  host -> waves   the area mean over the host cells inside the wave cell (eta,
                  U, V over the wet ones, WL over all), bilinear where the
                  host grid is the coarser one
The waves lag one interval: the forces the host gets at t_k are those
of the wave state at t_k, which was computed with the host fields of
t_(k-1); host and waves compute in parallel. The host ends the run
(stop message at ghostcell::final).

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

    // surfbeat runs (not with remote waves)
    bool surfbeat() const;

    // remote waves (A 750 2), host side: the stop message at the end of the run (rank 0, not collective)
    void remote_stop(lexer*);

    // sea state on the host grid (Hs [m], Tp [s]): in-process copies or the remote fields (the wave-dependent
    // drag W 207 7)
    slice4 Hsw, Tpw;

protected:
    // wave model set-up; the host then calls forces
    void ini_wave(lexer*, ghostcell*, const char *host);

    // remote waves (A 750 2), host side: exchange of the grids and the first fields; then every A 706 s
    void remote_ini(lexer*, ghostcell*, slice &eta, slice &WL, slice &U, slice &V);
    void remote_step(lexer*, ghostcell*, slice &eta, slice &WL, slice &U, slice &V, bool force);

    // remote waves: global cell centres of this model, the fields of this model gathered on its rank 0
    void centres(lexer*, ghostcell*, std::vector<double> &xs, std::vector<double> &ys);
    void gather(lexer*, ghostcell*, const std::vector<slice*> &f, std::vector<double> &out, double flag);
    void exchange(lexer*, ghostcell*, std::vector<double> &out, std::vector<double> &in, int tag_out, int tag_in);
    void bilinear(const std::vector<double> &xs, const std::vector<double> &ys, double x, double y,
                  int &a, int &b, double &wx, double &wy) const;
    bool covers(const std::vector<double> &xs, const std::vector<double> &ys, double x, double y) const;

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
    seastate_source *pnet;      // net source terms of the dissipation force (A 751 2)

    slice4 Fw,Gw,divM;          // wave force [m^2/s^2] and div(M) [m/s]
    slice4 Sxx,Sxy,Syy;         // radiation stress / (rho) [m^3/s^2]
    slice4 Mx,My;               // wave (and roller) mass flux [m^2/s]
    slice4 Fdx,Fdy,Q,B;         // vortex force: dissipation force, wave pressure, bottom-slope coefficient

    double tlast, tnext, tprint;
    int nstep;

    // remote waves (A 750 2)
    bool remote;                            // this object talks to the other app
    bool host_wind;                         // the host gets a wind stress (the wind input of the dissipation force)
    std::vector<double> rx, ry;             // global cell centres of the other model
    std::vector<double> sbuf, rbuf;         // send / receive buffers (header: flag, simtime)
    double t_wait;                          // rank 0: time spent waiting for the other app [s]

    met_forcing *pmetw = nullptr;           // wave-dependent drag W 207 7: told after every change of Hsw, Tpw
};

#endif
