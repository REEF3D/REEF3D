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

#ifndef SEASTATE_F_H_
#define SEASTATE_F_H_

#include"seastate.h"
#include"increment.h"
#include<fstream>
#include<vector>

class fdm_seastate;
class seastate_vtp;
class seastate_exchange;
class seastate_implicit;
class seastate_source;
class seastate_store;
class seastate_grid;
class seastate_surfbeat;
class seastate_roller;
class seastate_spc_series;
class seastate_wind_series;
class slice4;
class regression_dump;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - spectral (phase-averaged) wave model, A 10 7

2D horizontal grid of SFLOW (DIVEMesh), spectral grid (A 701-703),
block-sparse float32 action density (A 704), active cells from the
bathymetry (A 705), parametric initial spectrum (A 710), integrated wave
parameters, VTP output (P 181/P 182), integral log, regression dump.

Phase 1: implicit transport in x, y, sigma, theta (seastate_implicit),
nonstationary (A 700 1, time step A 706) or stationary (A 700 2),
iterations A 707, stationary convergence A 708, boundary spectrum
A 711 on the sides A 712, refraction A 713, frequency shift A 714,
halo exchange of the spectra (seastate_exchange).

Phase 2: source terms (seastate_source, A 730-745): wind input,
whitecapping and quadruplets (Komen, DIA), depth-induced breaking
(Battjes-Janssen), bottom friction (JONSWAP), triads (LTA); zero-
gradient sides (A 712 2). Stationary runs with the deep-water physics
iterate in pseudo time with the time step A 706.

Phase 3: coupling and handover.
  ini_coupled   set-up inside a host model (seastate_coupling: SFLOW
                seastate_sflow, NHFLOW seastate_nhflow): the
                host's time step, step counter and print counters are
                not touched, VTP files get their own numbering
  step_coupled  one wave step of length dt after the host has written
                the environment (depth, eta, U, V): cells of the initial
                active set dry and re-wet with the host's water depth
                (A 705), kinematics with dd/dt, transport, parameters,
                integral log; the host prints (printer, add_field)
  handover      2D spectra E(omega,theta) at the points A 760 for the
                irregular wave generation of FNPF/NHFLOW (B 85 11,
                spectrum-file-2d.dat format), sector A 761 around the
                mean direction, written every step

Phase 4: surfbeat (A 770 1), the wave-group action balance of XBeach
type (Roelvink 1993; Reniers et al. 2004) on a single representative
frequency (seastate_grid(f_rep, ndir)), nonstationary on the time step
of the host model.
  surfbeat_input     the boundary spectrum on the multi-frequency
                     grid A 701-703 (A 711 1 parametric, 2 SWAN file),
                     f_rep = 1/A 771 or 1/Tm-1,0 of that spectrum
  surfbeat_ini       wave-group generator (seastate_surfbeat): energy
                     envelope E(t,y) and bound long wave at the x- side
  surfbeat_boundary  x- boundary rows N(theta) = E(t,y) Dbar(theta)/sig
                     for the implicit solver at the new time level
  surfbeat_cap       H = sqrt(8E) <= A 747 h
  longwave           incoming long wave of a boundary row for the host
                     (zeta_b, q_b), A 774
  roller             roller energy balance (seastate_roller, A 748),
                     source: the breaking dissipation (A 740 2)

Phase 5a: external forcing (seastate_f_forcing.cpp).
  forcing_ini      opens the forcing files: seastate-boundary.spc with
                   spectra at many locations and times (A 711 3),
                   seastate-wind.dat with the wind field (A 730 2);
                   model time 0 = A 780 (YYYYMMDD.HHMMSS) or the first
                   time in the files
  boundary_series  boundary spectrum of every boundary cell of the sides
                   with A 712 1: the two nearest locations, weights by
                   inverse distance (= linear interpolation between two
                   locations on a straight side)
  forcing_update   spectra and wind at the new time level, linear in time
                   (clamped at the ends of the files); the wind of every
                   cell drives the source terms (seastate_implicit)
  Stationary runs with forcing are a sequence of stationary solves,
  one per A 706 interval (as SWAN's quasi-stationary mode).

  start   stand-alone run: ini, then the time loop calling step
  ini     set-up (environment, storage, initial spectrum)
  step    one time step of length A 706, so that a host model
          (SFLOW, NHFLOW) can drive the module later
--------------------------------------------------------------------*/

class seastate_f final : public seastate, public increment
{
public:
    seastate_f(lexer*, ghostcell*);
    virtual ~seastate_f();

    void start(lexer*, ghostcell*) override;

    void ini(lexer*, ghostcell*);
    void step(lexer*, ghostcell*);

    void ini_coupled(lexer*, ghostcell*);
    void step_coupled(lexer*, ghostcell*, double dt);

    fdm_seastate *e;

    const seastate_source *source() const {return psrc;}
    seastate_vtp *printer() {return pprint;}

    // surfbeat (A 770 1)
    bool surfbeat() const {return sb!=nullptr;}
    const seastate_surfbeat *generator() const {return sb;}
    seastate_roller *roller() {return proll;}
    bool longwave(lexer*, int j, double t, double &zeta, double &qx, double &qy) const;

private:
    void ini_common(lexer*, ghostcell*, bool coupled);
    void handover(lexer*, ghostcell*);
    void check_keys(lexer*, ghostcell*);
    void environment(lexer*, ghostcell*);
    void storage(lexer*, ghostcell*);
    void initial(lexer*, ghostcell*);
    void initial_parametric(lexer*, ghostcell*);
    void parametric_spectrum(lexer*, ghostcell*, const seastate_grid&, std::vector<float>&, const char*);
    void boundary(lexer*, ghostcell*);
    void boundary_spectrum(lexer*, ghostcell*, const seastate_grid&, std::vector<float>&);
    void sources(lexer*, ghostcell*);
    void kinematics(lexer*, ghostcell*, double);
    void transport(lexer*, ghostcell*);
    void parameters(lexer*, ghostcell*);

    void forcing_ini(lexer*, ghostcell*);
    void boundary_series(lexer*, ghostcell*);
    void forcing_update(lexer*, ghostcell*, double t);

    void surfbeat_input(lexer*, ghostcell*);
    void surfbeat_ini(lexer*, ghostcell*);
    void surfbeat_boundary(lexer*, double);
    void surfbeat_cap(lexer*);
    void surfbeat_step(lexer*, ghostcell*);

    void log_ini(lexer*);
    void log_step(lexer*);

    seastate_vtp *pprint;
    regression_dump *preg;
    seastate_exchange *pex;
    seastate_implicit *psolv;
    seastate_store *N0;
    seastate_source *psrc;

    std::vector<float> Nb;      // boundary spectrum
    int side[4];                // sides with the boundary spectrum: x-, x+, y-, y+
    int iter_max, iter_done;
    double dtw;                 // wave time step (A 706 stand-alone, the coupling interval in coupled runs)
    bool coupled;
    double conv;
    std::vector<double> hs_old;

    ofstream integral;
    double etot, hsmax, hsmean, nmin, cells_active;
    double starttime, endtime;

    // surfbeat
    seastate_surfbeat *sb;
    seastate_grid *gin;         // multi-frequency grid of the boundary spectrum
    seastate_roller *proll;
    std::vector<float> Nin;     // boundary spectrum on gin
    std::vector<float> Nbx;     // x- boundary rows (local j), ndir each
    std::vector<float> Nbx0;    // the rows of the previous time level
    double tnew;                // time of the new time level
    double trep;                // representative period [s]

    // external forcing (A 711 3, A 730 2)
    seastate_spc_series *bser;
    seastate_wind_series *wser;
    double tref;                            // file time [s since 1970] of model time 0
    std::vector<float> Nside[4];            // boundary spectra per side and boundary cell
    struct bweight {int a, b; float wa, wb;};
    std::vector<bweight> bw[4];
    slice4 *wU10, *wdir, *wUx, *wUy;        // wind of every cell: speed, direction [rad], components
};

#endif
