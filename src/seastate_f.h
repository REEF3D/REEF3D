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
#include<string>
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
class seastate_amr;
class seastate_fas2;
class seastate_bathy;
class slice4;
class seastate_obstacle;
class regression_dump;
class reefrun_ckpt;
class reefckpt;

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

Phase 5b: static mesh refinement on REEFAMR (seastate_amr, G 1, stand-
alone runs): refined patches from boxes (G 10, G 11) and criteria
(A 791 depth, A 792 coastline, A 793 bathymetry gradient), composite
solve over all levels per directional quadrant; bathymetry raster
seastate-bathy.dat for all grids (A 790 1, seastate_bathy). Also the
SWAN maximum energy with Battjes-Janssen breaking (A 737) and the
triad parameters (A 736).

Phase 7: performance (seastate_implicit, seastate_amr): restructured
cell kernel, spectral sparsity (A 795, stationary), second-order upwind
geographic fluxes (A 796 2), composite sweep across the refinement
levels (A 797 1); convergence history of stationary runs in
REEF3D_SEASTATE_Log/REEF3D_SEASTATE_convergence.dat.

Phase 18: checkpoints and the rebalance of stand-alone runs (seastate_f_ckpt.cpp,
reefrun_ckpt.h): a checkpoint every P 48 steps (stationary: every P 48
iterations of a step), the restart from the latest on any decomposition
(I 45 1): the spectra of the cells of all grids by their global index, the
active ranges of the spectral sparsity (A 795) and the last change of Hs of
the convergence test (A 799 1); a stationary step continues at the
iteration of the checkpoint. The measured work of every cell (the time of
its solves) gives the load of the ranks and the load map
REEF3D_SEASTATE_Log/REEF3D_SEASTATE_loadmap.dat (G 45) for the weighted
partition of DIVEMesh (M 20 5, M 35 1); with G 46 the run stops to
rebalance when the load max/mean of the ranks exceeds G 46 (checked every
G 47 steps, stationary every G 47 iterations) and continues on the
weighted partition (tools/reef3d_rebalance.sh, DIVEMesh M 50 1).

Phase 6: phase-decoupled diffraction (A 718, A 719,
seastate_f_diffraction.cpp), line obstacles (A 722) and, Phase 6b,
structures (A 725: d'Angremond, per frequency, porous; seastate_structure),
diffuse reflection (A 726), reflecting coasts (A 723, A 724;
seastate_obstacle), vegetation (A 754 - A 756, seastate_source), all of
them also on the patches of the mesh refinement.

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

    // Phase 18: the checkpoints (P 48, I 45) and the rebalance (G 46) of a stand-alone run (seastate_f_ckpt.cpp)
    reefrun_ckpt& ckpt(lexer*, ghostcell*);
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
    void auto_sector(lexer*, ghostcell*, const seastate_grid&);
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
    bool kba = false;                   // Phase 14: the sweeps pipelined across the ranks (A 781 1)
    double conv;
    std::vector<double> hs_old;
    std::ofstream convlog;              // stationary convergence history (rank 0)

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

    // mesh refinement (G 1), bathymetry raster (A 790 1)
    seastate_amr *pamr;
    seastate_bathy *bathy;

    // Phase 12: FAS on a uniform grid (A 758 without G 1): a V-cycle over coarse grids of 2 x 2 cells
    seastate_fas2 *pfas = nullptr;
    void transport_amr(lexer*, ghostcell*);

    // Phase 6: diffraction (A 718), obstacles (A 722)
    void diffraction(lexer*, ghostcell*);
    void obstacles(lexer*, ghostcell*);
    seastate_store *dca = nullptr, *dcax = nullptr, *dcay = nullptr;   // Ca and its gradient per frequency
    slice4 *dfS = nullptr, *dfT = nullptr, *dfK = nullptr, *dfC = nullptr;
    int dfsmooth = -1;                      // smoothing steps of the diffraction parameter (A 719)
    seastate_obstacle *pobs = nullptr;

    // Phase 6b: diffraction with mesh refinement (smoothed energy of level 0 for the patches, the mean
    // wavelength and the smallest level-0 cell), vegetation field (A 756 1)
    slice4 *dfE = nullptr;
    double dfL = 0.0, dfdmin = 0.0;
    slice4 *vegN = nullptr;
    seastate_bathy *vegr = nullptr;
    void vegetation_ini(lexer*, ghostcell*);

    // Phase 18: checkpoints and the rebalance (seastate_f_ckpt.cpp)
    reefrun_ckpt *rck = nullptr;
    bool restart = false;                   // I 45 1: from a checkpoint
    int ss_resume = 0;                      // stationary: the iterations done in the step of the checkpoint (continued)
    std::vector<double> dprev;              // stationary: the last change of Hs of every leaf cell (A 799 1)
    std::string ckpt_unsupported(lexer*) const;
    void ckpt_write(lexer*, reefckpt&);
    bool ckpt_read(lexer*, ghostcell*, reefckpt&, std::string&);
    void ckpt_iteration(lexer*, int k);     // stationary, after iteration k of the step that did not stop: P 48 and G 46
    void ckpt_step(lexer*, ghostcell*);     // after the outputs of a step: P 48 and G 46 (nonstationary), G 45
    // the measured work of the cells [s] (G 45, G 46): level 0 (the patches: seastate_amr), per level-0 column of the rank
    // the work of the last closed window and of the run, and the iterations in them
    slice4 *wk0 = nullptr;
    std::vector<double> wwin, wrun;
    std::vector<double> lmap0;              // rank 0: the last load map written, all level-0 columns (I*GNY+J)
    // whether a rebalance pays: the set-up time of the run, the work of the heaviest rank and of the mean rank in the
    // last window [s], the bytes of the state in a checkpoint, the iteration of the step (stationary)
    double t_ini = 0.0, busy_max = 0.0, busy_mean = 0.0, state_bytes = 0.0;
    int k_now = 0;

    // Phase 19: the time of the iterations (transport) and their number, the report of the pipelined sweeps
    double t_trans = 0.0;
    long n_trans = 0;
    void pipeline_report(lexer*, ghostcell*);
    std::vector<double> hconv, hperc;       // stationary: the largest change of Hs and the percentage within A 708 per iteration
    double nwin = 0.0, nrun = 0.0;
    int wit = 0;                            // iterations of the open window
    void work_columns(lexer*, std::vector<double>&, bool reset);
    double work_imbalance(lexer*, ghostcell*);
    void loadmap_write(lexer*, const std::vector<double>&, double iterations, const char *what);
};

#endif
