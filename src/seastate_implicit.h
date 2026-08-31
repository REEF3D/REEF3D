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

#ifndef SEASTATE_IMPLICIT_H_
#define SEASTATE_IMPLICIT_H_

#include"increment.h"
#include<vector>
#include<functional>
#include<memory>
#include<cstdint>

class lexer;
class ghostcell;
class fdm_seastate;
class seastate_store;
class seastate_exchange;
class seastate_source;
class slice;
class sliceint;

// spectra on the faces of the cells covered by a finer grid (REEFAMR, seastate_amr): face(i,j,s)
// is the spectrum that leaves covered cell (i,j) through its side s (0 x-, 1 x+, 2 y-, 3 y+),
// nullptr: the spectrum of the cell
struct seastate_faces
{
    virtual ~seastate_faces() = default;
    virtual const float *face(int i, int j, int s) const = 0;
};

using namespace std;

// neighbour of a cell: active cell, boundary with incoming spectrum, zero-gradient side, or nothing (land, open side)
struct seastate_neighbour
{
    const float *N = nullptr;     // spectrum (cell or boundary), nullptr: no inflow
    const float *cg = nullptr;    // group velocity, nullptr: no active cell
    double U = 0.0, V = 0.0;
    bool self = false;            // zero-gradient side: inflow of the cell's own spectrum
    const float *ca = nullptr;    // diffraction (A 718): Ca per frequency of the neighbour
    double tf = 1.0;              // obstacle on the face (A 722): energy transmission Kt^2
    const float *tff = nullptr;   // ... per frequency (structures, A 725), nullptr: tf
};

class seastate_grid;
class seastate_obstacle;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - implicit transport of the wave action density

  (N - N0)/dt + d(cx N)/dx + d(cy N)/dy + d(c_sig N)/dsig + d(c_theta N)/dtheta = 0

Finite volumes, first-order upwind in x, y, sigma and theta, implicit
(BDF1; stationary: 1/dt = 0). Face velocities: mean of the two cells
(geographic, sigma), analytic at the face angle (theta).

The matrix is an M-matrix: positive diagonal (1/dt + outflow
coefficients), non-positive off-diagonals (inflow coefficients), and
every column sums to 1/dt >= 0 (each outflow of a cell is the inflow of
its neighbour). Its inverse is non-negative, so N >= 0 without limiters
or clipping, for any time step.

Solver: Gauss-Seidel with four sweeps per iteration, one per directional
quadrant (theta in [q pi/2, (q+1) pi/2)), each in the downstream cell
order of its quadrant, so that geographic upwinding is exact within a
sweep. In each cell the directions of the quadrant are solved together
per frequency (tridiagonal in theta, Thomas algorithm); frequencies and
the directions of the other quadrants use the latest values. Halo
exchange after each sweep (block Gauss-Seidel across ranks).

Source terms (Phase 2, seastate_source): dN/dt = P - D N with P >= 0
added to the right-hand side and D >= 0 to the diagonal, computed per
cell from the latest spectrum in every sweep (Gauss-Seidel), so the
matrix stays an M-matrix. With the deep-water physics the change of N
per iteration is limited by the action density limiter of the source
terms (seastate_source::limit), which keeps N >= 0.

Boundaries: land (inactive cells) absorbs (outflow, no inflow). Domain
sides with side[s] = 1 (x-, x+, y-, y+) let the boundary spectrum Nb
enter; sides with side[s] = 2 are zero-gradient (the inflow spectrum
is the spectrum of the boundary cell itself, implicit as long as the
diagonal stays dominant, otherwise with the latest value); other sides
are open without incoming waves. Surfbeat (boundary_rows): the x- side
takes a spectrum per row j instead of Nb; with a boundary series (A 711 3,
boundary_sides) every side with side[s] = 1 takes a spectrum per boundary
cell. Frequency range ends: outflow only.

Surfbeat (second_order, A 775 2, cell_surfbeat): the wave groups travel
far compared with their length, so the transport is Crank-Nicolson in
time (theta = 1/2 with the old values N0 of the step) and the
geographic fluxes get a second-order correction
F2 - F1 = c/2 phi(r) (N_down - N_up) (van Leer limiter, r = (N_up -
N_upup)/(N_down - N_up)) between active cells, with the latest values
(deferred correction); the source terms stay implicit, N is clipped at
zero. The first-order path (cell) is unchanged.
Directions: periodic (full circle); with a computational sector (A 717,
Phase 10) the quadrant ranges m0..m1 hold only the active directions
(a quadrant without any is skipped), the directions outside hold zero.
Source iterations (A 738, Phase 7b; fixed in Phase 10): the action
density limiter bounds the change from the spectrum of the cell before
its first solve; a cell that does not contract keeps its first solve.

Mesh refinement (Phase 5b, seastate_amr): the same solver runs on every
grid of the hierarchy. range() limits the sweeps to the interior of a
patch (the cells around it hold the spectra of the neighbour grids and
act as inflow); skip() leaves out the cells covered by a finer grid,
whose spectra are restricted from it, and faces() gives the spectra
that leave the finer grid through the faces of the covered cells (the
mean of the two fine cells next to the face), so that the coarse grid
takes the outflow of the fine grid. Without these calls the solver
works on the whole grid as before. Composite sweep (A 797 1): visit()
calls back for each covered cell in the sweep order, solve() solves
one cell.

Phase 7 (performance):
- Cell kernel in three passes over the frequencies: A assembles the
  tridiagonal theta systems of all frequencies (coefficients per
  direction and frequency computed once per cell, the bin loop
  branch-free and vectorised), B computes the elimination factors of
  all frequencies together (the inner loop over the frequencies hides
  the latency of the divisions), C substitutes. Frequency shift only
  where c_sigma can be nonzero (depth change in time or along a current,
  current gradients); then C is Gauss-Seidel in sigma, else the
  frequencies are independent.
- Spectral sparsity (A 795 > 0, stationary): per cell, quadrant and
  frequency the range of directions holding energy; bins below A 795
  times the energy of the cell are zero. A sweep solves per frequency
  only the window of the bins that can receive energy (own range, one
  more direction with refraction, the next frequencies with frequency
  shift, the feeding frequencies of the triads, the ranges of the four
  neighbours, the ends of the quadrant next to energy in the neighbouring
  quadrants; the whole quadrant with wind input or DIA). Where the
  estimated inflow into the bin beyond an end of the window exceeds the
  threshold, the window grows and the frequency is solved again (with
  its source terms). Cells and frequencies without energy around them
  cost a few integer operations.
- Second-order upwind geographic fluxes (A 796 2, as SWAN SORDUP):
  F = c (3/2 N_up - 1/2 N_upup), implicit in the cell (the upwind cells
  are solved before it in the sweep), first order where a face has no
  two active upwind cells; not monotone, the solution is clipped at zero;
  two halo layers.

Phase 6 (seastate_obstacle): faces blocked by obstacles, structures or
reflecting coasts transmit Kt^2 of the inflow (per frequency for the
structures of A 725) and return Kr^2 of the outflow into the mirrored
direction (specular or diffuse); the spectral sparsity solves the whole
quadrant in a cell with a reflecting face, the second-order fluxes stay
first order next to a blocked face. Diffraction (A 718): face
velocities Ca c_g and the turning c_g dCa/dn.
--------------------------------------------------------------------*/

class seastate_implicit : public increment
{
public:
    seastate_implicit(lexer*, fdm_seastate*);

    void sources(seastate_source *s) {src = s;}

    // surfbeat: spectra of the x- boundary per row j (local, 0..knoy-1), nbin each, at the new
    // and the previous time level; nullptr: Nb
    void boundary_rows(const vector<float> *rows, const vector<float> *rows_old) {Nbs[0] = rows; Nbs0[0] = rows_old;}

    // boundary spectra per side and boundary cell (A 711 3): x-, x+ per local row j, y-, y+ per
    // local column i, nbin each; nullptr: Nb
    void boundary_sides(const vector<float> *const s[4]) {for(int k=0; k<4; ++k) Nbs[k] = s[k];}

    // wind field (A 730 2): U10 [m/s] and direction [rad] per cell for the source terms
    void wind_field(slice *U10, slice *dir) {wU = U10; wD = dir;}

    // vegetation field (A 756 1): stems per m^2 per cell for the source terms
    void vegetation_field(slice *nv) {vN = nv;}

    // surfbeat (A 775 2): Crank-Nicolson, second-order geographic fluxes (cell_surfbeat);
    // needs N0 (2 iterations) and 2 halo layers
    void second_order(bool s) {second = s;}

    void iterate(lexer*, ghostcell*, fdm_seastate*, seastate_exchange*, const seastate_store *N0,
                 double rdt, const vector<float> &Nb, const int side[4], bool refraction, bool fshift);

    // one sweep of the directional quadrant q (iterate: q = 0..3, halo exchange after each)
    void sweep(lexer*, fdm_seastate*, int q, const seastate_store *N0, double rdt,
               const vector<float> &Nb, const int side[4], bool refraction, bool fshift);

    // mesh refinement: sweeps over [i0,i1]x[j0,j1] only; cells with skip(i,j) == 1 not solved;
    // spectra on the faces of the skipped cells
    void range(int i0_, int i1_, int j0_, int j1_) {ri0 = i0_; ri1 = i1_; rj0 = j0_; rj1 = j1_; ranged = true;}
    void unrange() {ranged = false;}
    void skip(sliceint *s) {skp = s;}
    void faces(const seastate_faces *f) {fcs = f;}

    // Phase 7: spectral sparsity (A 795, threshold relative to the energy of the cell, 0 off; stationary)
    // and second-order geographic fluxes (A 796 2)
    void sparsity(double e_) {eps = e_;}
    void geographic_order(int o) {order2 = (o==2);}
    // FAS (A 758): residual mode (r = b - A N into R, no solve) and the coarse-grid source -tau (tau, the restricted
    // spectra Nt at the start of the coarse solve)
    void residual_out(seastate_store *R) {rout = R;}
    void fas_source(const seastate_store *tau, const seastate_store *Nt) {ftau = tau; ftil = Nt;}
    void source_iterations(int n, double tol) {nsrcit = std::max(n,1); srctol = tol;}

    // threads per rank (A 798): the cells of a quadrant sweep in wavefront order, the cells of one
    // diagonal i+j = const at the same time (level 0 without refinement, not the surfbeat model)
    void threads(int n, int nmin = 1) {nthreads = std::max(n,1); thread_min = nmin;}

    // Phase 6: diffraction (A 718): Ca and its gradient per cell and frequency; obstacles (A 722)
    void diffraction(const seastate_store *ca, const seastate_store *cax, const seastate_store *cay) {dca = ca; dcax = cax; dcay = cay;}
    void obstacles(const seastate_obstacle *o) {pob = o;}

    // mesh refinement, composite sweep: solve the one cell (i,j) for quadrant q
    void solve(lexer*, fdm_seastate*, int q, int ci, int cj, const seastate_store *N0, double rdt,
               const vector<float> &Nb, const int side[4], bool refraction, bool fshift);

    // mesh refinement, composite sweep: called for the covered cells in the sweep order (else skipped)
    void visit(std::function<void(int,int)> *f) {vis = f;}
    // ... with threads (A 798): the visit of a covered cell gets the thread index (0 .. A 798 - 1) first
    void visit_threads(std::function<void(int,int,int)> *f) {visT = f;}
    int thread_count() const {return nthreads;}
    // spectral sparsity: the active ranges of all managed cells of the grid, set before a threaded sweep
    void preset_ranges(lexer*, fdm_seastate*);

private:

    void cell(lexer*, fdm_seastate*, int q, int ci, int cj, const float *N0, double rdt,
              const vector<float> &Nb, const int side[4], bool refraction, bool fshift);

    void sweep_threads(lexer*, fdm_seastate*, int q, const seastate_store *N0, double rdt,
                       const vector<float> &Nb, const int side[4], bool refraction, bool fshift);
    int nthreads = 1;
    int thread_min = 1;                 // threads only on grids of at least this many cells
    const seastate_store *dca = nullptr, *dcax = nullptr, *dcay = nullptr;
    const seastate_obstacle *pob = nullptr;
    vector<double> cdp, cdm;            // diffraction: c_theta at the faces m+1/2, m-1/2 per frequency and direction

    void cell_surfbeat(lexer*, fdm_seastate*, int q, const seastate_store *N0, double rdt,
                       const vector<float> &Nb, const int side[4], bool refraction);

    int nsig, ndir;
    seastate_source *src;
    const vector<float> *Nbs[4], *Nbs0[4];
    slice *wU, *wD;
    slice *vN = nullptr;
    bool second;
    bool ranged = false;
    int ri0 = 0, ri1 = -1, rj0 = 0, rj1 = -1;
    sliceint *skp = nullptr;
    const seastate_faces *fcs = nullptr;
    vector<double> P, D;                // source terms of the cell (src != nullptr)
    vector<double> Ar;                  // sig/sinh(2kd) of the cell per frequency
    vector<double> qc, qs, cur;         // cos, sin and current shear term per direction of the quadrant
    vector<double> rth;                 // 1/dth per direction of the quadrant
    vector<double> tAp, tCp, tAm, tCm;  // c_theta = A tA + tC at the faces m+1/2 (p) and m-1/2 (m)
    vector<double> csg;                 // c_sigma per frequency and direction of the quadrant (nsig x ndir)
    vector<float> zero;                 // zero spectrum (no inflow)
    vector<double> zerod;               // zero source terms
    vector<double> tdi, tla, tup, trh, tsi, trd, tcp, tdp, tso;   // theta systems of all frequencies [l*ndir+n]
    vector<int> tok;                    // per frequency: 0 after a vanishing pivot
    void store(float *N, int l, int ma, int nq, double dmax, bool clip);

    // spectral sparsity
    double eps = 0.0;
    bool order2 = false;
    seastate_store *rout = nullptr;
    const seastate_store *ftau = nullptr, *ftil = nullptr;
    int nsrcit = 1;                     // source iterations per cell and sweep (A 738)
    double srctol = 0.0;                // ... until the relative energy change of the cell is below this
    vector<int> wlo, whi;               // window of directions solved per frequency (quadrant index)
    vector<double> thr;                 // threshold per frequency (action density)
    // the active ranges are shared with the copies of the solver that work the threads (A 798)
    struct shared_ranges
    {
    vector<uint16_t> rg, rv, rb;
    int rni = 0, rnj = 0, rimin = 0, rjmin = 0;
    };
    std::shared_ptr<shared_ranges> sr = std::make_shared<shared_ranges>();
    vector<uint16_t> &rg = sr->rg;      // active ranges per cell, quadrant, frequency
    vector<uint16_t> &rv = sr->rv;      // per cell: 1 ranges set, 2 triads active, 4<<q quadrant q holds energy,
                                        // 64<<2q / 128<<2q energy in the first / last direction of quadrant q
    void summary(int ci, int cj, int q, const uint16_t *r);
    vector<double> rsum;                // row sums of the active ranges
    vector<float> Nsi0, Nsi1;           // source iterations: the spectrum of the cell before and after the first solve
    const float *Nlim = nullptr;        // ... the spectrum the action density limiter refers to (before the first solve)
    vector<uint16_t> &rb = sr->rb;      // per cell and quadrant: band of frequencies with energy (lo<<8 | hi)
    int lmin = 0, lmax = -1;            // band of the windows of the current cell
    int &rni = sr->rni, &rnj = sr->rnj, &rimin = sr->rimin, &rjmin = sr->rjmin;
    std::function<void(int,int)> *vis = nullptr;
    std::function<void(int,int,int)> *visT = nullptr;
    double energy(const seastate_grid&, const float *N) const;
    bool managed(lexer*, fdm_seastate*, int ci, int cj) const;
    uint16_t *ranges(lexer*, fdm_seastate*, int ci, int cj);
    void windows(lexer*, fdm_seastate*, int q, int ci, int cj, const float *N, const seastate_neighbour &W,
                 const seastate_neighbour &E, const seastate_neighbour &S, const seastate_neighbour &Nn, bool rf, bool fs,
                 bool all=false);
    void keep(lexer*, fdm_seastate*, int q, int ci, int cj);
    const float *usable(lexer*, fdm_seastate*, int ci, int cj) const;
    vector<int> m0, m1;                 // first and last direction of each quadrant
    vector<double> la, di, up, rhs, sol, cp, dp;   // tridiagonal system of one frequency
};

#endif
