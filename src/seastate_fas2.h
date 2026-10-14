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

#ifndef SEASTATE_FAS2_H_
#define SEASTATE_FAS2_H_

#include"increment.h"
#include<vector>

class lexer;
class ghostcell;
class fdm_seastate;
class seastate_store;
class seastate_implicit;
class seastate_exchange;
class seastate_source;
class slice;
class slice4;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - FAS on a uniform grid (A 758 without mesh
refinement, Phase 12): a V-cycle over coarse grids of 2 x 2 cells

The slow mode of the stationary wind sea (the downshift of the peak
along the whole fetch, 29, 31, 46) is removed by a coarse grid that
sees the whole fetch (FAS, 45, 48). Without mesh refinement the
coarse grid is built here: every 2 x 2 cells of a rank form one coarse
cell (a last row or column of single cells where the rank has an odd
number of cells), so the coarse grid of a rank covers exactly its
cells and the decomposition needs nothing else.

Coarse cell: wet where a child is wet; depth, U and V the mean of the
wet children; k, cg and the gradients from its own depth
(seastate_kinematics); boundary spectra per side cell (A 711 3) the
mean of the two fine side cells; wind field (A 730 2) the first wet
child. The same implicit solver (seastate_implicit) and source terms
work the coarse grid, first-order geographic fluxes, no spectral
sparsity; halo exchange of its own (one layer).

The coarse grid is coarsened again in the same way (a V-cycle over up to
6 coarse grids, while every rank holds more than one cell): the slow
mode is smooth along the fetch, and the coarsest grids iterate it at a
small fraction of the cost of a fine iteration.

Every A 758 n iterations (A 758 m coarse iterations), as the FAS of the
mesh refinement (seastate_amr::fas), from the fine grid down:
  1  r = b - A N of the fine cells (residual mode of the solver), summed
     into the coarse cells directly (no fine store, residual_add)
  2  restriction: N_t and R r, the means over the wet children
  3  tau = r_c(N_t) - R r, the coarse residual at N_t added in place
Memory: three spectra per coarse cell (N, N_t, tau; nonstationary also
the restricted N0), i.e. about the spectra of the fine grid in all coarse
grids together.
  4  m iterations of the coarse problem with the source -tau (its loss
     part proportional to N), so that the fine solution restricted is
     a fixed point of it; with a further coarse grid m/2 iterations, the
     correction from that grid (the same steps one level down), m/2
     iterations; the coarsest grid 4 m iterations
  5  the children take the coarse change multiplicatively (N *= N_c/N_t
     within 1/4 to 4; additive where N_t vanishes); no correction under
     a coarse cell that covers a dry cell of the original grid (coasts:
     a coarse cell across a narrow spit of land carries energy across
     it, Phase 10 for the refinement)
Nonstationary: the time term of the step (N0 restricted), corrections
within the step but not after its last iteration.

Not with diffraction (A 718), obstacles and coasts (A 722, A 723),
vegetation (A 754), the surfbeat model (A 770) or a coupled host.
--------------------------------------------------------------------*/

class seastate_fas2 : public increment
{
public:
    // fine grid objects of seastate_f; Nside: boundary spectra per side cell (A 711 3) or nullptr; wU, wD: wind
    // field (A 730 2) or nullptr
    // depth: the number of this coarse grid (1 the first); further coarse grids are built while every rank holds
    // more than one cell and depth < maxdepth
    seastate_fas2(lexer*, ghostcell*, fdm_seastate*, seastate_implicit*, seastate_exchange*, seastate_source*,
                  const vector<float> *const Nside[4], slice *wU, slice *wD, int depth = 1, int maxdepth = 6,
                  const vector<char> *fine_coast = nullptr);
    virtual ~seastate_fas2();

    // one correction; N0 the spectra of the start of the step (nonstationary) or nullptr
    void correct(lexer*, ghostcell*, const seastate_store *N0, double rdt, const vector<float> &Nb, const int side[4],
                 bool refraction, bool fshift, int ncoarse);

    double last() const {return change;}
    long coarse_cells() const {return ncoarse_wet;}
    int levels() const {return child!=nullptr ? 1+child->levels() : 1;}
    long cells_of_level(int k) const {return k<=1 ? ncoarse_wet : (child!=nullptr ? child->cells_of_level(k-1) : 0);}
    long coarse_bytes() const;                  // all coarse grids of this rank

private:
    void build_lexer(lexer*, ghostcell*);
    void build_fields(lexer*, const vector<char> *fine_coast);
    void forcing();
    int child0(int I, int n, int nc) const;     // first fine index under coarse index I (n fine, nc coarse cells)
    int children(int I, int n, int nc, int c[2]) const;

    lexer *p0, *pc;
    fdm_seastate *ef, *ec;
    seastate_implicit *fsolv, *csolv;
    seastate_exchange *fex, *cex;
    seastate_source *src;
    seastate_store *Nt = nullptr, *Tc = nullptr, *N0c = nullptr;     // restricted N, tau, restricted N0 (nonstationary)
    const vector<float> *Nsf[4];
    vector<float> Nsc[4];
    slice *wUf, *wDf;
    slice4 *wUc = nullptr, *wDc = nullptr;
    vector<char> coast;                         // per coarse cell: a dry cell of the original grid under it
    int nbin, kx, ky, ms;
    int ci0, cj0, cgx, cgy;                     // global index of the coarse cell (0,0) and the coarse domain (boundary tests)
    long ncoarse_wet = 0;
    seastate_fas2 *child = nullptr;             // the next coarser grid (FAS V-cycle)
    int depth;
    double change = 0.0;
    int gcount[4] = {0,0,0,0};
};

#endif
