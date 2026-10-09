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


#ifndef SEASTATE_AMR_H_
#define SEASTATE_AMR_H_

#include"reefamr.h"
#include"seastate_implicit.h"
#include<fstream>
#include<string>
#include<unordered_map>
#include<vector>

class lexer;
class ghostcell;
class fdm_seastate;
class seastate_store;
class seastate_source;
class seastate_exchange;
class seastate_bathy;
class seastate_wind_series;
class slice4;
class sliceint4;
class seastate_obstacle;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - static mesh refinement on REEFAMR (Phase 5b, G 1)

Level 0 is the SEASTATE grid of seastate_f. Every refined patch is a
small domain of its own (reefamr: lexer, refinement ratio 2) with its
own field data (fdm_seastate), block-sparse spectra and implicit
solver, so the kernels of seastate_implicit run unchanged on every
grid.

Hierarchy: built once at the start (static), G 1 levels, tiles G 4,
buffer G 3, refinement boxes G 10, no-refinement boxes G 11, and the
criteria evaluated on each level for the next finer one:
  A 791  water depth below this value [m]
  A 792  within this many cells of a dry or land cell (coastline)
  A 793  relative depth change |d_nb - d|/max(d, d_nb) to a neighbour
         cell above this value (bathymetry gradient)
  A 762  water width below this many cells of the level (Phase 9): the
         shortest wet run through the level-0 cell along x, y and the
         diagonals, runs open to the domain edge do not count
No refinement within A 794 level-0 cells of the sides with boundary
spectra or zero gradient (A 712 1, 2): the forcing acts on level 0.

Bed of the patches: from the bathymetry raster (A 790 1,
seastate_bathy: mean of the raster nodes inside the cell, also on
level 0), else interpolated from level 0 (limited linear, as the
other REEFAMR modules). Active cells, k, cg and the gradients follow
from the depth of each grid (seastate_kinematics).

Composite solve, per directional quadrant q of every iteration
(block Gauss-Seidel over the levels, a V-cycle of sweeps):
  1  level 0: sweep q (cells covered by level 1 skipped)
  2  level l = 1..G 1: the cells next to the patch interiors take the
     spectra of the same level or of the next coarser level (constant
     prolongation: positive and conservative), then sweep q over the
     patch interiors (cells covered by level l+1 skipped); the patches
     of a level are swept in the downstream order of the quadrant and
     each takes the latest spectra of its neighbours just before its
     sweep (Gauss-Seidel across the patches of a level)
  3  restriction, finest first: a covered cell takes the mean of its
     four children; the faces of covered cells next to cells solved on
     the coarse grid take the mean of the two fine cells next to the
     face (the outflow of the fine grid, seastate_faces)
  4  the coarser levels again, G 1 - 1 down to 0, so that the cells
     downstream of a patch take its outflow within the same sweep
     (without this every pass through a patch costs one iteration);
     level 0 only in the window downstream of the level-1 patches of
     the rank; halo exchange of level 0
Composite sweep (A 797 1, Phase 7) instead of steps 1-4: one sweep of
level 0 per quadrant in which every covered cell, where the sweep
reaches it, sweeps its four children in the order of the quadrant
(recursively for finer levels), each fine cell after the ring cells
next to it took the latest spectra of their sources (copy or constant
prolongation), then takes the mean of the children and the spectra on
its faces. Every cell is solved after its upwind neighbours on all
levels, so one sweep carries the waves through the whole hierarchy and
the refined grids need the iterations of a uniform grid. Ring cells
held by other ranks: the exchange at the start of each quadrant sweep
(lagged, as the level-0 halo). With parents on other ranks (G 40) the
level-by-level sweeps are used.
Nonstationary runs keep the spectra of the start of the step on every
grid (N0); stationary runs iterate until the largest relative change
of Hs on all grids is below A 708.

Output: REEF3D_SEASTATE_AMR/ (VTR per patch with Hs, Tm01, Tp, dir,
spread, depth, level, at every VTP print of level 0; a .vtm per print
with all patches; the patch log). The handover points (A 760) and the
integral log use the finest grid.

FAS coarse-grid correction (A 758 n m, Phase 9; stationary, G 1 1):
every n iterations level 0 is the coarse grid of the level-1 patches.
The residuals r = b - A N of the patch cells (residual mode of the
solver) are restricted (mean of the children, as the spectra). Level 0
on its own (covered cells solved as coarse cells, no fine faces) gets
the source -tau: tau = r_c(R N) - R r on the covered cells, and on the
others the difference of their residual on their own and in the
composite problem (they differ next to the patches), so that the
converged composite solution is a fixed point of the coarse problem.
It is iterated m times (the loss part of -tau proportional to N), and
the children take the coarse change multiplicatively (N *= N_c/R N,
within 1/4 to 4; additive where R N vanishes), then the restriction
again.

Not with the coupling (A 750), the surfbeat model (A 770), regridding
or patches on several ranks (G 40): patches are cut at the rank boxes.
--------------------------------------------------------------------*/

struct seastate_amr_patch : public reefamr_patch
{
    fdm_seastate *e = nullptr;
    seastate_implicit *solv = nullptr;
    seastate_store *N0 = nullptr;           // spectra at the start of the step (nonstationary)
    seastate_store *res = nullptr;          // FAS (A 758): residual of the cell equations
    sliceint4 *cov = nullptr;               // 1: covered by a patch of the next finer level
    slice4 *wU = nullptr, *wD = nullptr;    // wind field (A 730 2)
    vector<int> ring;                       // fill entries next to the interior (inflow of the sweeps), local sources
    vector<int> ringr;                      // the same, held by other ranks

    // Phase 6b: obstacles and coasts, diffraction (Ca, its gradient; smoothed energy, work slices), vegetation
    seastate_obstacle *pobs = nullptr;
    seastate_store *dca = nullptr, *dcax = nullptr, *dcay = nullptr;
    slice4 *dS = nullptr, *dT = nullptr, *dK = nullptr, *dC = nullptr, *dE = nullptr;
    slice4 *vN = nullptr;
};

// spectra on the faces of the covered cells of one grid
struct seastate_amr_faces : public seastate_faces
{
    unordered_map<long,vector<float>> F;    // key(i,j,s) of the covered cell

    static long key(int i, int j, int s) {return ((long(i)+64L)*16777216L + (long(j)+64L))*4L + long(s);}
    const float *face(int i, int j, int s) const override
    {
        auto it = F.find(key(i,j,s));
        return it==F.end() ? nullptr : it->second.data();
    }
};

// a face spectrum and the two fine cells next to the face
struct seastate_amr_flink
{
    int g;          // coarse grid id
    long key;       // face key in the faces of g
    int id;         // fine patch
    int i0,j0,i1,j1;
};

class seastate_amr : public reefamr
{
public:
    seastate_amr(lexer*, ghostcell*);
    virtual ~seastate_amr();

    // level-0 objects of seastate_f
    struct level0
    {
        fdm_seastate *e;
        seastate_implicit *solv;
        seastate_exchange *pex;
        seastate_source *src;
        seastate_store *N0;
        const seastate_bathy *bathy;
        seastate_wind_series *wser;
        double tref;
        seastate_obstacle *pobs = nullptr;      // obstacles and coasts of level 0 (A 722 - A 726)
        const seastate_bathy *veg = nullptr;    // vegetation raster (A 756 1)
    };

    void ini(lexer*, ghostcell*, const level0&);
    bool active() const {return maxlev>0 && patches_total>0;}
    int level0_covered(int i, int j) const;

    // nonstationary step start: N0 of the patches
    void step_begin();

    // one composite iteration (four quadrant sweeps on all grids)
    void iterate(lexer*, ghostcell*, const seastate_store *N0, double rdt, const vector<float> &Nb,
                 const int side[4], bool refraction, bool fshift);

    // FAS coarse-grid correction (A 758, stationary, one refinement level): level 0 is the coarse grid of the
    // level-1 patches
    void fas(lexer*, ghostcell*, double rdt, const vector<float> &Nb, const int side[4], bool refraction, bool fshift, int ncoarse);
    double fas_last() const {return fas_change;}

    // integrated parameters of the patches; Hs of the patch interiors (leaf cells) for the
    // stationary convergence (appended to hs, the cell area relative to a level-0 cell to w),
    // the largest Hs and the smallest N of the leaf cells
    void parameters(vector<double> *hs, double &hmax, double &vmin, vector<double> *w=nullptr);

    // wind field on the patches (A 730 2), time t of the model
    void wind(double t);

    // Phase 6b: obstacles and coasts of the patches from their present spectra (before every iteration)
    void obstacles();

    // Phase 6b: diffraction parameter of the patches, level by level (A 718 mode 1: frequencies l0..l1
    // together, 2: one frequency), after level 0: the energy of every patch is smoothed with the values of
    // the next coarser grid in the ring around the interior (the smoothed energy and Ca, bilinear), the
    // number of steps 0.4 (L/dx)^2 of the patch (A 719 n: n 4^level), at most 400
    void diffraction(lexer*, int mode, int l0, int l1, double smax, double L, double dmin0, slice4 &E0,
                     slice4 &T0);

    // finest grid on this rank that holds the point: lexer, field data, cell
    bool locate(double x, double y, lexer *&q, fdm_seastate *&ee, int &ci, int &cj);

    void print(lexer*, ghostcell*);
    void convergence(lexer*, int iteration, double change, double percent);
    void log(lexer*, ghostcell*, int iterations);

protected:
    reefamr_patch* patch_new() override;
    void patch_objects(reefamr_patch*, ghostcell*) override;
    void patch_delete(reefamr_patch*) override;
    void tag(int, vector<unsigned char>&) override;
    void regrid_prepare(ghostcell*) override {}
    void regrid_static(ghostcell*) override;
    void regrid_state(ghostcell*, vector<reefamr_patch*>&) override;
    void regrid_finish(ghostcell*, int) override;

private:
    static seastate_amr_patch* SP(reefamr_patch *c) {return static_cast<seastate_amr_patch*>(c);}
    seastate_amr_patch* SP(int n) {return static_cast<seastate_amr_patch*>(P[n]);}

    // grid id -1: level 0
    fdm_seastate* gfd(int g);
    sliceint4* gcov(int g);
    seastate_amr_faces* gfaces(int g) {return &faces[g+1];}

    void environment(seastate_amr_patch&);
    bool coarse_value(int l, int I, int J, int what, slice4 &E0, slice4 &T0, double &v);
    double bed_at(int l, int I, int J);
    void covered();
    void restrict_all();
    void restrict_level(int l);
    void fill_level(int l);
    void fill_local(seastate_amr_patch*);
    vector<vector<int>> order[4];           // [q][l]: patches of level l in the downstream order of quadrant q
    void sweep_level(int l, int q, double rdt, bool refraction, bool fshift);

    // composite sweep (A 797 1): the cells of a patch are solved where the sweep of the coarser grid
    // reaches their parent cell, so that one sweep carries the waves through all levels
    static long long ckey(int g, int i, int j) {return ((long long)(g+1)<<42) | ((long long)(i+1048576)<<21) | (long long)(j+1048576);}
    unordered_map<long long,pair<int,int>> childof; // ckey(g,I,J) of a covered cell -> patch id, block
    vector<unordered_map<long long,int>> ghostof;   // [id]: ckey(-1,i,j) of a ring cell -> fill entry (this rank)
    unordered_map<long long,vector<pair<int,int>>> flinkof;   // ckey(g,I,J) -> faces of the cell: level, index in flinks[l]
    void composite_maps();
    void composite_sweep(lexer*, ghostcell*, int q, const seastate_store *N0, double rdt, const vector<float> &Nb,
                         const int side[4], bool refraction, bool fshift);
    void descend(int g, int I, int J, int q, double rdt, bool refraction, bool fshift);
    void ghost_fill(seastate_amr_patch*, int id, int fi, int fj);
    bool remote_blocks = false;             // a parent cell on another rank (restriction through the block plans)

    lexer *p0;
    level0 L0;
    int nbin, nsig;
    int fills, sweeps;
    vector<double> bed0;                    // level-0 bed of the whole domain (global index I*GNY+J)
    vector<double> width0;                  // water width of the level-0 cells [m] (A 762)
    seastate_store *fasT = nullptr, *fasN = nullptr, *fasR = nullptr, *fasC = nullptr;   // FAS: tau, restricted spectra, coarse and composite residual (level 0)
    double fas_change = 0.0;                // FAS: largest relative change of the coarse correction (log)
    void water_width(const vector<unsigned char> &wet0, double dx, double dy);
    unordered_map<unsigned long long,double> bedmemo;
    sliceint4 *cov0 = nullptr;
    vector<seastate_amr_faces> faces;       // [g+1]
    vector<vector<seastate_amr_flink>> flinks;  // [l]: the face spectra filled from level l
    vector<float> zero;
    double tm[4];
    ofstream logout, convout;
    int printcount_amr = 0;
};

#endif
