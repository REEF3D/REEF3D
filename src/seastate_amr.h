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
#include<functional>
#include<memory>
#include<mutex>
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
Faces across patches and ranks (Phase 16): a cell next to a covered cell
takes the face spectrum of the covered cell (the mean of the two fine
cells next to the face) also where the covered cell lies in another
patch of the level or on another rank: on level 0 the covered halo
cells are marked and their faces sent after every halo exchange; on the
patches the ring cells whose source cell is covered take its face, from
the fine cells on this rank or from the rank of the source patch with
the fills. The discrete problem no longer depends on where a level is
cut into patches, so a converged run on several ranks gives the result
of one rank.
Pipeline (Phase 17, A 781 1 with A 797 1): the composite sweep over the
level-0 tiles of all ranks as one dependency graph (seastate_implicit,
Phase 14). A boundary tile carries the level-0 halo cells, the faces of
its covered boundary cells and the ring cells (with their faces) of the
patches of the rank downstream whose sources lie in it; every cell then
sees its neighbours on all levels as on one rank, and the result is that
of one rank, bitwise. With the second-order fluxes (A 796 2) the patches
keep a second ring layer, used where it holds cells of the same level
not covered by a finer grid (as inside one patch). The level-by-level
sweeps (A 797 0) are not pipelined (their result depends on the order of
the patches of a level, which the rank boxes cut).
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
Phase 10: all levels (G 1 >= 1). The residuals of the uncovered cells
of every level, the covered cells taking the mean of the residuals of
their children from the finest level up; the coarse change of a level-0
cell goes into all cells under it on every level. No correction under a
coarse cell with a dry child (its restricted spectrum holds the dry
child as zero). Nonstationary runs: the residuals with the time term of
the step (N0 of every grid), a correction after every n iterations of
the step but not after the last one (a fine iteration must follow).

Threads (A 798, Phase 10): level by level, the level-0 sweeps (covered
cells skipped) and the sweeps of the patches of at least 64 cells run
over tiles of cells (seastate_implicit, Phase 13); Phase 12: on a level with at least
2 A 798 patches the patches themselves are shared out to the threads, one
thread per patch, in an order that keeps every patch to the values of the
serial order (a patch waits for the patches before it that it reads or
that read it through their rings), where this dependency graph leaves
enough parallelism (the cells of the level over the cells of its longest
chain at least 0.75 A 798; otherwise the patches are threaded inside, one
after the other, as before); in the composite sweep the level-0
tiles visit the covered cells with the index of its thread, which
solves the patch cells under them with its own copies of the patch
solvers and of the source terms (made once). The ring cells of the
patches are written under a lock. The result is that of one thread.

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
    vector<int> ring2, ring2r;              // Phase 17, second-order fluxes (A 796 2): the second layer, local / other ranks
    vector<int> ringsel;                    // ringr and ring2r: the entries the fills of the level store

    // Phase 6b: obstacles and coasts, diffraction (Ca, its gradient; smoothed energy, work slices), vegetation
    seastate_obstacle *pobs = nullptr;
    seastate_store *dca = nullptr, *dcax = nullptr, *dcay = nullptr;
    slice4 *dS = nullptr, *dT = nullptr, *dK = nullptr, *dC = nullptr, *dE = nullptr;
    slice4 *vN = nullptr;

    // Phase 18: the measured work of the cells [s] (G 45, G 46), nullptr: not measured
    slice4 *wk = nullptr;
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

class seastate_amr : public reefamr, public seastate_kba_extra
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

    // FAS coarse-grid correction (A 758): level 0 is the coarse grid of the patches of all levels; N0 the spectra
    // of the start of the step (nonstationary) or nullptr
    void fas(lexer*, ghostcell*, const seastate_store *N0, double rdt, const vector<float> &Nb, const int side[4], bool refraction, bool fshift, int ncoarse);
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

    // Phase 17: the composite sweep (A 797 1) pipelined across the ranks (A 781 1): the boundary tiles of level 0 also
    // carry the faces of the covered boundary cells and the ring cells (and their faces) of the patches of the rank
    // downstream whose sources lie in the tile. Collective; after ini, with the level-0 solver already pipelined
    void pipeline_setup(lexer*, ghostcell*);
    bool pipelined() const {return kba;}
    size_t kba_size(int q, int s, int idx, int tb, bool send) override;
    void kba_pack(int q, int s, int idx, int tb, float *buf) override;
    void kba_unpack(int q, int s, int idx, int tb, const float *buf) override;

    // finest grid on this rank that holds the point: lexer, field data, cell
    bool locate(double x, double y, lexer *&q, fdm_seastate *&ee, int &ci, int &cj);

    void print(lexer*, ghostcell*);
    void convergence(lexer*, int iteration, double change, double percent);
    void log(lexer*, ghostcell*, int iterations);

    // Phase 18: checkpoints (P 48, I 45; seastate_f_ckpt.cpp) and the measured work of the cells (G 45, G 46).
    // A patch: its level and box (global indices of the level), the offsets of the box in its arrays (local = global - I0 + ox),
    // its interior nx x ny, its lexer, field data, solver and covered marks
    struct gview {int lev, I0, I1, J0, J1, nx, ny, ox, oy; lexer *pp; fdm_seastate *e; seastate_implicit *solv; sliceint4 *cov;};
    int npatch() const {return (int)P.size();}
    gview patch_view(int id);
    // after a restart has set the spectra of the interiors of all grids: the restriction (the faces of the covered cells),
    // then again (the spectra of the checkpoint once more: its covered cells as they were, bitwise), the level-0 halo with
    // its faces and the fills of all levels (collective)
    void restored(lexer*, ghostcell*, const std::function<void()> &again);
    // the work fields of the patch solvers (before the first sweep, so that the copies for the threads take them)
    void work_fields();
    // the work of the patch cells added into the level-0 columns of this rank (col: i*knoy+j of level 0); reset: zero
    void work_columns(lexer*, vector<double> &col, bool reset);
    int print_count() const {return printcount_amr;}
    void print_count(int n) {printcount_amr = n;}
    void flush_logs() {logout.flush(); convout.flush();}

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
    void fill_local(seastate_amr_patch*, int id);
    vector<vector<int>> order[4];           // [q][l]: patches of level l in the downstream order of quadrant q
    void sweep_level(int l, int q, double rdt, bool refraction, bool fshift);

    // Phase 12: threads over the patches of a level (level by level, A 798, at least 2 A 798 patches on the level):
    // patch k of order[q][l] waits for the patches before it that it reads or that read it through their rings, so
    // that every patch sees the values of the serial order (bitwise the result of one thread); each patch is swept
    // by one thread. Only where the graph allows it: the cells of the level over the cells of its longest chain of
    // dependent patches at least 0.75 times the threads; otherwise the patches are threaded inside, one after the other
    vector<vector<int>> dag_npred[4];                   // [q][l][k]: number of predecessors of position k
    vector<double> dag_gain[4];                         // [q][l]: cells of the level / cells of the longest chain
    vector<vector<vector<int>>> dag_succ[4];            // [q][l][k]: successors of position k
    void dag_build();
    void dag_report();                                  // the choice per level in the log (rank 0)
    bool sweep_level_threads(int l, int q, double rdt, bool refraction, bool fshift);
    void thread_copies(int nt);

    // composite sweep (A 797 1): the cells of a patch are solved where the sweep of the coarser grid
    // reaches their parent cell, so that one sweep carries the waves through all levels
    static long long ckey(int g, int i, int j) {return ((long long)(g+1)<<42) | ((long long)(i+1048576)<<21) | (long long)(j+1048576);}
    unordered_map<long long,pair<int,int>> childof; // ckey(g,I,J) of a covered cell -> patch id, block
    vector<unordered_map<long long,int>> ghostof;   // [id]: ckey(-1,i,j) of a ring cell -> fill entry (this rank)
    unordered_map<long long,vector<pair<int,int>>> flinkof;   // ckey(g,I,J) -> faces of the cell: level, index in flinks[l]
    void composite_maps();
    void composite_sweep(lexer*, ghostcell*, int q, const seastate_store *N0, double rdt, const vector<float> &Nb,
                         const int side[4], bool refraction, bool fshift);
    void descend(int g, int I, int J, int q, double rdt, bool refraction, bool fshift, int t = 0);
    void ghost_fill(seastate_amr_patch*, int id, int fi, int fj);
    // threads in the composite sweep (A 798): per thread (1 .. A 798 - 1) a copy of every patch solver and of the
    // source terms; thread 0 uses the patch solvers themselves
    vector<vector<std::unique_ptr<seastate_implicit>>> tsolv;
    vector<std::unique_ptr<seastate_source>> tsrc;
    std::mutex ring_mutex;                  // writes into the ring cells of the patches (ghost_fill)
    seastate_implicit* solver_of(int id, int t) {return t==0 ? SP(id)->solv : tsolv[t][id].get();}
    bool remote_blocks = false;             // a parent cell on another rank (restriction through the block plans)

    // Phase 16: the faces of the covered level-0 cells at the rank boundaries. A covered cell next to a rank
    // boundary has a face spectrum towards the halo cell (the mean of the two fine cells next to the face, as the
    // faces inside the rank); the neighbour rank marks its halo cell as covered and takes the face spectrum for
    // its cell next to it, so that the cells next to a patch of another rank see the outflow of the fine grid as
    // on one rank. Sent after every level-0 halo exchange.
    vector<int> rface_send[4], rface_recv[4];   // [side]: interior cells (i,j,...) with a face, halo cells taking one
    bool rfaces = false;                        // any rank boundary with faces
    void remote_faces_plan(ghostcell*);
    void remote_faces(ghostcell*);
    void halo0(lexer*, ghostcell*);             // level-0 halo exchange, then the faces

    // Phase 17: the pipeline (pipeline_setup). Per side of the rank: the ring cells of the patches of this rank whose
    // source lies on the neighbour rank there (recv), and the sources here for the ring cells of the neighbour (send),
    // in the same order on both ranks; row: the level-0 row (side 0, 1) or column (side 2, 3) of the cell, local index
    struct kba_item
    {
        int row;
        bool face;
        int g = -2, si = 0, sj = 0;             // send: the source cell (grid, local index), g -2: none (zeros)
        seastate_amr_flink L;                   // send: the two fine cells of the face (L.g -1: the source cell itself)
        int id = -1, di = 0, dj = 0;            // recv: the ring cell (patch, local index)
        long key = 0;                           // recv: its face in the faces of patch id
    };
    vector<kba_item> ksend[4], krecv[4];
    vector<unsigned char> kface_send[4], kface_recv[4];   // [side][row]: the boundary (halo) cell of the row has a face
    bool kba = false;
    int kpos(int q, int s, int row) const;      // the position along side s in the order of quadrant q
    // Phase 17: the ring layers around the patch interiors that the sweeps read: 2 with the second-order fluxes (A 796 2)
    int ring_width() const;

    // Phase 16: the faces at the interfaces of the patches of one level (on this rank and across ranks). A ring cell
    // whose source cell (the same cell of a neighbouring patch of the level) is covered by the next finer level is
    // marked covered in its patch and takes the face spectrum of the source cell towards the patch (the mean of the
    // two fine cells next to the face), as a cell next to a covered cell inside one patch. Without this the result
    // depended on where the level is cut into patches, and so on the ranks (patches are cut at the rank boxes).
    vector<vector<seastate_amr_flink>> rlinks;  // [id]: g = id, key of the ring face, the fine cells (patch L.id)
    vector<vector<int>> rlink_of;               // [id][fill index]: index into rlinks[id], -1: none
    struct ring_plan
    {
        reefamr_xplan X;                        // the face spectra of covered source cells on other ranks
        vector<seastate_amr_flink> serve;       // served faces: the fine cells (L.id, i0..j1), L.g -1: none
        vector<pair<int,long>> dst;             // received faces: patch id, key in its faces
    };
    vector<ring_plan> rplan;                    // [l]
    void ring_faces_plan();
    void ring_face(const seastate_amr_flink&, float*);  // the mean of the two fine cells, nbin values
    void ring_faces_local(seastate_amr_patch*, int id);
    void ring_faces_remote(int l);

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
