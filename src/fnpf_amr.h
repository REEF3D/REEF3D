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

#ifndef FNPF_AMR_H_
#define FNPF_AMR_H_

#include"reefamr.h"
#include"amr_grid_lookup.h"
#include"fnpf_laplace.h"
#include<vector>
#include<fstream>

class lexer;
class fdm_fnpf;
class ghostcell;
class slice;
class slice4;
class solver;
class fnpf_fsf;
class fnpf_fsfbc;
class fnpf_sigma;
class fnpf_fsf_update;
class fnpf_bed_update;
class fnpf_laplace_cds2;
class reefmg_core;
class fnpf_body;
class solver;
class fnpf_fsf_update;

using namespace std;

//  Patch-based mesh refinement for REEF3D::FNPF, the FNPF module of REEFAMR (reefamr.h).
//
//  Level 0 is the native FNPF grid.  A refined patch refines x and y by 2 per level and,
//  with G 6 1, also the sigma layers (nested: every coarse node is a fine node).  Each
//  patch has its own lexer (horizontal geometry from the core, sigma grid and flags built
//  here), its own fdm_fnpf and its own instances of the FNPF classes (free-surface
//  discretisation, sigma transformation, Laplace assembly), so the kernels run unchanged.
//
//  Time stepping (fnpf_RK3 with the hooks stage_surface and the decorated Laplace solver):
//   - every RK stage, after level 0 has formed its stage values of eta and Fifsf, the patches
//     form theirs from their own tendencies (kinematic and dynamic FSBC), the fine values are
//     restricted into the coarser levels (averages of the 2x2 children), and the cells around
//     the patches are filled: from a patch of the same level, or interpolated from the next
//     coarser level (biquadratic; in the vertical linear in sigma)
//   - the Laplace equation is solved on all grids together: the unknowns are the nodes of
//     the leaf columns; a row takes its neighbours outside the grid from the finer level
//     (covered columns: average of the children) or the coarser level (the columns around a
//     patch, interpolated as above).  BiCGStab (reefamr_bicgstab), preconditioned by one FAC
//     sweep: a REEFMG V-cycle on level 0, then level by level the interpolated coarse
//     correction and patch-local REEFMG V-cycles on the residual it leaves
//   - one global time step, the stages of all grids in lockstep
//
//  Resolved body (fnpf_6DOF, X 10 1): every grid carries the body at its own resolution (ray
//  cast, footprint, body band); the phi and psi solves are composite solves over all grids,
//  every hull triangle is integrated on the finest grid that holds its centroid.  G 12 r
//  refines a margin r around the wetted hull at t = 0.
//
//  Several ranks: the patches are cut at the level-0 rank boxes, or with G 40 1 placed for the
//  load of the ranks; parents and old patches on other ranks are then reached through the block
//  plans and old_run of the core (restrict_sl/col, prolong_interior_sl/col, regrid_state), the
//  hull triangles are taken by the rank of the finest grid at their centroid (owns_point).  A
//  rectangle at a body is placed whole (place_whole_zones): the footprint extension and the
//  interior fill of the body are local to a patch.
//
//  Regridding (G 2 steps, default 4) only with the zone around the body (G 12): the zone is
//  aligned with x and y, and the layout is kept while it covers the flagged tiles (lazy layout,
//  reefamr_param::lazy), so that a moored body does not rebuild its patches every few steps.
//
//  Coastline (A 343 1-3): the patches are fully wet; no patch within 4 level-0 cells (A 343 2, 3:
//  8) of a dry or thin-film cell or of the coastline damping band (cell_unfit), checked at every
//  regrid; the coastline itself stays on level 0.  With A 343 1 and G 30 1 the patches may cover the
//  static coastline: the patches run fnpf_fsfbc_wd (wet mask and coastline distance from the fine
//  bed, clamp, coastline relaxation and damping as on level 0), the interpolation from a coarser
//  grid uses its wet cells only (a dry coarse cell: the mean of its wet neighbours) and the
//  restriction the wet children only.
//
//  Breaking (A 350 1): on every grid as fnpf_RK3: the implicit damping of the stage values with
//  the breaking viscosity (on a patch: its box, the cells around it filled with the coarse stage
//  values first), then the detection; the cells around a patch take the viscosity of their source
//  grid.
//
//  Scope of this version: refinement G 1 levels, G 10 boxes, G 11 boxes without
//  refinement, G 4 tile width, G 12 zone around the body; RK3 (A 310 3), wetting-drying
//  A 343 0-3 (patches off the coast; A 343 1 with G 30 1 also over it), breaking A 350 0 or 1, X 10 0 or 1, no ice (A 380 0),
//  A 324 0, A 328 0, 3D grids.
//  No refinement in the relaxation zones (B 96) and next to in- and outflow boundaries.

struct fnpf_amr_patch : public reefamr_patch
{
    fdm_fnpf *c = nullptr;
    fnpf_fsf *pf = nullptr;         // fnpf_fsfbc, or fnpf_fsfbc_wd with the coastline in the patches (wd)
    fnpf_sigma *psig = nullptr;
    fnpf_fsf_update *pfu = nullptr;
    fnpf_bed_update *pbu = nullptr;
    fnpf_laplace_cds2 *plap = nullptr;

    // RK stage values and tendencies, as fnpf_RK3
    slice4 *erk1 = nullptr, *erk2 = nullptr, *frk1 = nullptr, *frk2 = nullptr;
    slice4 *ek = nullptr, *fk = nullptr;
    slice4 *en = nullptr;           // eta at the start of the step (breaking detection, A 350 1)

    // columns with a solid neighbour (lexer index pairs): wall ghost nodes of Fi
    vector<int> wall;

    // composite Laplace
    vector<double*> kv;             // Krylov vectors, Fi layout
    vector<int> row;                // FIJK -> matrix row of the patch assembly, -1 none
    int serial = -1;                // unique number of the patch (body grids follow it)
    reefmg_core *mg = nullptr;      // patch-local multigrid of the preconditioner
    vector<int> fill1;              // the entries of fill next to the patch (the reach of the Laplace rows)
    vector<int> fill2;              // the entries of fill within two columns of the patch (the stencils of the prolongation into the next finer level)
};

class reefckpt;

class fnpf_amr : public reefamr, public amr_column_lookup
{
public:
    // checkpoints (P 48, I 45; fnpf_amr_ckpt.cpp): the fields of the patches of this rank; the restart
    // fills the patches of the hierarchy built by ini (any decomposition) from the patches of the checkpoint
    void ckpt_write(reefckpt&);
    bool ckpt_read(reefckpt&, std::string&);
    void flush_outputs() { gaugeout.flush(); logout.flush(); }   // a rebalance stop (G 46)
    // the output on the finest grid (amr_grid_lookup.h): the finest patch of this rank at (x,y)
    int finest_column(double x, double y, lexer *&q, void *&f) const override
    {
        int id;
        const int L = finest_point(x,y,id);
        if(L>0)
        {
            q = P[id]->pp;
            f = static_cast<fnpf_amr_patch*>(P[id])->c;
        }
        return L;
    }

    fnpf_amr(lexer*, fdm_fnpf*, ghostcell*);
    virtual ~fnpf_amr();

    void ini(lexer*, fdm_fnpf*, ghostcell*);
    void timestep(lexer*, fdm_fnpf*, ghostcell*);
    void print(lexer*, fdm_fnpf*, ghostcell*);

    // fnpf_RK3: stage s (0,1,2), the level-0 tendencies are formed: those of the patches
    void stage_tendency(lexer*, fdm_fnpf*, ghostcell*, int);
    // fnpf_RK3: stage s (0,1,2) values eta, Fifsf of level 0 are formed
    void stage_surface(lexer*, fdm_fnpf*, ghostcell*, slice&, slice&, int);
    // fnpf_RK3: end of the time step
    void step_end(lexer*, fdm_fnpf*, ghostcell*);

    // the Laplace solver for fnpf_RK3: level 0 assembled by plap0, then all grids together
    fnpf_laplace* laplace(fnpf_laplace *outer, fnpf_laplace_cds2 *plap0);
    bool active() const { return maxlev>0 && patches_total>0; }
    void lap_solve(lexer*, fdm_fnpf*, ghostcell*, solver*, fnpf_fsf*, double*, slice&);
    // a psi solve of the body loads on all grids: targets f[g+1] with free-surface data D[g+1];
    // boundary data are set by the caller
    void lap_solve_psi(lexer*, fdm_fnpf*, ghostcell*, solver*, fnpf_fsf*, double**, slice**);

    // resolved bodies (fnpf_6DOF) on the hierarchy
    void attach_body(fnpf_body*);
    int patches() const { return (int)P.size(); }
    lexer* patch_lexer(int n) { return P[n]->pp; }
    fdm_fnpf* patch_fdm(int n) { return FP(n)->c; }
    fnpf_fsf* patch_fsf(int n);
    slice& patch_tendency(int n, int m);
    int patch_level(int n) { return P[n]->lev; }
    int patch_serial(int n) { return FP(n)->serial; }
    void patch_walls_fi(int n, double *f) { walls_fi(*FP(n),f); }
    int finest_at(double, double);      // local grid id whose interior holds (x,y), -1: level 0
    // the hull triangle with centroid (x,y) belongs to local grid id: the finest grid that holds
    // the point is grid id of this rank (with placed patches the patch may be on another rank)
    bool owns_point(double, double, int);
    // Fi-layout arrays of the patches with need[n] from their coarser grids, coarse to fine
    // (initial guess of the fresh body grids); f[g+1] is the array of grid g.  Collective: all
    // ranks call it (the parent may be on another rank)
    void prolong_cols(double **f, const vector<char> &need);
    int layout() const { return layout_id; }    // changes when the patch set changes

    // vector space of the composite Laplace (fnpf_amr_lap.cpp, reefamr_bicgstab)
    void lap_apply(int, int, bool d1=false);
    void lap_prec(int, int);
    void lap_local(int, int, int);
    double lap_dot(int, int);
    void lap_start();
    void lap_restart();
    void lap_p(double, double);
    void lap_s(double);
    void lap_x(double, double);

    // subcycling (G 7 1, fnpf_amr_sub.cpp): level l takes 2^l steps per level-0 step
    bool sub_on() const { return sub==1 && active(); }
    // the level-0 objects of fnpf_RK3 (free surface, sigma grid, boundary conditions)
    void attach_level0(fnpf_fsf*, fnpf_sigma*, fnpf_fsf_update*, fnpf_bed_update*);
    // a psi solve of the body loads on the level-l patches (their parent columns fixed): edge
    // values of f from the parent phi_t, set by sub_psi_edges
    void lap_solve_psi_level(lexer*, ghostcell*, int, double**, slice**);
    void sub_psi_edges(int, double**);
    int sub_level_now() const { return slev; }
    bool sub_last_finest() const { return bfin==(1<<maxlev)-1; }
    double sub_time() const { return tsub; }

protected:
    // REEFAMR hooks
    reefamr_patch* patch_new() override;
    void patch_objects(reefamr_patch*, ghostcell*) override;
    void patch_delete(reefamr_patch*) override;
    void tag(int, vector<unsigned char>&) override;
    void regrid_prepare(ghostcell*) override;
    void regrid_static(ghostcell*) override;
    void regrid_state(ghostcell*, vector<reefamr_patch*>&) override;
    void regrid_finish(ghostcell*, int) override;
    void zone_bodies(vector<sixdof_obj*>&) override;
    bool cell_unfit(int, int) override;

private:
    fnpf_amr_patch* FP(int n) { return static_cast<fnpf_amr_patch*>(P[n]); }
    static fnpf_amr_patch* FP(reefamr_patch *c) { return static_cast<fnpf_amr_patch*>(c); }

    lexer* glx(int g) { return glex(g); }
    fdm_fnpf* gfd(int g) { return (g<0) ? c0 : FP(g)->c; }
    inline int fidx(lexer*, int, int, int) const;   // FIJK of lexer q (inline: fnpf_amr_fill.h)
    int gknoz(int g);                                // sigma layers of grid g

    // patch lexer: sigma grid, 3D flags, bed array
    void build_lexer3D(fnpf_amr_patch&);
    void free_lexer3D(fnpf_amr_patch&);
    void build_walls(fnpf_amr_patch&);
    void walls_fi(fnpf_amr_patch&, double*);

    // stage values of a grid: s 0,1: erk/frk of stage s, 2: eta/Fifsf; m 0 eta, 1 Fifsf
    slice& sval(int g, int s, int m);
    slice *Se0, *Sf0;               // level-0 stage values of the current stage
    int stg;

    // interpolation from the coarser grid g: cell (ic,jc), quadrant (ox,oy)
    double pq(slice&, lexer*, int, int, int, int);
    void pcol(int, int, int, int, int, const double*, int, double*);
    void pweights(lexer*, int, int, int, int, double*);

    // fill and restriction of slices and Fi-layout arrays
    template<class SEL> void fill_sl(int, int, int, SEL);
    template<class SEL> void fill_col(int, int, SEL, bool d1=false);
    template<class SEL> void restrict_sl(int, SEL);
    bool rcubic(lexer*, int, int);
    template<class F> double rc4(F);
    int rorder = 4;                 // restriction: 2 average of the children, 4 cubic
    int pord = 4;                   // prolongation: 3 biquadratic, 4 bicubic
    template<class SEL> void restrict_col(SEL);
    template<class ND, class SEL> void prolong_interior_sl(int, ND, int, SEL);
    template<class ND, class SRC, class DST> void prolong_interior_col(int, ND, SRC, DST);
    int klev(int l) const;          // sigma layers of level l
    void walls_sl(fnpf_amr_patch&, slice&, int);

    int layout_id = 0;
    int serial_next = 0;

    // composite Laplace (fnpf_amr_lap.cpp)
    double* lvec(int, int);         // grid, vector (-1: the target of the solve, ltgt)
    vector<double*> ltgt;           // [g+1]: unknowns of the current solve (Fi or psi)
    void lap_core(lexer*, ghostcell*);
    void lap_prepare(ghostcell*);
    void lap_sync(int, bool d1=false);  // d1: the columns next to the patches only (the operator)
    void lap_rows();
    int lap_layout = -1;
    int lap_wkey = 0;               // level window of the row lists (G 7 1)
    struct lgrid
    {
        vector<int> lq, lr;         // leaf unknowns: Fi index, matrix row
        vector<int> cq;             // covered unknowns: Fi index
        vector<int> aq, ar;         // all leaf rows (lq, lr: without the fixed rows)
    };
    vector<lgrid> lg;               // [g+1]
    vector<double*> kv0;            // level-0 Krylov vectors
    reefmg_core *mg0 = nullptr;     // level-0 multigrid of the preconditioner
    vector<int> rowmap0;            // level 0: FIJK -> matrix row
    long lap_it_total, lap_solves;
    int lap_it_last;
    double lap_res_last;
    int lap_kind = 0;               // solve of lap_core: 0 phi, 1 psi (body loads)
    int lap_it_phi_max = 0, lap_it_psi_max = 0, lap_solves_step = 0;   // this step (log)
    long lap_it_step = 0;
    int lap_capped = 0;             // solves that reached N 46
    int lap_restarts = 0;           // BiCGStab restarts after stagnation
    int lap_stalls = 0;             // solves stopped after a second stagnation (STAGTOL)

    // output
    void write_vtr(lexer*, fnpf_amr_patch&, int);
    void write_vtr0(lexer*, fdm_fnpf*);
    bool print_lagoon(lexer*, fdm_fnpf*, ghostcell*);
    void gauges(lexer*, fdm_fnpf*, ghostcell*);
    ofstream gaugeout, logout;

    // subcycling (G 7 1, fnpf_amr_sub.cpp)
    int sub = 0;
    double zr_user = 0.0;           // G 12 as given (G 7 1 with a body widens the margin)
    int wlo = 0, whi = -1;          // level window of the Laplace solve, restriction and fills (whi<0: maxlev)
    int wtop() const { return whi<0 ? maxlev : whi; }
    bool lap_edge = false;          // window above level 0: the parent columns of its lowest level are fixed
    bool lap_dir = false;           // fill_col of a Krylov vector (the fixed parent columns are 0)
    int slev = 0;                   // level of the current subcycled step
    int bfin = 0;                   // finest steps done in the level-0 step
    double tsub = 0.0;              // start time of the current step of slev
    vector<double> dtlev;           // step of each level in the current level-0 step
    struct tcol
    {
        int layout = -1;
        vector<int> ci, cj;         // the columns of the grid the fills of the next level read
        vector<double> old, live;   // their eta, Fifsf, Fi at the start of the step / current
        double *pt = nullptr;       // phi_t of the step at these columns (psi edges of the body)
        int ptn = 0;
    };
    vector<tcol> tc;                // [g+1]
    fnpf_fsf *pf0 = nullptr;
    fnpf_sigma *psig0 = nullptr;
    fnpf_fsf_update *pfu0 = nullptr;
    fnpf_bed_update *pbu0 = nullptr;
    int tc_nv(int g);
    void tc_build(int);
    void tc_pack(int, int, double*);
    void tc_unpack(int, int, const double*);
    void sub_snapshot(int);
    void sub_swap_in(int, double);
    void sub_swap_out(int);
    void sub_step(lexer*, fdm_fnpf*, ghostcell*);
    void sub_level(lexer*, ghostcell*, int, int, double, double);
    void sub_stage(lexer*, ghostcell*, int, int, int, double);
    void sub_sync(lexer*, ghostcell*, int, double);
    void sub_derived(lexer*, ghostcell*, int, double);
    void lap_prec_win(int, int);
    double rko(int s) const { return (s==1) ? 0.5 : 1.0; }

    // wave breaking on the patches (A 350 1): the implicit damping of the stage
    // values with the breaking viscosity on the patch box (the cells around it held), detection on
    // the patch, the viscosity of the cells around the patches from the source grid
    bool brk = false;
    void patch_damping(fnpf_amr_patch&, ghostcell*, slice&, int, double);
    void patch_breaking(int, ghostcell*, int, double);
    vector<double> dv[6];           // work vectors of the damping solve

    // coastline (A 343 1-3): the patches stay off the dry cells and the coastline band of level 0
    double dist_coast = 0.0;
    bool thin_warned = false;

    // A 343 1 with G 30 1: the patches may cover the static coastline: fnpf_fsfbc_wd on the patches
    // (wet mask and coastline distance from the fine bed, the clamp and the coastline relaxation as
    // on level 0); the interpolation from a coarser grid and the restriction use wet cells only (wfill)
    bool wd = false;
    bool wfill = false;
    void patch_coast(fnpf_amr_patch&, ghostcell*, slice&, slice&);
    void patch_wd_ini(fnpf_amr_patch&, ghostcell*);

    // the damping and the coastline of the stage values read the cells around a patch: these are
    // filled with the stage values of the coarser grid first (not the lagged ones of the last fill)
    bool stage_ring() const { return brk || wd; }
    void stage_damp_coast(fnpf_amr_patch&, ghostcell*, int, int);

    fdm_fnpf *c0;
    fnpf_laplace_cds2 *plap0 = nullptr;
    fnpf_body *body = nullptr;
    int vref;
    int gcval_eta, gcval_fifsf;
    double printtime_amr;
    int printcount_amr;
    double tm[6];
};

#endif
