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

#ifndef NHFLOW_AMR_H_
#define NHFLOW_AMR_H_

#include"reefamr.h"
#include"amr_grid_lookup.h"
#include"lexer.h"
#include"nhflow_momentum_func.h"
#include"nhflow_convection.h"
#include"nhflow_timestep.h"
#include"sediment_f.h"
#include<vector>
#include<fstream>
#include<unordered_map>

class lexer;
class fdm_nhf;
class ghostcell;
class slice;
class nhflow_signal_speed;
class nhflow_reconstruct;
class nhflow_convection;
class nhflow_diffusion;
class nhflow_pressure;
class nhflow_fsf;
class nhflow_forcing;
class nhflow_turbulence;
class vrans_nhflow;
class sixdof;
class sediment;
class sixdof_nhflow;
class solver;
class vec;
class matrix_diag;
class sixdof_obj;
class nhflow_momentum_SI;
class ioflow;
class ioflow_amr_block;
class patchBC_interface;
class reefmg_core;
class solver;

using namespace std;

//  Patch-based mesh refinement for REEF3D::NHFLOW, the NHFLOW module of REEFAMR (reefamr.h:
//  hierarchy, patch lexers, fill and exchange plans, restriction map, flux matching registry).
//
//  Level 0 is the native NHFLOW grid.  A refined patch refines x and y by 2 per level, with the
//  same sigma layers, or with G 6 1 the sigma layers of its parent halved (nested: the coarse
//  nodes are nodes of the fine grid).  Each patch has its own lexer (horizontal geometry from the core, 3D flags
//  and the sigma arrays built here), its own fdm_nhf and its own instances of the NHFLOW classes
//  (momentum RK2/RK3, reconstruction, HLL/HLLC, free surface, pressure), so the kernels run
//  unchanged on every grid.
//
//  Time stepping: the level-0 momentum object hands the time step to this class (stage runner).
//  Per RK stage:
//   - the cells around the patches are filled with the stage input (coarse to fine): a patch of
//     the same level, or interpolated from the next coarser level (bicubic in x,y, layer by
//     layer; eta, U, V, W and P, then WL = eta + depth of the patch and UH = WL U ...; with
//     G 6 U, V, W linear in sigma within each coarse layer, minmod-limited slope, mean
//     preserving, and P at the midpoints between the coarse nodes linear)
//   - continuity and momentum part of the stage (phase_F, phase_M), finest patches first: the
//     flux hook of a patch records the face fluxes on its box faces (FEx/FEy, Fx/Fy of UH, VH,
//     WH, per sigma layer, and dfx/dfy), a coarse face next to a finer patch takes the mean of
//     its two fine faces (G 6: and of both fine layers) -> mass and momentum conserved across
//     the interfaces
//   - restriction of WL, eta, UH, VH, WH to the covered coarse cells (2x2 averages per layer,
//     G 6: 2x2x2)
//   - pressure projection on all grids together: the Poisson rows of every grid (nhflow_poisson,
//     assembled per grid), the unknowns are the leaf nodes; BiCGStab with a FAC preconditioner
//     (REEFMG V-cycle on level 0, patch-local reefmg_core), as the composite Laplace of FNPF AMR
//     (nhflow_amr_press.cpp)
//   - restriction of the corrected UH, VH, WH and P, relaxation zones and ghost cells
//   - one global time step from the finest grid (nhflow_timestep with the patch hook)
//  With G 7 1 (subcycling, nhflow_amr_sub.cpp) level l takes 2^l steps per level-0 step: level 0
//  steps alone with its own pressure, every finer level after it with fills interpolated in time,
//  flux registers, a pressure solve per level and stage, refluxing and a synchronisation projection
//  over the finer levels after each of its steps.  Floating bodies need the zone (G 12): the finest
//  level advances the body with its stages and loads, the coarser levels a predicted copy.
//
//  Floating bodies (6DOF_nhflow, X 10 1/2): the body is advanced on level 0, before the patches
//  take their forcing (G 7 1: by the finest level); every patch casts the hull on its own sigma
//  grid and adds the direct forcing of the rigid-body velocity (nhflow_amr_6dof); the loads are
//  integrated once, every hull triangle on the finest grid at its centroid (pressure, free
//  surface, shear).  G 12 refines around the wetted hull (margin G 12, rectangle aligned with x
//  and y; with G 13 L a oriented along the motion plus a wake wedge of length L and half angle
//  a, as SFLOW); the hull triangles (X 185) are then sized for the finest level, as on a uniform
//  fine grid.  The zone follows the
//  body: a regrid every G 2 steps at the end of the step (G 2 0: static), the layout kept
//  while it covers the flagged tiles with at most 50 % excess, G 5 regrids of hysteresis.  A
//  fresh patch takes the state of the old patches of its level where they overlap; elsewhere it
//  is prolonged from its parent, conservatively (every 2x2 block keeps the water level and,
//  layer by layer, the momentum of its coarse cell).
//
//  Solution-adaptive refinement (as sflow_amr): a cell of level l flags level l+1 where its surface
//  jumps by more than G 20 to a neighbour or its second difference along x or y exceeds G 21;
//  G 3 buffer cells, regrid every G 2 steps with the machinery of the body zone (lazy layout,
//  G 5 hysteresis).  The hierarchy may be empty: the first patches appear with the first flags.
//  At t = 0 the patches take the initial water level boxes (F 72) on their own grid.
//
//  Scope of this version: static refinement boxes, the body zone and the adaptive flags (G 1
//  levels, G 10 boxes, G 11 boxes without refinement, G 4 tile width, G 3, G 20, G 21, G 12,
//  G 13, G 2, G 5, G 6, G 22, G 23, G 30, G 40), A 510 2/3, A 511 1/2, A 514 all, A 520 0/1/2,
//  A 512 0/1/2, A 560 0, A 550 0 (A 550 1 with A 510 2), B 200 0, X 10 0/1/2 (X 60 1, X 16 0,
//  A 516 0/1/3), S 10 0, no solids (A 580 1,
//  A 581-590), no membranes (X 330), nets (X 320), 3D grids.  Patches stay out of the relaxation
//  zones (B 96), the in- and outflow band and, checked at every regrid, 4 level-0 cells away from
//  dry or shallow cells (they are fully wet; level 0 keeps its wetting and drying).
//
//  Breaking (A 550 1, RK2): every grid detects breaking; the implicit diffusion with the breaking
//  viscosity (A 512 2) is one composite solve over the leaf cells of all grids (G 31 1, the
//  default: nhflow_amr_diff.cpp, a patch split in pieces gives the same result as the whole
//  patch), or with G 31 0 every grid solves its own (a bicgstab_ijk per patch, its global sums
//  local to the patch inside pscope).  The viscous flux across a coarse-fine interface is not
//  matched.  G 23 1 flags the cells with breaking viscosity for the adaptive mode.
//
//  Wetting and drying in the patches (G 30 1): the patches may cover dry and shallow cells, the
//  NHFLOW wetting and drying (A 540) runs on every grid.  The coupling: the cells around a patch
//  take the flags of their source cell (a fine face on the patch box carries mass only where the
//  coarse face can) and are interpolated from wet source cells only (pressure: wet and deep), a
//  dry source cell gives a dry cell on the fine bed; a fresh patch is well balanced at the
//  shoreline (a 2x2 block with a dry child keeps the surface of its wet parent, without the
//  conservative shift; dry children have no momentum) and takes its flags from its water level
//  (or from the old patch); a covered coarse cell is wet if all its children are, deep follows the rule of
//  nhflow_fsf_f::wetdry; the restriction of the pressure uses the wet and deep children.

struct nhflow_amr_patch : public reefamr_patch
{
    fdm_nhf *d = nullptr;
    nhflow_momentum_func *pmom = nullptr;
    nhflow_signal_speed *pss = nullptr;
    nhflow_reconstruct *precon = nullptr;
    nhflow_convection *pconv = nullptr;
    nhflow_diffusion *pdiff = nullptr;
    solver *psolv = nullptr;
    nhflow_pressure *ppress = nullptr;
    nhflow_fsf *pfsf = nullptr;
    nhflow_forcing *pdf = nullptr;
    nhflow_turbulence *pturb = nullptr;
    vrans_nhflow *pvrans = nullptr;
    sixdof *p6dof = nullptr;
    sediment *psed = nullptr;
    sediment_f *sed = nullptr;      // S 10: the sediment transport of the grid (nhflow_amr_sed.cpp)
    vector<double> sedrec[4];       // S 10: Exner face fluxes on the box faces, per side (as rec)
    ioflow *pflow = nullptr;
    ioflow_amr_block *pio = nullptr;    // multi-block: in- and outflow of a block 1.. through its open faces (G 54)
    patchBC_interface *pBC = nullptr;
    nhflow_stage_obj S;
    vector<int> wfix;               // G 30: flags of the filled cells, kept through wetdry (lexer::wetfix)
    vector<double> vbfill;          // A 550: breaking viscosity of the source of the filled cells (lexer::amrvb)

    // box faces recorded by the flux hook: rec[ipol][side][r*knoz+k], r fine face index along the
    // side, k layer of the patch; side 0 low x, 1 high x, 2 low y, 3 high y; ipol 0: dfx/dfy (r
    // only)
    vector<double> rec[5][4];

    // composite pressure
    vector<double*> kv;             // Krylov vectors, F layout
    vector<int> row;                // FIJK -> matrix row of the patch assembly, -1 none
    reefmg_core *mg = nullptr;      // patch-local multigrid of the preconditioner
    vector<int> fill1;              // the entries of fill next to the patch (the reach of the pressure rows)
    vector<int> fill2;              // the entries of fill within two columns of the patch (the stencils of the prolongation into the next finer level)
    // the overlap of the patch multigrid (pr_overlap): columns across the sides x-, x+, y-, y+ where
    // a patch of the level is next to it; mo: those of the multigrid as built; kc7: the matrix rows
    // of the patch (p n s w e t b), filled into the columns around it from the patches there
    int ov[4] = {0,0,0,0};
    int mo[4] = {-1,-1,-1,-1};
    vector<double*> kc7;
    double *ptgt = nullptr;         // unknowns of the current solve (P or PCORR)

    // composite implicit diffusion (G 31 1): Krylov vectors, cell layout
    vector<vector<double>> dkv;

    // G 7 1: fine face fluxes summed over the steps of the patch within one step of its parent
    // (dt * RK weight), freg[ipol][side][r*knoz+k], ipol 1-4
    vector<double> freg[5][4];
};

class reefckpt;

class nhflow_amr : public reefamr, public amr_column_lookup, public nhflow_stage_runner, public nhflow_flux_hook, public nhflow_timestep_hook, public sediment_grids
{
public:
    // checkpoints (P 48, I 45; nhflow_amr_ckpt.cpp): the arrays of nhflow_ckpt_fields and of the
    // turbulence model per patch (reefamr_ckpt.cpp)
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
            f = static_cast<nhflow_amr_patch*>(P[id])->d;
        }
        return L;
    }

    nhflow_amr(lexer*, fdm_nhf*, ghostcell*, nhflow_momentum*, nhflow_convection*, nhflow_timestep*, sixdof*);
    virtual ~nhflow_amr();

    void ini(lexer*, fdm_nhf*, ghostcell*);
    // S 10: sediment transport on all grids (nhflow_amr_sed.cpp): the level-0 sediment_f hands its step over
    void sed_attach(sediment*);
    void sediment_step_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*, sediment_f*) override;
    void print(lexer*, fdm_nhf*, ghostcell*);
    bool active() const { return (maxlev>0 && patches_total>0) || nblk>0; }

    // the time step of all grids (nhflow_momentum_RK2/RK3::start)
    void step(lexer*, fdm_nhf*, ghostcell*, nhflow_momentum_func*, nhflow_stage_obj&) override;

    // hydrostatic semi-implicit free surface on the multi-block grid (A 534 1, nhflow_amr_si.cpp): the
    // parts of the step of every block (nhflow_momentum_SI), the cells around the blocks filled in
    // between, one eta system over the cells of all blocks (BiCGStab with the Jacobi preconditioner)
    void si_step(lexer*, fdm_nhf*, ghostcell*, nhflow_momentum_func*, nhflow_stage_obj&);

    // vector space of the composite eta solve (reefamr_bicgstab)
    void si_apply(int, int);
    void si_prec(int, int);
    double si_dot(int, int);
    void si_opstart();
    void si_op_p(double, double);
    void si_op_s(double);
    void si_op_x(double, double);
    void si_restart();

    // face fluxes of grid id (0: level 0, n: patch n-1) before the divergence
    void flux_hook(lexer*, fdm_nhf*, int, int, double*, double*) override;

    // time step: patch maxima and cell sizes
    double dt_local_max(int) override;
    void dt_cell_size(int, double&, double&) override;
    void dt_si_cmin(double, double&, double&) override;

    // vector space of the composite pressure solve (nhflow_amr_press.cpp, reefamr_bicgstab)
    void pr_apply(int, int, bool d1=false);
    void pr_prec(int, int);
    double pr_dot(int, int);
    void pr_start();
    void pr_p(double, double);
    void pr_s(double);
    void pr_x(double, double);

    // vector space of the composite implicit diffusion (nhflow_amr_diff.cpp, reefamr_bicgstab)
    void df_apply(int, int);
    void df_prec(int, int);
    double df_dot(int, int);
    void df_start();
    void df_p(double, double);
    void df_s(double);
    void df_x(double, double);

    // the matrix of a grid as nhflow_idiff assembles it (capture solver of phase_D)
    void df_capture(double*, vec&, matrix_diag&, int);

protected:
    // REEFAMR hooks
    reefamr_patch* patch_new() override;
    void patch_objects(reefamr_patch*, ghostcell*) override;
    void patch_delete(reefamr_patch*) override;
    void tag(int, vector<unsigned char>&) override;
    void regrid_prepare(ghostcell*) override;
    void regrid_ids() override;
    void regrid_static(ghostcell*) override;
    void regrid_state(ghostcell*, vector<reefamr_patch*>&) override;
    void regrid_finish(ghostcell*, int) override;
    void zone_bodies(vector<sixdof_obj*>&) override;
    bool cell_unfit(int, int) override;
    bool cell_unfit_patch(reefamr_patch*, int, int) override;

private:
    // a patch kernel runs: ghostcell exchange off and the fdm of the ghostcell on the patch (the
    // V-type wall conditions of ghostcell read eta, dfx, UH, U ... of the fdm it holds)
    struct pscope
    {
        ghostcell *g;
        fdm_nhf *dl0;
        bool old, oldlocal;
        pscope(ghostcell*, fdm_nhf*, fdm_nhf*);
        ~pscope();
    };

    nhflow_amr_patch* NP(int n) { return static_cast<nhflow_amr_patch*>(P[n]); }
    static nhflow_amr_patch* NP(reefamr_patch *c) { return static_cast<nhflow_amr_patch*>(c); }

    fdm_nhf* gfd(int g) { return (g<0) ? d0 : NP(g)->d; }
    nhflow_momentum_func* gmom(int g) { return (g<0) ? mom0 : NP(g)->pmom; }
    int fidx(lexer *q, int ii, int jj, int kk) const    // FIJK of lexer q
    { return (ii-q->imin)*q->jmax*q->kmaxF + (jj-q->jmin)*q->kmaxF + kk - q->kmin; }
    int cidx(lexer *q, int ii, int jj, int kk) const    // IJK of lexer q
    { return (ii-q->imin)*q->jmax*q->kmax + (jj-q->jmin)*q->kmax + kk - q->kmin; }

    // patch lexer: 3D flags, boundary list, sigma arrays
    void build_lexer3D(nhflow_amr_patch&);
    void open_lists(nhflow_amr_patch&, bool, bool);     // G 54: the 3D in- and outflow cells of the open faces of a block
    void free_lexer3D(nhflow_amr_patch&);

    // interpolation from the coarser grid g: cell (ic,jc), quadrant (ox,oy); the last argument
    // (G 30) selects the coarse cells the stencil may use: 0 fluid, 1 wet (surface, velocities),
    // 2 wet and deep (pressure)
    void pweights(lexer*, int, int, int, int, double*, int=0);
    double pq(slice&, lexer*, int, int, int, int, int=0);
    double pq3(const double*, lexer*, int, int, int, int, int, bool, int=0);
    double pqw(slice&, lexer*, int, int, const double*);
    double pq3w(const double*, lexer*, int, int, int, bool, const double*);
    double plin(slice&, lexer*, int, int, int, int);
    void pcol(int, int, int, int, int, const double*, int, double*, int=0);
    void vcell(lexer*, const double*, int, int, double*);

    // vertical refinement (G 6): the sigma layers doubled on every level (vr 2), nested
    int vr = 1;
    int klev(int l) const { return p0->knoz*((vr==2) ? (1<<l) : 1); }

    // stage input arrays of grid g at stage s (s<0: end of the step)
    struct stg { slice *WL; double *UH,*VH,*WH; };
    stg stage_in(int g, int s);
    stg stage_out(int g, int s);

    // fill, restriction
    void fill_stage(ghostcell*, int, int);
    void fill_bed(ghostcell*, int);
    void bc_patch(ghostcell*, nhflow_amr_patch&, int);
    void restrict_surface(int);
    void restrict_momentum(int, bool);
    void halo0(ghostcell*, int);
    void exchange_fluxes(int);
    vector<vector<double>> rval;    // [target grid id+1][rmatch index * NF]: fine face values from other ranks
    int NF = 0;                     // values per face entry: 4 variables x layers + dfx
    void prolong_patch(ghostcell*, int);         // fresh patches of level l
    void ini_boxes(int);
    template<class SEL> void fill_col(int, int, SEL);
    template<class SEL> void restrict_col(SEL);
    void rcol_block(reefamr_patch*, int, const double*, int, double*);
    int pord = 4;                   // prolongation: 3 biquadratic, 4 bicubic

    // composite pressure (nhflow_amr_press.cpp)
    void press_solve(lexer*, ghostcell*, int);
    double* pvec(int, int);         // grid, vector (-1: the unknowns of the solve)
    void pr_core(lexer*, ghostcell*);
    void pr_prepare(ghostcell*);
    void pr_sync(int, bool d1=false);   // d1: the columns next to the patches only (the operator)
    void pr_rows();
    int pr_layout = -1;
    struct pgrid
    {
        vector<int> lq, lr;         // leaf unknowns: F index, matrix row (without the fixed rows)
        vector<int> cq;             // covered unknowns: F index
        vector<int> aq, ar;         // all leaf rows
        vector<int> aw;             // their column (slice index), for the wet flag
        vector<int> fq;             // identity rows of wet columns (shallow): P = 0 after the solve
    };
    vector<pgrid> pg;               // [g+1]
    vector<double*> kv0;            // level-0 Krylov vectors
    reefmg_core *mg0 = nullptr;     // level-0 multigrid of the preconditioner
    vector<int> rowmap0;            // level 0: FIJK -> matrix row
    double *ptgt0 = nullptr;
    long pr_it_total = 0, pr_solves = 0;
    int pr_it_last = 0;
    double pr_res_last = 0.0;
    int layout_id = 0;

    // stencils and index lists of the composite solve, rebuilt in pr_prepare: they depend on the
    // layout and on the flags, wet and deep, which do not change during a solve.  The same
    // arithmetic in the same order as pcol and restrict_col.
    struct pstencil { int nw; int off[25]; double w[25]; };     // pcol: source offsets, weights
    void pst_make(int, int, int, int, int, pstencil&);
    void pst_col(const pstencil&, const double*, int, int, double*) const;
    struct prblock { int mode, src; double wa[4]; };                // restriction of a 2x2 block
    vector<vector<prblock>> pr_rbk;         // [patch id][block]
    vector<vector<pstencil>> pr_pik;        // [l][4 key + child]: interior prolongation from the parent
    vector<vector<char>> pr_pik_ok;         // made
    std::unordered_map<const reefamr_fill*, pstencil> pr_fs;    // fills from the coarser grid
    vector<long> pr_l0a, pr_l0z;            // level 0: active MG index, F index; inactive F index
    vector<int> pr_l0f;
    struct pact { long lq; int qq, rw; };
    vector<vector<pact>> pr_pa;             // [patch id]: active unknowns of the patch MG
    vector<vector<pact>> pr_po;             // [patch id]: active overlap cells of the patch MG (qq in the patch arrays)
    bool pr_ovl = false;                    // some patch of some rank has overlap columns
    void pr_overlap();                      // the overlap columns of the patches (layout)
    void pr_kc7_fill();                     // the rows of the patches into the columns around them (siblings)
    template<class SEL> void pr_fill_sib(int, int, int, SEL);   // a fill from the patches of the level only
    bool pr_sib(const reefamr_patch*, int, int);                 // column (di,dj) of the patch lies in a patch of its level
    void pr_stencils();
    template<class SEL> void pr_restrict(SEL);
    template<class SEL> void pr_fill(int, int, SEL, bool d1=false);
    template<class SEL> void pr_fill0(SEL);  // multi-block: the columns around the blocks and the ring
    int npass0 = 2;                          // multi-block: passes of the block-local V-cycles (G 51)
    // multi-block: coarse space across the blocks (G 53 w): aggregates of w x w columns of every
    // level-0 grid, the coarse operator summed from the rows, factorized on every rank; the
    // preconditioner corrects with it before the block-local V-cycles
    struct cs_impl;
    cs_impl *cs = nullptr;
    int csw = 0;
    int cs_n = 0;
    void cs_build(ghostcell*);
    void cs_correct(int, int);
    void cs_free();
    void pr_prolong(int, int);               // level, vector

    // composite implicit diffusion (G 31 1 with A 512 2, nhflow_amr_diff.cpp): phase_D of every
    // grid assembles the rows of one component (a capture solver records them), one BiCGStab
    // over the leaf cells of all grids, then phase_M with a diffusion object that keeps the result
    bool cdiff = false;
    void diff_solve(lexer*, ghostcell*, int);

    // breaking with G 31 1 (A 550 1): the patches run the detection, the flags and the viscosity
    // of nhflow_breaking as three parts after phase_F of all grids; between them the cells around
    // a patch that come from a sibling patch take its detection and flags of this stage (a patch
    // split in pieces breaks as the whole patch)
    bool brk_split = false;
    void breaking_patches(ghostcell*);
    void brk_exchange(int, int);
    void df_free();
    void df_core(lexer*);
    void df_sync(int);
    double* dvec(int, int);         // grid, vector (-1: the unknowns of the solve)
    struct dgrid
    {
        double *x = nullptr;        // UHDIFF, VHDIFF or WHDIFF of the grid
        vec *rhs = nullptr;
        matrix_diag *M = nullptr;
        int var = 0;
        vector<int> lq, lr;         // leaf cells: IJK, matrix row
        vector<int> aq;             // all computed cells (LOOP): IJK
    };
    vector<dgrid> dg;               // [g+1]
    int dcur = -1;                  // grid whose phase_D runs
    vector<vector<double>> dkv0;    // level-0 vectors
    struct dstencil { int nw; int off[25]; double w[25]; };     // fills from the coarser grid
    std::unordered_map<const reefamr_fill*, dstencil> df_fs;
    int df_stage = 0;
    solver *dcap = nullptr;         // capture solver
    nhflow_diffusion *dkeep = nullptr;  // diffusion object of phase_M: keeps UHDIFF ...
    long df_it_total = 0, df_solves = 0;
    int df_it_last = 0;
    double tdiff = 0.0;

    // wetting and drying in the patches (G 30 1): wet-aware interpolation, the flags of fresh
    // patches and of the covered coarse cells
    bool shore = false;
    int nshore = 0;
    bool flagbreak = false;         // breaking flag (G 23): cells with breaking viscosity                 // shoreline flag (G 22): cells within nshore cells of the other wet state
    bool wet_at(lexer*, int, int, int) const;
    void dry_cell(lexer*, fdm_nhf*, slice&, double*, double*, double*, int, int);
    void patch_flags(nhflow_amr_patch&);
    void restrict_flags(ghostcell*, int);
    void deep_rule(lexer*, slice&);

    // solution-adaptive flags (G 20 surface jump, G 21 second difference) and the regrid of a
    // step
    double tol_eta = 0.0, tol_curv = 0.0;
    bool adaptive = false;
    void regrid_step(lexer*, ghostcell*);

    // floating bodies (X 10 1/2): the level-0 bodies, the finest grid at (x,y) (local grid id, -1:
    // level 0), the hull on a patch, the loads from the finest grids
    sixdof_nhflow *b6 = nullptr;
    int finest_at(double, double);
    void body_patch(ghostcell*, nhflow_amr_patch&);
    void body_loads(lexer*, ghostcell*);

    // G 7 1: subcycling (nhflow_amr_sub.cpp)
    int sub = 0;                    // 1: subcycled
    int hstage = 0;                 // RK stage of the grids running (flux register weights)
    int wlo = 0, whi = -1;          // level window of the restriction and of the composite solve (whi<0: maxlev)
    int wtop() const { return whi<0 ? maxlev : whi; }
    bool pr_edge = false;           // composite solve with its lowest level > 0: the parent columns are fixed
    bool pr_dir = false;            // pr_fill of a Krylov vector (the fixed parent columns are 0)
    double zr_user = 0.0;           // G 12 as given (G 7 1 with bodies widens the zone margin)
    int pr_wkey = -1;               // window of the current row lists
    struct tcol                     // parent columns read by the fills of the next level, their start-of-step values
    {
        vector<int> ci, cj;
        vector<double> old, live;
        int layout = -1;
    };
    vector<tcol> tc;                // [g+1]
    int tc_nv(int g);
    void tc_build(int);
    void tc_pack(int, int, double*);
    void tc_unpack(int, int, const double*);
    void sub_snapshot(int);
    void sub_swap_in(int, double);
    void sub_swap_out(int);
    void sub_fill(ghostcell*, int, int, double);
    vector<vector<double>> cregL, cregR, rflux;   // [id+1]: coarse face fluxes summed over the step (4 K per entry), received fine sums
    vector<double> dtlev;           // step of each level in the current level-0 step
    void sub_creg_reset(int);
    void sub_dfx();
    void sub_step(lexer*, fdm_nhf*, ghostcell*, nhflow_momentum_func*, nhflow_stage_obj&);
    void sub_level(lexer*, ghostcell*, int, int, double, double);
    void sub_sync(lexer*, ghostcell*, int);
    void sub_reflux(int);
    void sub_press_level(lexer*, ghostcell*, int, int, int);
    void sub_project(lexer*, ghostcell*, int);
    int bfin = 0;                   // finest steps of the level-0 step done (the bodies)
    void sub_body_stage(lexer*, ghostcell*, int, int, double, double);
    void sub_body_end(lexer*, ghostcell*);
    void sub_body_mask(lexer*, fdm_nhf*);
    long sub_lv_it = 0, sub_lv_n = 0, sub_sy_it = 0, sub_sy_n = 0;   // level and synchronisation solves: iterations, count
    double tsync = 0.0;
    double rkw(int) const;          // weight of stage s in the step
    double rkc(int) const;          // time of the input of stage s in the step
    double rko(int) const;          // time of the output of stage s in the step

    // multi-block level 0 (G 50 1 with DIVEMesh R 20 1, nhflow_amr_block.cpp): blocks 1.. are level-0
    // patches (lev[0]), block 0 the native grid with its ring; before every stage the cells around the
    // blocks and the ring take the stage input of the cells they stand for
    void blocks_check(lexer*);
    void blocks_ini(lexer*, ghostcell*);
    void fill_level0(ghostcell*, int);
    void fill_bed0(ghostcell*);
    void blocks_sigma(ghostcell*);
    void ring_halo(ghostcell*, int);
    void blocks_phase(ghostcell*, int, int);
    void dump_grids(lexer*);
    void potential(lexer*, ghostcell*);     // I 11: potential flow initialisation over the blocks
    bool ring0(int ii, int jj) const { return nblk>0 && nring[(ii-p0->imin)*p0->jmax + (jj-p0->jmin)]==1; }

    // output
    void write_vtr(lexer*, nhflow_amr_patch&, int);
    void write_vtr0(lexer*, fdm_nhf*);
    bool print_lagoon(lexer*, fdm_nhf*, ghostcell*);
    void gauges(lexer*, fdm_nhf*, ghostcell*);
    double mass(lexer*, fdm_nhf*, ghostcell*);
    ofstream gaugeout, logout;

    fdm_nhf *d0;
    nhflow_momentum_func *mom0;
    nhflow_stage_obj *S0p = nullptr;

    // A 534 1 (nhflow_amr_si.cpp): the grids of the composite eta system (-1 the native grid, then
    // lev[0]), their own rows (i, j, row of the si_eta2D assembly) and the Krylov vectors per grid
    std::vector<int> sigrid;
    std::vector<std::vector<int>> sirow;
    std::vector<std::vector<slice*>> sivec;
    nhflow_momentum_SI* gsi(int);
    bool si_own(int, int, int);
    void si_setup();
    template<class F> void si_each(ghostcell*, F);
    template<class A> void si_fill3(ghostcell*, int, int, A);
    void si_fillv(ghostcell*, int);
    void si_fillF(ghostcell*);
    void si_dirichlet(ghostcell*);
    void si_solve(lexer*, ghostcell*);
    ghostcell *sipgc = nullptr;
    ioflow *pflowv;                 // objects shared by all patches (their NHFLOW calls do nothing)
    patchBC_interface *pBCv;
    sediment *psedv;
    nhflow_convection *pconv0;
    nhflow_timestep *pstep0;
    int cur_stage = -1;             // RK stage of the time step, -1: outside
    int gcval_eta;
    int bc5, bc6;                   // boundary types of the bed and the free surface (gcb4)
    double printtime_amr;
    int printcount_amr;
    double m0 = 0.0;
    double tm[8];
    
    // S 10 (nhflow_amr_sed.cpp)
    sediment_f *sed0 = nullptr;     // the sediment of level 0
    sediment_f* gsed(int g) { return (g<0) ? sed0 : NP(g)->sed; }
    void sed_patch_ini(ghostcell*, nhflow_amr_patch&);
    void sed_bed_ghosts(ghostcell*);
    vector<vector<double>> sedr;    // remote matched faces: the fine flux of rmatch[g][n]
};

#endif
