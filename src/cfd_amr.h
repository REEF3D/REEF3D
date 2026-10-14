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

#ifndef CFD_AMR_H_
#define CFD_AMR_H_

#include"reefamr3d.h"
#include"amr_grid_lookup.h"
#include<vector>
#include<functional>

class lexer;
class fdm;
class ghostcell;
class field;
class field4;
class convection;
class diffusion;
class solver;
class pressure;
class pjm_corr;
class poisson;
class turbulence;
class heat;
class concentration;
class reini;
class ioflow;
class patchBC_interface;
class vrans;
class fsi;
class sixdof;
class sixdof_cfd;
class cfd_amr_6dof;
class momentum;
class momentum_rk;
class reefmg_core;
class initialize;

using namespace std;

// subcycling (G 7 1): a grid that is the coarser level of a patch keeps its state at the start of
// its step (u, v, w, phi, press, eddyv) and gets the values at the time of a stage of the finer level
// (linear in time between the start and the end of its step)
struct cfd_amr_tlev
{
    cfd_amr_tlev(lexer*, int n);
    ~cfd_amr_tlev();
    int n;
    field *o[8];        // start of the step: u, v, w, phi, press, eddyv (and kin, eps with turbulence)
    field *t[8];        // at the time of the finer stage
};

//  REEF3D::CFD with mesh refinement (A 10 6, G 1 > 0, static boxes G 15).
//
//  The patches of reefamr3d (3D boxes, refined by 2 in every active direction) each carry an fdm
//  and the objects of a CFD step (convection, diffusion, level set, pressure, momentum_rk).  All
//  levels take the time step of the finest one (G 7 0) and go through the SSP-RK stages together
//  (N 40 2, 3):
//
//    stage s, level by level from level 0 up:
//      - the cells around the patches filled: velocities, level set, pressure, eddy viscosity
//        (from a patch of the same level, or prolonged from the next coarser level), density and
//        viscosity from the level set
//      - level-set transport; the cells around the patch get the level set of the coarser level
//        at the end of its stage, then relaxation and reinitialisation on the patch (those cells
//        fixed)
//      - momentum: convection, sources, diffusion (the cells around the patch are Dirichlet values
//        of the coarser level's stage velocity)
//    the composite projection of all levels (cfd_amr_press.cpp): the faces between a patch and the
//    coarser level take the coarse face velocity (prolonged, conservative), the pressure
//    correction is solved on the leaf cells with the fluxes through those faces matched, so the
//    velocities of all leaf cells are divergence free together
//    the end of the stage on every grid (level set, density, viscosity), then the restriction of
//    velocities, level set, pressure onto the cells and faces under a patch
//
//  The step size is the minimum of all levels.  Scope of this version: two-phase or single-phase
//  flow (F 30, N 40 2/3), RANS turbulence, solids and topography, floating bodies (X 10 1); no
//  sediment, porous media, heat or concentration; wave generation and absorption on level 0 (the
//  patches keep out of the zones).
//
//  G 7 1, subcycling in time (cfd_amr_sub.cpp): every level its own step, level l+1 two steps per
//  step of level l; the cells around a patch linear in time between the start and the end of the
//  step of the coarser level; the projection of a level alone, with the pressure of the coarser
//  level at the patch edges (pr_level), and at the end of the level-0 step the synchronisation
//  projection over all levels (pr_sync_project).  F 46 2/3: the volume correction of the level set
//  over the leaf cells of all grids at the end of the step.

struct cfd_amr_patch : public r3patch
{
    ~cfd_amr_patch();

    fdm *a = nullptr;

    convection *pconvec = nullptr, *pfsfdisc = nullptr;
    diffusion *pdiff = nullptr;
    solver *psolv_in = nullptr, *psolv = nullptr, *ppoissonsolv = nullptr;
    pjm_corr *ppress = nullptr;
    poisson *ppois = nullptr;
    turbulence *pturb = nullptr;
    convection *pturbdisc = nullptr;
    diffusion *pturbdiff = nullptr;
    heat *pheat = nullptr;
    concentration *pconc = nullptr;
    reini *preini = nullptr;
    patchBC_interface *pBC = nullptr;
    ioflow *pflow = nullptr;
    vrans *pvrans = nullptr;
    momentum_rk *pmom = nullptr;
    cfd_amr_6dof *p6dof = nullptr;  // the floating bodies on the patch (X 10 1)

    vector<signed char> mkind;      // per cell of the array: -1 interior, else the fill kind of the cell
    vector<double> edge[3];         // per face of u, v, w: the Dirichlet value of the implicit diffusion at the high patch edges (NaN: none)
    vector<int> row;                // per cell: row of the Poisson matrix (LOOP order), -1

    // composite pressure
    vector<field4*> kv;
    reefmg_core *mg = nullptr;

    // subcycling (G 7 1)
    cfd_amr_tlev *ts = nullptr;
    int comp = -1;      // connected component of the patches of its level (level projection)
    double dt_l = 0.0;  // its step
};

// a face between a patch cell and a cell of the coarser level
struct cfd_amr_cf
{
    int pid;            // patch
    int d, dir;         // direction 0,1,2 and side (+1: the coarse cell at +d)
    int qf, nf;         // fine cell: index in the patch arrays, row
    int qm;             // the cell around the patch across the face
    int gc, qc, nc;     // coarse cell (leaf): grid, index, row
    int gp, qp;         // covered parent of the fine cell: grid, index
    int gface;          // grid of the coarse face value (the grid that computes it)
    int cf[3];          // coarse face: cell index (local on gface) that holds it in field d
    int ff[3];          // fine face: cell index in the patch (fi or fi - e_d)
    int o[3];           // offset of the fine cell in its parent
    int fi[3], ci[3], pi[3];   // fine cell, coarse cell, parent (local indices)
    double dist;        // distance of the cell centres along d
    double area;        // fine face area / coarse face area
    double dxnf, dxnc;  // cell widths along d
    double ustar;       // predictor of the stage
    double kf, kc;      // flux coefficients of the fine and the coarse row
    double rof;         // face density of the stage (pr_weights)
    int ns;             // the coarse value at the position of the fine cell: sum of sw x[sg][sq]
    int sg[5], sq[5];   // (the coarse cell and its leaf neighbours in the transverse directions,
    double sw[5];       // weights of the stage: pr_weights)
    int sc[5][3];       // (their local cell indices)
    double ufin;        // G 7 1 level projection: the corrected face velocity
    double hf, hs;      // hydrostatic pressure of the fine cell and at its position on the coarse side
    bool hok;           // (computed)
    int tq[3][2];       // transverse neighbours of the coarse cell (index, -1: not a leaf cell)
    double tdel[3];     // offset of the fine cell from the coarse centre (coarse cells), 0: none
};

class cfd_amr : public reefamr3d, public amr_grid_lookup
{
public:
    cfd_amr(lexer*, fdm*, ghostcell*, momentum*, pressure*, poisson*, ioflow*, vrans*, fsi*, initialize*, turbulence*, sixdof*);
    virtual ~cfd_amr();

    // false: the run is outside the scope of the refinement (a message has been printed)
    static bool scope(lexer*, ghostcell*);

    void ini(lexer*, fdm*, ghostcell*);

    // the finest patch whose interior holds the point two cells inside its edges (the loads of
    // the floating bodies, X 60 1)
    bool finest_grid(double, double, double, lexer*&, fdm*&) const override;
    void step(lexer*, fdm*, ghostcell*, vrans*, sixdof*);
    void timestep(lexer*, fdm*, ghostcell*);
    void print(lexer*, bool after=false);

    // the composite projection (cfd_amr_press.cpp), public for the vector space of the solver
    void pr_apply(int kx, int ky);
    void pr_prec(int kr, int kz);
    void pr_patch_solve(int id, int kt, int kz, int kc);
    void pr_residual(int kr, int kz, int kt, int lmin);
    double pr_dot(int ka, int kb);
    void pr_start();
    void pr_p(double beta, double om);
    void pr_s(double alp);
    void pr_x(double alp, double om);

    // subcycling (cfd_amr_sub.cpp)
    void step_sub(lexer*, fdm*, ghostcell*, vrans*, sixdof*);

    double tm[8] = {0,0,0,0,0,0,0,0};
    double tp[4] = {0,0,0,0};

protected:
    r3patch* patch_new() override;
    void patch_objects(r3patch*) override;

private:
    cfd_amr_patch* CP(int id) { return static_cast<cfd_amr_patch*>(P[id]); }
    fdm* gfd(int g) { return g<0 ? a0 : CP(g)->a; }
    momentum_rk* gmom(int g) { return g<0 ? mom0 : CP(g)->pmom; }

    // fills: a field of every grid, its kind (cell centred, or the faces of a direction) and the
    // prolongation (0: linear, 1: linear with the MC limiter, 2: injection)
    enum { CELL=-1 };
    struct fspec
    {
        int type;       // CELL or direction 0, 1, 2
        int mode;
        std::function<field&(int)> f;      // source (and destination)
        std::function<field&(int)> fd;     // destination, if not f
    };
    void fill(int l, vector<fspec> &fs, bool siblings_only=false);   // siblings_only: copies from the patches of the same level (kind 0, 2) alone
    void fill_all(vector<fspec> &fs);
    double prolong(const fspec &F, int g, const int *s, const int *o);
    double cell_slope(field &f, int g, const int *s, int d, int mode);
    void margin_rovisc(int id);
    // the boundary values of the state of patch id (before it is the parent of a new patch)
    void bc_state(int id);

    // restriction onto the cells and faces under the patches of level l (from level l)
    void restrict_level(int l, const vector<std::function<field&(int)>> &cells, const vector<std::function<field&(int)>> &faces);
    void restrict_vel(int s);
    void restrict_scalars(int s);

    // stage parts
    void fill_stage_in(int l, int s);
    void fill_phi_out(int l, int s);
    void fill_velout(int l, int s);
    void fill_vel_after(int s);
    void edge_faces(int l, int s);

    // the hierarchy after setup or a regrid: checks, cell kinds, composite projection, subcycling state
    void hierarchy(bool first);

    // adaptive refinement (cfd_amr_regrid.cpp): flags G 24-26, regrid every G 2 steps
    int adapt = 0;
    long nregrid = 0;
    vector<vector<int>> tage;       // [l][tile] regrids since the last flag (hysteresis G 5)
    void flags(vector<vector<unsigned char>> &seed, bool init);
    bool flagged(int g, int i, int j, int k, int lt, bool init);
    void regrid(lexer*, ghostcell*);
    void patch_fields_ini(int id);
    void free_patch(cfd_amr_patch*);

    // subcycling (cfd_amr_sub.cpp)
    cfd_amr_tlev* gts(int g) { return g<0 ? ts0 : CP(g)->ts; }
    field& tsf(int g, int m) { return *gts(g)->t[m]; }
    field& stf(int g, int m);
    void advance(int l, double dtc, double tc);
    void patch_step(int l, double dt, double th0, double t);
    void save_old(int l);
    void interp(int l, double th);
    void fill_level(int l, double th);
    void restrict_sub(int l);
    double tin[3], tout[3];         // stage input and output times in units of the step
    vector<double> th_end;          // per level: end of its last step in the step of the coarser level

    // composite volume of the level set (F 46 2/3): leaf cells of all grids
    void volume(double &V, double &dV);
    void volume_correct(lexer*, ghostcell*);
    void rovisc_all();
    double vol_step = 0.0, vol_net = 0.0;
    bool vol_ini = false;

    // composite projection (cfd_amr_press.cpp)
    void pr_setup();
    void pr_stencils();
    void pr_weights();
    void pr_hydro(int lsub=-1);         // lsub >= 0 (G 7 1, the projection of level lsub): the coarser grids at the time of the stage
    double pr_xs(const cfd_amr_cf&, int k);
    double pr_ps(const cfd_amr_cf&);
    void pr_matrix(int s);
    void pr_predict(int s);
    void pr_correct(int s);
    void pr_project(lexer*, ghostcell*, int s);
    void pr_sync(int k);
    void pr_level(lexer*, ghostcell*, int l, int s);
    void pr_sync_project(lexer*, ghostcell*);   // G 7 1: the synchronisation projection at the end of the level-0 step
    void pr_components();
    const vector<int>& UQ(int g) { return prl<0 ? lq[g+1] : aq[g+1]; }
    const vector<int>& UR(int g) { return prl<0 ? lr[g+1] : ar[g+1]; }
    double* pvec(int g, int k);
    field& pfield(int g, int k);

    lexer *p0;
    fdm *a0;
    momentum_rk *mom0;
    pjm_corr *press0;
    poisson *pois0;
    ioflow *pflow0;
    vrans *pvrans0;
    fsi *pfsi0;
    initialize *pini0;

    int gcval_phi;
    int sub = 0;                        // G 7 1
    turbulence *pturb0 = nullptr;
    int turb = 0;                       // RANS on the patches (T 10 1, 2)
    int solids = 0;                     // solids or topography (direct forcing) on the patches
    sixdof_cfd *b6 = nullptr;           // the floating bodies (X 10 1): moved on level 0 (G 7 0), by the finest level (G 7 1)
    int bodies = 0;
    // G 7 1: the finest level advances the bodies (amr_stage_cfd, the loads after each of its
    // stages on the finest grid); the coarser levels step with a predicted copy (amr_save,
    // amr_restore); bfin counts the finest steps within the level-0 step
    void body_stage(int l, int s, double t, double dt);
    void body_loads(int s, double t, bool fin);
    void body_end();
    int bfin = 0;
    int patch_forcing() const { return bodies ? 1 : (solids ? 2 : 0); }
    int nts = 6;                        // fields of the subcycling state (8 with turbulence)
    field& kfield(int g, int m);        // kin (m 0), eps (m 1) of the turbulence model of grid g
    // the state of a grid transferred to new patches (ini, regrid): 0 phi, 1 press, 2 eddyv, 3-5 u,
    // v, w; with turbulence 6 kin, 7 eps, 10 walld; with solids 8 solid, 9 topo
    vector<int> state_ids();
    static field& sfield(fdm *a, turbulence *t, int m);
    vector<fspec> state_specs();
    void turb_step(int l);              // the turbulence model on the patches of level l (cells around filled)
    cfd_amr_tlev *ts0 = nullptr;

    // composite pressure
    vector<cfd_amr_cf> cf;
    vector<int> row0;                   // level-0 rows
    vector<vector<int>> lq, lr, cq;     // [g+1] leaf cells: index, row; covered cells
    vector<vector<int>> aq, ar;         // [g+1] all cells with a row: index, row (level projection)
    int prl = -1;                       // the projection: -1 composite, l level l alone (G 7 1)
    vector<int> prg;                    // its grids
    int ncomp = 0;                      // connected components of the patches of a level
    vector<field4*> kv0;
    reefmg_core *mg0 = nullptr;
    vector<long> l0a, l0f;              // level-0 V-cycle: active rows of the multigrid, cells
    int pr_it_last = 0;
    int pr_it_sync = 0;                 // G 7 1: iterations of the synchronisation projection
    vector<vector<double>> sync_x;      // G 7 1: the correction of the last synchronisation (leaf cells per grid), its start
    bool sync_x_ok = false;             // (false after a regrid)
    double t_sync = 0.0;                // G 7 1: time of the synchronisation projections
    double sync_bref = 0.0;             // G 7 1: the largest right-hand side norm of the level projections of the level-0 step
    bool press_step = false;            // G 7 1: false until the end of the first level-0 step (the initial pressure is not that of a step)
    double pr_res_last = 0.0;
    long pr_it_total = 0, pr_solves = 0;
    bool pr_ready = false;
};

#endif
