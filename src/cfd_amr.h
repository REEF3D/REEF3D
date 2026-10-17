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
#include<map>
#include<tuple>

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
class wave_interface;
class bc_zone_set;
class sixdof_cfd;
class cfd_amr_6dof;
class momentum;
class momentum_rk;
class reefmg_core;
class initialize;
class field1;
class field2;
class field3;
class weno_flux_nug;
class fluid_update_fsf_heat;

using namespace std;

// subcycling (G 7 1): a grid that is the coarser level of a patch keeps its state at the start of
// its step (the fields of the registry with AF_TIME: u, v, w, phi, press, eddyv; kin, eps) and gets
// the values at the time of a stage of the finer level (linear in time between the start and the
// end of its step)
struct cfd_amr_tlev
{
    cfd_amr_tlev(lexer*, const vector<int> &types);     // per field: CELL (-1) or the face direction
    ~cfd_amr_tlev();
    int n;
    vector<field*> o;   // start of the step
    vector<field*> t;   // at the time of the finer stage
};

// The registry of the patch state (cfd_amr_fields.cpp): every field of a grid that the refinement
// moves between the grids, once, with what it takes part in.  The fills (cells around the patches),
// the restriction, the time levels of the subcycling, the transfer at a regrid and the boundary
// values of the grids are built from it; new physics (heat, concentration) adds its fields there.
enum : unsigned
{
    AF_STATE   = 1u<<0,     // the state of new patches (ini, regrid: from the parent or an old patch) and its fill
    AF_TIME    = 1u<<1,     // the time levels of the subcycling (G 7 1): kept at the start of a step, in time for the finer level
    AF_TMID    = 1u<<2,     //   its value belongs to the middle of the step (the pressure; held at the free surface from a patch)
    AF_STAGEIN = 1u<<3,     // filled around the patches at the start of a stage
    AF_TURBIN  = 1u<<4,     // filled around the patches before the turbulence step (G 7 0)
    AF_RSTEP   = 1u<<5,     // restricted at the end of a step of the finer level (G 7 1)
    AF_RSTAGE  = 1u<<6,     // restricted at the end of a stage (G 7 0)
    AF_BCSTATE = 1u<<7,     // boundary values of a patch before it is a parent (bc_state)
    AF_BCINI   = 1u<<8,     // boundary values of the initial patches
    AF_BCSTEP  = 1u<<9,     // boundary values of level 0 after the restriction of a step (G 7 1)
    AF_BCSTAGE = 1u<<10,    // boundary values of level 0 after the restriction of a stage (G 7 0)
    AF_CONS    = 1u<<11     // conserved quantity: its fluxes at the patch edges matched (flux registers, cfd_amr_scalar.cpp)
};
struct cfd_amr_patch;
struct cfd_amr_field
{
    const char *name;
    int type;           // CELL (-1) or the face direction 0, 1, 2
    int prolong;        // 0 linear, 1 linear with the MC limiter, 2 injection
    unsigned roles;     // AF_*
    int gc, gcval;      // its boundary values: ghostcell start1-4 (1-4), start4a (5), none (0); the code
    std::function<field&(fdm*,turbulence*)> get;      // the field of a grid (its fdm and turbulence model)
    int tl = -1;        // AF_TIME: its index in the time levels
    // a field of an object of the grid (heat T, concentration C): of a patch, of level 0
    std::function<field&(cfd_amr_patch*)> pget;
    std::function<field&()> get0;
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
//  flow (F 30, N 40 2/3), RANS turbulence, solids and topography, floating bodies (X 10 1), the
//  porous structures B 200 1 (vrans_f on every patch, the zones at its cells; CPOR n at the faces of
//  the patch edges); heat or concentration (H 10 3, C 10 3: transported by the refinement on all
//  grids in flux form, the patch edges with flux registers, cfd_amr_scalar.cpp); no sediment or
//  vegetation; the convection schemes D 10
//  0-8, 10-16, 60 and F 35 1-8, 10-16, the diffusion D 20 0-3; the relaxation zones of the wave
//  generation and the beach (B 98 2, B 99 1/2) with a wave theory on the patches as on level 0
//  (cfd_amr_relax.cpp, the eddy viscosity as well), the other wave generation and absorption on
//  level 0; the tiles at an in- or outflow boundary stay on level 0 (reefamr3d::mark_tiles); the
//  pressure reference B 30 1, 2, 4 over all grids at the end of the step (cfd_amr::reference), B 30
//  3 (the atmosphere) on all grids.  The probes and wave gauges of the printer (P 51, P 61, P 64)
//  read the finest grid at their position (finest_grid).
//
//  G 40 1, the placement of the patches for the load of the ranks (reefamr3d, cfd_amr_place.cpp):
//  a level-1 piece with the finer patches inside it (a unit) on any rank, its coarser level a
//  shadow of the level-0 cells under it and around it (filled from level 0, giving back what the
//  patches write onto level 0); the outputs and the loads of the bodies on the rank of the finest
//  grid at the point (amr_cfd_mine).  G 40 2: the test of it (small units, placed anew).
//
//  Periodic boundaries in x and y (C 21, C 22 of DIVEMesh; on one rank or across ranks): a level is
//  refined alike on both sides of a periodic edge (reefamr3d::mark_tiles), so no patch edge lies on
//  it; the margin of a patch across it comes from the patches at the other end (the fill plans
//  wrap the index), the composite solve has the level-0 multigrid periodic, the patches next to
//  each other across the edge in one joint multigrid (shifted by the length; periodic when it goes
//  around), a patch over the whole length periodic over itself, the overlap of the multigrids at
//  the sides to patches of another rank.  Not with G 40 (it falls back to G 40 0) nor in z.
//
//  G 7 1, subcycling in time (cfd_amr_sub.cpp): every level its own step, level l+1 two steps per
//  step of level l; the cells around a patch linear in time between the start and the end of the
//  step of the coarser level; the projection of a level alone, with the pressure of the coarser
//  level at the patch edges (pr_level), and at the end of the level-0 step the synchronisation
//  projection over all levels (pr_sync_project).  F 46 2/3: the volume correction of the level set
//  over the leaf cells of all grids at the end of the step.

// heat, concentration (cfd_amr_scalar.cpp): the work fields of a grid
struct cfd_amr_scgrid
{
    field *s0 = nullptr, *kap = nullptr, *L = nullptr;
    field *fx = nullptr, *fy = nullptr, *fz = nullptr;
    weno_flux_nug *weno = nullptr;
};

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
    cfd_amr_scgrid *sc = nullptr;   // heat, concentration: the work fields

    vector<signed char> mkind;      // per cell of the array: -1 interior, else the fill kind of the cell
    vector<double> edge[3];         // per face of u, v, w: the Dirichlet value of the implicit diffusion at the high patch edges (NaN: none)
    vector<int> row;                // per cell: row of the Poisson matrix (LOOP order), -1

    // composite pressure
    vector<field4*> kv;
    reefmg_core *mg = nullptr;
    // G 40: the overlap of the patch multigrid into the patches of the level next to it (cells below
    // and above per direction; a restricted additive Schwarz across the cuts between the units), the
    // overlap its multigrid was made with, the matrix rows of the patch (filled into the overlap of
    // the others)
    int olo[3] = {0,0,0}, ohi[3] = {0,0,0};
    int mlo[3] = {0,0,0}, mhi[3] = {0,0,0};
    vector<field4*> kc7;

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
    bool finest_grid(double, double, double, lexer*&, fdm*&, turbulence*&) const override;
    bool placed() const override { return place>0; }
    int finest_level(double, double, double) const override;
    int finest_id(double, double, double) const;
    void step(lexer*, fdm*, ghostcell*, vrans*, sixdof*);
    void timestep(lexer*, fdm*, ghostcell*);
    void print(lexer*, ghostcell*, bool after=false);
    bool print_lagoon(lexer*, ghostcell*, int num);    // P 18: the blocks into the LAGOON store

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

    // heat (H 10) or concentration (C 10) on all grids (cfd_amr_scalar.cpp): the driver's objects
    // of level 0 (before ini); true: the refinement transports the scalar on all grids (the driver
    // leaves its step out)
    void scalar_objects(heat*, concentration*);
    bool scalars() const { return scal>0; }

    double tm[8] = {0,0,0,0,0,0,0,0};
    double tp[4] = {0,0,0,0};
    double t0l[3] = {0,0,0};        // G 7 0: the level-0 parts of tm[0], tm[1], tm[3]

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
    void fill(int l, vector<fspec> &fs, bool siblings_only=false, int maxdepth=1<<20);   // siblings_only: copies from the patches of the same level (kind 0, 2) alone;
                                                                                        // maxdepth: the cells around the patch up to that distance
    void fill_all(vector<fspec> &fs, int maxdepth=1<<20);
    double prolong(const fspec &F, int g, const int *s, const int *o);
    double prolong_f(field &f, int type, int mode, int g, const int *s, const int *o);
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
    void regrid_old(const vector<int> &fresh, vector<vector<char>> &got);   // G 40: from the old patches of all ranks
    void patch_fields_ini(int id);
    void free_patch(cfd_amr_patch*);

    // the registry of the patch state (cfd_amr_fields.cpp)
    vector<cfd_amr_field> reg;
    void register_fields();
    turbulence* gturb(int g) { return g<0 ? pturb0 : CP(g)->pturb; }
    field& rf(const cfd_amr_field &e, int g)
    {
        if(e.pget)
        return g<0 ? e.get0() : e.pget(CP(g));
        return e.get(gfd(g),gturb(g));
    }
    field& rfp(const cfd_amr_field &e, cfd_amr_patch *c) { return e.pget ? e.pget(c) : e.get(c->a,c->pturb); }
    vector<fspec> specs(unsigned role);                         // fills
    vector<std::function<field&(int)>> rcells(unsigned role);    // restriction: cell fields
    vector<std::function<field&(int)>> rfaces(unsigned role);    // restriction: u, v, w (by direction)
    void bc_grid(int g, unsigned role);                         // boundary values of grid g (level 0: -1)
    vector<int> tl_types();
    int tl_phi = -1, tl_press = -1;     // the time levels of the level set and the pressure

    // subcycling (cfd_amr_sub.cpp)
    cfd_amr_tlev* gts(int g) { return g<0 ? ts0 : CP(g)->ts; }
    field& tsf(int g, int m) { return *gts(g)->t[m]; }
    field& stf(int g, int m);
    int nts() const { return tl_n; }
    int tl_n = 0;
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
    void pr_covered_coef();                              // the coarse couplings under the finer patches from the fine faces
    void pr_halo0(field &f);                             // the halo of a level-0 vector (ranks, periodic edges)
    bool pr_wrap(int id, int d) const;                   // patch id over the whole length of periodic direction d
    // the overlap of the patch multigrids at the sides where a patch of the level is next to it (a
    // restricted additive Schwarz): G 40 (instead of the joint multigrids of the level); periodic (the
    // patches not in a joint multigrid)
    bool pr_ras() const { return place>0 || per[0]>0 || per[1]>0; }
    void pr_level(lexer*, ghostcell*, int l, int s);
    void pr_sync_project(lexer*, ghostcell*);   // G 7 1: the synchronisation projection at the end of the level-0 step
    double psi_band(int g);                     // the half width of the interface band of grid g
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
    field& kfield(int g, int m);        // kin (m 0), eps (m 1) of the turbulence model of grid g
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
    // the joint multigrid of the patches of a level on this rank that touch each other (their
    // bounding box, the cells outside the patches identity rows): the V-cycle of a level in the
    // preconditioners (the level projection of G 7 1, the FAC sweep of the composite projection
    // and the synchronisation) couples the patches next to each other; a patch alone its own
    struct cfd_amr_lmg
    {
        reefmg_core *mg = nullptr;
        vector<int> ids;                        // its patches
        int lo[3] = {0,0,0}, n[3] = {0,0,0};    // box: first cell (index of the level), extent
        vector<int> own;                        // per cell of the box: the patch, -1
        std::map<int,vector<int>> off;          // per patch: the box position of its cell 0 (3); a patch
                                                // joined across a periodic edge lies past the end of the box
        bool wrap[3] = {false,false,false};     // the box over the whole length of a periodic direction:
                                                // the multigrid periodic there
        int o(int id, int d) const { return off.at(id)[d]; }
        // periodic: the overlap across the sides where a member meets a patch of another rank (OV
        // cells of its margin, own -2): the rows and the right-hand side of that patch, from the margin
        // of the member ovm at index ovq; the correction kept in the members
        vector<int> ovm, ovq;
    };
    vector<cfd_amr_lmg> lmg;            // the joint multigrids (all levels)
    vector<int> lmg_of;                 // per patch: its joint multigrid, -1 (its own)
    void pr_level_mg_setup();
    void pr_level_mg_fill(int m);
    void pr_level_mg_fill_comp(int m);
    void pr_level_mg_ovrows(int m, double dshift);   // periodic: the rows of the overlap cells of a joint multigrid
    void pr_level_vcycle(int m, int kr, int kz);
    void pr_level_solve(int l, int kr, int kz);     // the V-cycles of level l: joint or patch-local
    // the fine level of the multigrid of patch id from its rows (diagonal times dshift), the couplings
    // to the cells around it dropped, G 40: over its box with the overlap (pr_overlap, pr_coef_fill)
    void pr_patch_mg(int id, double dshift);
    void pr_overlap();
    void pr_coef_fill(int l);

    // wave generation and the numerical beach on the patches (cfd_amr_relax.cpp): the relaxation
    // zones (B 98 2, B 99 1/2) with a wave theory of level 0, evaluated at the cells of the patch
    // and the time of its step (G 7 1) or of the step (G 7 0)
    wave_interface *pwave0 = nullptr;
    bc_zone_set *zones0 = nullptr;
    int zgen = 0, zbeach = 0;
    void relax_setup(lexer*, ghostcell*);
    bool patch_in_zone(int id);
    void patch_relax_vel(int id, field &u, field &v, field &w);
    void patch_relax_phi(int id, field &f);
    bool patch_relax_turb(int id);      // the eddy viscosity in the generation zone (as iowave::turb_relax)

    // the pressure reference (B 30 1, 2, 4) over all grids, once at the end of the step
    // (cfd_amr_relax.cpp): the patches have none of their own, level 0 none during the step
    int b30 = 0;
    void reference(lexer*, ghostcell*);
    void step_comp(lexer*, fdm*, ghostcell*, vrans*, sixdof*);
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

    // G 40, the placement of the units (cfd_amr_place.cpp): the copy plans between level 0 and the
    // shadows of the level-1 patches (reefamr3d: sh_down, all cells of the shadow arrays)
    r3cplan pu_cell;                    // up: the restriction onto the cells under the level-1 patches (set)
    r3cplan pu_face[3];                 //     the faces under them (set; the - face at the low edge too)
    r3cplan pu_cf[3];                   //     the coarse faces at the edges of the level-1 patches (set)
    r3cplan pu_cc;                      //     the coarse cells at those edges (added: pr_apply, refluxing)
    r3cplan pu_key;                     //     the coarse cell and side of the faces: the folds of the level-0
                                        //     rows (pjm_corr), the couplings of the level-0 multigrid (pr_matrix)
    vector<vector<double>> kval;        // per shadow: the sum of kc of its keys (pu_key order)
    std::map<std::tuple<int,int,int,int,int,int>,int> kidx;    // (shadow, coarse cell, d, dir) -> its key
    void sh_plans();
    void sh_run(r3cplan &X, const vector<std::function<field&(int)>> &f, bool add);
    void sh_fill_down(const vector<fspec> &fs);    // the sources of a fill of level 1 into the shadows
    void sh_down_f(const vector<std::function<field&(int)>> &f) { sh_run(sh_down,f,false); }
    void sh_zero(r3cplan &X, const std::function<field&(int)> &f);  // the sent items of X in f (the shadows)
    void sh_folds();
    void sh_ks();
    bool is_shadow(int g) const { return g>=npat; }

    // heat or concentration (cfd_amr_scalar.cpp): a conservative transport of the scalar on all
    // grids (SSP-RK3 in flux form; the diffusion explicit in the stages, D 20 1, or implicit after
    // them, D 20 2), the fluxes through the faces at the patch edges matched with flux registers
    int scal = 0;                       // 0 none, 1 heat T (H 10), 2 concentration C (C 10)
    int sc_conv = 0;                    // convection H 15 / C 15: 0 none, 1 first-order upwind, 4 WENO (weno_flux_nug)
    int sc_diff = 0;                    // D 20: 0 none, 1 explicit, 2 implicit
    double sc_sig = 1.0;                // the turbulent Prandtl (Schmidt) number of the eddy viscosity
    heat *pheat0 = nullptr;
    concentration *pconc0 = nullptr;
    fluid_update_fsf_heat *smat = nullptr;      // heat: the material tables of the density and viscosity
    cfd_amr_scgrid sc0;                 // level 0 (the patches: cfd_amr_patch::sc)
    cfd_amr_scgrid& SG(int g) { return g<0 ? sc0 : *CP(g)->sc; }
    void sc_make(int g);
    void sc_release(cfd_amr_scgrid&);
    vector<double> fr_c, fr_f;          // flux registers per face of cf: the coarse and the fine flux (time integral x fine face area)
    int tl_s = -1;                      // the time level of the scalar
    field& sfield(int g);               // the scalar of grid g
    double sval(int g, int i, int j, int k);
    double sval(int g, int q);
    double rmix(const lexer*, double H, double s);   // density of the mixture (as fluid_update_fsf_heat / _concentration)
    double vmix(const lexer*, double H, double s);   // viscosity
    void sc_alloc();
    void sc_free();
    void sc_kappa(int g);
    void sc_save(int g);
    void sc_flux(int g);
    void sc_record(int g, double w, bool coarse, int lfine);
    void sc_record_diff(int g, double w, bool coarse, int lfine);
    void sc_update(int g, int s, double dt);
    void sc_bc(int g);
    void sc_diffuse(int g);
    void sc_fill(int l, int mode);      // the cells around the patches of level l: 0 the current scalar, 1 the time level of the coarser grid
    void sc_step_comp();                // G 7 0: all grids
    void sc_step_level(int l, double dt, double th0);  // G 7 1: the grids of level l
    void sc_reflux(int l);              // the flux registers of the patch edges of level l onto level l-1
    void sc_reset(int lfine);
    void sc_record_shadows(double w, bool diff);    // G 40: the coarse fluxes of level 0 at the level-1 edges
    double sc_total(double &smin, double &smax);
};

#endif
