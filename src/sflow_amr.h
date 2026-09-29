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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef SFLOW_AMR_H_
#define SFLOW_AMR_H_

#include"increment.h"
#include<vector>
#include<fstream>
#include<unordered_map>
#include<cstdint>

class lexer;
class fdm2D;
class ghostcell;
class slice;
class ioflow;
class sixdof;
class patchBC_interface;
class sflow_HLL;
class sflow_signal_speed;
class sflow_reconstruct;
class sflow_diffusion;
class sflow_pressure;
class sflow_fsf;
class sflow_forcing;
class sflow_momentum_RK3;
class sflow_pressure_nh;
class reefmg_core;
class reefmg2D;
class vec2D;
class sixdof_sflow;
class sflow_amr_ship;
class solver2D;

using namespace std;

//  Patch-based mesh refinement for REEF3D::SFLOW (hydrostatic HLL, A 220 0).
//
//  Level 0 is the native SFLOW grid.  Each refined patch is a small SFLOW domain
//  of its own: its own lexer (control keys copied, patch geometry, solid flags and
//  wall lists built from flagslice4 as for level 0), its own fdm2D and its own
//  instances of the SFLOW classes, so the kernels run unchanged on every grid.
//
//   - a patch computes EXT cells beyond its box on every side (redundant work,
//     results discarded): the wall, wet-dry and limiter steps inside a stage that
//     need the new values of neighbour cells then see the same data as on one
//     large grid, so the patch layout does not change the solution
//   - all other cells of the patch arrays are filled before each stage: from a
//     patch of the same level, or prolonged from the next coarser level (well
//     balanced: the surface is interpolated, the depth follows the fine bed);
//     cells owned by another rank are computed there and sent
//   - after every RK3 stage a patch is restricted into the coarser level
//     (conservative averages of WL, UH, VH)
//   - a coarse face on a coarse-fine interface takes the mean of the two fine
//     face fluxes and face depths (sflow_HLL calls hll_hook between flux_bc and
//     the divergence); across a partition edge the fine values are sent
//   - one global time step, the stages of all grids in lockstep
//   - every A 271 steps the patches are rebuilt from refinement flags (A 273
//     surface jump, A 274 shoreline, A 276 boxes) with a buffer of A 272 cells.
//     Flags mark tiles of A 275 cells on the global index space of each level;
//     the tile maps are global, so the refined region does not depend on the
//     domain decomposition.  Marked tiles are merged into rectangles and cut at
//     the partition edges.
//   - moving body (X 10 2/3): the body stays on level 0 (sixdof_sflow); the patches
//     evaluate it on their own cells (sflow_amr_ship.cpp): the level set is interpolated
//     from level 0, the draft is ray-cast from the hull triangles at the patch cell centres,
//     and the pressure (X 10 3) or the direct forcing (X 10 2) is applied in the patch
//     kernels.  A 278 refines a margin around the hull, A 279 a wake wedge behind the bow,
//     both moving with the body.
//
//  Index conventions: every level has a global cell index space, level l+1
//  refines level l by 2.  A patch covers the global box [I0,I1]x[J0,J1] of its
//  level; its lexer index is i = I-I0+EXT.

struct sflow_amr_fill
{
    int di,dj;          // destination cell (patch lexer index)
    int kind;           // 0: copy from grid g, 1: prolong from grid g (coarser level), 2: remote
    int g;              // grid id (-1: level 0)
    int si,sj;          // source cell on grid g
    int ox,oy;          // prolongation: quadrant of the fine cell (-1/+1)
    int slot;           // remote: position in the receive buffer of the level
    double depth;       // prolongation: still water depth of the fine cell
};

struct sflow_amr_patch
{
    int lev;
    int I0,I1,J0,J1;                // global box of the patch on its level
    int nx,ny;

    lexer *pp;
    fdm2D *b;
    sflow_HLL *phll;
    sflow_signal_speed *pss;
    sflow_reconstruct *precon;
    sflow_diffusion *pdiff;
    sflow_pressure *ppress;
    sflow_fsf *pfsf;
    sflow_forcing *psfdf;
    sflow_momentum_RK3 *pmom;

    // boundary face values recorded by hll_hook: rec[ipol][side][fine index]
    // side 0: low x, 1: high x, 2: low y, 3: high y; ipol 0 holds dfx/dfy
    vector<double> rec[5][4];

    // cells filled before each stage
    vector<sflow_amr_fill> fill;

    // restriction: coarse target (grid id, i, j) for every 2x2 block (-2: wall)
    vector<int> rgrid, ric, rjc;

    bool fresh;                     // created by the current regrid

    // non-hydrostatic pressure (A 220 1), solved on all grids together
    sflow_pressure_nh *pnh = nullptr;
    vector<slice*> nv;              // Krylov vectors
    vector<signed char> act;        // -2 no row, -1 covered by a finer patch, 0 q = 0 row, 1 active
    vector<int> row;                // matrix row of a cell (SLICELOOP4 order), -1 none
    reefmg_core *mg = nullptr;      // patch-local multigrid of the preconditioner
    vector<double> qrec[4];         // boundary face gradients of q

    // moving body on the patch (X 10 2/3)
    sflow_amr_ship *pship = nullptr;

    // line solver of the Boussinesq u_a inversion (A 220 4)
    solver2D *psolv = nullptr;
};

// a coarse face overridden with fine values
struct sflow_amr_match
{
    int dir;          // 0: x face, 1: y face
    int fi,fj;        // local face index on the target grid
    int child;        // patch that recorded the fine faces (-1: remote)
    int side,r;       // recorded side and first fine index
    double val[6];    // remote: averaged fine values (0: face depth, 1-4: fluxes, 5: gradient of q)
};

// point-to-point exchange with a fixed set of peers
struct sflow_amr_xplan
{
    vector<int> speer, rpeer;            // ranks
    vector<vector<int>> sitem;           // per send peer: item indices (meaning depends on the plan)
    vector<int> rcount;                  // per receive peer: number of items
    vector<vector<double>> sbuf, rbuf;
};

class sflow_amr : public increment
{
public:
    sflow_amr(lexer*, fdm2D*, ghostcell*, patchBC_interface*, sixdof*, sflow_HLL*, sflow_momentum_RK3*);
    virtual ~sflow_amr();

    void ini(lexer*, fdm2D*, ghostcell*);
    void step_begin(lexer*, fdm2D*, ghostcell*);
    void stage_begin(lexer*, fdm2D*, ghostcell*, int);
    void stage_end(lexer*, fdm2D*, ghostcell*, int);
    void step_end(lexer*, fdm2D*, ghostcell*);
    void timestep(lexer*, fdm2D*, ghostcell*);
    void print(lexer*, fdm2D*, ghostcell*);

    // called by sflow_HLL between flux_bc and the divergence (ipol 1-4)
    void hll_hook(lexer*, fdm2D*, int, int);

    int patches_total;

    // composite non-hydrostatic pressure, called by the level-0 sflow_pjm_lin after its assembly
    bool nh_patches() const { return maxlev>0 && patches_total>0 && nh==1; }
    void nh_solve(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&, double);

private:
    // composite non-hydrostatic solve (sflow_amr_nh.cpp)
    int nh;
    struct nhg { lexer *q; fdm2D *b; vector<signed char> *act; vector<int> *row; vector<slice*> *v; };
    nhg nh_grid(int);
    slice& nh_vec(int, int);        // grid, vector (-1: press)
    void nh_prepare(ghostcell*);
    void nh_restrict_vec(int);
    void nh_sync(int);
    void nh_qfill(int, int);
    double nh_eval(const sflow_amr_fill&, int);
    void nh_apply(int, int);
    double nh_dot(int, int);
    void nh_prec(int, int);
    vector<signed char> nh0_act;
    vector<int> nh0_row;
    vector<slice*> nh0_v;
    reefmg2D *nhmg0 = nullptr;
    vec2D *nhr0 = nullptr;
    bool nh_rebuild0;
    long nh_it_total, nh_solves;
    int nh_it_last;
    vector<int> nlevg;              // patches per level, all ranks

    // moving ship (X 10 2/3): body fields on the patches, refinement zone A 278/A 279
    int shipmode;                   // X 10 of the body, 0: none
    sixdof_sflow *ship6;
    struct shipzone { double cx,cy,ex,ey,smin,smax,nmin,nmax,sfront,wake,bx0,bx1,by0,by1; };
    vector<shipzone> zones;
    vector<double> ship_x0, ship_y0;    // initial position of every body
    vector<vector<int>> shiptri;    // per body: triangles reaching below the still water level
    void ship_setup(lexer*);
    bool ship_zone(double, double);
    void ship_fields(sflow_amr_patch&, bool);
    void ship_patches(bool);
    double fs0_at(sflow_amr_patch&, int, int);

    // grid handles: id -1 is level 0, otherwise an index into P
    struct gh
    {
        lexer *q;
        fdm2D *b;
        sflow_momentum_RK3 *m;
        int oi,oj;    // local index = global index - oi
    };
    gh grid(int);
    int patch_at(int, int, int);            // level, I, J: interior patch on this rank, -1: level 0 (l==0), -2: none
    int owner(int, int);                    // rank of the level-0 cell (I,J), -1: outside the domain
    int flag0(int, int);                    // level-0 flagslice4 of the global cell (I,J), from the rank box + halo
    void build_flags(lexer*);
    vector<int> fl0;
    static const int FH = 8;

    // regridding
    void regrid(lexer*, fdm2D*, ghostcell*, bool);
    void tag_level(int, vector<unsigned char>&);
    void global_or(vector<unsigned char>&);
    void build_tiles();
    sflow_amr_patch* make_patch(lexer*, ghostcell*, int, int, int, int, int);
    void free_patch(sflow_amr_patch*);
    void build_lexer(lexer*, sflow_amr_patch&);
    void build_bc(ghostcell*, sflow_amr_patch&);
    double bed_at(int, int, int);
    unordered_map<uint64_t,double> bedmemo;     // the bed is static (S 10 0)
    void build_plans(ghostcell*);
    void ini_patch_state(lexer*, ghostcell*, sflow_amr_patch&, vector<sflow_amr_patch*>&);

    // coupling
    void cache_stage(int);
    void fill_level(ghostcell*, int, int);
    void eval_fill(const sflow_amr_fill&, double*);
    void prolong(int, int, int, int, int, double, double*);
    void apply_bc(ghostcell*, sflow_amr_patch&, int);
    void restrict_patch(lexer*, sflow_amr_patch&, int);
    void exchange_level0(lexer*, fdm2D*, ghostcell*, int);
    void exchange_fluxes(int);

    // MPI helpers
    void xsetup(vector<vector<int>>&, vector<vector<int>>&, int);
    void xrun(sflow_amr_xplan&, int, int);

    double mass(lexer*, fdm2D*, ghostcell*);
    void write_vtr(lexer*, sflow_amr_patch&, int);
    void write_vtr0(lexer*, fdm2D*);
    void gauges(lexer*, fdm2D*, ghostcell*);
    ofstream gaugeout;

    vector<sflow_amr_patch*> P;
    vector<vector<int>> lev;        // patch ids per level (index 0 unused)
    int maxlev, nest, tile, nbuf, regrid_int;
    double tol_eta;
    int shore;

    // rank boxes of level 0 (global cell indices), all ranks
    int O0i,O0j,NX0,NY0,GNX,GNY;
    vector<int> rbx0,rbx1,rby0,rby1;
    int last_owner;
    int bxlo(int l) const { return O0i<<l; }
    int bxhi(int l) const { return ((O0i+NX0)<<l)-1; }
    int bylo(int l) const { return O0j<<l; }
    int byhi(int l) const { return ((O0j+NY0)<<l)-1; }

    // global tile maps per level (1 = refined), and the local tile -> patch map
    vector<vector<unsigned char>> gtile;
    vector<vector<unsigned char>> tage;     // regrids since a tile was last flagged (A 280)
    int keep;
    vector<int> gtnx,gtny;
    vector<vector<int>> tmap;
    vector<int> tti0,ttj0,tnx,tny;

    // flux matching entries per target grid (index id+1): local and remote
    vector<vector<sflow_amr_match>> match;
    vector<vector<sflow_amr_match>> rmatch;

    // ghost service per level: items are sflow_amr_fill descriptors of the cells served
    vector<sflow_amr_xplan> gplan;
    vector<vector<sflow_amr_fill>> gserve;      // [l]: served cells (kind 0/1)
    vector<vector<sflow_amr_fill*>> grecv;      // [l]: receive slot -> fill entry

    // fine face records sent across partition edges, per level
    vector<sflow_amr_xplan> fplan;
    vector<vector<int>> fsend;                  // [l]: (patch, side, r) triplets
    vector<vector<int>> frecv;                  // [l]: target grid id+1 and index into rmatch, -1: skip

    // cached stage input arrays per grid (index id+1)
    vector<slice*> sWL,sUH,sVH,sWH;

    // no refinement near in- and outflow boundaries and in A 277 boxes: global level-0 cells
    vector<unsigned char> forbid0;

    lexer *p0;
    fdm2D *b0;
    ghostcell *pgc0;
    sflow_HLL *phll0;
    sflow_momentum_RK3 *pmom0;
    patchBC_interface *pBC;
    sixdof *p6dof;
    ioflow *pflow_void;

    double m0, printtime_amr;
    double tm[9];
    int printcount_amr, regrids;
    long cells_total;
    ofstream logout;
    const double eps;
    int EXT;                        // extra computed cells on each side of a patch (2, 3 with B 60)
    static const int NVMAX = 14;
    int NV;                         // values per filled cell: WL UH VH eta U V wet deep WH W (+ Boussinesq: UA VA MX MY)
    int bous;                       // Boussinesq (A 220 4): UH,VH hold V, u_a and M are filled as well
};

#endif
