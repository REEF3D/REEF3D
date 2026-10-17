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

#ifndef REEFAMR3D_H_
#define REEFAMR3D_H_

#include"increment.h"
#include<vector>
#include<unordered_map>
#include<mpi.h>

class lexer;
class ghostcell;

using namespace std;

//  REEFAMR3D: patch-based mesh refinement with 3D boxes on a Cartesian grid (the CFD module,
//  cfd_amr).  The other modules (SFLOW, NHFLOW, FNPF) use the column core reefamr.h.
//
//  Level 0 is the native grid of the run.  Level l+1 refines level l by 2 in every active
//  direction (i_dir, j_dir, k_dir; an inactive direction, such as y of a 2D x-z case, keeps one
//  cell).  The refined region of a level is a union of tiles of G 4 cells on the global index
//  space of the level, the same on every rank; tiles are merged into boxes and the boxes are cut
//  at the rank boxes of level 0, so every patch lies on the rank of the level-0 cells below it,
//  and its parent cells are on that rank too (restriction is local).
//
//  A patch has its own lexer with the CFD margin (3 cells) around its box.  The margin cells
//  inside the domain are filled before they are used: from a patch of the same level (on this
//  rank, or sent from the rank that holds it), or prolonged from the next coarser level.  Margin
//  cells outside the domain are the physical boundary of the patch: the patch lexer has the
//  boundary lists of level 0 below it, and ghostcell sets them.
//
//  Placement (G 40 1): the level-1 patches are not cut at the rank boxes; a level-1 piece with the
//  patches of the finer levels inside it (a unit) is placed on a rank for the load of the ranks.
//  The coarser level of a level-1 patch is then a shadow: a copy of the level-0 cells under it and
//  around it (two cells, the array margin three more) on the rank of the unit, filled from the
//  ranks that hold those cells (down) and giving back what the patch writes onto level 0 (up):
//  every access of the patches to level 0 is local, to the shadow.  The finer levels of a unit are
//  cut at its edges, so their parents are in the unit; a unit is cut only where the finer levels
//  continue across the cut (no edge of a finer level on it).
//
//  Static refinement (step 15b): the boxes G 15 xs xe ys ye zs ze are refined to level G 1, and
//  every coarser level covers the next finer one with a nesting band of 3 cells.  Adaptive
//  refinement: the module gives seed tiles per level (setup, rebuild), the boxes G 15 are added on
//  the finest level; rebuild keeps the patches whose box is unchanged and makes the others anew.

struct r3fill
{
    int d[3];           // destination cell (patch index, in the margin)
    int kind;           // 0: copy from grid g (same level), 1: prolong from grid g (next coarser level),
                        // 2: remote (copy from a patch of another rank: slot in the receive buffer), 3: boundary
    int g;              // grid id (-1: level 0)
    int s[3];           // source cell on grid g (kind 1: the parent cell)
    int o[3];           // kind 1: position of the child in the parent (0, 1; 0 in an inactive direction)
    int slot;           // kind 2
    int depth;          // the distance from the patch (1: next to it, across a face, edge or corner)
};

// the level-0 cell data of a box (G 40, the shadows): what a patch lexer takes from level 0 below it
struct r3raw
{
    int lo[3] = {0,0,0}, n[3] = {0,0,0};    // global level-0 cells of the box (with the array margin)
    vector<int> flag4;                      // per cell
    vector<int> fs1, fs2, fs4, wet, wetn, deep;         // per column (i,j)
    vector<double> sbed, tbed, gbed, bed, depth, WL;
    std::unordered_map<long long,int> grp;  // the boundary group of a global cell and side
    vector<double> xn2d, yn2d;              // the 2D nodes from lo, (n0+1)(n1+1)
    size_t c3(int I, int J, int K) const { return ((size_t)(I-lo[0])*n[1] + (J-lo[1]))*n[2] + (K-lo[2]); }
    size_t c2(int I, int J) const { return (size_t)(I-lo[0])*n[1] + (J-lo[1]); }
    size_t n2(int I, int J) const { return (size_t)(I-lo[0])*(n[1]+1) + (J-lo[1]); }
    static long long gkey(int I, int J, int K, int side) { return (((long long)(I+16)*65536 + (J+16))*65536 + (K+16))*8 + side; }
};

struct r3patch
{
    virtual ~r3patch(){ delete raw; }

    int lev = 0;
    int gid = -1;       // index into the global list
    bool shadow = false;    // G 40: a copy of level-0 cells (the coarser level of a level-1 patch), not a patch
    int sh = -1;            // G 40: the shadow of the unit of the patch (its coarser level for level 1)
    r3raw *raw = nullptr;   // a shadow: the level-0 data of its array
    int lo[3], hi[3];   // box on the global index space of the level (inclusive)
    int n[3];           // cells

    lexer *pp = nullptr;

    vector<r3fill> fill;

    // the coarse cells of the patch (index of level lev-1) inside the local grids of that level
    struct pblock { int g; int lo[3], hi[3]; };
    vector<pblock> par;
};

struct r3gpatch
{
    int lev, lo[3], hi[3], rank, lid;
    int unit = -1;      // G 40: the level-1 patch (global id) of its unit
};

// items served to other ranks (copy from local grid g at cell s)
struct r3serve
{
    int g;
    int s[3];
};

struct r3xplan
{
    vector<int> scnt, sdsp, rcnt, rdsp;     // per rank, in items
    vector<r3serve> serve;
    int nsend = 0, nrecv = 0;
    vector<double> sbuf, rbuf;
};

// G 40: a copy plan between the level-0 cells (on the ranks that hold them) and grids of this rank
// (the shadows): down the level-0 values into the grids, up the values of the grids onto level 0
// (set or added).  Every item: a grid and cell index here, a level-0 cell index there; a tag per
// item (the receiving side keeps it).
struct r3cplan
{
    vector<int> sp, rp;                 // peers: sent to, received from
    vector<vector<int>> sg, sq;         // per send peer: grid, cell index of the items
    vector<vector<int>> rg, rq, rt;     // per receive peer: grid, cell index, tag
    vector<vector<double>> sb, rb;
    bool empty() const { return sp.empty() && rp.empty(); }
};

class reefamr3d : public increment
{
public:
    reefamr3d(lexer*, ghostcell*);
    virtual ~reefamr3d();

    int maxlev;
    lexer *p0;
    ghostcell *pgc0;
    int myrank, nranks;

    int rr[3];                      // refinement ratio per direction (2, 1 inactive)
    int gn0[3];                     // global cells of level 0
    int org[3], kn0[3];             // level-0 box of this rank
    vector<int> rbox;               // level-0 boxes of all ranks: origin, cells (6 per rank)
    int tile[3];                    // tile size in cells of the refined level
    int per[3] = {0,0,0};           // periodic directions of level 0 (1 on one rank, 2 across ranks; no z)
    vector<vector<int>> tn;         // [l][3] tiles per direction
    vector<vector<unsigned char>> tmask;   // [l] refined tiles of level l (l >= 1)

    vector<r3patch*> P;             // local patches
    vector<vector<int>> lev;        // [l] local patch ids
    vector<r3gpatch> GP;            // all patches
    vector<vector<int>> glev;       // [l] global ids
    vector<r3xplan> xp;             // [l] remote fills
    vector<double> gx[3];           // global level-0 nodes, with margin nodes (index + marge)

    // G 40: the placement of the units (level-1 pieces with their finer patches) and the shadows
    int place = 0;                  // 0 cut at the rank boxes; 1 placed; 2 test of 1 (small pieces, placed anew)
    double rebalance = 0.1;         // a new placement when it lowers the predicted maximum load by more
    int npat = 0;                   // P[0..npat-1]: the patches, P[npat..]: the shadows
    vector<int> col0;               // the rank of every level-0 column (I,J)
    int owner0(const int *I) const; // the rank that holds level-0 cell I (clamped into the domain)
    r3cplan sh_down;                // level 0 down into the arrays of the shadows (all their cells)
    double load_w(int l) const;     // the work of a cell of level l (G 7 1: the steps of the level)
    // items[peer]: grid, cell index here, global level-0 cell I0 I1 I2, tag (6 ints each)
    void cplan_build(r3cplan&, bool down, vector<vector<int>> &items);
    // the transfer: the items packed from src(m, grid) (an array, at the cell index), received into
    // X.rb (per receive peer: field m of item q at m*n + q); collective
    template<class SRC>
    void cplan_xfer(r3cplan &X, int nv, SRC &&src)
    {
        const int ns = (int)X.sp.size(), nr = (int)X.rp.size();
        X.sb.resize(ns);
        X.rb.resize(nr);
        vector<MPI_Request> rq(ns+nr);
        for(int r=0; r<nr; ++r)
        {
            X.rb[r].resize(X.rq[r].size()*(size_t)nv);
            MPI_Irecv(X.rb[r].data(),(int)X.rb[r].size(),MPI_DOUBLE,X.rp[r],3170,MPI_COMM_WORLD,&rq[r]);
        }
        for(int s=0; s<ns; ++s)
        {
            vector<double> &b = X.sb[s];
            const size_t n = X.sq[s].size();
            b.resize(n*nv);
            for(int m=0; m<nv; ++m)
            for(size_t q=0; q<n; ++q)
            b[(size_t)m*n+q] = src(m,X.sg[s][q])[X.sq[s][q]];
            MPI_Isend(b.data(),(int)b.size(),MPI_DOUBLE,X.sp[s],3170,MPI_COMM_WORLD,&rq[nr+s]);
        }
        MPI_Waitall(ns+nr,rq.data(),MPI_STATUSES_IGNORE);
    }
    // down: sent from level 0 (grid -1), received into the grids; up: the other way; set or added
    template<class SRC, class DST>
    void cplan_run(r3cplan &X, int nv, SRC &&src, DST &&dst, bool add)
    {
        cplan_xfer(X,nv,src);
        const int nr = (int)X.rp.size();
        for(int r=0; r<nr; ++r)
        {
            const vector<double> &b = X.rb[r];
            const size_t n = X.rq[r].size();
            for(int m=0; m<nv; ++m)
            for(size_t q=0; q<n; ++q)
            {
                double &v = dst(m,X.rg[r][q])[X.rq[r][q]];
                v = add ? v + b[(size_t)m*n+q] : b[(size_t)m*n+q];
            }
        }
    }
    // the relaxation zones of the wave generation (B 98 2) and the beach (B 99 1/2) may be refined
    // (the model relaxes its patches there; CFD with a wave theory), else they stay on level 0
    bool zone_gen_refine = false, zone_beach_refine = false;

    // index helpers
    int gn(int l, int d) const;                         // global cells of level l
    int ratio(int l, int d) const;                      // rr^l
    bool in_domain(int l, const int *I) const;
    bool wrap(int l, const int *I, int *J) const;       // J: I with the periodic directions wrapped; false if J is
                                                        // outside the domain
    bool wrap_local0(const int *I, int *L) const;       // level-0 cell I (global, may be past a periodic
                                                        // edge) in the box of this rank: L local; else false
    bool refined(int l, const int *I) const;            // level-l cell in the region of level l (l >= 1)
    bool covered(int l, const int *I) const;            // level-l cell covered by level l+1
    int gpatch_at(int l, const int *I) const;           // global patch whose box holds I, -1
    lexer* glex(int g);
    void goff(int g, int *o);                           // local index = global index - o
    int glevel(int g) const { return g<0 ? 0 : P[g]->lev; }
    double node(int l, int d, int I) const;             // coordinate of node I of level l (direction d)

    // fill: eval(kind, g, s, o, double *v) gives nv values for a local source (kind 0 copy,
    // kind 1 prolong; served items are kind 0), store(r3patch*, int id, const r3fill&, const
    // double*) puts them into the patch.  Collective: every rank calls it for every level.
    template<class EV, class ST>
    void fill_run(int l, int nv, EV &&eval, ST &&store, int maxdepth=1<<20)
    {
        static const int o0[3] = {0,0,0};
        r3xplan &X = xp[l];

        if(nranks>1)
        {
            X.sbuf.resize((size_t)X.nsend*nv);
            for(int m=0; m<X.nsend; ++m)
            eval(0,X.serve[m].g,X.serve[m].s,o0,&X.sbuf[(size_t)m*nv]);

            X.rbuf.resize((size_t)X.nrecv*nv);
            xrun(X,nv);
        }

        vector<double> v(nv);
        for(int id : lev[l])
        {
            r3patch *c = P[id];
            for(const r3fill &f : c->fill)
            {
                if(f.kind==3 || f.depth>maxdepth)
                continue;
                if(f.kind==2)
                {
                    store(c,id,f,&X.rbuf[(size_t)f.slot*nv]);
                    continue;
                }
                eval(f.kind,f.g,f.s,f.o,&v[0]);
                store(c,id,f,&v[0]);
            }
        }
    }

protected:
    // setup: hierarchy, patch lexers, plans; seed (optional): refined tiles per level [l][tile]
    void setup_base(lexer*, ghostcell*);
    void setup(lexer*, ghostcell*, const vector<vector<unsigned char>> *seed=nullptr);
    virtual r3patch* patch_new()=0;
    virtual void patch_objects(r3patch*)=0;

    // new hierarchy from the seeds: false if the boxes are unchanged (nothing done); else P, lev, GP,
    // the plans are new, gone holds the old patches not kept (the caller deletes them), fresh the
    // local ids of the new patches (lexer and module objects made, no data)
    bool rebuild(const vector<vector<unsigned char>> &seed, vector<r3patch*> &gone, vector<int> &fresh);
    vector<r3gpatch> GPold;         // the global list before the last rebuild (its patches: Pold, by lid)
    vector<r3patch*> Pold;          // (valid until the caller deletes gone)
    size_t tid(int l, int ti, int tj, int tk) const { return ((size_t)ti*tn[l][1] + tj)*tn[l][2] + tk; }
    bool verbose = true;

    void build_lexer(r3patch&);
    void xrun(r3xplan&, int);

private:
    void global_nodes();
    void mark_tiles(const vector<vector<unsigned char>> *seed);
    void make_boxes();
    void place_units(vector<vector<int>> &box);     // G 40: the level-1 pieces and their ranks
    bool conforming(int d, int X, const int *lo, const int *hi) const;   // a cut of level 1 at X along d
    double unit_work(const int *lo, const int *hi) const;
    void make_shadows(vector<r3patch*> *old, vector<int> &fresh_sh);   // G 40: the shadows of the level-1 patches
    void fetch_raw(const vector<int> &ids);
    void plans();
    bool tile_marked(int l, int ti, int tj, int tk) const;
};

// ghostcell without MPI exchange and with rank-local reductions (patch work)
struct reefamr3d_local
{
    ghostcell *g;
    bool oldc, oldl;
    reefamr3d_local(ghostcell*);
    ~reefamr3d_local();
};

#endif
