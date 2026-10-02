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

#ifndef REEFAMR_H_
#define REEFAMR_H_

#include"increment.h"
#include<vector>
#include<cstddef>

class lexer;
class ghostcell;
class sixdof_obj;

using namespace std;

//  REEFAMR: the module-independent part of the patch-based mesh refinement of REEF3D.
//
//  Level 0 is the native grid of the module.  Every refined patch is a small domain of
//  its own: its own lexer (control keys copied, patch geometry, solid flags from level 0),
//  and, owned by the module, its own fdm and solver objects, so the kernels of the module
//  run unchanged on every grid.  A patch computes EXT cells beyond its box on every side;
//  these cells are filled before each stage from a patch of the same level, prolonged from
//  the next coarser level, or received from the rank that holds them.
//
//  The core owns everything that does not depend on the equations:
//   - the hierarchy: global tile maps per level, tagging (tag_cell), buffer, hysteresis
//     (keep), proper nesting, merging the marked tiles into rectangles cut at the
//     partition edges, regridding with reuse of unchanged patches
//   - the patch lexers: horizontal geometry subdividing the level-0 nodes, solid flags,
//     2D wall lists, the vertical grid (vref: vertical refinement per level)
//   - index conventions, owners, the fill descriptors of the cells around the patches
//     and their exchange plans, the restriction map (2x2 blocks to their coarse cell),
//     the registry of coarse faces next to finer patches (flux matching) and the
//     exchange of the fine face values across partition edges
//   - refinement boxes, no-refinement boxes and in/outflow margins, and the refinement
//     zone that follows moving bodies (margin around the hull, wake wedge)
//   - the point-to-point MPI plans (xsetup/xrun) and a BiCGStab driver for composite
//     elliptic solves on the leaf cells (reefamr_krylov.h)
//
//  The module derives from reefamr and supplies what depends on its equations through the
//  hooks below: its patch objects, the refinement criteria, the state of new patches, the
//  values of the filled cells (prolongation), the restriction, the flux matching of its
//  fluxes, output.  sflow_amr (REEF3D::SFLOW) is the first module.
//
//  Index conventions: every level has a global cell index space, level l+1 refines level l
//  by 2 in x and y.  A patch covers the global box [I0,I1]x[J0,J1] of its level; its lexer
//  index is i = I-I0+EXT.  Grid id -1 is level 0, otherwise an index into P.  Per grid
//  arrays of the core (match, rmatch) are indexed with id+1.
//
//  Ownership: holder(l,I,J) is the rank that holds cell (I,J) of level l and its coarser
//  ancestors.  Patches are cut at the level-0 rank boxes, so this is the owner of the
//  level-0 cell below.  All plans ask holder(), so that patches can later be placed on
//  other ranks (load balancing) by changing holder() and the patch placement only.

// a cell outside the patch interior, filled before each stage
struct reefamr_fill
{
    int di,dj;          // destination cell (patch lexer index)
    int kind;           // 0: copy from grid g, 1: prolong from grid g (coarser level), 2: remote, 3: served without data
    int g;              // grid id (-1: level 0)
    int si,sj;          // source cell on grid g; remote: receive peer index in si
    int ox,oy;          // prolongation: quadrant of the fine cell (-1/+1)
    int slot;           // remote: position in the receive buffer of the level
    double aux;         // module value of the destination cell (SFLOW: still water depth)
};

// base of the module patches
struct reefamr_patch
{
    virtual ~reefamr_patch(){}

    int lev;
    int I0,I1,J0,J1;                // global box of the patch on its level
    int nx,ny;

    lexer *pp = nullptr;

    // cells filled before each stage
    vector<reefamr_fill> fill;

    // restriction: coarse target (grid id, i, j) for every 2x2 block (-2: wall)
    vector<int> rgrid, ric, rjc;

    bool fresh;                     // created by the current regrid
    bool bc2D = false;              // 2D wall lists built (build_bc2D)
    bool zown = false;              // own vertical grid arrays (vref)
};

// a coarse face next to a finer patch, overridden with fine values
struct reefamr_match
{
    static const int NVAL = 6;
    int dir;          // 0: x face, 1: y face
    int fi,fj;        // local face index on the target grid
    int child;        // patch that recorded the fine faces (-1: remote)
    int side,r;       // recorded side (0: low x, 1: high x, 2: low y, 3: high y) and first fine index
    double val[NVAL]; // remote: averaged fine values (module layout)
};

// point-to-point exchange with a fixed set of peers
struct reefamr_xplan
{
    vector<int> speer, rpeer;            // ranks
    vector<vector<int>> sitem;           // per send peer: item indices (meaning depends on the plan)
    vector<int> rcount;                  // per receive peer: number of items
    vector<vector<double>> sbuf, rbuf;
};

// refinement zone of a moving body
struct reefamr_zone
{
    double cx,cy,ex,ey,smin,smax,nmin,nmax,sfront,wake,bx0,bx1,by0,by1;
};

// parameters set by the module (configure)
struct reefamr_param
{
    const char *name = "AMR";   // prefix of the messages
    int maxlev = 0;             // levels above level 0
    int regrid = 0;             // steps between regrids, 0: static
    int nbuf = 0;               // buffer cells around flagged cells
    int tile = 8;               // tile width (even, >= 4)
    int nest = 3;               // proper nesting width in coarse cells
    int keep = 0;               // regrids a refined tile is kept after its last flag (<= 250)
    int ext = 2;                // computed cells beyond the patch box on every side
    vector<int> vref;           // vertical refinement (1 or 2) of level l over l-1, index l; empty: 1
    vector<double> rbox;        // refinement boxes: xs,xe,ys,ye
    vector<double> fbox;        // no-refinement boxes: xs,xe,ys,ye
    int ioband = 4;             // no refinement within this many level-0 cells of in- and outflow
    bool zones = false;         // refinement zone around the moving bodies (zone_bodies)
    double zr = 0.0;            // zone: margin around the hull
    double zL = 0.0;            // zone: length of the wake wedge (at most)
    double za = 0.0;            // zone: half angle of the wake wedge in degrees
};

// switches the MPI exchange of the ghostcell class off while patch kernels run
struct reefamr_comms_off
{
    ghostcell *g;
    bool old;
    reefamr_comms_off(ghostcell *gg);
    ~reefamr_comms_off();
};

class reefamr : public increment
{
public:
    reefamr(lexer*, ghostcell*);
    virtual ~reefamr();

    int patches_total;

protected:
    void configure(const reefamr_param&);

    // ---- module hooks
    // allocate an (empty) module patch
    virtual reefamr_patch* patch_new() = 0;
    // module objects of a new patch (lexer set up, comms off); usually build_bc2D first
    virtual void patch_objects(reefamr_patch*, ghostcell*) = 0;
    // delete the module objects of a patch
    virtual void patch_delete(reefamr_patch*) = 0;
    // cells of level l on this rank that need level l+1: tag_cell(l+1,...)
    virtual void tag(int, vector<unsigned char>&) = 0;
    // regrid: final state of the current hierarchy with filled cells, before the tagging
    virtual void regrid_prepare(ghostcell*) = 0;
    // regrid: the patch list is new (grid ids changed)
    virtual void regrid_ids() {}
    // regrid: time-independent data of the fresh patches, before the plans
    virtual void regrid_static(ghostcell*) = 0;
    // regrid: state of the fresh patches (coarse to fine), old patches still alive
    virtual void regrid_state(ghostcell*, vector<reefamr_patch*>&) = 0;
    // regrid: end, old patches freed, counts updated (restriction to level 0, halo)
    virtual void regrid_finish(ghostcell*, int) = 0;
    // aux value of a filled cell of a patch (reefamr_fill::aux) and of a cell served to another rank
    virtual double fill_aux(reefamr_patch*, int, int) {return 0.0;}
    virtual double serve_aux(int, int, int) {return 0.0;}
    // moving bodies of the refinement zone
    virtual void zone_bodies(vector<sixdof_obj*>&) {}

    // ---- hierarchy
    void setup(lexer*, ghostcell*);
    void regrid(lexer*, ghostcell*, bool);
    void free_patches();                    // all patches (module destructor)
    int patch_at(int, int, int);            // level, I, J: interior patch on this rank, -1: level 0 (l==0), -2: none, -3: outside the rank box
    int owner(int, int);                    // rank of the level-0 cell (I,J), -1: outside the domain
    int holder(int, int, int);              // rank that holds cell (I,J) of level l and its ancestors
    int flag0(int, int);                    // level-0 flagslice4 of the global cell (I,J), from the rank box + halo
    void goff(int, int&, int&);             // grid id: local index = global index - offset
    lexer* glex(int);                       // grid id: lexer
    void tag_cell(int, int, int, int, vector<unsigned char>&);   // level lf, coarse cell I,J of level lf-1, buffer
    void global_or(vector<unsigned char>&);
    void build_bc2D(ghostcell*, reefamr_patch&);

    // ---- exchange
    void xsetup(vector<vector<int>>&, vector<vector<int>>&, int);
    void xrun(reefamr_xplan&, int, int);

    // cells around the level-l patches: eval(const reefamr_fill&, double*) gives nv values of a
    // local source (kind 0, 1, 3), store(reefamr_patch*, int id, const reefamr_fill&, const double*)
    // puts them into the patch
    template<class EV, class ST>
    void fill_run(int l, int nv, int tag, EV &&eval, ST &&store)
    {
        reefamr_xplan &X = gplan[l];

        for(size_t k=0; k<X.speer.size(); ++k)
        {
            vector<double> &sb = X.sbuf[k];
            sb.resize(X.sitem[k].size()*nv);
            for(size_t m=0; m<X.sitem[k].size(); ++m)
            eval(gserve[l][X.sitem[k][m]],&sb[m*nv]);
        }

        xrun(X,nv,tag);

        fillv.resize(nv);
        for(int id : lev[l])
        {
            reefamr_patch *c = P[id];
            for(auto &f : c->fill)
            {
                if(f.kind==2)
                {
                    store(c,id,f,&X.rbuf[f.si][(size_t)f.slot*nv]);
                    continue;
                }
                eval(f,&fillv[0]);
                store(c,id,f,&fillv[0]);
            }
        }
    }

    // as fill_run, but stores only the fill entries sel(reefamr_patch*) lists (a const vector<int>*
    // of indices into its fill, nullptr: none); all served cells are still evaluated and sent
    template<class EV, class ST, class SEL>
    void fill_run_sub(int l, int nv, int tag, EV &&eval, ST &&store, SEL &&sel)
    {
        reefamr_xplan &X = gplan[l];

        for(size_t k=0; k<X.speer.size(); ++k)
        {
            vector<double> &sb = X.sbuf[k];
            sb.resize(X.sitem[k].size()*nv);
            for(size_t m=0; m<X.sitem[k].size(); ++m)
            eval(gserve[l][X.sitem[k][m]],&sb[m*nv]);
        }

        xrun(X,nv,tag);

        fillv.resize(nv);
        for(int id : lev[l])
        {
            reefamr_patch *c = P[id];
            const vector<int> *E = sel(c);
            if(E==nullptr)
            continue;
            for(int e : *E)
            {
                const reefamr_fill &f = c->fill[e];
                if(f.kind==2)
                {
                    store(c,id,f,&X.rbuf[f.si][(size_t)f.slot*nv]);
                    continue;
                }
                eval(f,&fillv[0]);
                store(c,id,f,&fillv[0]);
            }
        }
    }

    // fine face values of the level-l patches on the partition edges, to the rank of the coarse
    // cell: pack(reefamr_patch*, side, r, double*) gives nv values of the fine faces r, r+1,
    // unpack(reefamr_match&, const double*) stores them in the remote match entry
    template<class PK, class UP>
    void face_run(int l, int nv, int tag, PK &&pack, UP &&unpack)
    {
        reefamr_xplan &F = fplan[l];

        for(size_t k=0; k<F.speer.size(); ++k)
        {
            vector<double> &sb = F.sbuf[k];
            sb.resize(F.sitem[k].size()*nv);
            for(size_t m=0; m<F.sitem[k].size(); ++m)
            {
                int it = F.sitem[k][m];
                reefamr_patch *c = P[fsend[l][3*it]];
                pack(c,fsend[l][3*it+1],fsend[l][3*it+2],&sb[m*nv]);
            }
        }

        xrun(F,nv,tag);

        size_t it=0;
        for(size_t k=0; k<F.rpeer.size(); ++k)
        for(int m=0; m<F.rcount[k]; ++m)
        {
            int tgt = frecv[l][2*it], idx = frecv[l][2*it+1];
            if(tgt>=0)
            unpack(rmatch[tgt][idx],&F.rbuf[k][(size_t)m*nv]);
            ++it;
        }
    }

    // ---- moving bodies
    void zone_setup(lexer*);
    bool zone_test(double, double);
    vector<reefamr_zone> zones;
    vector<vector<int>> ztri;       // per body: triangles reaching below the still water level

    // ---- data
    lexer *p0;
    ghostcell *pgc0;
    reefamr_param par;
    vector<reefamr_patch*> P;
    vector<vector<int>> lev;        // patch ids per level (index 0 unused)
    vector<int> nlevg;              // patches per level, all ranks
    int maxlev, nest, tile, nbuf, regrid_int, keep;
    int EXT;                        // extra computed cells on each side of a patch
    int regrids;
    long cells_total;

    // rank boxes of level 0 (global cell indices), all ranks
    int O0i,O0j,NX0,NY0,GNX,GNY;
    vector<int> rbx0,rbx1,rby0,rby1;
    int bxlo(int l) const { return O0i<<l; }
    int bxhi(int l) const { return ((O0i+NX0)<<l)-1; }
    int bylo(int l) const { return O0j<<l; }
    int byhi(int l) const { return ((O0j+NY0)<<l)-1; }

    // flux matching entries per target grid (index id+1): local and remote
    vector<vector<reefamr_match>> match;
    vector<vector<reefamr_match>> rmatch;

    // ghost service per level: items are reefamr_fill descriptors of the cells served
    vector<reefamr_xplan> gplan;
    vector<vector<reefamr_fill>> gserve;        // [l]: served cells (kind 0/1/3)
    vector<vector<reefamr_fill*>> grecv;        // [l]: receive slot -> fill entry

    // fine face records sent across partition edges, per level
    vector<reefamr_xplan> fplan;
    vector<vector<int>> fsend;                  // [l]: (patch, side, r) triplets
    vector<vector<int>> frecv;                  // [l]: target grid id+1 and index into rmatch, -1: skip

private:
    reefamr_patch* make_patch(lexer*, ghostcell*, int, int, int, int, int);
    void free_patch(reefamr_patch*);
    void build_lexer(lexer*, reefamr_patch&);
    void build_vertical(lexer*, reefamr_patch&);
    void build_flags(lexer*);
    void build_tiles();
    void build_plans(ghostcell*);

    vector<int> fl0;
    static const int FH = 8;
    int last_owner;

    // global tile maps per level (1 = refined), and the local tile -> patch map
    vector<vector<unsigned char>> gtile;
    vector<vector<unsigned char>> tage;     // regrids since a tile was last flagged (keep)
    vector<int> gtnx,gtny;
    vector<vector<int>> tmap;
    vector<int> tti0,ttj0,tnx,tny;

    // vertical refinement factor of each level over level 0
    vector<int> vfac;

    // no refinement: global level-0 cells
    vector<unsigned char> forbid0;

    // initial position of every body of the refinement zone
    vector<double> zx0, zy0;

    vector<double> fillv;
};

#endif
