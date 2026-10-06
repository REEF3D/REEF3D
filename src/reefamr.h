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
#include<algorithm>

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
//  Ownership: holder(l,I,J) is the rank that holds cell (I,J) of level l: the rank of the
//  level-l patch there (the global patch table GP), else the holder of its parent cell, down to
//  the owner of the level-0 cell.  Patches are cut at the level-0 rank boxes, so this is the
//  owner of the level-0 cell below.  All plans ask holder(), and the parent cells of the 2x2
//  blocks of a patch are reached through the block plans (block_up, block_down), so that patches
//  can be placed on other ranks (load balancing) by changing the patch placement only.

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

// the parent cell of a 2x2 block of a patch, on the rank that holds it: grid id, local cell
struct reefamr_block
{
    int g, ic, jc;
};

// a patch of the hierarchy, known to all ranks: level, global box, rank, index into P on that rank
struct reefamr_gpatch
{
    int lev, I0, I1, J0, J1, rank, lid;
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

    // restriction: coarse target (grid id, i, j) for every 2x2 block (-2: wall, -3: on another
    // rank, served through the block plans: block_up, block_down)
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
    double sback,nlo,nhi;       // hull rectangle in (s,n): [sback,sfront] x [nlo,nhi]
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
    bool zalign = false;        // zone: rectangle aligned with x and y around the wetted hull, grown by
                                // the distance travelled until the next regrid (moored and
                                // oscillating bodies), instead of aligned with the direction of motion
    double lazy = 0.0;          // regrid: the refined tiles of the current layout are kept as long as
                                // the union with the flagged tiles has at most lazy times the flagged
                                // tiles, an unchanged layout skips the regrid (0: off)
    int dryband = 0;            // regrid: no refinement within this many level-0 cells of a cell the
                                // module reports as unfit for a patch (cell_unfit, e.g. dry), 0: off
    int place = 0;              // patches on several ranks: 0 cut at the rank boxes, on the rank of the
                                // level-0 cells below; 1 not cut, split to the work per rank and placed
                                // for the load of the ranks (the module must reach parents through the
                                // block plans and old patches through old_run); 2 test of 1: small
                                // pieces, placed anew at every regrid without the preference for the
                                // rank below, so that parents and old patches are mostly on other ranks
    double rebalance = 0.1;     // place 1: a new placement when it lowers the predicted maximum load by
                                // more than this fraction, else the patches keep their ranks
    bool place_whole_zones = false; // place 1/2: a rectangle that reaches the hull of a body (its box
                                // grown by zr) is placed whole, not split (FNPF: the footprint
                                // extension of the body is local to a patch)
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
    // level-0 cell (ii,jj) of this rank that no patch may cover at this regrid (par.dryband > 0)
    virtual bool cell_unfit(int, int) {return false;}

    // ---- hierarchy
    void setup(lexer*, ghostcell*);
    void regrid(lexer*, ghostcell*, bool);
    void free_patches();                    // all patches (module destructor)
    int patch_at(int, int, int);            // level, I, J: interior patch on this rank, -1: level 0 (l==0), -2: none, -3: outside the rank box
    int owner(int, int);                    // rank of the level-0 cell (I,J), -1: outside the domain
    int holder(int, int, int);              // rank that holds cell (I,J) of level l (the patch of level l there, else its parent)
    int gpatch_at(int, int, int);           // level, I, J: index into GP of the patch of that level there, -1: none
    bool covered(int l, int I, int J) { return l>=1 && gpatch_at(l,I,J)>=0; }   // a level-l patch (any rank) holds (I,J)
    double gnode(int, int, int);            // place 1: x (dir 0) or y (dir 1) of global node N of level l
    int finest_rank(double, double);        // place 1: rank of the finest grid at (x,y), -1 outside
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

    // as face_run, unpack(int tgt, int idx, const double*) gets the remote match entry by its
    // target grid (id+1) and its index in rmatch[tgt] (module arrays parallel to rmatch)
    template<class PK, class UP>
    void face_run_at(int l, int nv, int tag, PK &&pack, UP &&unpack)
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
            unpack(tgt,idx,&F.rbuf[k][(size_t)m*nv]);
            ++it;
        }
    }

    // parent cells of the 2x2 blocks of the level-l patches.  Every block has a key, unique on the
    // rank for the level and fixed for the layout: the blocks of the local patches (lev[l], block
    // order) first, then the blocks this rank serves to other ranks; block_keys(l) is their number.
    //
    // down: eval(const reefamr_block&, int key, double*) gives nv values of the parent on its rank,
    // store(reefamr_patch*, int id, int k, const double*) puts them into block k of patch id
    template<class EV, class ST>
    void block_down(int l, int nv, int tag, EV &&eval, ST &&store)
    {
        block_down_if(l,nv,tag,[](reefamr_patch*) { return true; },eval,store);
    }

    // as block_down, for the local patches with need(reefamr_patch*) only (the served parents are
    // all evaluated and sent)
    template<class ND, class EV, class ST>
    void block_down_if(int l, int nv, int tag, ND &&need, EV &&eval, ST &&store)
    {
        reefamr_xplan &X = bdn[l];

        for(size_t k=0; k<X.speer.size(); ++k)
        {
            vector<double> &sb = X.sbuf[k];
            sb.resize(X.sitem[k].size()*nv);
            for(size_t m=0; m<X.sitem[k].size(); ++m)
            {
                const int it = X.sitem[k][m];
                if(bsrv[l][it].g<-1)
                {
                    for(int v=0; v<nv; ++v)
                    sb[m*nv+v] = 0.0;
                    continue;
                }
                eval(bsrv[l][it],bnloc[l]+it,&sb[m*nv]);
            }
        }

        xrun(X,nv,tag);

        blockv.resize(nv);
        int key=0;
        for(int id : lev[l])
        {
            reefamr_patch *c = P[id];
            const bool nd = need(c);
            for(size_t k=0; k<c->rgrid.size(); ++k)
            if(c->rgrid[k]>=-1)
            {
                if(nd)
                {
                    eval(reefamr_block{c->rgrid[k],c->ric[k],c->rjc[k]},key,&blockv[0]);
                    store(c,id,(int)k,&blockv[0]);
                }
                ++key;
            }
        }

        size_t it=0;
        for(size_t k=0; k<X.rpeer.size(); ++k)
        for(int m=0; m<X.rcount[k]; ++m)
        {
            const int id = bloc[l][2*it], kb = bloc[l][2*it+1];
            if(need(P[id]))
            store(P[id],id,kb,&X.rbuf[k][(size_t)m*nv]);
            ++it;
        }
    }

    // up: compute(reefamr_patch*, int id, int k, double*) gives nv values of block k of patch id,
    // store(const reefamr_block&, int key, const double*) puts them into the parent on its rank
    template<class CP, class ST>
    void block_up(int l, int nv, int tag, CP &&compute, ST &&store)
    {
        reefamr_xplan &X = bup[l];

        for(size_t k=0; k<X.speer.size(); ++k)
        {
            vector<double> &sb = X.sbuf[k];
            sb.resize(X.sitem[k].size()*nv);
            for(size_t m=0; m<X.sitem[k].size(); ++m)
            {
                const int it = X.sitem[k][m];
                const int id = bloc[l][2*it], kb = bloc[l][2*it+1];
                compute(P[id],id,kb,&sb[m*nv]);
            }
        }

        xrun(X,nv,tag);

        blockv.resize(nv);
        int key=0;
        for(int id : lev[l])
        {
            reefamr_patch *c = P[id];
            for(size_t k=0; k<c->rgrid.size(); ++k)
            if(c->rgrid[k]>=-1)
            {
                compute(c,id,(int)k,&blockv[0]);
                store(reefamr_block{c->rgrid[k],c->ric[k],c->rjc[k]},key++,&blockv[0]);
            }
        }

        size_t it=0;
        for(size_t k=0; k<X.rpeer.size(); ++k)
        for(int m=0; m<X.rcount[k]; ++m)
        {
            if(bsrv[l][it].g>=-1)
            store(bsrv[l][it],bnloc[l]+(int)it,&X.rbuf[k][(size_t)m*nv]);
            ++it;
        }
    }

    int block_keys(int l) const { return bnloc[l] + (int)bsrv[l].size(); }

    // rank and number of ranks (lexer is incomplete in this header: the templates use these)
    int my_rank() const;
    int n_ranks() const;

    // the cells of the fresh level-l patches that an old patch of the same level held (regrid_state,
    // old patches still alive): pack(reefamr_patch *old, int io, int jo, double*) gives nv values of
    // the old patch on its rank, unpack(reefamr_patch*, int id, int ii, int jj, const double*) puts
    // them into the fresh patch (lexer indices)
    template<class PK, class UP>
    void old_run(int l, int nv, int tag, vector<reefamr_patch*> &oldP, PK &&pack, UP &&unpack)
    {
        const int me = my_rank();
        const int np = n_ranks();
        vector<vector<int>> req(np), srv, dst(np);
        vector<double> v(nv);

        for(int id : lev[l])
        {
            reefamr_patch *c = P[id];
            if(!c->fresh)
            continue;

            for(const reefamr_gpatch &G : GPold)
            {
                if(G.lev!=l)
                continue;
                const int I0 = std::max(c->I0,G.I0), I1 = std::min(c->I1,G.I1);
                const int J0 = std::max(c->J0,G.J0), J1 = std::min(c->J1,G.J1);
                if(I0>I1 || J0>J1)
                continue;

                if(G.rank==me)
                {
                    reefamr_patch *o = oldP[G.lid];
                    if(o==c)
                    continue;
                    for(int I=I0; I<=I1; ++I)
                    for(int J=J0; J<=J1; ++J)
                    {
                        pack(o,I-o->I0+EXT,J-o->J0+EXT,&v[0]);
                        unpack(c,id,I-c->I0+EXT,J-c->J0+EXT,&v[0]);
                    }
                    continue;
                }

                for(int I=I0; I<=I1; ++I)
                for(int J=J0; J<=J1; ++J)
                {
                    req[G.rank].push_back(G.lid); req[G.rank].push_back(I); req[G.rank].push_back(J);
                    dst[G.rank].push_back(id); dst[G.rank].push_back(I-c->I0+EXT); dst[G.rank].push_back(J-c->J0+EXT);
                }
            }
        }

        xsetup(req,srv,3);

        reefamr_xplan X;
        for(int r=0; r<np; ++r)
        {
            if(!srv[r].empty())
            {
                X.speer.push_back(r);
                X.sitem.push_back(vector<int>());
                vector<double> sb(srv[r].size()/3*nv);
                for(size_t k=0; k<srv[r].size(); k+=3)
                {
                    reefamr_patch *o = oldP[srv[r][k]];
                    pack(o,srv[r][k+1]-o->I0+EXT,srv[r][k+2]-o->J0+EXT,&sb[k/3*nv]);
                }
                X.sbuf.push_back(sb);
            }
            if(!req[r].empty())
            {
                X.rpeer.push_back(r);
                X.rcount.push_back((int)req[r].size()/3);
            }
        }

        xrun(X,nv,tag);

        for(size_t k=0; k<X.rpeer.size(); ++k)
        {
            const vector<int> &D = dst[X.rpeer[k]];
            for(size_t m=0; m<D.size(); m+=3)
            unpack(P[D[m]],D[m],D[m+1],D[m+2],&X.rbuf[k][m/3*nv]);
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
    int regrids_skipped = 0;        // regrids with an unchanged layout (par.lazy)
    long cells_total;
    long cells_local = 0;           // refined columns of the patches on this rank

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

    // patches of all ranks (built at every regrid, the same on every rank) and, per level and
    // global tile, the patches that overlap the tile; the table of the previous layout
    vector<reefamr_gpatch> GP, GPold;
    vector<vector<int>> gtoff, gtlist;          // [l]: per global tile the range gtoff[t]..gtoff[t+1]-1 of gtlist

    // block plans per level: on the rank of a patch the blocks whose parent is on another rank
    // ((id,k) pairs, in the order of the peers), on the rank of the parent the blocks it serves
    vector<reefamr_xplan> bup, bdn;
    vector<vector<int>> bloc;
    vector<vector<reefamr_block>> bsrv;
    vector<int> bnloc;                          // [l]: blocks of the local patches with a local parent

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
    void build_gtable();
    void place_patches(lexer*, ghostcell*, vector<vector<unsigned char>>&, vector<reefamr_patch*>&, vector<char>&,
                       vector<reefamr_patch*>&, vector<vector<int>>&);
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

    // no refinement: global level-0 cells, static and (dryband) at the current regrid
    vector<unsigned char> forbid0;
    vector<unsigned char> forbidd;

    // initial position of every body of the refinement zone
    vector<double> zx0, zy0;

    vector<double> fillv, blockv;

    // place 1: global level-0 nodes (index K+marge) and solid flags
    vector<double> gxn, gyn;
    vector<short> gfl0;
    double x0g(int) const;
    double y0g(int) const;
};

#endif
