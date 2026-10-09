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
//  Static refinement (step 15b): the boxes G 15 xs xe ys ye zs ze are refined to level G 1, and
//  every coarser level covers the next finer one with a nesting band of 3 cells.

struct r3fill
{
    int d[3];           // destination cell (patch index, in the margin)
    int kind;           // 0: copy from grid g (same level), 1: prolong from grid g (next coarser level),
                        // 2: remote (copy from a patch of another rank: slot in the receive buffer), 3: boundary
    int g;              // grid id (-1: level 0)
    int s[3];           // source cell on grid g (kind 1: the parent cell)
    int o[3];           // kind 1: position of the child in the parent (0, 1; 0 in an inactive direction)
    int slot;           // kind 2
};

struct r3patch
{
    virtual ~r3patch(){}

    int lev = 0;
    int gid = -1;       // index into the global list
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
    vector<vector<int>> tn;         // [l][3] tiles per direction
    vector<vector<unsigned char>> tmask;   // [l] refined tiles of level l (l >= 1)

    vector<r3patch*> P;             // local patches
    vector<vector<int>> lev;        // [l] local patch ids
    vector<r3gpatch> GP;            // all patches
    vector<vector<int>> glev;       // [l] global ids
    vector<r3xplan> xp;             // [l] remote fills
    vector<double> gx[3];           // global level-0 nodes, with margin nodes (index + marge)

    // index helpers
    int gn(int l, int d) const;                         // global cells of level l
    int ratio(int l, int d) const;                      // rr^l
    bool in_domain(int l, const int *I) const;
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
    void fill_run(int l, int nv, EV &&eval, ST &&store)
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
                if(f.kind==3)
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
    // setup: hierarchy, patch lexers, plans
    void setup(lexer*, ghostcell*);
    virtual r3patch* patch_new()=0;
    virtual void patch_objects(r3patch*)=0;

    void build_lexer(r3patch&);
    void xrun(r3xplan&, int);

private:
    void global_nodes();
    void mark_tiles();
    void make_boxes();
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
