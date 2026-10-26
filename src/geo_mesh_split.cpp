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

#include"geo_mesh.h"
#include"lexer.h"
#include<iostream>
#include<cstdlib>
#include<algorithm>
#include<cmath>

// G 55 m: the blocks 1.. of a multi-block grid (G 50 1, DIVEMesh R 20 1) cut into pieces, so that the
// blocks 1.. spread over the ranks (every block 1.. is a level-0 grid of REEFAMR on one rank,
// reefamr::blocks_setup places them for the load).
//
// The pieces of a block are slices along its i direction (along j for a block with more cells in j and
// without open faces, G 54, B 522: the open -i and +i faces stay whole on the first and the last piece),
// of at most (cells of the blocks 1..) / (m ranks) cells and at least 16 cells deep (the cells beyond a face
// reach EXT + margin into the neighbour).  Every piece is a block of its own: the first piece keeps the
// number of its block, the others follow the blocks of DIVEMesh.  The interfaces are rebuilt cell by cell
// along every face of every piece (block 0 included): a cell beyond the face inside the same block of
// DIVEMesh is in the next piece (identity), else the interface of DIVEMesh maps it into the neighbour block,
// into the piece that holds it there; the transforms compose (the index maps are linear).  Consecutive
// cells with the same piece and transform form a segment.  A check maps every segment back.
//
// Run as for the blocks of DIVEMesh; the input (G 54 b, B 522 b) keeps the numbers of DIVEMesh (face_piece).
//
// G 56 a dz: the cells count with their load (cell_load): a wet cell 1, a cell whose bed lies above the still
// water level + dz the load a.  The size limit is then (load of the blocks 1..) / (m ranks), and a block is cut
// where its slices reach equal parts of its load (at least 16 cells deep), so that the pieces carry about the
// same work, not the same number of cells.

namespace
{
struct piece { int id, oi, oj, ni, nj; };

struct segkey
{
    int nbr=-1, nface=0, perm[2]={0,1}, sign[2]={1,1}, off[2]={0,0};
    bool operator==(const segkey &o) const
    {
        return nbr==o.nbr && nface==o.nface && perm[0]==o.perm[0] && perm[1]==o.perm[1]
            && sign[0]==o.sign[0] && sign[1]==o.sign[1] && off[0]==o.off[0] && off[1]==o.off[1];
    }
};
}

int geo_mesh::face_piece(int b, int f) const
{
    if(b<0 || b>=nblk_orig)
    return -1;

    int found=-1, n=0;
    for(int q=0; q<int(blk.size()); ++q)
    {
        const blk_block &B = blk[q];
        if(B.orig!=b)
        continue;
        const bool on = (f==0 && B.oi==0) || (f==1 && B.oi+B.ni==orig_ni[b])
                     || (f==2 && B.oj==0) || (f==3 && B.oj+B.nj==orig_nj[b]);
        const bool whole = (f<2) ? (B.oj==0 && B.nj==orig_nj[b]) : (B.oi==0 && B.ni==orig_ni[b]);
        if(on)
        {
            ++n;
            if(whole)
            found = q;
        }
    }
    return (n==1) ? found : -1;
}

// G 56 a dz: the still water level the dry cells are measured from: F 60, else F 62 with F 65 (river
// initialisation), else none (every cell counts 1)
double geo_mesh::load_level(const lexer *p)
{
    if(p->F60>-1.0e19)
    return p->F60;
    if(p->F65>0 && p->F62>-1.0e19)
    return p->F62;
    return -1.0e20;
}

bool geo_mesh::load_weighted(const lexer *p)
{
    return p->G56==1 && p->G56_a!=1.0 && load_level(p)>-1.0e19;
}

// the load of a cell with bed level zb, relative to a wet cell
double geo_mesh::cell_load(const lexer *p, double zb)
{
    if(!load_weighted(p))
    return 1.0;
    return (zb > load_level(p) + p->G56_z) ? p->G56_a : 1.0;
}

// the load of cell (ci,cj) of block b (1..): its bed is the mean of its four nodes
double geo_mesh::blk_cell_load(const lexer *p, int b, int ci, int cj) const
{
    const blk_block &B = blk[b];
    const size_t w = size_t(B.ni+1);
    const double zb = 0.25*(B.zb[ci+w*cj] + B.zb[ci+1+w*cj] + B.zb[ci+w*(cj+1)] + B.zb[ci+1+w*(cj+1)]);
    return cell_load(p,zb);
}

// the load of block b (1..)
double geo_mesh::blk_load(const lexer *p, int b) const
{
    const blk_block &B = blk[b];
    if(!load_weighted(p))
    return double(B.ni)*double(B.nj);
    double w=0.0;
    for(int cj=0; cj<B.nj; ++cj)
    for(int ci=0; ci<B.ni; ++ci)
    w += blk_cell_load(p,b,ci,cj);
    return w;
}

void geo_mesh::split_blocks(lexer *p)
{
    const int nb = int(blk.size());
    const int m = p->G55;
    const int np = std::max(p->mpi_size,1);
    const int MIN = 16;

    if(m<=0 || nb<=1)
    return;

    auto fail = [&](const string &msg)
    {
        if(p->mpirank==0)
        cout<<"!!! G 55: "<<msg<<" !!!"<<endl;
        exit(1);
    };

    // G 56: the load of the cells (1 each without it)
    const bool wgt = load_weighted(p);
    long long tot=0;
    double wtot=0.0;
    vector<double> wb(nb,0.0);
    for(int b=1; b<nb; ++b)
    {
        tot += (long long)blk[b].ni*blk[b].nj;
        wb[b] = blk_load(p,b);
        wtot += wb[b];
    }
    const long long S = std::max(1LL,(tot + (long long)m*np - 1)/((long long)m*np));
    const double SW = wtot/double((long long)m*np);

    // a block with open faces is cut along i only (its -i and +i faces stay whole)
    auto open = [&](int b)
    {
        for(int n=0; n<p->G54; ++n)
        if(p->G54_b[n]==b)
        return true;
        for(int n=0; n<p->B522; ++n)
        if(p->B522_b[n]==b)
        return true;
        return false;
    };

    // the pieces of every block
    vector<vector<piece>> pc(nb);
    pc[0].push_back({0,0,0,blk[0].ni,blk[0].nj});
    int next = nb;
    bool any = false;
    for(int b=1; b<nb; ++b)
    {
        const blk_block &B = blk[b];
        const int dir = (B.nj>B.ni && !open(b)) ? 1 : 0;
        const int L = dir ? B.nj : B.ni;
        const long long c = (long long)B.ni*B.nj;
        int k = wgt ? int(std::ceil(wb[b]/std::max(SW,1.0e-12) - 1.0e-9)) : int((c + S - 1)/S);
        k = std::min(k,std::max(1,L/MIN));
        k = std::max(k,1);

        // G 56: the cuts where the slices along the cut reach q/k of the load of the block
        vector<int> cut(k+1,0);
        cut[k] = L;
        for(int q=1; q<k; ++q)
        cut[q] = int((long long)q*L/k);
        if(wgt && k>1)
        {
            vector<double> cum(L+1,0.0);
            for(int x=0; x<L; ++x)
            {
                double w=0.0;
                const int M = dir ? B.ni : B.nj;
                for(int y=0; y<M; ++y)
                w += dir ? blk_cell_load(p,b,y,x) : blk_cell_load(p,b,x,y);
                cum[x+1] = cum[x] + w;
            }
            for(int q=1; q<k; ++q)
            {
                const double target = cum[L]*double(q)/double(k);
                int x = 0;
                while(x<L && cum[x+1]<=target)
                ++x;
                if(x<L && target-cum[x] > cum[x+1]-target)
                ++x;
                cut[q] = std::max(cut[q-1]+MIN, std::min(x, L-(k-q)*MIN));
            }
        }

        for(int q=0; q<k; ++q)
        {
            const int a = cut[q], e = cut[q+1];
            piece P;
            P.id = (q==0) ? b : next++;
            P.oi = dir ? 0 : a;
            P.oj = dir ? a : 0;
            P.ni = dir ? B.ni : e-a;
            P.nj = dir ? e-a : B.nj;
            pc[b].push_back(P);
        }
        if(k>1)
        any = true;
    }

    if(!any)
    {
        if(p->mpirank==0)
        {
        if(wgt)
        cout<<"G 55 "<<m<<": no block 1.. has a load above "<<SW<<" (load of the blocks 1.. / "<<m*np<<", G 56 "<<p->G56_a<<" "<<p->G56_z<<"), the blocks stay whole"<<endl;
        else
        cout<<"G 55 "<<m<<": no block 1.. has more than "<<S<<" cells (cells of the blocks 1.. / "<<m*np<<"), the blocks stay whole"<<endl;
        }
        return;
    }

    const vector<blk_block> ob = blk;
    blk.assign(next,blk_block());

    // geometry of the pieces: size, place in the block of DIVEMesh, nodes
    for(int b=0; b<nb; ++b)
    for(const piece &P : pc[b])
    {
        blk_block &N = blk[P.id];
        const blk_block &O = ob[b];
        N.ni = P.ni; N.nj = P.nj;
        N.orig = b; N.oi = P.oi; N.oj = P.oj;
        for(int f=0; f<4; ++f)
        N.ring[f] = (b==0) ? O.ring[f] : 0;
        if(b>0)
        {
            const size_t nn = size_t(P.ni+1)*size_t(P.nj+1);
            N.x.resize(nn); N.y.resize(nn); N.zb.resize(nn);
            for(int L=0; L<=P.nj; ++L)
            for(int K=0; K<=P.ni; ++K)
            {
                const size_t s = size_t(K+P.oi) + size_t(O.ni+1)*size_t(L+P.oj);
                const size_t d = size_t(K) + size_t(P.ni+1)*size_t(L);
                N.x[d] = O.x[s]; N.y[d] = O.y[s]; N.zb[d] = O.zb[s];
            }
        }
    }

    // the piece of block b that holds its cell (ci,cj)
    auto find = [&](int b, int ci, int cj) -> const piece*
    {
        for(const piece &P : pc[b])
        if(ci>=P.oi && ci<P.oi+P.ni && cj>=P.oj && cj<P.oj+P.nj)
        return &P;
        return nullptr;
    };

    // the interfaces of every piece, cell by cell along its faces
    int lost = 0;
    for(int b=0; b<nb; ++b)
    for(const piece &P : pc[b])
    {
        blk_block &N = blk[P.id];
        const blk_block &O = ob[b];
        const int o[2] = {P.oi,P.oj};
        for(int f=0; f<4; ++f)
        {
            const int len = (f<2) ? P.nj : P.ni;
            segkey cur;
            int start = -1;
            for(int k=0; k<=len; ++k)
            {
                segkey key;
                bool has = false;
                if(k<len)
                {
                    // the cell beyond the face, in the indices of the block of DIVEMesh
                    int c[2];
                    c[0] = (f==0) ? -1 : (f==1) ? P.ni : k;
                    c[1] = (f==2) ? -1 : (f==3) ? P.nj : k;
                    c[0] += P.oi; c[1] += P.oj;

                    if(c[0]>=0 && c[0]<O.ni && c[1]>=0 && c[1]<O.nj)
                    {
                        // the next piece of the same block
                        const piece *Q = find(b,c[0],c[1]);
                        if(Q!=nullptr)
                        {
                            key.nbr = Q->id; key.nface = f^1;
                            key.off[0] = P.oi - Q->oi; key.off[1] = P.oj - Q->oj;
                            has = true;
                        }
                    }
                    else
                    {
                        // through the interface of DIVEMesh
                        const int ko = (f<2) ? c[1] : c[0];
                        for(const blk_seg &S : O.seg)
                        {
                            if(S.face!=f || ko<S.lo || ko>S.hi)
                            continue;
                            const int a = S.sign[0]*c[S.perm[0]] + S.off[0];
                            const int d = S.sign[1]*c[S.perm[1]] + S.off[1];
                            const piece *Q = find(S.nbr,a,d);
                            if(Q==nullptr)
                            {
                                ++lost;
                                break;
                            }
                            key.nbr = Q->id; key.nface = S.nface;
                            for(int q=0; q<2; ++q)
                            {
                                key.perm[q] = S.perm[q];
                                key.sign[q] = S.sign[q];
                            }
                            key.off[0] = S.off[0] + S.sign[0]*o[S.perm[0]] - Q->oi;
                            key.off[1] = S.off[1] + S.sign[1]*o[S.perm[1]] - Q->oj;
                            has = true;
                            break;
                        }
                    }
                }

                // the run of equal cells ends
                if(start>=0 && (!has || !(key==cur)))
                {
                    blk_seg s;
                    s.face = f; s.lo = start; s.hi = k-1;
                    s.nbr = cur.nbr; s.nface = cur.nface;
                    for(int q=0; q<2; ++q)
                    {
                        s.perm[q] = cur.perm[q];
                        s.sign[q] = cur.sign[q];
                        s.off[q] = cur.off[q];
                    }
                    N.seg.push_back(s);
                    start = -1;
                }
                if(has && start<0)
                {
                    start = k;
                    cur = key;
                }
            }
        }
    }

    if(lost>0)
    fail(std::to_string(lost)+" cells beyond the faces map outside their neighbour block");

    // check: every interface cell maps back to the cell it came from
    int bad = 0;
    for(int q=0; q<int(blk.size()); ++q)
    {
        const blk_block &B = blk[q];
        for(const blk_seg &S : B.seg)
        for(int k=S.lo; k<=S.hi; ++k)
        {
            int cb[2], ci[2];
            cb[0] = (S.face==0) ? -1 : (S.face==1) ? B.ni : k;
            cb[1] = (S.face==2) ? -1 : (S.face==3) ? B.nj : k;
            ci[0] = (S.face==0) ? 0 : (S.face==1) ? B.ni-1 : k;
            ci[1] = (S.face==2) ? 0 : (S.face==3) ? B.nj-1 : k;
            const int a = S.sign[0]*cb[S.perm[0]] + S.off[0];
            const int d = S.sign[1]*cb[S.perm[1]] + S.off[1];
            const blk_block &Q = blk[S.nbr];
            const bool at = (S.nface==0 && a==0) || (S.nface==1 && a==Q.ni-1) || (S.nface==2 && d==0) || (S.nface==3 && d==Q.nj-1);
            if(a<0 || a>=Q.ni || d<0 || d>=Q.nj || !at)
            {
                ++bad;
                continue;
            }
            int qb[2];
            qb[0] = (S.nface==0) ? -1 : (S.nface==1) ? Q.ni : a;
            qb[1] = (S.nface==2) ? -1 : (S.nface==3) ? Q.nj : d;
            const int kq = (S.nface<2) ? d : a;
            bool ok = false;
            for(const blk_seg &R : Q.seg)
            {
                if(R.face!=S.nface || kq<R.lo || kq>R.hi)
                continue;
                const int a2 = R.sign[0]*qb[R.perm[0]] + R.off[0];
                const int d2 = R.sign[1]*qb[R.perm[1]] + R.off[1];
                ok = (R.nbr==q && a2==ci[0] && d2==ci[1]);
                break;
            }
            if(!ok)
            ++bad;
        }
    }
    if(bad>0)
    fail(std::to_string(bad)+" interface cells of the pieces do not map back (grid-geometry.dat, BLKS)");

    if(p->mpirank==0)
    {
        if(wgt)
        cout<<"G 55 "<<m<<": the blocks 1.. cut into pieces of about "<<SW<<" load (G 56 "<<p->G56_a<<" "<<p->G56_z<<": load "
            <<wtot<<" of "<<tot<<" cells / "<<m*np<<"), "<<nb<<" blocks -> "<<blk.size()<<endl;
        else
        cout<<"G 55 "<<m<<": the blocks 1.. cut into pieces of at most "<<S<<" cells ("<<tot<<" cells / "<<m*np<<"), "
            <<nb<<" blocks -> "<<blk.size()<<endl;
        for(int b=1; b<nb; ++b)
        if(pc[b].size()>1)
        {
            cout<<"  block "<<b<<" ("<<ob[b].ni<<" x "<<ob[b].nj<<"): pieces";
            for(const piece &P : pc[b])
            {
            cout<<" "<<P.id<<" ("<<P.ni<<" x "<<P.nj;
            if(wgt)
            cout<<", load "<<blk_load(p,P.id);
            cout<<")";
            }
            cout<<endl;
        }
    }
}
