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
#include<string>

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
// The first and the last piece of a block reach at least the width of the relaxation and beach zones on the -i and
// +i face of the block (B 522), as far as the other pieces stay 16 cells deep (zone_cells).
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

// interface cells of the blocks that do not map back to the cell they came from (through the segment of the
// neighbour block at the face they arrive)
int geo_mesh::seg_check() const
{
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
    return bad;
}

// G 55 with the zones of B 522: the cells of block b (1..) within the width of its relaxation and beach zones from
// its face f (0 -i, 1 +i), the most over its rows (the distance from the centre of the face along the cell centres)
int geo_mesh::zone_cells(const lexer *p, int b, int f) const
{
    int nz = 0;
    const blk_block &B = blk[b];
    const size_t w = size_t(B.ni+1);
    auto ctr = [&](int ci, int cj, double &x, double &y)
    {
        x = 0.25*(B.x[ci+w*cj] + B.x[ci+1+w*cj] + B.x[ci+w*(cj+1)] + B.x[ci+1+w*(cj+1)]);
        y = 0.25*(B.y[ci+w*cj] + B.y[ci+1+w*cj] + B.y[ci+w*(cj+1)] + B.y[ci+1+w*(cj+1)]);
    };
    for(int m=0; m<p->B522; ++m)
    {
        if(p->B522_b[m]!=b)
        continue;
        int e=0, meth=0;
        double wz=0.0;
        for(int q=0; q<p->B521; ++q)
        if(p->B521_id[q]==p->B522_id[m])
        {
            e = p->B521_edge[q];
            wz = p->B521_w[q];
        }
        for(int q=0; q<p->B520; ++q)
        if(p->B520_id[q]==p->B522_id[m])
        meth = p->B520_method[q];
        if((meth!=1 && meth!=2) || e<1 || e>4 || face_turned(b,e-1)!=f)
        continue;

        for(int cj=0; cj<B.nj; ++cj)
        {
            const int K = (f==0) ? 0 : B.ni;
            double xo = 0.5*(B.x[K+w*cj] + B.x[K+w*(cj+1)]);
            double yo = 0.5*(B.y[K+w*cj] + B.y[K+w*(cj+1)]);
            double s = 0.0;
            for(int n=0; n<B.ni; ++n)
            {
                double xc,yc;
                ctr((f==0) ? n : B.ni-1-n,cj,xc,yc);
                s += sqrt((xc-xo)*(xc-xo) + (yc-yo)*(yc-yo));
                xo = xc;
                yo = yc;
                if(s>=wz)
                break;
                nz = std::max(nz,n+1);
            }
        }
    }
    return nz;
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
        const int kwant = std::max(1, wgt ? int(std::ceil(wb[b]/std::max(SW,1.0e-12) - 1.0e-9)) : int((c + S - 1)/S));

        // the cuts of a direction of length Ln into kn parts: at equal cells, with G 56 where the slices across
        // that direction reach equal parts of the load of the block; at least MIN cells deep
        auto cuts = [&](int dn, int Ln, int kn) -> vector<int>
        {
            vector<int> cut(kn+1,0);
            cut[kn] = Ln;
            for(int q=1; q<kn; ++q)
            cut[q] = int((long long)q*Ln/kn);
            if(wgt && kn>1)
            {
                vector<double> cum(Ln+1,0.0);
                const int Mn = dn ? B.ni : B.nj;
                for(int x=0; x<Ln; ++x)
                {
                    double w=0.0;
                    for(int y=0; y<Mn; ++y)
                    w += dn ? blk_cell_load(p,b,y,x) : blk_cell_load(p,b,x,y);
                    cum[x+1] = cum[x] + w;
                }
                for(int q=1; q<kn; ++q)
                {
                    const double target = cum[Ln]*double(q)/double(kn);
                    int x = 0;
                    while(x<Ln && cum[x+1]<=target)
                    ++x;
                    if(x<Ln && target-cum[x] > cum[x+1]-target)
                    ++x;
                    cut[q] = std::max(cut[q-1]+MIN, std::min(x, Ln-(kn-q)*MIN));
                }
            }
            // the relaxation and beach zones of the -i and +i face (B 522) on the first and the last piece whole,
            // as far as the pieces stay MIN cells deep
            if(dn==0 && kn>1)
            {
                const int zlo = zone_cells(p,b,0), zhi = zone_cells(p,b,1);
                if(zlo>cut[1])
                {
                    cut[1] = std::min(zlo, Ln-(kn-1)*MIN);
                    for(int q=2; q<kn; ++q)
                    cut[q] = std::max(cut[q], cut[q-1]+MIN);
                }
                if(Ln-zhi<cut[kn-1])
                cut[kn-1] = std::max(Ln-zhi, cut[kn-2]+MIN);
            }
            return cut;
        };

        // strips along dir; G 58 1 (blocks without open faces): tiles in both directions, about square in cells
        int ki, kj;
        if(p->G58==1 && !open(b))
        {
            const int mi = std::max(1,B.ni/MIN), mj = std::max(1,B.nj/MIN);
            ki = std::max(1, std::min(mi, int(std::lround(std::sqrt(double(kwant)*B.ni/std::max(B.nj,1))))));
            kj = std::max(1, std::min(mj, int(std::ceil(double(kwant)/ki - 1.0e-9))));
            ki = std::max(1, std::min(mi, int(std::ceil(double(kwant)/kj - 1.0e-9))));
        }
        else
        {
            const int k1 = std::max(1, std::min(kwant, std::max(1,L/MIN)));
            ki = dir ? 1 : k1;
            kj = dir ? k1 : 1;
        }
        const int k = ki*kj;
        const vector<int> ci = cuts(0,B.ni,ki), cj = cuts(1,B.nj,kj);

        int q = 0;
        for(int a=0; a<ki; ++a)
        for(int e=0; e<kj; ++e)
        {
            piece P;
            P.id = (q==0) ? b : next++;
            P.oi = ci[a];
            P.oj = cj[e];
            P.ni = ci[a+1]-ci[a];
            P.nj = cj[e+1]-cj[e];
            pc[b].push_back(P);
            ++q;
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
    const int bad = seg_check();
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

// ---------------------------------------------------------------------------------------------------------------
// G 54 and B 522 on any face of a block 1..: the blocks turned by quarter turns of their indices
//
// The open faces of a block are its -i face (G 54 discharge in, a relaxation zone, an edge or beach zone) and its
// +i face (G 54 water level out, an edge or beach zone): ioflow_amr_block and the boundary lists of the solvers
// know these two (a relaxation zone generates waves on the -i face as on the x- side of the native grid, B 98 2, and
// a Riemann edge with the waves of B 524 sources has their non-hydrostatic pressure there).  An open face on a j
// face, or a discharge / relaxation zone / wave edge at the +i or a level at the -i face, turns the block so that its
// open faces become the -i and the +i face.  A quarter turn takes cell (i,j) of a block of ni x nj cells to
// (j, ni-1-i) of the turned block of nj x ni cells: the faces -i, +i, -j, +j become +j, -j, -i, +i, the
// orientation of the indices is kept (no mirror), and the interface transforms of the block and of its
// neighbours compose with the turn.  The physical grid is the same; the cells are numbered the other way.
// Two open faces of a block must be opposite faces (-i/+i or -j/+j) with the discharge (or the relaxation zone)
// on the one and the level on the other, else the block cannot be turned.

namespace
{
const int turn_face[4] = {3,2,0,1};     // the face a quarter turn takes face 0..3 to
}

int geo_mesh::g54_kind(int t)
{
    return (t>=1 && t<=16) ? (t-1)%4 + 1 : 0;
}

// t 1-4 the -i (discharge) and +i (level) faces, 5-8 -j and +j, 9-12 +i and -i, 13-16 +j and -j
int geo_mesh::g54_face(int t)
{
    if(t<1 || t>16)
    return -1;
    const int g = (t-1)/4, k = (t-1)%4 + 1;
    const int fq[4] = {0,2,1,3};
    return (k==1 || k==3) ? fq[g] : (fq[g]^1);
}

int geo_mesh::face_turned(int b, int f) const
{
    if(f<0 || f>3)
    return f;
    const int r = (b>=0 && b<int(blk_turn.size())) ? blk_turn[b] : 0;
    for(int q=0; q<r; ++q)
    f = turn_face[f];
    return f;
}

void geo_mesh::turn_block(int b)
{
    blk_block &B = blk[b];
    const int ni = B.ni, nj = B.nj;

    // the nodes: (K,L) -> (L, ni-K), the turned block has nj+1 nodes along its i
    {
    vector<double> x(B.x.size()), y(B.y.size()), z(B.zb.size());
    for(int L=0; L<=nj; ++L)
    for(int K=0; K<=ni; ++K)
    {
        const size_t s = size_t(K) + size_t(ni+1)*size_t(L);
        const size_t d = size_t(L) + size_t(nj+1)*size_t(ni-K);
        x[d] = B.x[s]; y[d] = B.y[s]; z[d] = B.zb[s];
    }
    B.x.swap(x); B.y.swap(y); B.zb.swap(z);
    }

    // the segments of the block: the face, the range along it (the cells beyond the face turned), the
    // transform from the turned indices (c = (ni-1-c'[1], c'[0]))
    for(blk_seg &S : B.seg)
    {
        int lo2=0, hi2=0;
        for(int e=0; e<2; ++e)
        {
            const int k = (e==0) ? S.lo : S.hi;
            const int c0 = (S.face==0) ? -1 : (S.face==1) ? ni : k;
            const int c1 = (S.face==2) ? -1 : (S.face==3) ? nj : k;
            const int n0 = c1, n1 = ni-1-c0;
            const int kn = (turn_face[S.face]<2) ? n1 : n0;
            if(e==0)
            lo2 = hi2 = kn;
            else
            {
            lo2 = std::min(lo2,kn);
            hi2 = std::max(hi2,kn);
            }
        }
        S.face = turn_face[S.face];
        S.lo = lo2;
        S.hi = hi2;
        for(int q=0; q<2; ++q)
        {
            if(S.perm[q]==0)
            {
                S.off[q] += S.sign[q]*(ni-1);
                S.sign[q] = -S.sign[q];
                S.perm[q] = 1;
            }
            else
            S.perm[q] = 0;
        }
    }

    // the segments into the block (its own included): its cells turned, (a,d) -> (d, ni-1-a)
    for(blk_block &Q : blk)
    for(blk_seg &S : Q.seg)
    if(S.nbr==b)
    {
        const int p0 = S.perm[0], s0 = S.sign[0], o0 = S.off[0];
        S.perm[0] = S.perm[1];
        S.sign[0] = S.sign[1];
        S.off[0] = S.off[1];
        S.perm[1] = p0;
        S.sign[1] = -s0;
        S.off[1] = ni-1-o0;
        S.nface = turn_face[S.nface];
    }

    B.ni = nj;
    B.nj = ni;
}

void geo_mesh::turn_blocks(lexer *p)
{
    const int nb = int(blk.size());
    blk_turn.assign(nb,0);
    if(nb<=1)
    return;

    // the open faces of every block of DIVEMesh and the faces they may become (bit 0: -i, bit 1: +i, bit 2: -j,
    // bit 3: +j; jwant: in the second search, when no turn puts all of them on -i / +i)
    vector<vector<int>> face(nb), want(nb), jwant(nb);
    vector<vector<string>> what(nb);
    for(int n=0; n<p->G54; ++n)
    {
        const int b = p->G54_b[n], f = g54_face(p->G54_t[n]), k = g54_kind(p->G54_t[n]);
        if(b<1 || b>=nb || f<0)
        continue;
        face[b].push_back(f);
        want[b].push_back((k==1 || k==3) ? 1 : 2);
        jwant[b].push_back(want[b].back());
        what[b].push_back((k==1 || k==3) ? "G 54 discharge" : "G 54 water level");
    }
    for(int m=0; m<p->B522; ++m)
    {
        const int b = p->B522_b[m];
        int e = 0, meth = 0;
        for(int q=0; q<p->B521; ++q)
        if(p->B521_id[q]==p->B522_id[m])
        e = p->B521_edge[q];
        for(int q=0; q<p->B520; ++q)
        if(p->B520_id[q]==p->B522_id[m])
        meth = p->B520_method[q];
        if(b<1 || b>=nb || e<1 || e>4)
        continue;
        face[b].push_back(e-1);
        // a relaxation zone (wave generation) and a Riemann edge with the waves of sources (B 524) on the -i face, as
        // the generation zone and the wave edge of the native grid at x- (the non-hydrostatic pressure of the waves)
        bool src = false;
        for(int q=0; q<p->B524; ++q)
        if(p->B524_id[q]==p->B522_id[m])
        src = true;
        want[b].push_back((meth==1 || (meth==3 && src)) ? 1 : 3);
        // an edge (B 520 methods 3-6) without waves and a beach zone may also lie on a j face (-j: side 3, +j: side
        // 2 of the ghost cell lists, as the y edges of the native grid)
        jwant[b].push_back((meth==1 || (meth==3 && src)) ? 1 : 15);
        what[b].push_back(meth==1 ? "B 522 relaxation zone" : (meth==3 && src) ? "B 522 Riemann edge with waves" : "B 522 zone");
    }

    int bad = 0;
    for(int b=1; b<nb; ++b)
    {
        if(face[b].empty())
        continue;

        // the first turn that puts all open faces on -i / +i; else the first one with the edges and beach zones on any
        // face (two open faces that are not opposite: one of them on a j face)
        int r = -1;
        for(int pass=0; pass<2 && r<0; ++pass)
        for(int t=0; t<4 && r<0; ++t)
        {
            bool ok = true;
            for(size_t q=0; q<face[b].size(); ++q)
            {
                int f = face[b][q];
                for(int u=0; u<t; ++u)
                f = turn_face[f];
                const int w = (pass==0) ? want[b][q] : jwant[b][q];
                if(!((w>>f)&1))
                ok = false;
            }
            // the second search: one open face per face
            if(pass==1)
            for(size_t q=0; q<face[b].size() && ok; ++q)
            for(size_t u=q+1; u<face[b].size(); ++u)
            if(face[b][q]==face[b][u])
            ok = false;
            if(ok)
            r = t;
        }

        if(r<0)
        {
            if(p->mpirank==0)
            {
                const char *fn[4] = {"-i","+i","-j","+j"};
                cout<<"!!! block "<<b<<": its open faces (G 54, B 522)";
                for(size_t q=0; q<face[b].size(); ++q)
                cout<<(q>0 ? ", " : " ")<<fn[face[b][q]]<<" ("<<what[b][q]<<")";
                cout<<": the block cannot be turned to take them (a discharge, a relaxation zone or a Riemann edge with waves on its -i face, a water level on its +i face; edges without waves and beach zones on any face; one open face per face) !!!"<<endl;
            }
            ++bad;
            continue;
        }

        blk_turn[b] = r;
        for(int t=0; t<r; ++t)
        turn_block(b);
    }
    if(bad>0)
    exit(1);

    bool any = false;
    for(int b=1; b<nb; ++b)
    if(blk_turn[b]>0)
    any = true;
    if(!any)
    return;

    // the sizes of the blocks of DIVEMesh as REEF3D has them (G 55, face_piece)
    for(int b=1; b<nb; ++b)
    {
        orig_ni[b] = blk[b].ni;
        orig_nj[b] = blk[b].nj;
    }

    const int nbad = seg_check();
    if(nbad>0)
    {
        if(p->mpirank==0)
        cout<<"!!! turning the blocks with open faces: "<<nbad<<" interface cells do not map back !!!"<<endl;
        exit(1);
    }

    if(p->mpirank==0)
    for(int b=1; b<nb; ++b)
    if(blk_turn[b]>0)
    cout<<"block "<<b<<": turned by "<<90*blk_turn[b]<<" deg in its indices (now "<<blk[b].ni<<" x "<<blk[b].nj
        <<" cells), so that its open faces lie where the solvers take them"<<endl;
}
