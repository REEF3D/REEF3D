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

#include"reefamr.h"
#include"lexer.h"
#include"ghostcell.h"
#include"geo_mesh.h"
#include"ioflow_amr_block.h"
#include<algorithm>
#include<numeric>

// Multi-block level 0 (G 50 1 with the BLKS section of DIVEMesh R 20 1), see reefamr.h.
//
// Interfaces are point matching: a segment of a face of block b maps the cells next to it on both
// sides with a transform (8 orientations); on level l the transform acts on the 2^l x 2^l children
// of the level-0 cells.  A cell beyond a face is found in the neighbour block (two hops at most:
// corners of four blocks); cells no interface reaches are outside the domain (walls).
// Velocities are Cartesian components, so the filled values need no rotation.

namespace
{
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }
// floor(a/2^l) for negative a as well
inline int fsh(int a, int l) { return (a>=0) ? (a>>l) : -(((-a)-1)>>l)-1; }
}

void reefamr::block_own(int l, int b, int I, int J, int &i, int &j) const
{
    const rblock &B = blk[b];
    if(b==0)
    {
        i = I - (B.ring[0]<<l);
        j = J - (B.ring[2]<<l);
        return;
    }
    i = I - (B.vi<<l);
    j = J - (B.vj<<l);
}

void reefamr::block_virt(int l, int b, int i, int j, int &I, int &J) const
{
    const rblock &B = blk[b];
    if(b==0)
    {
        I = i + (B.ring[0]<<l);
        J = j + (B.ring[2]<<l);
        return;
    }
    I = i + (B.vi<<l);
    J = j + (B.vj<<l);
}

int reefamr::block_near(int I, int J) const
{
    if(nblk==0)
    return -1;
    if(I>=0 && I<GNX && J>=0 && J<GNY)
    return 0;
    const int M = BGAP/2;
    for(int b=1; b<nblk; ++b)
    {
        const rblock &B = blk[b];
        if(I>=B.vi-M && I<B.vi+B.ni+M && J>=B.vj-M && J<B.vj+B.nj+M)
        return b;
    }
    return -1;
}

bool reefamr::native_ring(int I, int J) const
{
    if(nblk==0)
    return false;
    int i,j;
    block_own(0,0,I,J,i,j);
    return i<0 || i>=blk[0].ni || j<0 || j>=blk[0].nj;
}

// cell (i,j) of block b on level l (own indices, possibly beyond the faces) -> the block and own
// cell that hold it
bool reefamr::block_resolve(int l, int b, int i, int j, int &b2, int &i2, int &j2, int hops) const
{
    const rblock &B = blk[b];
    const int s = 1<<l;
    const int ni = B.ni*s, nj = B.nj*s;

    if(i>=0 && i<ni && j>=0 && j<nj)
    {
        b2=b; i2=i; j2=j;
        return true;
    }

    if(hops>=2)
    return false;

    int fc[2], nf=0;
    if(i<0) fc[nf++]=0;
    if(i>=ni) fc[nf++]=1;
    if(j<0) fc[nf++]=2;
    if(j>=nj) fc[nf++]=3;

    for(int q=0; q<nf; ++q)
    {
        const int f = fc[q];
        // position along the face, clamped (a corner continues the segment at the end)
        const int k = (f<2) ? std::min(std::max(j,0),nj-1) : std::min(std::max(i,0),ni-1);
        for(const rseg &S : B.seg)
        {
            if(S.face!=f || k<S.lo*s || k>(S.hi+1)*s-1)
            continue;
            const int c[2] = {i,j};
            const int a = S.sign[0]*c[S.perm[0]] + S.off[0]*s + (S.sign[0]<0 ? s-1 : 0);
            const int d = S.sign[1]*c[S.perm[1]] + S.off[1]*s + (S.sign[1]<0 ? s-1 : 0);
            if(block_resolve(l,S.nbr,a,d,b2,i2,j2,hops+1))
            return true;
        }
    }
    return false;
}

// cell (i,j) of block b beyond its face f, through that face first (corners: the other face is
// crossed in the neighbour block)
bool reefamr::block_cross(int l, int b, int f, int i, int j, int &b2, int &i2, int &j2) const
{
    const rblock &B = blk[b];
    const int s = 1<<l;
    const int ni = B.ni*s, nj = B.nj*s;
    const int k = (f<2) ? std::min(std::max(j,0),nj-1) : std::min(std::max(i,0),ni-1);
    for(const rseg &S : B.seg)
    {
        if(S.face!=f || k<S.lo*s || k>(S.hi+1)*s-1)
        continue;
        const int c[2] = {i,j};
        const int a = S.sign[0]*c[S.perm[0]] + S.off[0]*s + (S.sign[0]<0 ? s-1 : 0);
        const int d = S.sign[1]*c[S.perm[1]] + S.off[1]*s + (S.sign[1]<0 ? s-1 : 0);
        return block_resolve(l,S.nbr,a,d,b2,i2,j2,1);
    }
    return false;
}

// flux matching: patch of block b, coarse level L, patch side, outside coarse cell (Ic,Jc) and inside
// cell (Ii,Ji) (index space of block b).  Beyond the block face the outside cell becomes the cell of
// the neighbour block, side2 the side of the coarse face seen there (1/3: the inside cell lower) and
// sgn -1 where the face normal is turned against the fine faces.  false: beyond a wall
bool reefamr::block_face_map(int L, int b, int side, int &Ic, int &Jc, int Ii, int Ji, int &side2, int &sgn) const
{
    side2 = side;
    sgn = 1;
    const rblock &B = blk[b];
    const int s = 1<<L;
    const int ni = B.ni*s, nj = B.nj*s;
    int io,jo,ii,ji;
    block_own(L,b,Ic,Jc,io,jo);
    if(io>=0 && io<ni && jo>=0 && jo<nj)
    return true;
    block_own(L,b,Ii,Ji,ii,ji);

    const int k = (side<2) ? jo : io;
    for(const rseg &S : B.seg)
    {
        if(S.face!=side || k<S.lo*s || k>(S.hi+1)*s-1)
        continue;
        auto map = [&](int i, int j, int &a, int &d)
        {
            const int c[2] = {i,j};
            a = S.sign[0]*c[S.perm[0]] + S.off[0]*s + (S.sign[0]<0 ? s-1 : 0);
            d = S.sign[1]*c[S.perm[1]] + S.off[1]*s + (S.sign[1]<0 ? s-1 : 0);
        };
        int ao,dd,ai,di;
        map(io,jo,ao,dd);
        map(ii,ji,ai,di);
        const rblock &N = blk[S.nbr];
        if(ao<0 || ao>=N.ni*s || dd<0 || dd>=N.nj*s)
        return false;
        const int dir2 = (ao!=ai) ? 0 : 1;
        const bool lower = (dir2==0) ? (ai<ao) : (di<dd);
        side2 = 2*dir2 + (lower ? 1 : 0);
        sgn = ((side&1)==(side2&1)) ? 1 : -1;
        block_virt(L,S.nbr,ao,dd,Ic,Jc);
        return true;
    }
    return false;
}

// own cells of block b on level l in the index space (inclusive)
void reefamr::block_region(int l, int b, int &X0, int &X1, int &Y0, int &Y1) const
{
    const rblock &B = blk[b];
    const int a = (b==0) ? B.ring[0] : B.vi;
    const int d = (b==0) ? B.ring[2] : B.vj;
    X0 = a<<l; X1 = ((a+B.ni)<<l)-1;
    Y0 = d<<l; Y1 = ((d+B.nj)<<l)-1;
}

int reefamr::block_at(int l, int I, int J) const
{
    if(nblk==0)
    return 0;
    const int b = block_near(fsh(I,l),fsh(J,l));
    if(b<0)
    return -1;
    int X0,X1,Y0,Y1;
    block_region(l,b,X0,X1,Y0,Y1);
    return (I>=X0 && I<=X1 && J>=Y0 && J<=Y1) ? b : -1;
}

// the index of the cell an index around a block stands for (one block: identity); false: beyond a wall
bool reefamr::vcanon(int l, int I, int J, int &I2, int &J2) const
{
    I2 = I;
    J2 = J;
    if(nblk==0)
    return true;
    const int b = block_near(fsh(I,l),fsh(J,l));
    if(b<0)
    return false;
    return canon(l,b,I,J,I2,J2);
}

bool reefamr::canon(int l, int b, int I, int J, int &I2, int &J2) const
{
    int i,j,b2,i2,j2;
    block_own(l,b,I,J,i,j);
    if(!block_resolve(l,b,i,j,b2,i2,j2))
    return false;
    block_virt(l,b2,i2,j2,I2,J2);
    return true;
}

// level-0 node (K,L) of block b (own node indices), beyond the faces through the interfaces
bool reefamr::block_node(int b, int K, int L, double &x, double &y, double &z, int hops) const
{
    const rblock &B = blk[b];

    if(K>=0 && K<=B.ni && L>=0 && L<=B.nj)
    {
        if(b==0)
        {
            p0->curv_node(K+B.ring[0],L+B.ring[2],x,y,z);
            return true;
        }
        const geo_mesh::blk_block &G = p0->gridgeo->blk[b];
        const size_t a = K + size_t(B.ni+1)*L;
        x = G.x[a]; y = G.y[a]; z = G.zb[a];
        return true;
    }

    if(hops>=2)
    return false;

    int fc[2], nf=0;
    if(K<0) fc[nf++]=0;
    if(K>B.ni) fc[nf++]=1;
    if(L<0) fc[nf++]=2;
    if(L>B.nj) fc[nf++]=3;

    for(int q=0; q<nf; ++q)
    {
        const int f = fc[q];
        const int k = (f<2) ? std::min(std::max(L,0),B.nj) : std::min(std::max(K,0),B.ni);
        for(const rseg &S : B.seg)
        {
            if(S.face!=f || k<S.lo || k>S.hi+1)
            continue;
            const int N[2] = {K,L};
            const int K2 = S.sign[0]*N[S.perm[0]] + S.off[0] + (S.sign[0]<0 ? 1 : 0);
            const int L2 = S.sign[1]*N[S.perm[1]] + S.off[1] + (S.sign[1]<0 ? 1 : 0);
            if(block_node(S.nbr,K2,L2,x,y,z,hops+1))
            return true;
        }
    }
    return false;
}

// as block_node, and beyond the walls by point reflection at the boundary nodes (as lexer::curv_node);
// a node beyond a wall and an interface goes through the interface first and is reflected at the wall of
// the neighbour block, where the wall continues (as the one-block grid reflects it)
void reefamr::block_node_any(int b, int K, int L, double &x, double &y, double &z) const
{
    if(block_node(b,K,L,x,y,z))
    return;

    const rblock &B = blk[b];
    int fc[2], nf=0;
    if(K<0) fc[nf++]=0;
    if(K>B.ni) fc[nf++]=1;
    if(L<0) fc[nf++]=2;
    if(L>B.nj) fc[nf++]=3;

    if(nf==2)
    for(int q=0; q<nf; ++q)
    {
        const int f = fc[q];
        const int k = (f<2) ? std::min(std::max(L,0),B.nj) : std::min(std::max(K,0),B.ni);
        for(const rseg &S : B.seg)
        {
            if(S.face!=f || k<S.lo || k>S.hi+1)
            continue;
            const int N[2] = {K,L};
            const int K2 = S.sign[0]*N[S.perm[0]] + S.off[0] + (S.sign[0]<0 ? 1 : 0);
            const int L2 = S.sign[1]*N[S.perm[1]] + S.off[1] + (S.sign[1]<0 ? 1 : 0);
            const rblock &M = blk[S.nbr];
            // beyond one face of the neighbour only: reflected there
            const int out = (K2<0 || K2>M.ni ? 1 : 0) + (L2<0 || L2>M.nj ? 1 : 0);
            if(out==1)
            {
                block_node_reflect(S.nbr,K2,L2,x,y,z);
                return;
            }
        }
    }

    block_node_reflect(b,K,L,x,y,z);
}

void reefamr::block_node_reflect(int b, int K, int L, double &x, double &y, double &z) const
{
    if(block_node(b,K,L,x,y,z))
    return;

    const rblock &B = blk[b];
    auto get = [&](int k, int q, double &xx, double &yy, double &zz)
    {
        if(block_node(b,k,q,xx,yy,zz))
        return;
        block_node(b,std::min(std::max(k,0),B.ni),std::min(std::max(q,0),B.nj),xx,yy,zz);
    };

    int Kb=K, Kr=K, Lb=L, Lr=L;
    if(K<0)    { Kb=0;    Kr=-K; }
    if(K>B.ni) { Kb=B.ni; Kr=2*B.ni-K; }
    if(L<0)    { Lb=0;    Lr=-L; }
    if(L>B.nj) { Lb=B.nj; Lr=2*B.nj-L; }

    double xr,yr,zr,xb,yb,zb;
    get(Kr,Lr,xr,yr,zr);
    z = zr;

    // reflect in i at the boundary node (Kb,Lr), then in j at (K,Lb)
    if(Kb!=K)
    {
        get(Kb,Lr,xb,yb,zb);
        xr = 2.0*xb - xr;
        yr = 2.0*yb - yr;
        z = zb;
    }
    if(Lb!=L)
    {
        double xc,yc,zc,xd,yd,zd;
        get(K,Lb,xc,yc,zc);
        if(Kb!=K)
        {
            get(Kb,Lb,xc,yc,zc);
            get(Kr,Lb,xd,yd,zd);
            xc = 2.0*xc - xd;
            yc = 2.0*yc - yd;
        }
        xr = 2.0*xc - xr;
        yr = 2.0*yc - yr;
        z = zc;
    }
    x = xr;
    y = yr;
}

// blocks from the geometry file, their place in the index space, their ranks and patches
void reefamr::blocks_setup(lexer *p, ghostcell *pgc)
{
    const geo_mesh *G = p->gridgeo;
    nblk = int(G->blk.size());
    blk.assign(nblk,rblock());

    for(int b=0; b<nblk; ++b)
    {
        const geo_mesh::blk_block &S = G->blk[b];
        rblock &B = blk[b];
        B.ni = S.ni; B.nj = S.nj;
        for(int f=0; f<4; ++f)
        B.ring[f] = S.ring[f];
        for(const geo_mesh::blk_seg &s : S.seg)
        {
            rseg r;
            r.face=s.face; r.lo=s.lo; r.hi=s.hi; r.nbr=s.nbr; r.nface=s.nface;
            for(int c=0; c<2; ++c)
            {
                r.perm[c]=s.perm[c]; r.sign[c]=s.sign[c]; r.off[c]=s.off[c];
            }
            B.seg.push_back(r);
        }
    }

    blocks_open(p);

    // block 0 with its ring is the native grid
    const rblock &B0 = blk[0];
    if(B0.ni+B0.ring[0]+B0.ring[1]!=GNX || B0.nj+B0.ring[2]+B0.ring[3]!=GNY)
    {
        if(p->mpirank==0)
        cout<<"!!! "<<par.name<<": block 0 with its ring has "<<B0.ni+B0.ring[0]+B0.ring[1]<<" x "<<B0.nj+B0.ring[2]+B0.ring[3]
            <<" cells, the grid "<<GNX<<" x "<<GNY<<" !!!"<<endl;
        exit(1);
    }

    // the cells beyond a face reach EXT + margin into the neighbour (fills), the ring of block 0 has
    // to hold the cells the native grid computes redundantly (EXT) and its ghost cells (margin)
    const int need = EXT + p->margin;
    int bad=0;
    for(int b=0; b<nblk; ++b)
    for(const rseg &S : blk[b].seg)
    {
        const rblock &N = blk[S.nbr];
        const int depth = (S.nface<2) ? N.ni : N.nj;
        if(depth<need)
        {
            if(p->mpirank==0)
            cout<<"!!! "<<par.name<<": block "<<S.nbr<<" is "<<depth<<" cells deep behind its face "<<S.nface
                <<", the cells beyond the faces reach "<<need<<" (EXT + margin) !!!"<<endl;
            ++bad;
        }
        if(b==0 && B0.ring[S.face]<need)
        {
            if(p->mpirank==0)
            cout<<"!!! "<<par.name<<": the ring of block 0 has "<<B0.ring[S.face]<<" cells, at least "<<need
                <<" needed (DIVEMesh R 21) !!!"<<endl;
            ++bad;
        }
    }
    if(bad>0)
    exit(1);

    // index space: the blocks to the right of the native box, with gaps on all sides, their first cell
    // on a multiple of 8 (tiles of 8 cells start at the block corner on every level)
    int I = GNX + BGAP;
    GNYV = GNY;
    for(int b=1; b<nblk; ++b)
    {
        I = ((I+7)/8)*8;
        blk[b].vi = I;
        blk[b].vj = BGAP;
        I += blk[b].ni + BGAP;
        GNYV = std::max(GNYV,blk[b].nj+2*BGAP);
    }
    GNXV = I;

    // ranks: the largest block first onto the rank with the least cells (native box + blocks)
    const int np = p->mpi_size;
    vector<long> load(np,0);
    for(int r=0; r<np; ++r)
    load[r] = long(rbx1[r]-rbx0[r]+1)*long(rby1[r]-rby0[r]+1);
    vector<int> ord;
    for(int b=1; b<nblk; ++b)
    ord.push_back(b);
    std::stable_sort(ord.begin(),ord.end(),[&](int a, int c){ return long(blk[a].ni)*blk[a].nj > long(blk[c].ni)*blk[c].nj; });
    for(int b : ord)
    {
        int r = int(std::min_element(load.begin(),load.end())-load.begin());
        blk[b].rank = r;
        load[r] += long(blk[b].ni)*blk[b].nj;
    }

    // their patches (level 0), in block order on every rank: the local index is known everywhere
    vector<int> nloc(np,0);
    for(int b=1; b<nblk; ++b)
    {
        rblock &B = blk[b];
        B.lid = nloc[B.rank]++;
        if(B.rank!=p->mpirank)
        continue;
        reefamr_patch *c = make_patch(p,pgc,0,B.vi,B.vi+B.ni-1,B.vj,B.vj+B.nj-1,b);
        c->fresh = false;
        P.push_back(c);
        lev[0].push_back(int(P.size())-1);
    }

    // the ring of the native grid on this rank
    nring.assign(p->imax*p->jmax,0);
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    if(native_ring(ii+O0i,jj+O0j))
    nring[lij(p,ii,jj)] = 1;

    // refinement over the index space with the blocks
    blocks_refine_setup(p);

    if(p->mpirank==0)
    {
        cout<<par.name<<": multi-block grid, "<<nblk<<" blocks (block 0 the native grid, ring "<<B0.ring[0]<<" "<<B0.ring[1]
            <<" "<<B0.ring[2]<<" "<<B0.ring[3]<<")";
        for(int b=1; b<nblk; ++b)
        {
            cout<<"; block "<<b<<": "<<blk[b].ni<<" x "<<blk[b].nj;
            if(G->blk[b].orig!=b || G->nblk_orig<nblk)
            cout<<" (of block "<<G->blk[b].orig<<" at "<<G->blk[b].oi<<" "<<G->blk[b].oj<<")";
            cout<<" on rank "<<blk[b].rank;
        }
        cout<<endl;
        if(maxlev>=1 && nsing>0)
        cout<<par.name<<": "<<nsing<<" block corner(s) at singular vertices (not four blocks around the node), no refinement around them"<<endl;
    }
}

// tile maps and no-refinement cells over the index space with the blocks; the singular vertices; the
// open faces of the blocks 1.. (G 54 b t v): t 1 discharge v through the -i face of block b, t 2 water level v at
// its +i face, t 3 and 4 the same from the time series file discharge-<v>.dat or waterlevel-<v>.dat (lines
// "t value", ioflow_amr_block); the face has no interface, one line per face.  b is the block of DIVEMesh: with
// G 55 the face belongs to the piece that holds it (geo_mesh::face_piece)
void reefamr::blocks_open(lexer *p)
{
    bopen[0].assign(nblk,-1);
    bopen[1].assign(nblk,-1);
    const geo_mesh *G = p->gridgeo;
    int bad=0;
    for(int n=0; n<p->G54; ++n)
    {
        const int bo = p->G54_b[n], t = p->G54_t[n];
        const int f = (t==1 || t==3) ? 0 : 1;
        const int b = G->face_piece(bo,f);
        const char *why = nullptr;
        string file;
        if(bo<1 || bo>=G->nblk_orig)
        why = "no block 1.. of the grid";
        else
        if(t<1 || t>4)
        why = "type 1 (discharge in, -i face), 2 (water level out, +i face), 3 (discharge time series in) or 4 (water level time series out)";
        else
        if(b<0)
        why = "the face is cut into pieces (G 55)";
        else
        if(bopen[f][b]>=0)
        why = "the face is open already";
        else
        for(const rseg &S : blk[b].seg)
        if(S.face==f)
        why = "the face has an interface";

        // time series: every rank reads the file (the block's rank uses it)
        vector<double> st, sv;
        if(why==nullptr && t>=3)
        {
            file = ioflow_amr_block::series_file(t,int(p->G54_v[n]));
            if(p->G54_v[n]<0.0 || p->G54_v[n]!=double(int(p->G54_v[n])))
            why = "the id of the time series file is a whole number >= 0";
            else
            if(!ioflow_amr_block::series_read(file,st,sv))
            why = "the time series file is missing, has no line \"t value\" or its times do not increase";
        }

        if(why!=nullptr)
        {
            if(p->mpirank==0)
            cout<<"!!! "<<par.name<<": G 54 "<<bo<<" "<<t<<" "<<p->G54_v[n]<<(file.empty() ? "" : " ("+file+")")<<": "<<why<<" !!!"<<endl;
            ++bad;
            continue;
        }
        bopen[f][b] = n;
        if(p->mpirank==0)
        {
            cout<<par.name<<": block "<<bo;
            if(b!=bo)
            cout<<" (piece "<<b<<")";
            cout<<(f==0?" inflow ":" outflow ");
            if(t<=2)
            cout<<(f==0?"Q = ":"water level ")<<p->G54_v[n]<<(f==0?" m3/s":" m");
            else
            cout<<(f==0?"Q":"water level")<<" from "<<file<<" ("<<st.size()<<" values, t "<<st.front()<<" .. "<<st.back()<<" s)";
            cout<<(f==0?" through its -i face":" at its +i face")<<endl;
        }
    }

    // edge zones on the faces of the blocks (B 520 methods 3-6 with B 522 b, NHFLOW): the -i or +i face of a
    // block 1.. without interface and without G 54, the background of the zone given (B 510)
    bzone[0].assign(nblk,-1);
    bzone[1].assign(nblk,-1);
    bzones = bc_zone_set::from_input(p).bedges;
    for(int q=0; q<(int)bzones.size(); ++q)
    {
        const bc_zone &z = bzones[q];
        const int bo = z.block, f = z.edge-1;
        const int b = G->face_piece(bo,f);
        const char *why = nullptr;
        if(bo>=G->nblk_orig)
        why = "no block 1.. of the grid";
        else
        if(p->A10!=5)
        why = "edge zones on the blocks are available for NHFLOW (A 10 5) only, so far";
        else
        if(b<0)
        why = "the face is cut into pieces (G 55)";
        else
        if(bopen[f][b]>=0)
        why = "the face is open by G 54 already";
        else
        for(const rseg &S : blk[b].seg)
        if(S.face==f)
        why = "the face has an interface";

        if(why==nullptr && z.bg>0)
        {
            bool found = false;
            for(int m=0; m<p->B510; ++m)
            if(p->B510_id[m]==z.bg)
            found = true;
            if(!found)
            why = "its background (B 523) has no B 510";
        }

        if(why!=nullptr)
        {
            if(p->mpirank==0)
            cout<<"!!! "<<par.name<<": zone "<<z.id<<" on face "<<z.edge<<" of block "<<bo<<" (B 522): "<<why<<" !!!"<<endl;
            ++bad;
            continue;
        }
        bzone[f][b] = q;
        if(p->mpirank==0)
        {
            const char *m = z.method==bc_method::riemann ? "Riemann edge" : z.method==bc_method::flather ? "Flather edge"
                          : z.method==bc_method::clamp_level ? "clamped level edge" : "clamped discharge edge";
            cout<<par.name<<": block "<<bo;
            if(b!=bo)
            cout<<" (piece "<<b<<")";
            cout<<": zone "<<z.id<<", "<<m;
            if(z.bg>0)
            cout<<" with background "<<z.bg;
            if(z.has_Q)
            cout<<", Q = "<<z.Q<<" m3/s (ramp "<<z.Q_tramp<<" s)";
            cout<<" on its "<<(f==0 ? "-i" : "+i")<<" face"<<endl;
        }
    }
    if(bad>0)
    exit(1);

    // the backgrounds of the edges, read on every rank (a missing or broken time series stops all ranks here)
    if(!bzones.empty())
    bzbg.read(p);
}

// tiles around every block and the tiles of the cells they stand for (refinement across interfaces)
void reefamr::blocks_refine_setup(lexer *p)
{
    for(int l=1; l<=maxlev; ++l)
    {
        gtnx[l] = ((GNXV<<l)+tile-1)/tile;
        gtny[l] = ((GNYV<<l)+tile-1)/tile;
    }

    // no refinement: the native cells as set up, the G 11 boxes over the blocks (cell centres)
    vector<unsigned char> F((size_t)GNXV*GNYV,0);
    for(int I=0; I<GNX; ++I)
    for(int J=0; J<GNY; ++J)
    F[(size_t)I*GNYV+J] = forbid0[(size_t)I*GNY+J];

    auto inbox = [&](double x, double y)
    {
        for(size_t q=0; q+3<par.fbox.size(); q+=4)
        if(x>=par.fbox[q] && x<=par.fbox[q+1] && y>=par.fbox[q+2] && y<=par.fbox[q+3])
        return true;
        return false;
    };

    if(!par.fbox.empty())
    {
        for(int b=1; b<nblk; ++b)
        for(int i=0; i<blk[b].ni; ++i)
        for(int j=0; j<blk[b].nj; ++j)
        {
            double xc=0.0, yc=0.0, x,y,z;
            for(int a=0; a<2; ++a)
            for(int d=0; d<2; ++d)
            {
                block_node(b,i+a,j+d,x,y,z);
                xc += 0.25*x;
                yc += 0.25*y;
            }
            if(inbox(xc,yc))
            F[(size_t)(blk[b].vi+i)*GNYV+(blk[b].vj+j)] = 1;
        }
    }

    // no refinement next to the open faces (G 54, the edge zones B 522), as next to the in- and outflow of the native grid
    for(int b=1; b<nblk; ++b)
    for(int t=0; t<2; ++t)
    if(bopen[t][b]>=0 || bzone[t][b]>=0)
    for(int i=0; i<=std::min(par.ioband,blk[b].ni-1); ++i)
    for(int j=0; j<blk[b].nj; ++j)
    F[(size_t)(blk[b].vi+(t==0 ? i : blk[b].ni-1-i))*GNYV+(blk[b].vj+j)] = 1;

    // singular vertices: the corner cell beyond both faces is another cell through the one face
    // than through the other (three or five blocks around the node), or only one way reaches it
    const int R = std::max(4,(EXT+p->margin+1)/2+nest+1);
    nsing = 0;
    for(int b=0; b<nblk; ++b)
    {
        const rblock &B = blk[b];
        for(int cx=0; cx<2; ++cx)
        for(int cy=0; cy<2; ++cy)
        {
            const int ci = cx ? B.ni : -1, cj = cy ? B.nj : -1;
            int bA,iA,jA,bB,iB,jB;
            const bool okA = block_cross(0,b,cx ? 1 : 0,ci,cj,bA,iA,jA);
            const bool okB = block_cross(0,b,cy ? 3 : 2,ci,cj,bB,iB,jB);
            const bool sing = (okA!=okB) || (okA && okB && (bA!=bB || iA!=iB || jA!=jB));
            if(!sing)
            continue;
            ++nsing;
            const int i0 = cx ? std::max(B.ni-R,0) : 0, i1 = cx ? B.ni-1 : std::min(R,B.ni)-1;
            const int j0 = cy ? std::max(B.nj-R,0) : 0, j1 = cy ? B.nj-1 : std::min(R,B.nj)-1;
            for(int i=i0; i<=i1; ++i)
            for(int j=j0; j<=j1; ++j)
            {
                int I,J;
                block_virt(0,b,i,j,I,J);
                F[(size_t)I*GNYV+J] = 1;
            }
        }
    }

    global_or(F);
    forbid0.swap(F);

    // tiles around the blocks (in the ring of block 0, in the gaps of the others) and the tiles of
    // the cells they stand for
    hmap.assign(maxlev+1,vector<int>());
    for(int l=1; l<=maxlev; ++l)
    {
        const int T = tile;
        const int H = std::min(T + 2*(nest+nbuf+2),(BGAP/2<<l)-1);
        vector<long long> pr;
        for(int b=0; b<nblk; ++b)
        {
            int X0,X1,Y0,Y1;
            block_region(l,b,X0,X1,Y0,Y1);
            int hx0 = X0-H, hx1 = X1+H, hy0 = Y0-H, hy1 = Y1+H;
            if(b==0)
            {
                hx0 = std::max(hx0,0); hx1 = std::min(hx1,(GNX<<l)-1);
                hy0 = std::max(hy0,0); hy1 = std::min(hy1,(GNY<<l)-1);
            }
            for(int I=hx0; I<=hx1; ++I)
            for(int J=hy0; J<=hy1; ++J)
            {
                if(I>=X0 && I<=X1 && J>=Y0 && J<=Y1)
                continue;
                int I2,J2;
                if(!canon(l,b,I,J,I2,J2))
                continue;
                const long long s = (long long)(I/T)*gtny[l] + J/T;
                const long long d = (long long)(I2/T)*gtny[l] + J2/T;
                if(s!=d)
                pr.push_back(s*(1LL<<31) + d);
            }
        }
        std::sort(pr.begin(),pr.end());
        pr.erase(std::unique(pr.begin(),pr.end()),pr.end());
        for(long long v : pr)
        {
            hmap[l].push_back(int(v>>31));
            hmap[l].push_back(int(v&((1LL<<31)-1)));
        }
    }
}

// level-l tile map: forward, a marked tile around a block marks the tiles of the cells it stands for
// (refinement continues across the interfaces); back, a tile around a block stays marked only where
// all those tiles are marked (proper nesting across the interfaces)
void reefamr::blocks_propagate(int l, vector<unsigned char> &M, bool back)
{
    if(nblk==0 || l>=(int)hmap.size())
    return;
    const vector<int> &H = hmap[l];
    for(int pass=0; pass<16; ++pass)
    {
        bool changed = false;
        for(size_t k=0; k<H.size(); k+=2)
        {
            const int s = H[k], d = H[k+1];
            if(!back && M[s] && !M[d])
            {
                M[d] = 1;
                changed = true;
            }
            if(back && M[s] && !M[d])
            {
                M[s] = 0;
                changed = true;
            }
        }
        if(!changed)
        break;
    }
}

// fill plans of level 0: the cells around the blocks and the ring of the native grid
void reefamr::blocks_plans(ghostcell *pgc)
{
    const int me = p0->mpirank;
    const int np = p0->mpi_size;

    vector<vector<int>> req(np), srv;
    vector<vector<int>> dst(np);
    int lost = 0;

    // holder of cell (i,j) of block b: rank, grid id on that rank, local index there
    auto src = [&](int b, int i, int j, int &r, int &g, int &si, int &sj)
    {
        if(b==0)
        {
            int I,J;
            block_virt(0,0,i,j,I,J);
            r = owner(I,J);
            g = -1;
            si = I - rbx0[r];
            sj = J - rby0[r];
            return;
        }
        r = blk[b].rank;
        g = blk[b].lid;
        si = i + EXT;
        sj = j + EXT;
    };

    auto add = [&](int b, int i, int j, int di, int dj, int dg, vector<reefamr_fill> &L)
    {
        int b2,i2,j2;
        if(!block_resolve(0,b,i,j,b2,i2,j2))
        {
            ++lost;
            return;
        }
        int r,g,si,sj;
        src(b2,i2,j2,r,g,si,sj);

        reefamr_fill f;
        f.di=di; f.dj=dj; f.ox=f.oy=0; f.slot=-1; f.aux=0.0;
        if(r==me)
        {
            f.kind=0; f.g=g; f.si=si; f.sj=sj;
            L.push_back(f);
            return;
        }
        f.kind=2; f.g=-1; f.si=f.sj=0;
        L.push_back(f);
        req[r].push_back(g); req[r].push_back(si); req[r].push_back(sj);
        dst[r].push_back(dg); dst[r].push_back(int(L.size())-1);
    };

    for(int id : lev[0])
    {
        reefamr_patch *c = P[id];
        lexer *pp = c->pp;
        c->fill.clear();
        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            if(ii>=EXT && ii<EXT+c->nx && jj>=EXT && jj<EXT+c->ny)
            continue;
            if(pp->flagslice4[(ii-pp->imin)*pp->jmax + (jj-pp->jmin)]<0)
            continue;
            add(c->blk,ii-EXT,jj-EXT,ii,jj,id,c->fill);
        }
    }

    nfill.clear();
    for(int ii=0; ii<p0->knox; ++ii)
    for(int jj=0; jj<p0->knoy; ++jj)
    if(nring[lij(p0,ii,jj)]==1 && p0->flagslice4[lij(p0,ii,jj)]>0)
    {
        int i,j;
        block_own(0,0,ii+O0i,jj+O0j,i,j);
        add(0,i,j,ii,jj,-1,nfill);
    }

    xsetup(req,srv,3);

    reefamr_xplan &X = gplan[0];
    X = reefamr_xplan();
    gserve[0].clear();
    grecv[0].clear();

    for(int r=0; r<np; ++r)
    if(!srv[r].empty())
    {
        X.speer.push_back(r);
        vector<int> items;
        for(size_t k=0; k<srv[r].size(); k+=3)
        {
            reefamr_fill f;
            f.di=f.dj=0; f.ox=f.oy=0; f.slot=-1; f.aux=0.0;
            f.kind=0; f.g=srv[r][k]; f.si=srv[r][k+1]; f.sj=srv[r][k+2];
            items.push_back(int(gserve[0].size()));
            gserve[0].push_back(f);
        }
        X.sitem.push_back(items);
        X.sbuf.push_back(vector<double>());
    }

    for(int r=0; r<np; ++r)
    if(!req[r].empty())
    {
        const int kp = int(X.rpeer.size());
        X.rpeer.push_back(r);
        X.rcount.push_back(int(req[r].size())/3);
        for(size_t k=0; k<dst[r].size(); k+=2)
        {
            reefamr_fill *f = (dst[r][k]<0) ? &nfill[dst[r][k+1]] : &P[dst[r][k]]->fill[dst[r][k+1]];
            f->si = kp;
            f->slot = int(k/2);
            grecv[0].push_back(f);
        }
    }

    lost = pgc->globalisum(lost);
    if(lost>0 && me==0)
    cout<<par.name<<": "<<lost<<" cells beyond the block faces without a source (corners at singular vertices); they keep their values"<<endl;
}
