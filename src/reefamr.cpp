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
#include"mgcslice1.h"
#include"mgcslice2.h"
#include"mgcslice4.h"
#include<mpi.h>

namespace
{
// floor(a/2^l) for negative a as well
inline int fsh(int a, int l) { return (a>=0) ? (a>>l) : -(((-a)-1)>>l)-1; }

// index of cell (i,j) in the lexer 2D arrays (wet, flagslice4)
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

inline bool inarr(const lexer *q, int ii, int jj)
{
    return ii>=q->imin && ii<q->imin+q->imax && jj>=q->jmin && jj<q->jmin+q->jmax;
}
}

reefamr_comms_off::reefamr_comms_off(ghostcell *gg) : g(gg)
{
    old = g->set_comms(false);
}

reefamr_comms_off::~reefamr_comms_off()
{
    g->set_comms(old);
}

reefamr::reefamr(lexer *p, ghostcell *pgc) : patches_total(0), p0(p), pgc0(pgc)
{
    maxlev = 0;
    nest = 3;
    tile = 8;
    nbuf = 0;
    regrid_int = 0;
    keep = 0;
    EXT = 2;
    regrids = 0;
    cells_total = 0;
    last_owner = 0;
    O0i = O0j = NX0 = NY0 = GNX = GNY = 0;

    match.resize(1);
    rmatch.resize(1);
}

reefamr::~reefamr()
{
}

void reefamr::configure(const reefamr_param &q)
{
    par = q;
    maxlev = q.maxlev;
    regrid_int = q.regrid;
    nbuf = q.nbuf;
    tile = q.tile;
    nest = q.nest;
    keep = q.keep;
    EXT = q.ext;
}

// --------------------------------------------------------------------- grids
void reefamr::goff(int g, int &oi, int &oj)
{
    if(g<0)
    {
        oi = O0i; oj = O0j;
        return;
    }
    oi = P[g]->I0-EXT;
    oj = P[g]->J0-EXT;
}

lexer* reefamr::glex(int g)
{
    return (g<0) ? p0 : P[g]->pp;
}

int reefamr::patch_at(int l, int I, int J)
{
    if(I<bxlo(l) || I>bxhi(l) || J<bylo(l) || J>byhi(l))
    return -3;

    if(l==0)
    return -1;

    int ti = I/tile - tti0[l];
    int tj = J/tile - ttj0[l];
    return tmap[l][ti*tny[l]+tj];
}

int reefamr::owner(int I, int J)
{
    if(I<0 || I>=GNX || J<0 || J>=GNY)
    return -1;

    int r = last_owner;
    if(I>=rbx0[r] && I<=rbx1[r] && J>=rby0[r] && J<=rby1[r])
    return r;

    for(r=0; r<(int)rbx0.size(); ++r)
    if(I>=rbx0[r] && I<=rbx1[r] && J>=rby0[r] && J<=rby1[r])
    {
        last_owner = r;
        return r;
    }

    return -1;
}

// the patch of level l at (I,J), from the global patch table: index into GP, -1 none
int reefamr::gpatch_at(int l, int I, int J)
{
    if(l<1 || l>maxlev || I<0 || J<0)
    return -1;

    const int ti = I/tile, tj = J/tile;
    if(ti>=gtnx[l] || tj>=gtny[l] || gtoff[l].empty())
    return -1;

    const size_t t = (size_t)ti*gtny[l]+tj;
    for(int m=gtoff[l][t]; m<gtoff[l][t+1]; ++m)
    {
        const reefamr_gpatch &G = GP[gtlist[l][m]];
        if(I>=G.I0 && I<=G.I1 && J>=G.J0 && J<=G.J1)
        return gtlist[l][m];
    }
    return -1;
}

// the rank of the level-l patch at (I,J), else the holder of the parent cell, down to the owner of
// the level-0 cell (with the patches cut at the rank boxes: the owner of the level-0 cell below)
int reefamr::holder(int l, int I, int J)
{
    for(int m=l; m>=1; --m)
    {
        const int q = gpatch_at(m,I,J);
        if(q>=0)
        return GP[q].rank;
        I = fsh(I,1);
        J = fsh(J,1);
    }
    return owner(I,J);
}

// global patch table: the patches of all ranks, and per level the patches of every global tile
void reefamr::build_gtable()
{
    const int np = p0->mpi_size;

    vector<int> mine;
    for(int id=0; id<(int)P.size(); ++id)
    {
        reefamr_patch *c = P[id];
        mine.push_back(c->lev); mine.push_back(c->I0); mine.push_back(c->I1);
        mine.push_back(c->J0); mine.push_back(c->J1); mine.push_back(id);
    }

    int n = (int)mine.size();
    vector<int> cnt(np), off(np,0);
    MPI_Allgather(&n,1,MPI_INT,&cnt[0],1,MPI_INT,MPI_COMM_WORLD);
    int tot=0;
    for(int r=0; r<np; ++r)
    {
        off[r]=tot;
        tot+=cnt[r];
    }
    vector<int> all(MAX(tot,1));
    MPI_Allgatherv(mine.empty() ? nullptr : &mine[0],n,MPI_INT,&all[0],&cnt[0],&off[0],MPI_INT,MPI_COMM_WORLD);

    GPold = GP;
    GP.clear();
    for(int r=0; r<np; ++r)
    for(int k=off[r]; k<off[r]+cnt[r]; k+=6)
    GP.push_back(reefamr_gpatch{all[k],all[k+1],all[k+2],all[k+3],all[k+4],r,all[k+5]});

    gtoff.assign(maxlev+1,vector<int>());
    gtlist.assign(maxlev+1,vector<int>());
    for(int l=1; l<=maxlev; ++l)
    {
        const size_t nt = (size_t)gtnx[l]*gtny[l];
        vector<int> &O = gtoff[l];
        O.assign(nt+1,0);
        for(const reefamr_gpatch &G : GP)
        if(G.lev==l)
        for(int ti=G.I0/tile; ti<=G.I1/tile; ++ti)
        for(int tj=G.J0/tile; tj<=G.J1/tile; ++tj)
        ++O[(size_t)ti*gtny[l]+tj+1];
        for(size_t t=0; t<nt; ++t)
        O[t+1] += O[t];

        vector<int> pos(O.begin(),O.end()-1);
        gtlist[l].assign(O[nt],0);
        for(int q=0; q<(int)GP.size(); ++q)
        {
            const reefamr_gpatch &G = GP[q];
            if(G.lev!=l)
            continue;
            for(int ti=G.I0/tile; ti<=G.I1/tile; ++ti)
            for(int tj=G.J0/tile; tj<=G.J1/tile; ++tj)
            gtlist[l][pos[(size_t)ti*gtny[l]+tj]++] = q;
        }
    }
}

// level-0 solid flags of the rank box and a halo of FH cells (the ghostcell halo
// of flagslice4 is only one cell deep)
int reefamr::flag0(int I, int J)
{
    if(I<0 || I>=GNX || J<0 || J>=GNY)
    return -10;

    int ii = I-O0i+FH, jj = J-O0j+FH;
    if(ii<0 || jj<0 || ii>=NX0+2*FH || jj>=NY0+2*FH)
    return -10;

    return fl0[(size_t)ii*(NY0+2*FH)+jj];
}

void reefamr::build_flags(lexer *p)
{
    const int me = p->mpirank;
    const int ny = NY0+2*FH;
    fl0.assign((size_t)(NX0+2*FH)*ny,-10);

    vector<vector<int>> req(p->mpi_size), srv;
    vector<vector<int>> pos(p->mpi_size);

    for(int ii=0; ii<NX0+2*FH; ++ii)
    for(int jj=0; jj<ny; ++jj)
    {
        int I = ii-FH+O0i, J = jj-FH+O0j;
        int o = owner(I,J);
        if(o<0)
        continue;

        if(o==me)
        {
            fl0[(size_t)ii*ny+jj] = p->flagslice4[lij(p,I-O0i,J-O0j)];
            continue;
        }

        req[o].push_back(I); req[o].push_back(J);
        pos[o].push_back(ii*ny+jj);
    }

    xsetup(req,srv,2);

    reefamr_xplan X;
    for(int r=0; r<p->mpi_size; ++r)
    {
        if(!srv[r].empty())
        {
            X.speer.push_back(r);
            vector<double> v;
            for(size_t k=0; k<srv[r].size(); k+=2)
            v.push_back(double(p->flagslice4[lij(p,srv[r][k]-O0i,srv[r][k+1]-O0j)]));
            X.sbuf.push_back(v);
        }
        if(!req[r].empty())
        {
            X.rpeer.push_back(r);
            X.rcount.push_back((int)req[r].size()/2);
        }
    }
    xrun(X,1,7002);

    for(size_t k=0; k<X.rpeer.size(); ++k)
    {
        int r = X.rpeer[k];
        for(size_t m=0; m<pos[r].size(); ++m)
        fl0[pos[r][m]] = (int)X.rbuf[k][m];
    }
}

// --------------------------------------------------------------------- setup
// rank boxes, level arrays, tiles, level-0 flags and the no-refinement cells
void reefamr::setup(lexer *p, ghostcell *pgc)
{
    O0i = p->origin_i;
    O0j = p->origin_j;
    NX0 = p->knox;
    NY0 = p->knoy;
    GNX = p->gknox;
    GNY = p->gknoy;

    // rank boxes
    int mine[4] = {O0i, O0i+NX0-1, O0j, O0j+NY0-1};
    vector<int> all(4*p->mpi_size);
    MPI_Allgather(mine,4,MPI_INT,&all[0],4,MPI_INT,MPI_COMM_WORLD);
    rbx0.resize(p->mpi_size); rbx1.resize(p->mpi_size);
    rby0.resize(p->mpi_size); rby1.resize(p->mpi_size);
    for(int r=0; r<p->mpi_size; ++r)
    {
        rbx0[r]=all[4*r]; rbx1[r]=all[4*r+1];
        rby0[r]=all[4*r+2]; rby1[r]=all[4*r+3];
    }

    lev.assign(maxlev+1,vector<int>());
    gtile.assign(maxlev+1,vector<unsigned char>());
    tage.assign(maxlev+1,vector<unsigned char>());
    gtnx.assign(maxlev+1,0);
    gtny.assign(maxlev+1,0);
    tmap.assign(maxlev+1,vector<int>());
    tti0.assign(maxlev+1,0); ttj0.assign(maxlev+1,0);
    tnx.assign(maxlev+1,0); tny.assign(maxlev+1,0);
    gplan.assign(maxlev+1,reefamr_xplan());
    gserve.assign(maxlev+1,vector<reefamr_fill>());
    grecv.assign(maxlev+1,vector<reefamr_fill*>());
    fplan.assign(maxlev+1,reefamr_xplan());
    fsend.assign(maxlev+1,vector<int>());
    frecv.assign(maxlev+1,vector<int>());
    bup.assign(maxlev+1,reefamr_xplan());
    bdn.assign(maxlev+1,reefamr_xplan());
    bloc.assign(maxlev+1,vector<int>());
    bsrv.assign(maxlev+1,vector<reefamr_block>());
    bnloc.assign(maxlev+1,0);
    gtoff.assign(maxlev+1,vector<int>());
    gtlist.assign(maxlev+1,vector<int>());

    vfac.assign(maxlev+1,1);
    for(int l=1; l<=maxlev; ++l)
    {
        int vr = (l<(int)par.vref.size()) ? par.vref[l] : 1;
        vfac[l] = vfac[l-1]*(vr==2 ? 2 : 1);
    }

    for(int l=1; l<=maxlev; ++l)
    {
        gtnx[l] = ((GNX<<l)+tile-1)/tile;
        gtny[l] = ((GNY<<l)+tile-1)/tile;
    }
    build_tiles();
    build_flags(p);

    // no refinement next to in- and outflow boundaries and in the no-refinement boxes
    forbid0.assign((size_t)GNX*GNY,0);
    const int B = par.ioband;
    auto forbid_around = [&](int ii, int jj)
    {
        int I=ii+O0i, J=jj+O0j;
        for(int a=MAX(I-B,0); a<=MIN(I+B,GNX-1); ++a)
        for(int d=MAX(J-B,0); d<=MIN(J+B,GNY-1); ++d)
        forbid0[(size_t)a*GNY+d]=1;
    };
    for(n=0; n<p->gcslin_count; ++n)
    forbid_around(p->gcslin[n][0],p->gcslin[n][1]);
    for(n=0; n<p->gcslout_count; ++n)
    forbid_around(p->gcslout[n][0],p->gcslout[n][1]);

    for(size_t q=0; q+3<par.fbox.size(); q+=4)
    SLICELOOP4
    if(p->XP[IP]>=par.fbox[q] && p->XP[IP]<=par.fbox[q+1] && p->YP[JP]>=par.fbox[q+2] && p->YP[JP]<=par.fbox[q+3])
    forbid0[(size_t)(i+O0i)*GNY+(j+O0j)]=1;

    global_or(forbid0);
}

void reefamr::build_tiles()
{
    for(int l=1; l<=maxlev; ++l)
    {
        tti0[l] = bxlo(l)/tile;
        ttj0[l] = bylo(l)/tile;
        tnx[l] = bxhi(l)/tile - tti0[l] + 1;
        tny[l] = byhi(l)/tile - ttj0[l] + 1;
        tmap[l].assign(tnx[l]*tny[l],-2);

        for(int id : lev[l])
        {
            reefamr_patch *c = P[id];
            for(int ti=c->I0/tile; ti<=c->I1/tile; ++ti)
            for(int tj=c->J0/tile; tj<=c->J1/tile; ++tj)
            tmap[l][(ti-tti0[l])*tny[l]+(tj-ttj0[l])] = id;
        }
    }
}

void reefamr::global_or(vector<unsigned char> &v)
{
    if(p0->mpi_size>1 && !v.empty())
    MPI_Allreduce(MPI_IN_PLACE,&v[0],(int)v.size(),MPI_UNSIGNED_CHAR,MPI_BOR,MPI_COMM_WORLD);
}

// the level-lf tiles covering the children of the level lf-1 cell (I,J) and buf cells around it
void reefamr::tag_cell(int lf, int I, int J, int buf, vector<unsigned char> &M)
{
    const int T = tile;
    int a0 = MAX((2*(I-buf))/T,0), a1 = MIN((2*(I+buf)+1)/T,gtnx[lf]-1);
    int d0 = MAX((2*(J-buf))/T,0), d1 = MIN((2*(J+buf)+1)/T,gtny[lf]-1);
    if(2*(I-buf)<0) a0=0;
    if(2*(J-buf)<0) d0=0;
    for(int a=a0; a<=a1; ++a)
    for(int d=d0; d<=d1; ++d)
    M[(size_t)a*gtny[lf]+d]=1;
}

// --------------------------------------------------------------------- patches
// allocate the patch, its lexer and (through the module) its objects; the state is set later
reefamr_patch* reefamr::make_patch(lexer *p, ghostcell *pgc, int l, int I0, int I1, int J0, int J1)
{
    reefamr_patch *c = patch_new();
    c->lev = l;
    c->I0=I0; c->I1=I1; c->J0=J0; c->J1=J1;
    c->nx = I1-I0+1;
    c->ny = J1-J0+1;
    c->fresh = true;

    c->pp = new lexer(*p,1);
    build_lexer(p,*c);

    {
    reefamr_comms_off guard(pgc);
    patch_objects(c,pgc);
    }

    return c;
}

void reefamr::free_patch(reefamr_patch *c)
{
    lexer *pp = c->pp;

    patch_delete(c);

    const int nxa = pp->knox+1+4*marge;
    const int nya = pp->knoy+1+4*marge;
    const int nsl = pp->imax*pp->jmax;

    pp->del_Darray(pp->XN,nxa); pp->del_Darray(pp->XP,nxa); pp->del_Darray(pp->DXN,nxa); pp->del_Darray(pp->DXP,nxa);
    pp->del_Darray(pp->YN,nya); pp->del_Darray(pp->YP,nya); pp->del_Darray(pp->DYN,nya); pp->del_Darray(pp->DYP,nya);

    if(c->zown)
    {
        const int nza = pp->knoz+1+4*marge;
        pp->del_Darray(pp->ZN,nza); pp->del_Darray(pp->ZP,nza); pp->del_Darray(pp->DZN,nza); pp->del_Darray(pp->DZP,nza);
    }

    pp->del_Iarray(pp->flagslice1,nsl);
    pp->del_Iarray(pp->flagslice2,nsl);
    pp->del_Iarray(pp->flagslice4,nsl);
    pp->del_Iarray(pp->wet,nsl);
    pp->del_Iarray(pp->wet_n,nsl);
    pp->del_Iarray(pp->deep,nsl);

    if(c->bc2D)
    {
    pp->del_Iarray(pp->IOSL,nsl);
    pp->del_Iarray(pp->sizeS1,5);
    pp->del_Iarray(pp->sizeS2,5);
    pp->del_Iarray(pp->sizeS4,5);
    pp->del_Iarray(pp->gcbsl1,pp->gcbsl1_count,6);
    pp->del_Iarray(pp->gcbsl2,pp->gcbsl2_count,6);
    pp->del_Iarray(pp->gcbsl4,pp->gcbsl4_count,6);
    pp->del_Iarray(pp->dgcsl1,pp->dgcsl1_count,3);
    pp->del_Iarray(pp->dgcsl2,pp->dgcsl2_count,3);
    pp->del_Iarray(pp->dgcsl4,pp->dgcsl4_count,3);
    }

    delete pp;
    delete c;
}

void reefamr::free_patches()
{
    for(auto c : P)
    free_patch(c);
    P.clear();
}

// patch geometry: nodes of the global level-l grid, subdividing the level-0 nodes
void reefamr::build_lexer(lexer *p, reefamr_patch &c)
{
    lexer *pp = c.pp;
    const int ms = pp->margin;   // ghost layers of the slices
    const int m = marge;         // offset of the coordinate arrays (IP = i+marge)
    const int l = c.lev;
    const double rl = double(1<<l);

    pp->knox = c.nx+2*EXT;
    pp->knoy = c.ny+2*EXT;
    pp->knoz = p->knoz;
    pp->imin = -ms;
    pp->jmin = -ms;
    pp->imax = pp->knox+2*ms;
    pp->jmax = pp->knoy+2*ms;
    pp->kmin = p->kmin;
    pp->kmax = p->kmax;
    pp->kmaxF = p->kmaxF;

    pp->i_dir = p->i_dir; pp->j_dir = p->j_dir; pp->k_dir = p->k_dir;
    pp->x_dir = p->x_dir; pp->y_dir = p->y_dir; pp->z_dir = p->z_dir;

    pp->ulast = pp->vlast = pp->wlast = pp->flast = 0;
    pp->ulastsflow = 0;

    // coordinates
    const int nxa = pp->knox+1+4*m;
    const int nya = pp->knoy+1+4*m;
    const int pnxa = p->knox+1+4*m;
    const int pnya = p->knoy+1+4*m;

    pp->Darray(pp->XN,nxa); pp->Darray(pp->XP,nxa); pp->Darray(pp->DXN,nxa); pp->Darray(pp->DXP,nxa);
    pp->Darray(pp->YN,nya); pp->Darray(pp->YP,nya); pp->Darray(pp->DYN,nya); pp->Darray(pp->DYP,nya);
    pp->ZN = p->ZN; pp->ZP = p->ZP; pp->DZN = p->DZN; pp->DZP = p->DZP;
    build_vertical(p,c);

    auto x0 = [&](int K)   // global level-0 node K, linear extrapolation outside the rank arrays
    {
        int a = K-O0i+m;
        if(a<0)      return p->XN[0] + a*(p->XN[1]-p->XN[0]);
        if(a>=pnxa)  return p->XN[pnxa-1] + (a-pnxa+1)*(p->XN[pnxa-1]-p->XN[pnxa-2]);
        return p->XN[a];
    };
    auto y0 = [&](int K)
    {
        int a = K-O0j+m;
        if(a<0)      return p->YN[0] + a*(p->YN[1]-p->YN[0]);
        if(a>=pnya)  return p->YN[pnya-1] + (a-pnya+1)*(p->YN[pnya-1]-p->YN[pnya-2]);
        return p->YN[a];
    };

    for(int a=0; a<nxa; ++a)
    {
        int N = a-m-EXT+c.I0;      // global level-l node
        int K = fsh(N,l);
        int r = N-(K<<l);
        pp->XN[a] = (r==0) ? x0(K) : x0(K) + (x0(K+1)-x0(K))*(double(r)/rl);
    }
    for(int a=0; a<nya; ++a)
    {
        int N = a-m-EXT+c.J0;
        int K = fsh(N,l);
        int r = N-(K<<l);
        pp->YN[a] = (r==0) ? y0(K) : y0(K) + (y0(K+1)-y0(K))*(double(r)/rl);
    }

    for(int a=0; a<nxa-1; ++a) { pp->XP[a] = 0.5*(pp->XN[a]+pp->XN[a+1]); pp->DXN[a] = pp->XN[a+1]-pp->XN[a]; }
    for(int a=0; a<nya-1; ++a) { pp->YP[a] = 0.5*(pp->YN[a]+pp->YN[a+1]); pp->DYN[a] = pp->YN[a+1]-pp->YN[a]; }
    pp->XP[nxa-1] = pp->XP[nxa-2] + pp->DXN[nxa-2]; pp->DXN[nxa-1] = pp->DXN[nxa-2];
    pp->YP[nya-1] = pp->YP[nya-2] + pp->DYN[nya-2]; pp->DYN[nya-1] = pp->DYN[nya-2];
    for(int a=0; a<nxa-1; ++a) pp->DXP[a] = pp->XP[a+1]-pp->XP[a];
    for(int a=0; a<nya-1; ++a) pp->DYP[a] = pp->YP[a+1]-pp->YP[a];
    pp->DXP[nxa-1] = pp->DXP[nxa-2];
    pp->DYP[nya-1] = pp->DYP[nya-2];

    // mean spacing: exact fraction of level 0, independent of the patch size
    pp->DXM = p->DXM/rl;
    pp->DYM = pp->DXM;
    pp->DXD = p->DXD;
    pp->DYD = p->DYD;

    // flags from level 0, wet arrays
    const int nsl = pp->imax*pp->jmax;
    pp->Iarray(pp->flagslice1,nsl);
    pp->Iarray(pp->flagslice2,nsl);
    pp->Iarray(pp->flagslice4,nsl);
    pp->Iarray(pp->wet,nsl);
    pp->Iarray(pp->wet_n,nsl);
    pp->Iarray(pp->deep,nsl);

    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        int I = ii-EXT+c.I0, J = jj-EXT+c.J0;
        pp->flagslice4[lij(pp,ii,jj)] = flag0(fsh(I,l),fsh(J,l));
    }

    // face flags as mgcslice1/2::makemgc, over the whole array
    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        int q = lij(pp,ii,jj);
        int f4 = pp->flagslice4[q];
        pp->flagslice1[q] = f4;
        pp->flagslice2[q] = f4;

        if(f4>0 && inarr(pp,ii+1,jj) && pp->flagslice4[lij(pp,ii+1,jj)]<0)
        pp->flagslice1[q] = pp->flagslice4[lij(pp,ii+1,jj)];

        if(f4>0 && inarr(pp,ii,jj+1) && pp->flagslice4[lij(pp,ii,jj+1)]<0)
        pp->flagslice2[q] = pp->flagslice4[lij(pp,ii,jj+1)];
    }

    for(int q=0; q<nsl; ++q)
    {
        pp->wet[q] = pp->wet_n[q] = pp->deep[q] = 0;
    }

    // boundary lists are built in build_bc2D; no partition neighbours
    pp->gcbsl1 = pp->gcbsl2 = pp->gcbsl3 = pp->gcbsl4 = pp->gcbsl4a = nullptr;
    pp->dgcsl1 = pp->dgcsl2 = pp->dgcsl3 = pp->dgcsl4 = nullptr;
    pp->gcslin = pp->gcslout = pp->gcslawa1 = pp->gcslawa2 = nullptr;
    pp->gcbsl1_count = pp->gcbsl2_count = pp->gcbsl3_count = pp->gcbsl4_count = pp->gcbsl4a_count = 0;
    pp->gcslin_count = pp->gcslout_count = 0;
    pp->gcslawa1_count = pp->gcslawa2_count = 0;
    pp->dgcsl1_count = pp->dgcsl2_count = pp->dgcsl3_count = pp->dgcsl4_count = 0;
    pp->gcsldfeta4_count = pp->gcsldfbed4_count = 0;
    pp->gcslpara1_count = pp->gcslpara2_count = pp->gcslpara3_count = pp->gcslpara4_count = 0;
    pp->gcslparaco1_count = pp->gcslparaco2_count = pp->gcslparaco3_count = pp->gcslparaco4_count = 0;
    pp->nb1 = pp->nb2 = pp->nb3 = pp->nb4 = pp->nb5 = pp->nb6 = -2;
    pp->periodic1 = pp->periodic2 = pp->periodic3 = 0;   // no periodic ghost copies on a patch (the lexer copy leaves them unset)
    pp->mpi_size = 1;
    pp->mpirank = p->mpirank;

    // the domain boundary never lies on a patch: all boundary faces are walls (bc 21)
    pp->origin_i = pp->origin_j = pp->origin_k = -100000000;
    pp->gknox = p->gknox; pp->gknoy = p->gknoy; pp->gknoz = p->gknoz;

    pp->vec2Dlength = nsl;
    pp->veclength = 0;
    pp->cellnum2D = pp->knox*pp->knoy;
    pp->slicenum = nsl;

    pp->global_xmin = p->global_xmin; pp->global_ymin = p->global_ymin;
    pp->global_xmax = p->global_xmax; pp->global_ymax = p->global_ymax;
    pp->originx = p->originx; pp->originy = p->originy;

    pp->wd = p->wd;
    pp->phimean = p->phimean;
    pp->phiout = p->phiout;
    pp->dt = p->dt;
    pp->dt_old = p->dt_old;
    pp->simtime = p->simtime;
    pp->count = p->count;
}

// vertical grid of the patch: the level-0 arrays, or (vref) the level-0 nodes subdivided by
// the vertical factor of the level.  The coarse nodes are nodes of the fine grid (nested),
// so injection and interpolation in the vertical need no search.
void reefamr::build_vertical(lexer *p, reefamr_patch &c)
{
    lexer *pp = c.pp;
    const int f = vfac[c.lev];
    c.zown = false;

    if(f==1)
    return;

    const int m = marge;
    const int ms = pp->margin;
    const double rf = double(f);

    pp->knoz = p->knoz*f;
    pp->kmin = -ms;
    pp->kmax = pp->knoz+2*ms;
    pp->kmaxF = pp->knoz+1+2*ms;
    pp->gknoz = p->gknoz*f;

    const int nza = pp->knoz+1+4*m;
    const int pnza = p->knoz+1+4*m;

    pp->ZN = pp->ZP = pp->DZN = pp->DZP = nullptr;
    pp->Darray(pp->ZN,nza); pp->Darray(pp->ZP,nza); pp->Darray(pp->DZN,nza); pp->Darray(pp->DZP,nza);
    c.zown = true;

    auto z0 = [&](int K)   // level-0 node K, linear extrapolation outside the arrays
    {
        int a = K+m;
        if(a<0)      return p->ZN[0] + a*(p->ZN[1]-p->ZN[0]);
        if(a>=pnza)  return p->ZN[pnza-1] + (a-pnza+1)*(p->ZN[pnza-1]-p->ZN[pnza-2]);
        return p->ZN[a];
    };

    for(int a=0; a<nza; ++a)
    {
        int N = a-m;
        int K = (N>=0) ? N/f : -((-N+f-1)/f);
        int r = N-K*f;
        pp->ZN[a] = (r==0) ? z0(K) : z0(K) + (z0(K+1)-z0(K))*(double(r)/rf);
    }

    for(int a=0; a<nza-1; ++a) { pp->ZP[a] = 0.5*(pp->ZN[a]+pp->ZN[a+1]); pp->DZN[a] = pp->ZN[a+1]-pp->ZN[a]; }
    pp->ZP[nza-1] = pp->ZP[nza-2] + pp->DZN[nza-2]; pp->DZN[nza-1] = pp->DZN[nza-2];
    for(int a=0; a<nza-1; ++a) pp->DZP[a] = pp->ZP[a+1]-pp->ZP[a];
    pp->DZP[nza-1] = pp->DZP[nza-2];
}

// 2D wall lists of the patch, as driver::makegrid2D for level 0
void reefamr::build_bc2D(ghostcell *pgc, reefamr_patch &c)
{
    lexer *pp = c.pp;

    pp->gridini2D();

    mgcslice1 m1(pp);
    mgcslice2 m2(pp);
    mgcslice4 m4(pp);

    m1.gcb_seed(pp);
    m2.gcb_seed(pp);
    m4.gcb_seed(pp);

    pgc->gcsl_setbc1(pp);
    pgc->gcsl_setbc2(pp);
    pgc->gcsl_setbc4(pp);

    pgc->gcsl_setbcio(pp);

    pgc->dgcslini1(pp);
    pgc->dgcslini2(pp);
    pgc->dgcslini4(pp);

    c.bc2D = true;
}

// --------------------------------------------------------------------- MPI helpers
// req[r]: keys (K ints per item) this rank wants from rank r; srv[r]: keys rank r wants from this rank
void reefamr::xsetup(vector<vector<int>> &req, vector<vector<int>> &srv, int K)
{
    const int np = p0->mpi_size;
    srv.assign(np,vector<int>());

    if(np==1)
    return;

    vector<int> sc(np), rc(np), sd(np,0), rd(np,0);
    for(int r=0; r<np; ++r)
    sc[r] = (int)req[r].size();

    MPI_Alltoall(&sc[0],1,MPI_INT,&rc[0],1,MPI_INT,MPI_COMM_WORLD);

    int ns=0, nr=0;
    for(int r=0; r<np; ++r)
    {
        sd[r]=ns; ns+=sc[r];
        rd[r]=nr; nr+=rc[r];
    }

    vector<int> sb(MAX(ns,1)), rb(MAX(nr,1));
    for(int r=0; r<np; ++r)
    for(int k=0; k<sc[r]; ++k)
    sb[sd[r]+k] = req[r][k];

    MPI_Alltoallv(&sb[0],&sc[0],&sd[0],MPI_INT,&rb[0],&rc[0],&rd[0],MPI_INT,MPI_COMM_WORLD);

    for(int r=0; r<np; ++r)
    srv[r].assign(rb.begin()+rd[r],rb.begin()+rd[r]+rc[r]);
}

void reefamr::xrun(reefamr_xplan &X, int nv, int tag)
{
    vector<MPI_Request> rq;
    X.rbuf.resize(X.rpeer.size());

    for(size_t k=0; k<X.rpeer.size(); ++k)
    {
        X.rbuf[k].resize((size_t)X.rcount[k]*nv);
        if(!X.rbuf[k].empty())
        {
            rq.push_back(MPI_REQUEST_NULL);
            MPI_Irecv(&X.rbuf[k][0],(int)X.rbuf[k].size(),MPI_DOUBLE,X.rpeer[k],tag,MPI_COMM_WORLD,&rq.back());
        }
    }

    for(size_t k=0; k<X.speer.size(); ++k)
    if(!X.sbuf[k].empty())
    {
        rq.push_back(MPI_REQUEST_NULL);
        MPI_Isend(&X.sbuf[k][0],(int)X.sbuf[k].size(),MPI_DOUBLE,X.speer[k],tag,MPI_COMM_WORLD,&rq.back());
    }

    if(!rq.empty())
    MPI_Waitall((int)rq.size(),&rq[0],MPI_STATUSES_IGNORE);
}
