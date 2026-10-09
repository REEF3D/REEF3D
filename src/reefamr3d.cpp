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

#include"reefamr3d.h"
#include"lexer.h"
#include"ghostcell.h"
#include"grid_helper.h"
#include"mgcslice1.h"
#include"mgcslice2.h"
#include"mgcslice4.h"
#include<unordered_map>
#include<algorithm>
#include<cmath>
#include<iostream>

namespace
{
// floor(a/r) for negative a as well
inline int fdiv(int a, int r) { return (a>=0) ? a/r : -((-a+r-1)/r); }
}

reefamr3d_local::reefamr3d_local(ghostcell *gg) : g(gg)
{
    oldc = g->set_comms(false);
    oldl = g->set_local(true);
}

reefamr3d_local::~reefamr3d_local()
{
    g->set_comms(oldc);
    g->set_local(oldl);
}

reefamr3d::reefamr3d(lexer *p, ghostcell *pgc) : p0(p), pgc0(pgc)
{
    maxlev = MAX(p->G1,0);
    myrank = p->mpirank;
    nranks = p->mpi_size;
}

reefamr3d::~reefamr3d()
{
}

// --------------------------------------------------------------------- index helpers
int reefamr3d::ratio(int l, int d) const
{
    int r=1;
    for(int q=0; q<l; ++q)
    r *= rr[d];
    return r;
}

int reefamr3d::gn(int l, int d) const
{
    return gn0[d]*ratio(l,d);
}

bool reefamr3d::in_domain(int l, const int *I) const
{
    for(int d=0; d<3; ++d)
    if(I[d]<0 || I[d]>=gn(l,d))
    return false;
    return true;
}

bool reefamr3d::tile_marked(int l, int ti, int tj, int tk) const
{
    if(ti<0 || tj<0 || tk<0 || ti>=tn[l][0] || tj>=tn[l][1] || tk>=tn[l][2])
    return false;
    return tmask[l][((size_t)ti*tn[l][1] + tj)*tn[l][2] + tk]!=0;
}

bool reefamr3d::refined(int l, const int *I) const
{
    if(l<1 || l>maxlev || !in_domain(l,I))
    return false;
    return tile_marked(l,I[0]/tile[0],I[1]/tile[1],I[2]/tile[2]);
}

bool reefamr3d::covered(int l, const int *I) const
{
    if(l>=maxlev)
    return false;
    int C[3];
    for(int d=0; d<3; ++d)
    C[d] = I[d]*rr[d];
    return refined(l+1,C);
}

int reefamr3d::gpatch_at(int l, const int *I) const
{
    for(int g : glev[l])
    {
        const r3gpatch &q = GP[g];
        if(I[0]>=q.lo[0] && I[0]<=q.hi[0] && I[1]>=q.lo[1] && I[1]<=q.hi[1] && I[2]>=q.lo[2] && I[2]<=q.hi[2])
        return g;
    }
    return -1;
}

lexer* reefamr3d::glex(int g)
{
    return g<0 ? p0 : P[g]->pp;
}

void reefamr3d::goff(int g, int *o)
{
    for(int d=0; d<3; ++d)
    o[d] = (g<0) ? org[d] : P[g]->lo[d];
}

double reefamr3d::node(int l, int d, int I) const
{
    const int R = ratio(l,d);
    const int I0 = fdiv(I,R);
    const double f = double(I - I0*R)/double(R);
    const vector<double> &x = gx[d];
    const int a = I0+marge;
    const int nmax = (int)x.size()-1;
    const int ac = MAX(MIN(a,nmax-1),0);
    return x[ac] + f*(x[ac+1]-x[ac]);
}

// --------------------------------------------------------------------- setup
void reefamr3d::setup(lexer *p, ghostcell *pgc)
{
    rr[0] = (p->i_dir==1 && p->gknox>1) ? 2 : 1;
    rr[1] = (p->j_dir==1 && p->gknoy>1) ? 2 : 1;
    rr[2] = (p->k_dir==1 && p->gknoz>1) ? 2 : 1;

    gn0[0] = p->gknox; gn0[1] = p->gknoy; gn0[2] = p->gknoz;
    org[0] = p->origin_i; org[1] = p->origin_j; org[2] = p->origin_k;
    kn0[0] = p->knox; kn0[1] = p->knoy; kn0[2] = p->knoz;

    // level-0 boxes of all ranks
    rbox.assign(6*nranks,0);
    {
        int mine[6] = {org[0],org[1],org[2],kn0[0],kn0[1],kn0[2]};
        MPI_Allgather(mine,6,MPI_INT,&rbox[0],6,MPI_INT,MPI_COMM_WORLD);
    }

    for(int d=0; d<3; ++d)
    {
        tile[d] = 1;
        if(rr[d]==2)
        {
            tile[d] = MAX(p->G4,2);
            tile[d] += tile[d]%2;
        }
    }

    global_nodes();
    mark_tiles();
    make_boxes();

    // local patches, lexers, module objects
    lev.assign(maxlev+1,vector<int>());
    for(size_t g=0; g<GP.size(); ++g)
    if(GP[g].rank==myrank)
    {
        r3patch *c = patch_new();
        c->gid = (int)g;
        c->lev = GP[g].lev;
        for(int d=0; d<3; ++d)
        {
            c->lo[d] = GP[g].lo[d];
            c->hi[d] = GP[g].hi[d];
            c->n[d] = c->hi[d]-c->lo[d]+1;
        }
        P.push_back(c);
        lev[c->lev].push_back((int)P.size()-1);
    }

    for(auto c : P)
    build_lexer(*c);

    for(auto c : P)
    patch_objects(c);

    plans();
}

// level-0 nodes of the whole domain (every rank), with marge nodes beyond each end
void reefamr3d::global_nodes()
{
    const double *xn[3] = {p0->XN,p0->YN,p0->ZN};

    for(int d=0; d<3; ++d)
    {
        const int nn = gn0[d]+1+2*marge;
        vector<double> v(nn,-1.0e300);

        int i0 = (org[d]==0) ? -marge : 0;
        int i1 = (org[d]+kn0[d]==gn0[d]) ? kn0[d]+marge : kn0[d];
        for(int i=i0; i<=i1; ++i)
        v[org[d]+i+marge] = xn[d][i+marge];

        gx[d].assign(nn,0.0);
        MPI_Allreduce(&v[0],&gx[d][0],nn,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    }
}

// refined tiles of every level: the boxes G 15 on level G 1, every coarser level around the next
// finer one with a nesting band; no refinement in the wave generation and absorption zones
void reefamr3d::mark_tiles()
{
    tn.assign(maxlev+1,vector<int>(3,0));
    tmask.assign(maxlev+1,vector<unsigned char>());

    for(int l=1; l<=maxlev; ++l)
    {
        for(int d=0; d<3; ++d)
        tn[l][d] = (gn(l,d)+tile[d]-1)/tile[d];
        tmask[l].assign((size_t)tn[l][0]*tn[l][1]*tn[l][2],0);
    }

    auto mark_cells = [&](int l, const int *lo, const int *hi)
    {
        int a[3], b[3];
        for(int d=0; d<3; ++d)
        {
            a[d] = MAX(lo[d],0)/tile[d];
            b[d] = MIN(hi[d],gn(l,d)-1)/tile[d];
            if(lo[d]>hi[d] || hi[d]<0 || lo[d]>gn(l,d)-1)
            return;
        }
        for(int ti=a[0]; ti<=b[0]; ++ti)
        for(int tj=a[1]; tj<=b[1]; ++tj)
        for(int tk=a[2]; tk<=b[2]; ++tk)
        tmask[l][((size_t)ti*tn[l][1] + tj)*tn[l][2] + tk] = 1;
    };

    // cells of level l whose centre lies in [xs,xe] along d
    auto range = [&](int l, int d, double xs, double xe, int &lo, int &hi)
    {
        lo = gn(l,d); hi = -1;
        if(rr[d]==1)
        {
            lo = 0; hi = gn(l,d)-1;
            return;
        }
        for(int I=0; I<gn(l,d); ++I)
        {
            double xc = 0.5*(node(l,d,I)+node(l,d,I+1));
            if(xc>=xs && xc<=xe)
            {
                lo = MIN(lo,I);
                hi = MAX(hi,I);
            }
        }
    };

    if(maxlev>=1)
    for(int b=0; b<p0->G15; ++b)
    {
        int lo[3], hi[3];
        const double xs[3] = {p0->G15_xs[b],p0->G15_ys[b],p0->G15_zs[b]};
        const double xe[3] = {p0->G15_xe[b],p0->G15_ye[b],p0->G15_ze[b]};
        for(int d=0; d<3; ++d)
        range(maxlev,d,xs[d],xe[d],lo[d],hi[d]);
        mark_cells(maxlev,lo,hi);
    }

    // nesting: every level covers the parents of the next finer one with a band of nest cells
    const int nest = 3;
    for(int l=maxlev-1; l>=1; --l)
    for(int ti=0; ti<tn[l+1][0]; ++ti)
    for(int tj=0; tj<tn[l+1][1]; ++tj)
    for(int tk=0; tk<tn[l+1][2]; ++tk)
    if(tile_marked(l+1,ti,tj,tk))
    {
        int lo[3], hi[3];
        const int t[3] = {ti,tj,tk};
        for(int d=0; d<3; ++d)
        {
            lo[d] = fdiv(t[d]*tile[d],rr[d]) - (rr[d]==2 ? nest : 0);
            hi[d] = fdiv((t[d]+1)*tile[d]-1,rr[d]) + (rr[d]==2 ? nest : 0);
        }
        mark_cells(l,lo,hi);
    }

    // wave generation and absorption zones (relaxation, B 98 2, B 99 1/2) stay on level 0
    double zx0 = 1.0e20, zx1 = -1.0e20;     // forbidden: x < zx0 and x > zx1
    const double xmin = gx[0][marge], xmax = gx[0][gn0[0]+marge];
    zx0 = -1.0e20; zx1 = 1.0e20;
    if(p0->B98==2 && p0->B96_1>0.0)
    zx0 = xmin + p0->B96_1;
    if((p0->B99==1 || p0->B99==2) && p0->B96_2>0.0)
    zx1 = xmax - p0->B96_2;

    int cleared = 0;
    for(int l=1; l<=maxlev; ++l)
    for(int ti=0; ti<tn[l][0]; ++ti)
    {
        // x extent of the tile, with the margin of 3 fine cells the patch reaches
        const int I0 = ti*tile[0] - 3, I1 = MIN((ti+1)*tile[0],gn(l,0)) + 3;
        const double xa = node(l,0,MAX(I0,0)), xb = node(l,0,MIN(I1,gn(l,0)));
        if(!(xa<zx0 || xb>zx1))
        continue;
        for(int tj=0; tj<tn[l][1]; ++tj)
        for(int tk=0; tk<tn[l][2]; ++tk)
        {
            unsigned char &m = tmask[l][((size_t)ti*tn[l][1] + tj)*tn[l][2] + tk];
            if(m)
            {
                m = 0;
                ++cleared;
            }
        }
    }

    // proper nesting after the zones: a tile of level l needs its parents and one more cell on
    // level l-1
    int dropped = 0;
    for(int l=2; l<=maxlev; ++l)
    for(int ti=0; ti<tn[l][0]; ++ti)
    for(int tj=0; tj<tn[l][1]; ++tj)
    for(int tk=0; tk<tn[l][2]; ++tk)
    if(tile_marked(l,ti,tj,tk))
    {
        const int t[3] = {ti,tj,tk};
        int lo[3], hi[3];
        for(int d=0; d<3; ++d)
        {
            lo[d] = MAX(fdiv(t[d]*tile[d],rr[d]) - (rr[d]==2 ? 1 : 0),0);
            hi[d] = MIN(fdiv((t[d]+1)*tile[d]-1,rr[d]) + (rr[d]==2 ? 1 : 0),gn(l-1,d)-1);
        }
        bool ok = true;
        int I[3];
        for(I[0]=lo[0]; I[0]<=hi[0] && ok; ++I[0])
        for(I[1]=lo[1]; I[1]<=hi[1] && ok; ++I[1])
        for(I[2]=lo[2]; I[2]<=hi[2] && ok; ++I[2])
        if(!refined(l-1,I))
        ok = false;

        if(!ok)
        {
            tmask[l][((size_t)ti*tn[l][1] + tj)*tn[l][2] + tk] = 0;
            ++dropped;
        }
    }

    if(myrank==0 && (cleared>0 || dropped>0))
    cout<<"CFD AMR: "<<cleared<<" tiles in the wave zones and "<<dropped<<" tiles without nesting not refined"<<endl;

    // the edges of a level must not lie on a partition edge of level 0 (the coarse cell across a
    // patch edge on the rank of the patch): such a level is extended by the tiles across, finest
    // level first, and the next coarser level nests the extension
    auto rank_of = [&](const int *I0) -> int
    {
        for(int r=0; r<nranks; ++r)
        {
            bool in = true;
            for(int d=0; d<3; ++d)
            if(I0[d]<rbox[6*r+d] || I0[d]>=rbox[6*r+d]+rbox[6*r+3+d])
            in = false;
            if(in)
            return r;
        }
        return -1;
    };

    // the tiles of level l to extend; returns their number
    auto edges = [&](int l, bool fix) -> int
    {
        vector<size_t> add;
        int I[3];
        for(I[0]=0; I[0]<gn(l,0); ++I[0])
        for(I[1]=0; I[1]<gn(l,1); ++I[1])
        for(I[2]=0; I[2]<gn(l,2); ++I[2])
        {
            if(!refined(l,I))
            continue;
            for(int d=0; d<3; ++d)
            {
                if(rr[d]!=2)
                continue;
                const int R = ratio(l,d);
                for(int dir=-1; dir<=1; dir+=2)
                {
                    int J[3] = {I[0],I[1],I[2]};
                    J[d] += dir;
                    if(!in_domain(l,J) || refined(l,J))
                    continue;
                    const int N = (dir>0) ? I[d]+1 : I[d];
                    if(N%R!=0)
                    continue;
                    int a[3], b[3];
                    for(int e=0; e<3; ++e)
                    a[e] = b[e] = fdiv(I[e],ratio(l,e));
                    a[d] = N/R-1;
                    b[d] = N/R;
                    if(rank_of(a)==rank_of(b))
                    continue;
                    add.push_back(((size_t)(J[0]/tile[0])*tn[l][1] + J[1]/tile[1])*tn[l][2] + J[2]/tile[2]);
                }
            }
        }
        if(fix)
        for(size_t t : add)
        tmask[l][t] = 1;
        return (int)add.size();
    };

    int extended = 0;
    for(int l=maxlev; l>=1; --l)
    {
        for(int pass=0; pass<64; ++pass)
        {
            const int n = edges(l,true);
            extended += n;
            if(n==0)
            break;
        }

        if(l>=2)
        for(int ti=0; ti<tn[l][0]; ++ti)
        for(int tj=0; tj<tn[l][1]; ++tj)
        for(int tk=0; tk<tn[l][2]; ++tk)
        if(tile_marked(l,ti,tj,tk))
        {
            int lo[3], hi[3];
            const int t[3] = {ti,tj,tk};
            for(int d=0; d<3; ++d)
            {
                lo[d] = fdiv(t[d]*tile[d],rr[d]) - (rr[d]==2 ? nest : 0);
                hi[d] = fdiv((t[d]+1)*tile[d]-1,rr[d]) + (rr[d]==2 ? nest : 0);
            }
            mark_cells(l-1,lo,hi);
        }
    }

    // checks after the extensions: the wave zones, the edges
    int bad_zone = 0, bad_edge = 0;
    for(int l=1; l<=maxlev; ++l)
    {
        for(int ti=0; ti<tn[l][0]; ++ti)
        {
            const int I0 = ti*tile[0] - 3, I1 = MIN((ti+1)*tile[0],gn(l,0)) + 3;
            const double xa = node(l,0,MAX(I0,0)), xb = node(l,0,MIN(I1,gn(l,0)));
            if(!(xa<zx0 || xb>zx1))
            continue;
            for(int tj=0; tj<tn[l][1]; ++tj)
            for(int tk=0; tk<tn[l][2]; ++tk)
            if(tile_marked(l,ti,tj,tk))
            ++bad_zone;
        }
        bad_edge += edges(l,false);
    }

    if(myrank==0 && extended>0)
    cout<<"CFD AMR: "<<extended<<" tiles added: the edges of the refined region kept off the partition edges"<<endl;

    if(bad_zone>0 || bad_edge>0)
    {
        if(myrank==0)
        cout<<"CFD AMR: the refined region reaches into a wave zone or keeps an edge on a partition edge; "
              "move or resize the boxes G 15"<<endl;
        MPI_Abort(MPI_COMM_WORLD,-3102);
    }
}

// tiles to boxes (greedy: runs along x, extended along y and z), cut at the rank boxes
void reefamr3d::make_boxes()
{
    GP.clear();
    glev.assign(maxlev+1,vector<int>());

    for(int l=1; l<=maxlev; ++l)
    {
        const int nx = tn[l][0], ny = tn[l][1], nz = tn[l][2];
        vector<unsigned char> used((size_t)nx*ny*nz,0);
        auto id = [&](int a, int b, int c) { return ((size_t)a*ny + b)*nz + c; };
        auto free_ = [&](int a, int b, int c) { return tmask[l][id(a,b,c)] && !used[id(a,b,c)]; };

        for(int tk=0; tk<nz; ++tk)
        for(int tj=0; tj<ny; ++tj)
        for(int ti=0; ti<nx; ++ti)
        {
            if(!free_(ti,tj,tk))
            continue;

            int ti2 = ti;
            while(ti2+1<nx && free_(ti2+1,tj,tk))
            ++ti2;

            int tj2 = tj;
            for(;;)
            {
                if(tj2+1>=ny)
                break;
                bool ok = true;
                for(int a=ti; a<=ti2; ++a)
                if(!free_(a,tj2+1,tk)) { ok=false; break; }
                if(!ok)
                break;
                ++tj2;
            }

            int tk2 = tk;
            for(;;)
            {
                if(tk2+1>=nz)
                break;
                bool ok = true;
                for(int a=ti; a<=ti2 && ok; ++a)
                for(int b=tj; b<=tj2; ++b)
                if(!free_(a,b,tk2+1)) { ok=false; break; }
                if(!ok)
                break;
                ++tk2;
            }

            for(int a=ti; a<=ti2; ++a)
            for(int b=tj; b<=tj2; ++b)
            for(int c=tk; c<=tk2; ++c)
            used[id(a,b,c)] = 1;

            int blo[3] = {ti*tile[0], tj*tile[1], tk*tile[2]};
            int bhi[3] = {MIN((ti2+1)*tile[0],gn(l,0))-1, MIN((tj2+1)*tile[1],gn(l,1))-1, MIN((tk2+1)*tile[2],gn(l,2))-1};

            // cut at the rank boxes
            for(int r=0; r<nranks; ++r)
            {
                r3gpatch q;
                q.lev = l;
                q.rank = r;
                q.lid = -1;
                bool empty = false;
                for(int d=0; d<3; ++d)
                {
                    const int R = ratio(l,d);
                    const int a = rbox[6*r+d]*R, b = (rbox[6*r+d]+rbox[6*r+3+d])*R-1;
                    q.lo[d] = MAX(blo[d],a);
                    q.hi[d] = MIN(bhi[d],b);
                    if(q.lo[d]>q.hi[d])
                    empty = true;
                }

                if(empty)
                continue;
                glev[l].push_back((int)GP.size());
                GP.push_back(q);
            }
        }
    }

    // local ids on the owning rank, in the order of the global list
    vector<int> cnt(nranks,0);
    for(auto &q : GP)
    q.lid = cnt[q.rank]++;
}

// --------------------------------------------------------------------- patch lexer
// geometry, flags and boundary lists of a patch from level 0 below it; the setup of the
// level-0 lexer (read_grid, gridini, flagini, flagfield, makegrid, makegrid2D) on the patch
void reefamr3d::build_lexer(r3patch &c)
{
    reefamr3d_local guard(pgc0);

    lexer *p = p0;
    c.pp = new lexer(*p,1);
    lexer *pp = c.pp;
    const int l = c.lev;
    int R[3];
    for(int d=0; d<3; ++d)
    R[d] = ratio(l,d);

    pp->amr_patch = 1;
    pp->mpirank = p->mpirank;
    pp->mpi_size = 1;

    pp->knox = c.n[0]; pp->knoy = c.n[1]; pp->knoz = c.n[2];
    pp->gknox = gn(l,0); pp->gknoy = gn(l,1); pp->gknoz = gn(l,2);
    pp->origin_i = c.lo[0]; pp->origin_j = c.lo[1]; pp->origin_k = c.lo[2];
    pp->i_dir = p->i_dir; pp->j_dir = p->j_dir; pp->k_dir = p->k_dir;
    pp->nb1 = pp->nb2 = pp->nb3 = pp->nb4 = pp->nb5 = pp->nb6 = -2;
    pp->periodic1 = pp->periodic2 = pp->periodic3 = 0;
    pp->periodicX1 = pp->periodicX2 = pp->periodicX3 = pp->periodicX4 = pp->periodicX5 = pp->periodicX6 = 0;
    pp->mx = pp->my = pp->mz = 1;
    pp->P150 = 0;
    pp->cms_flag = p->cms_flag;
    pp->dx = p->dx/double(R[0]);

    pp->originx = node(l,0,c.lo[0]); pp->endx = node(l,0,c.hi[0]+1);
    pp->originy = node(l,1,c.lo[1]); pp->endy = node(l,1,c.hi[1]+1);
    pp->originz = node(l,2,c.lo[2]); pp->endz = node(l,2,c.hi[2]+1);
    pp->global_xmin = p->global_xmin; pp->global_ymin = p->global_ymin; pp->global_zmin = p->global_zmin;
    pp->global_xmax = p->global_xmax; pp->global_ymax = p->global_ymax; pp->global_zmax = p->global_zmax;
    pp->global_orig_x = p->global_orig_x; pp->global_orig_y = p->global_orig_y; pp->alpha_grid = p->alpha_grid;
    pp->gridgeo = p->gridgeo;
    pp->solidread = 0;
    pp->toporead = 0;
    pp->maxlength = p->maxlength;
    pp->wd = p->wd;
    pp->phimean = p->phimean;
    pp->phiout = p->phiout;

    pp->assign_margin();

    // ---- arrays of read_grid
    const int m = pp->margin;
    auto idx = [&](lexer *q, int ii, int jj, int kk) { return (ii-q->imin)*q->jmax*q->kmax + (jj-q->jmin)*q->kmax + kk-q->kmin; };
    auto lij = [&](lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); };

    pp->Iarray(pp->flag4,pp->imax*pp->jmax*pp->kmax);
    pp->Darray(pp->flag_solid,pp->imax*pp->jmax*pp->kmax);
    pp->Darray(pp->flag_topo,pp->imax*pp->jmax*pp->kmax);
    pp->Darray(pp->solidbed,pp->imax*pp->jmax);
    pp->Darray(pp->topobed,pp->imax*pp->jmax);
    pp->Darray(pp->geobed,pp->imax*pp->jmax);
    pp->Darray(pp->bed,pp->imax*pp->jmax);
    pp->Iarray(pp->wet,pp->imax*pp->jmax);
    pp->Iarray(pp->wet_n,pp->imax*pp->jmax);
    pp->Iarray(pp->deep,pp->imax*pp->jmax);
    pp->Darray(pp->depth,pp->imax*pp->jmax);
    pp->Darray(pp->WL,pp->imax*pp->jmax);
    pp->Darray(pp->data,pp->imax*pp->jmax);
    pp->Iarray(pp->flagslice1,pp->imax*pp->jmax);
    pp->Iarray(pp->flagslice2,pp->imax*pp->jmax);
    pp->Iarray(pp->flagslice4,pp->imax*pp->jmax);

    // level-0 cell below a patch cell (local index of level 0)
    auto below = [&](int ii, int jj, int kk, int *L)
    {
        const int a[3] = {ii,jj,kk};
        for(int d=0; d<3; ++d)
        L[d] = fdiv(c.lo[d]+a[d],R[d]) - org[d];
    };

    for(int ii=-m; ii<pp->knox+m; ++ii)
    for(int jj=-m; jj<pp->knoy+m; ++jj)
    for(int kk=-m; kk<pp->knoz+m; ++kk)
    {
        int L[3];
        below(ii,jj,kk,L);
        const int q = idx(pp,ii,jj,kk);
        pp->flag4[q] = p->flag4[idx(p,L[0],L[1],L[2])];
        pp->flag_solid[q] = 1.0e8;
        pp->flag_topo[q] = 1.0e8;
    }

    for(int ii=-m; ii<pp->knox+m; ++ii)
    for(int jj=-m; jj<pp->knoy+m; ++jj)
    {
        int L[3];
        below(ii,jj,0,L);
        const int q = lij(pp,ii,jj), q0 = lij(p,L[0],L[1]);
        pp->flagslice1[q] = p->flagslice1[q0];
        pp->flagslice2[q] = p->flagslice2[q0];
        pp->flagslice4[q] = p->flagslice4[q0];
        pp->solidbed[q] = p->solidbed[q0];
        pp->topobed[q] = p->topobed[q0];
        pp->geobed[q] = p->geobed[q0];
        pp->bed[q] = p->bed[q0];
        pp->wet[q] = p->wet[q0];
        pp->wet_n[q] = p->wet_n[q0];
        pp->deep[q] = p->deep[q0];
        pp->depth[q] = p->depth[q0];
        pp->WL[q] = p->WL[q0];
        pp->data[q] = 0.0;
    }

    // boundary surfaces: the fluid cells next to a non-fluid cell, the group of level 0 below
    unordered_map<long long,int> grp;
    auto key = [&](int ii, int jj, int kk, int side) { return (((long long)(ii+16)*65536 + (jj+16))*65536 + (kk+16))*8 + side; };
    for(int q=0; q<p->gcb4_count; ++q)
    grp[key(p->gcb4[q][0],p->gcb4[q][1],p->gcb4[q][2],p->gcb4[q][3])] = p->gcb4[q][4];

    const int off[7][3] = {{0,0,0},{-1,0,0},{0,1,0},{0,-1,0},{1,0,0},{0,0,-1},{0,0,1}};
    vector<int> surf;
    pp->gcin_count = pp->gcout_count = 0;

    for(int ii=0; ii<pp->knox; ++ii)
    for(int jj=0; jj<pp->knoy; ++jj)
    for(int kk=0; kk<pp->knoz; ++kk)
    {
        if(pp->flag4[idx(pp,ii,jj,kk)]<0)
        continue;
        for(int cs=1; cs<=6; ++cs)
        {
            const int a = ii+off[cs][0], b = jj+off[cs][1], e = kk+off[cs][2];
            if(pp->flag4[idx(pp,a,b,e)]>0)
            continue;
            int L[3];
            below(ii,jj,kk,L);
            auto it = grp.find(key(L[0],L[1],L[2],cs));
            int group = (it!=grp.end()) ? it->second : 21;
            surf.push_back(ii); surf.push_back(jj); surf.push_back(kk); surf.push_back(cs); surf.push_back(group);
            if(group==1 || group==6)
            ++pp->gcin_count;
            if(group==2 || group==7 || group==8)
            ++pp->gcout_count;
        }
    }

    const int ns = (int)surf.size()/5;
    pp->gcwall_count = ns;
    pp->gcb1_count = pp->gcb2_count = pp->gcb3_count = pp->gcb4_count = pp->gcb4a_count = ns;
    pp->gcb_fix = pp->gcb_solid = pp->gcb_topo = pp->gcb_fb = ns;
    pp->solid_gcb_est = pp->topo_gcb_est = 0;
    pp->solid_gcbextra_est = pp->topo_gcbextra_est = pp->tot_gcbextra_est = 0;
    pp->gcextra4 = 0;

    pp->gcpara1_count=pp->gcpara2_count=pp->gcpara3_count=pp->gcpara4_count=pp->gcpara5_count=pp->gcpara6_count=0;
    pp->gcparaco1_count=pp->gcparaco2_count=pp->gcparaco3_count=pp->gcparaco4_count=pp->gcparaco5_count=pp->gcparaco6_count=0;
    pp->gcslpara1_count=pp->gcslpara2_count=pp->gcslpara3_count=pp->gcslpara4_count=0;
    pp->gcslparaco1_count=pp->gcslparaco2_count=pp->gcslparaco3_count=pp->gcslparaco4_count=0;
    pp->gcpara_sum = 0;
    pp->gcparaco_sum = 0;

    if(ns>0)
    {
    pp->Iarray(pp->gcb1,ns,6);
    pp->Iarray(pp->gcb2,ns,6);
    pp->Iarray(pp->gcb3,ns,6);
    pp->Iarray(pp->gcb4,ns,6);
    pp->Iarray(pp->gcb4a,ns,6);
    pp->Darray(pp->gcd1,ns);
    pp->Darray(pp->gcd2,ns);
    pp->Darray(pp->gcd3,ns);
    pp->Darray(pp->gcd4,ns);
    pp->Darray(pp->gcd4a,ns);
    }

    for(int n=0; n<ns; ++n)
    for(int q=0; q<5; ++q)
    pp->gcb4[n][q] = surf[5*n+q];

    pp->Iarray(pp->gcin,pp->gcin_count,6);
    pp->Iarray(pp->gcout,pp->gcout_count,6);

    pp->Iarray(pp->gcpara1,0,16); pp->Iarray(pp->gcpara2,0,16); pp->Iarray(pp->gcpara3,0,16);
    pp->Iarray(pp->gcpara4,0,16); pp->Iarray(pp->gcpara5,0,16); pp->Iarray(pp->gcpara6,0,16);
    pp->Iarray(pp->gcparaco1,0,3); pp->Iarray(pp->gcparaco2,0,3); pp->Iarray(pp->gcparaco3,0,3);
    pp->Iarray(pp->gcparaco4,0,3); pp->Iarray(pp->gcparaco5,0,3); pp->Iarray(pp->gcparaco6,0,3);

    pp->gcbsl1_count=pp->gcbsl2_count=pp->gcbsl3_count=pp->gcbsl4_count=pp->gcbsl4a_count=1;
    pp->Iarray(pp->gcbsl1,1,6); pp->Iarray(pp->gcbsl2,1,6); pp->Iarray(pp->gcbsl3,1,6);
    pp->Iarray(pp->gcbsl4,1,6); pp->Iarray(pp->gcbsl4a,1,6);
    pp->Iarray(pp->gcslin,pp->gcin_count,6);
    pp->Iarray(pp->gcslout,pp->gcout_count,6);
    pp->Iarray(pp->gcslpara1,0,2); pp->Iarray(pp->gcslpara2,0,2); pp->Iarray(pp->gcslpara3,0,2); pp->Iarray(pp->gcslpara4,0,2);
    pp->Iarray(pp->gcslparaco1,0,4); pp->Iarray(pp->gcslparaco2,0,4); pp->Iarray(pp->gcslparaco3,0,4); pp->Iarray(pp->gcslparaco4,0,4);

    // nodes
    pp->Darray(pp->XN,pp->knox+1+2*marge);
    pp->Darray(pp->YN,pp->knoy+1+2*marge);
    pp->Darray(pp->ZN,pp->knoz+1+2*marge);
    for(int ii=-marge; ii<pp->knox+1+marge; ++ii)
    pp->XN[ii+marge] = node(l,0,c.lo[0]+ii);
    for(int jj=-marge; jj<pp->knoy+1+marge; ++jj)
    pp->YN[jj+marge] = node(l,1,c.lo[1]+jj);
    for(int kk=-marge; kk<pp->knoz+1+marge; ++kk)
    pp->ZN[kk+marge] = node(l,2,c.lo[2]+kk);

    // ---- gridini (without the solids of the grid file)
    pp->gridspacing(pgc0);
    pp->DXM = p->DXM/double(R[0]);
    pp->DXD = p->DXD/double(R[0]);
    pp->DYD = p->DYD/double(R[1]);
    pp->DYM = p->DYM/double(R[1]);
    pp->DZM = p->DZM/double(R[2]);
    pp->gcd_ini(pgc0);

    // velocity loops: the last face of a direction is a boundary face only at the domain end
    pp->ulast = (c.hi[0]==gn(l,0)-1) ? p->ulast : 0;
    pp->vlast = (c.hi[1]==gn(l,1)-1) ? p->vlast : 0;
    pp->wlast = (c.hi[2]==gn(l,2)-1) ? p->wlast : 0;
    pp->flast = p->flast;
    pp->ulastsflow = p->ulastsflow;

    // ---- flagini, flagfield, makegrid, makegrid2D
    pp->flagini();
    pgc0->flagfield(pp);

    pgc0->flagx(pp,pp->flag1);
    pgc0->flagx(pp,pp->flag2);
    pgc0->flagx(pp,pp->flag3);
    pgc0->flagx(pp,pp->flag4);
    pgc0->gcxupdate(pp);

    pp->vecsize(pgc0);

    {
    grid_helper gridgen(pp);
    gridgen.fillgcb1(pp);
    gridgen.fillgcb2(pp);
    gridgen.fillgcb3(pp);
    gridgen.fillgcb4a(pp);
    gridgen.fillgcb4_wall(pp);
    gridgen.make_dgc(pp);
    gridgen.fill_dgc1(pp);
    gridgen.fill_dgc2(pp);
    gridgen.fill_dgc3(pp);
    gridgen.fill_dgc4(pp);
    }

    pgc0->gcslflagx(pp,pp->flagslice4);
    {
    mgcslice1 m1(pp);
    mgcslice2 m2(pp);
    mgcslice4 m4(pp);
    m1.makemgc(pp);
    pgc0->gcslflagx(pp,pp->flagslice1);
    m1.gcb_seed(pp);
    m2.makemgc(pp);
    pgc0->gcslflagx(pp,pp->flagslice2);
    m2.gcb_seed(pp);
    m4.makemgc(pp);
    m4.gcb_seed(pp);
    }
    pgc0->gcsl_setbc1(pp);
    pgc0->gcsl_setbc2(pp);
    pgc0->gcsl_setbc4(pp);
    pgc0->gcsl_setbcio(pp);
    pgc0->dgcslini1(pp);
    pgc0->dgcslini2(pp);
    pgc0->dgcslini4(pp);

    // interface thickness of the level (F 45 x the spacing of the level), as initialize::inipsi;
    // the fluid properties take the local width of the grid (interface_width.h), so per level
    if(p->j_dir==0)
    pp->psi = p->F45*0.5*(pp->DXM+pp->DZM);
    else
    pp->psi = p->F45*(1.0/3.0)*(pp->DXM+pp->DYM+pp->DZM);
    pp->psi0 = pp->psi;

    pp->dt = p->dt;
    pp->dt_old = p->dt_old;
    pp->simtime = p->simtime;
    pp->count = p->count;
}

// --------------------------------------------------------------------- plans
void reefamr3d::plans()
{
    xp.assign(maxlev+1,r3xplan());
    const int m = marge>=3 ? 3 : marge;     // the CFD margin: 3 cells
    (void)m;

    for(int l=1; l<=maxlev; ++l)
    {
        vector<vector<int>> req(nranks);    // per rank: lid, s0, s1, s2
        vector<vector<int>> reqpos(nranks); // the fill entries waiting for them: patch id, entry

        for(int id : lev[l])
        {
            r3patch &c = *P[id];
            lexer *pp = c.pp;
            const int mg = pp->margin;
            c.fill.clear();

            // parent blocks
            c.par.clear();
            {
                int clo[3], chi[3];
                for(int d=0; d<3; ++d)
                {
                    clo[d] = fdiv(c.lo[d],rr[d]);
                    chi[d] = fdiv(c.hi[d],rr[d]);
                }
                auto add = [&](int g, const int *a, const int *b)
                {
                    r3patch::pblock B;
                    B.g = g;
                    for(int d=0; d<3; ++d)
                    {
                        B.lo[d] = MAX(clo[d],a[d]);
                        B.hi[d] = MIN(chi[d],b[d]);
                        if(B.lo[d]>B.hi[d])
                        return;
                    }
                    c.par.push_back(B);
                };
                if(l==1)
                {
                    int a[3], b[3];
                    for(int d=0; d<3; ++d)
                    {
                        a[d] = org[d];
                        b[d] = org[d]+kn0[d]-1;
                    }
                    add(-1,a,b);
                }
                else
                for(int g : lev[l-1])
                add(g,P[g]->lo,P[g]->hi);
            }

            for(int ii=-mg; ii<pp->knox+mg; ++ii)
            for(int jj=-mg; jj<pp->knoy+mg; ++jj)
            for(int kk=-mg; kk<pp->knoz+mg; ++kk)
            {
                if(ii>=0 && ii<pp->knox && jj>=0 && jj<pp->knoy && kk>=0 && kk<pp->knoz)
                continue;

                r3fill f;
                f.d[0]=ii; f.d[1]=jj; f.d[2]=kk;
                f.kind = 3;
                f.g = -1;
                f.slot = -1;
                for(int d=0; d<3; ++d)
                f.s[d] = f.o[d] = 0;

                const int G[3] = {c.lo[0]+ii, c.lo[1]+jj, c.lo[2]+kk};
                if(!in_domain(l,G))
                {
                    c.fill.push_back(f);
                    continue;
                }

                const int gp = gpatch_at(l,G);
                if(gp>=0)
                {
                    const r3gpatch &q = GP[gp];
                    for(int d=0; d<3; ++d)
                    f.s[d] = G[d]-q.lo[d];
                    if(q.rank==myrank)
                    {
                        f.kind = 0;
                        f.g = q.lid;
                    }
                    else
                    {
                        f.kind = 2;
                        req[q.rank].push_back(q.lid);
                        req[q.rank].push_back(f.s[0]);
                        req[q.rank].push_back(f.s[1]);
                        req[q.rank].push_back(f.s[2]);
                        reqpos[q.rank].push_back(id);
                        reqpos[q.rank].push_back((int)c.fill.size());
                    }
                    c.fill.push_back(f);
                    continue;
                }

                // the parent, on a local grid of level l-1 whose array holds it with one more cell
                int Gp[3];
                for(int d=0; d<3; ++d)
                {
                    Gp[d] = fdiv(G[d],rr[d]);
                    f.o[d] = G[d]-Gp[d]*rr[d];
                }
                f.kind = 1;

                if(l==1)
                {
                    f.g = -1;
                    for(int d=0; d<3; ++d)
                    f.s[d] = Gp[d]-org[d];
                }
                else
                {
                    int best = -1, bestscore = -1;
                    for(int g : lev[l-1])
                    {
                        r3patch &q = *P[g];
                        int score = 2;
                        for(int d=0; d<3; ++d)
                        {
                            if(Gp[d]<q.lo[d]-q.pp->margin+1 || Gp[d]>q.hi[d]+q.pp->margin-1)
                            score = -1;
                            else if(score>=0 && (Gp[d]<q.lo[d] || Gp[d]>q.hi[d]))
                            score = 1;
                        }
                        if(score>bestscore)
                        {
                            bestscore = score;
                            best = g;
                        }
                    }
                    if(best<0)
                    {
                        cout<<"CFD AMR: no parent for a cell of a level-"<<l<<" patch on rank "<<myrank<<endl;
                        MPI_Abort(MPI_COMM_WORLD,-3101);
                    }
                    f.g = best;
                    for(int d=0; d<3; ++d)
                    f.s[d] = Gp[d]-P[best]->lo[d];
                }
                c.fill.push_back(f);
            }
        }

        // remote plan of the level
        r3xplan &X = xp[l];
        X.rcnt.assign(nranks,0);
        X.scnt.assign(nranks,0);
        for(int r=0; r<nranks; ++r)
        X.rcnt[r] = (int)req[r].size()/4;

        MPI_Alltoall(&X.rcnt[0],1,MPI_INT,&X.scnt[0],1,MPI_INT,MPI_COMM_WORLD);

        X.rdsp.assign(nranks,0);
        X.sdsp.assign(nranks,0);
        for(int r=1; r<nranks; ++r)
        {
            X.rdsp[r] = X.rdsp[r-1]+X.rcnt[r-1];
            X.sdsp[r] = X.sdsp[r-1]+X.scnt[r-1];
        }
        X.nrecv = X.rdsp[nranks-1]+X.rcnt[nranks-1];
        X.nsend = X.sdsp[nranks-1]+X.scnt[nranks-1];

        vector<int> sreq(4*X.nrecv), rreq(4*X.nsend+4);
        vector<int> c4(nranks), d4(nranks), c4r(nranks), d4r(nranks);
        for(int r=0; r<nranks; ++r)
        {
            for(size_t q=0; q<req[r].size(); ++q)
            sreq[4*X.rdsp[r]+q] = req[r][q];
            c4[r] = 4*X.rcnt[r]; d4[r] = 4*X.rdsp[r];
            c4r[r] = 4*X.scnt[r]; d4r[r] = 4*X.sdsp[r];
        }
        MPI_Alltoallv(sreq.empty()?nullptr:&sreq[0],&c4[0],&d4[0],MPI_INT,&rreq[0],&c4r[0],&d4r[0],MPI_INT,MPI_COMM_WORLD);

        X.serve.resize(X.nsend);
        for(int q=0; q<X.nsend; ++q)
        {
            X.serve[q].g = rreq[4*q];
            X.serve[q].s[0] = rreq[4*q+1];
            X.serve[q].s[1] = rreq[4*q+2];
            X.serve[q].s[2] = rreq[4*q+3];
        }

        for(int r=0; r<nranks; ++r)
        for(size_t q=0; q<reqpos[r].size()/2; ++q)
        P[reqpos[r][2*q]]->fill[reqpos[r][2*q+1]].slot = X.rdsp[r]+(int)q;
    }
}

void reefamr3d::xrun(r3xplan &X, int nv)
{
    vector<int> sc(nranks), sd(nranks), rc(nranks), rd(nranks);
    for(int r=0; r<nranks; ++r)
    {
        sc[r] = X.scnt[r]*nv; sd[r] = X.sdsp[r]*nv;
        rc[r] = X.rcnt[r]*nv; rd[r] = X.rdsp[r]*nv;
    }
    double dummy = 0.0;
    MPI_Alltoallv(X.sbuf.empty()?&dummy:&X.sbuf[0],&sc[0],&sd[0],MPI_DOUBLE,
                  X.rbuf.empty()?&dummy:&X.rbuf[0],&rc[0],&rd[0],MPI_DOUBLE,MPI_COMM_WORLD);
}
