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
#include"patchBC_codes.h"
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
#include<functional>
#include<map>
#include<cstdlib>

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

bool reefamr3d::wrap(int l, const int *I, int *J) const
{
    bool in = true;
    for(int d=0; d<3; ++d)
    {
        const int n = gn(l,d);
        J[d] = I[d];
        if(per[d] && (J[d]<0 || J[d]>=n))
        J[d] = ((J[d]%n)+n)%n;
        if(J[d]<0 || J[d]>=n)
        in = false;
    }
    return in;
}

// a level-0 cell index wrapped across the periodic edges (the other directions unchanged, also
// outside the domain): the cell that holds the data of a shadow cell (G 40)
void reefamr3d::wrap0(int *I) const
{
    for(int d=0; d<3; ++d)
    if(per[d])
    {
        const int n = gn(0,d);
        I[d] = ((I[d]%n)+n)%n;
    }
}

bool reefamr3d::wrap_local0(const int *I, int *L) const
{
    int J[3];
    if(!wrap(0,I,J))
    return false;
    for(int d=0; d<3; ++d)
    {
        L[d] = J[d]-org[d];
        if(L[d]<0 || L[d]>=kn0[d])
        return false;
    }
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

int reefamr3d::owner0(const int *I) const
{
    const int a = MIN(MAX(I[0],0),gn0[0]-1), b = MIN(MAX(I[1],0),gn0[1]-1);
    return col0[(size_t)a*gn0[1]+b];
}

double reefamr3d::load_w(int l) const
{
    // all levels take the step of the finest one (G 7 0), or level l its 2^l steps (G 7 1); a patch
    // cell the work wpatch of a level-0 cell (G 41)
    return ((p0->G7==1) ? double(1<<l) : 1.0)*(l>0 ? wpatch : 1.0);
}

// a copy plan from the items of this rank (per peer: grid, index here; the level-0 cell there; tag)
void reefamr3d::cplan_build(r3cplan &X, bool down, vector<vector<int>> &items)
{
    X = r3cplan();
    vector<int> scnt(nranks,0), rcnt(nranks,0);
    for(int r=0; r<nranks; ++r)
    scnt[r] = (int)items[r].size()/6;
    MPI_Alltoall(scnt.data(),1,MPI_INT,rcnt.data(),1,MPI_INT,MPI_COMM_WORLD);

    // the far side: the level-0 cells and tags of the items of every peer
    vector<int> sc4(nranks), sd4(nranks), rc4(nranks), rd4(nranks);
    int ns=0, nr=0;
    for(int r=0; r<nranks; ++r)
    {
        sc4[r] = 4*scnt[r]; sd4[r] = ns; ns += sc4[r];
        rc4[r] = 4*rcnt[r]; rd4[r] = nr; nr += rc4[r];
    }
    vector<int> sb(ns+1), rb(nr+1);
    for(int r=0; r<nranks; ++r)
    for(int q=0; q<scnt[r]; ++q)
    for(int e=0; e<4; ++e)
    sb[sd4[r]+4*q+e] = items[r][6*q+2+e];
    MPI_Alltoallv(sb.data(),sc4.data(),sd4.data(),MPI_INT,rb.data(),rc4.data(),rd4.data(),MPI_INT,MPI_COMM_WORLD);

    for(int r=0; r<nranks; ++r)
    {
        // my items with peer r
        if(scnt[r]>0)
        {
            vector<int> g(scnt[r]), q(scnt[r]), t(scnt[r]);
            for(int m=0; m<scnt[r]; ++m)
            {
                g[m] = items[r][6*m];
                q[m] = items[r][6*m+1];
                t[m] = items[r][6*m+5];
            }
            if(down)
            {
                X.rp.push_back(r);
                X.rg.push_back(g);
                X.rq.push_back(q);
                X.rt.push_back(t);
            }
            else
            {
                X.sp.push_back(r);
                X.sg.push_back(g);
                X.sq.push_back(q);
            }
        }
        // the items of peer r with my level-0 cells
        if(rcnt[r]>0)
        {
            vector<int> g(rcnt[r],-1), q(rcnt[r]), t(rcnt[r]);
            for(int m=0; m<rcnt[r]; ++m)
            {
                const int *I = &rb[rd4[r]+4*m];
                const int i = I[0]-org[0], j = I[1]-org[1], k = I[2]-org[2];
                q[m] = (i-p0->imin)*p0->jmax*p0->kmax + (j-p0->jmin)*p0->kmax + k-p0->kmin;
                t[m] = I[3];
            }
            if(down)
            {
                X.sp.push_back(r);
                X.sg.push_back(g);
                X.sq.push_back(q);
            }
            else
            {
                X.rp.push_back(r);
                X.rg.push_back(g);
                X.rq.push_back(q);
                X.rt.push_back(t);
            }
        }
    }
}

// --------------------------------------------------------------------- setup
// refinement ratios, level-0 boxes of the ranks, tiles, global nodes
void reefamr3d::setup_base(lexer *p, ghostcell *pgc)
{
    rr[0] = (p->i_dir==1 && p->gknox>1) ? 2 : 1;
    rr[1] = (p->j_dir==1 && p->gknoy>1) ? 2 : 1;
    rr[2] = (p->k_dir==1 && p->gknoz>1) ? 2 : 1;

    gn0[0] = p->gknox; gn0[1] = p->gknoy; gn0[2] = p->gknoz;
    // periodic boundaries in x and y (1: the direction on one rank, 2: across ranks); not in z
    per[0] = (rr[0]==2 && p->periodic1>0) ? p->periodic1 : 0;
    per[1] = (rr[1]==2 && p->periodic2>0) ? p->periodic2 : 0;
    per[2] = 0;
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

    // the rank of every level-0 column (no partition in z)
    col0.assign((size_t)gn0[0]*gn0[1],0);
    for(int r=0; r<nranks; ++r)
    for(int I=rbox[6*r]; I<rbox[6*r]+rbox[6*r+3]; ++I)
    for(int J=rbox[6*r+1]; J<rbox[6*r+1]+rbox[6*r+4]; ++J)
    col0[(size_t)I*gn0[1]+J] = r;

    tn.assign(maxlev+1,vector<int>(3,0));
    for(int l=1; l<=maxlev; ++l)
    for(int d=0; d<3; ++d)
    tn[l][d] = (gn(l,d)+tile[d]-1)/tile[d];
}

void reefamr3d::setup(lexer *p, ghostcell *pgc, const vector<vector<unsigned char>> *seed)
{
    if(tn.empty())
    setup_base(p,pgc);

    mark_tiles(seed);
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
    npat = (int)P.size();

    // G 40: the shadows of the level-1 patches (the level-0 data of their units)
    vector<int> fresh_sh;
    if(place)
    make_shadows(nullptr,fresh_sh);

    for(auto c : P)
    build_lexer(*c);

    for(auto c : P)
    patch_objects(c);

    plans();
}

// a new hierarchy from the seeds (adaptive refinement)
bool reefamr3d::rebuild(const vector<vector<unsigned char>> &seed, vector<r3patch*> &gone, vector<int> &fresh)
{
    gone.clear();
    fresh.clear();

    const vector<r3gpatch> GPprev = GP;
    const vector<vector<unsigned char>> told = tmask;

    mark_tiles(&seed);

    if(tmask==told)
    return false;

    return rebuild_from(GPprev,gone,fresh);
}

bool reefamr3d::replace_units(vector<r3patch*> &gone, vector<int> &fresh)
{
    gone.clear();
    fresh.clear();
    const vector<r3gpatch> GPprev = GP;
    return rebuild_from(GPprev,gone,fresh);
}

// the patches of the tiles as they are now marked (make_boxes), the old ones of the same box on
// the same rank kept
bool reefamr3d::rebuild_from(const vector<r3gpatch> &GPprev, vector<r3patch*> &gone, vector<int> &fresh)
{
    GPold = GPprev;
    make_boxes();

    bool same = (GP.size()==GPold.size());
    for(size_t g=0; same && g<GP.size(); ++g)
    {
        const r3gpatch &a = GP[g], &b = GPold[g];
        if(a.lev!=b.lev || a.rank!=b.rank)
        same = false;
        for(int d=0; d<3; ++d)
        if(a.lo[d]!=b.lo[d] || a.hi[d]!=b.hi[d])
        same = false;
    }
    if(same)
    {
        GP = GPprev;
        GPold.clear();
        return false;
    }

    // the old patches (by their local id) and the old shadows
    vector<r3patch*> Sold(P.begin()+npat,P.end());
    Pold.assign(P.begin(),P.begin()+npat);
    vector<char> kept(Pold.size(),0);

    P.clear();
    lev.assign(maxlev+1,vector<int>());
    for(size_t g=0; g<GP.size(); ++g)
    if(GP[g].rank==myrank)
    {
        r3patch *c = nullptr;
        for(size_t q=0; q<Pold.size(); ++q)
        {
            r3patch *o = Pold[q];
            if(kept[q] || o->lev!=GP[g].lev)
            continue;
            bool eq = true;
            for(int d=0; d<3; ++d)
            if(o->lo[d]!=GP[g].lo[d] || o->hi[d]!=GP[g].hi[d])
            eq = false;
            if(eq)
            {
                c = o;
                kept[q] = 1;
                break;
            }
        }

        if(c==nullptr)
        {
            c = patch_new();
            c->lev = GP[g].lev;
            for(int d=0; d<3; ++d)
            {
                c->lo[d] = GP[g].lo[d];
                c->hi[d] = GP[g].hi[d];
                c->n[d] = c->hi[d]-c->lo[d]+1;
            }
            fresh.push_back((int)P.size());
        }
        c->gid = (int)g;
        P.push_back(c);
        lev[c->lev].push_back((int)P.size()-1);
    }
    npat = (int)P.size();

    // G 40: the shadows (those of an unchanged level-1 box kept)
    vector<int> fresh_sh;
    if(place)
    make_shadows(&Sold,fresh_sh);
    for(r3patch *o : Sold)
    if(o!=nullptr)
    gone.push_back(o);

    for(int id : fresh)
    build_lexer(*P[id]);
    for(int id : fresh_sh)
    build_lexer(*P[id]);
    for(int id : fresh)
    patch_objects(P[id]);
    for(int id : fresh_sh)
    patch_objects(P[id]);

    for(size_t q=0; q<Pold.size(); ++q)
    if(!kept[q])
    gone.push_back(Pold[q]);

    plans();
    return true;
}

// G 40: the shadow of every level-1 patch: the level-0 cells under it and two around (inside the
// domain), the array margin three more; an old shadow of the same box is kept (old: set to null).
// The finer patches take the shadow of their unit (the source of their lexer).
void reefamr3d::make_shadows(vector<r3patch*> *old, vector<int> &fresh_sh)
{
    fresh_sh.clear();
    vector<int> unit_sh(GP.size(),-1);
    for(int id : lev[1])
    {
        r3patch *c = P[id];
        int lo[3], hi[3];
        for(int d=0; d<3; ++d)
        {
            const int b = (rr[d]==2) ? 2 : 0;
            lo[d] = MAX(fdiv(c->lo[d],rr[d])-b,0);
            hi[d] = MIN(fdiv(c->hi[d],rr[d])+b,gn0[d]-1);
        }
        r3patch *q = nullptr;
        if(old!=nullptr)
        for(r3patch *&o : *old)
        if(o!=nullptr && o->lo[0]==lo[0] && o->lo[1]==lo[1] && o->lo[2]==lo[2] && o->hi[0]==hi[0] && o->hi[1]==hi[1] && o->hi[2]==hi[2])
        {
            q = o;
            o = nullptr;
            break;
        }
        if(q==nullptr)
        {
            q = patch_new();
            q->shadow = true;
            q->lev = 0;
            q->gid = -1;
            for(int d=0; d<3; ++d)
            {
                q->lo[d] = lo[d];
                q->hi[d] = hi[d];
                q->n[d] = hi[d]-lo[d]+1;
            }
            fresh_sh.push_back((int)P.size());
        }
        c->sh = (int)P.size();
        unit_sh[c->gid] = c->sh;
        P.push_back(q);
    }
    for(int l=2; l<=maxlev; ++l)
    for(int id : lev[l])
    P[id]->sh = unit_sh[GP[P[id]->gid].unit];

    fetch_raw(fresh_sh);
}

// G 40: the level-0 data of the arrays of the shadows ids from the ranks that hold the cells
// (collective): flag4 and the boundary groups per cell, the slices per column, the 2D nodes
void reefamr3d::fetch_raw(const vector<int> &ids)
{
    const int m = 3;
    auto lidx = [&](int i, int j, int k) { return (i-p0->imin)*p0->jmax*p0->kmax + (j-p0->jmin)*p0->kmax + k-p0->kmin; };
    auto lij = [&](int i, int j) { return (i-p0->imin)*p0->jmax + (j-p0->jmin); };

    for(int id : ids)
    {
        r3patch *c = P[id];
        c->raw = new r3raw;
        r3raw &R = *c->raw;
        for(int d=0; d<3; ++d)
        {
            R.lo[d] = c->lo[d]-m;
            R.n[d] = c->n[d]+2*m;
        }
        R.flag4.assign((size_t)R.n[0]*R.n[1]*R.n[2],-1);
        const size_t nc = (size_t)R.n[0]*R.n[1];
        R.fs1.assign(nc,0); R.fs2.assign(nc,0); R.fs4.assign(nc,0);
        R.wet.assign(nc,0); R.wetn.assign(nc,0); R.deep.assign(nc,0);
        R.sbed.assign(nc,0.0); R.tbed.assign(nc,0.0); R.gbed.assign(nc,0.0);
        R.bed.assign(nc,0.0); R.depth.assign(nc,0.0); R.WL.assign(nc,0.0);
        R.xn2d.assign((size_t)(R.n[0]+1)*(R.n[1]+1),0.0);
        R.yn2d.assign((size_t)(R.n[0]+1)*(R.n[1]+1),0.0);
    }

    // the boundary groups of the level-0 cells of this rank
    std::unordered_map<long long,int> grp;
    for(int q=0; q<p0->gcb4_count; ++q)
    grp[r3raw::gkey(p0->gcb4[q][0],p0->gcb4[q][1],p0->gcb4[q][2],p0->gcb4[q][3])] = p0->gcb4[q][4];

    // one exchange: requests of global indices (sin ints), answers of sout doubles each
    auto exchange = [&](vector<vector<int>> &req, int sin, int sout, std::function<void(const int*, double*)> serve,
                        vector<vector<double>> &ans)
    {
        vector<int> sc(nranks), rc(nranks), sd(nranks), rd(nranks);
        for(int r=0; r<nranks; ++r)
        sc[r] = (int)req[r].size();
        MPI_Alltoall(sc.data(),1,MPI_INT,rc.data(),1,MPI_INT,MPI_COMM_WORLD);
        int ns=0, nr=0;
        for(int r=0; r<nranks; ++r)
        {
            sd[r] = ns; ns += sc[r];
            rd[r] = nr; nr += rc[r];
        }
        vector<int> sb(ns+1), rb(nr+1);
        for(int r=0; r<nranks; ++r)
        std::copy(req[r].begin(),req[r].end(),sb.begin()+sd[r]);
        MPI_Alltoallv(sb.data(),sc.data(),sd.data(),MPI_INT,rb.data(),rc.data(),rd.data(),MPI_INT,MPI_COMM_WORLD);

        vector<double> ab((size_t)(nr/sin)*sout+1);
        for(int q=0; q<nr/sin; ++q)
        serve(&rb[(size_t)q*sin],&ab[(size_t)q*sout]);

        vector<int> ac(nranks), ad(nranks), bc(nranks), bd(nranks);
        for(int r=0; r<nranks; ++r)
        {
            ac[r] = rc[r]/sin*sout; ad[r] = rd[r]/sin*sout;
            bc[r] = sc[r]/sin*sout; bd[r] = sd[r]/sin*sout;
        }
        vector<double> bb((size_t)(ns/sin)*sout+1);
        MPI_Alltoallv(ab.data(),ac.data(),ad.data(),MPI_DOUBLE,bb.data(),bc.data(),bd.data(),MPI_DOUBLE,MPI_COMM_WORLD);
        ans.assign(nranks,vector<double>());
        for(int r=0; r<nranks; ++r)
        ans[r].assign(bb.begin()+bd[r],bb.begin()+bd[r]+bc[r]);
    };

    // cells: flag4, the groups of the sides 1-6 (-1 none)
    {
        vector<vector<int>> req(nranks);
        vector<vector<std::pair<int,size_t>>> back(nranks);
        for(int id : ids)
        {
            r3raw &R = *P[id]->raw;
            for(int a=0; a<R.n[0]; ++a)
            for(int b=0; b<R.n[1]; ++b)
            for(int e=0; e<R.n[2]; ++e)
            {
                int I[3] = {R.lo[0]+a,R.lo[1]+b,R.lo[2]+e};
                wrap0(I);
                const int r = owner0(I);
                req[r].insert(req[r].end(),I,I+3);
                back[r].push_back({id,((size_t)a*R.n[1]+b)*R.n[2]+e});
            }
        }
        vector<vector<double>> ans;
        exchange(req,3,7,[&](const int *I, double *v)
                 {
                     const int i = I[0]-org[0], j = I[1]-org[1], k = I[2]-org[2];
                     v[0] = p0->flag4[lidx(i,j,k)];
                     for(int s=1; s<=6; ++s)
                     {
                         auto it = grp.find(r3raw::gkey(i,j,k,s));
                         v[s] = (it!=grp.end()) ? it->second : -1;
                     }
                 },ans);
        for(int r=0; r<nranks; ++r)
        for(size_t q=0; q<back[r].size(); ++q)
        {
            r3raw &R = *P[back[r][q].first]->raw;
            const size_t c = back[r][q].second;
            const double *v = &ans[r][7*q];
            R.flag4[c] = (int)v[0];
            const int a = (int)(c/((size_t)R.n[1]*R.n[2])), b = (int)((c/R.n[2])%R.n[1]), e = (int)(c%R.n[2]);
            for(int s=1; s<=6; ++s)
            if(v[s]>=0.0)
            R.grp[r3raw::gkey(R.lo[0]+a,R.lo[1]+b,R.lo[2]+e,s)] = (int)v[s];
        }
    }

    // columns: the slices
    {
        vector<vector<int>> req(nranks);
        vector<vector<std::pair<int,size_t>>> back(nranks);
        for(int id : ids)
        {
            r3raw &R = *P[id]->raw;
            for(int a=0; a<R.n[0]; ++a)
            for(int b=0; b<R.n[1]; ++b)
            {
                int I[3] = {R.lo[0]+a,R.lo[1]+b,0};
                wrap0(I);
                const int r = owner0(I);
                req[r].insert(req[r].end(),I,I+2);
                back[r].push_back({id,(size_t)a*R.n[1]+b});
            }
        }
        vector<vector<double>> ans;
        exchange(req,2,12,[&](const int *I, double *v)
                 {
                     const int q = lij(I[0]-org[0],I[1]-org[1]);
                     v[0] = p0->flagslice1[q]; v[1] = p0->flagslice2[q]; v[2] = p0->flagslice4[q];
                     v[3] = p0->wet[q]; v[4] = p0->wet_n[q]; v[5] = p0->deep[q];
                     v[6] = p0->solidbed[q]; v[7] = p0->topobed[q]; v[8] = p0->geobed[q];
                     v[9] = p0->bed[q]; v[10] = p0->depth[q]; v[11] = p0->WL[q];
                 },ans);
        for(int r=0; r<nranks; ++r)
        for(size_t q=0; q<back[r].size(); ++q)
        {
            r3raw &R = *P[back[r][q].first]->raw;
            const size_t c = back[r][q].second;
            const double *v = &ans[r][12*q];
            R.fs1[c] = (int)v[0]; R.fs2[c] = (int)v[1]; R.fs4[c] = (int)v[2];
            R.wet[c] = (int)v[3]; R.wetn[c] = (int)v[4]; R.deep[c] = (int)v[5];
            R.sbed[c] = v[6]; R.tbed[c] = v[7]; R.gbed[c] = v[8];
            R.bed[c] = v[9]; R.depth[c] = v[10]; R.WL[c] = v[11];
        }
    }

    // the 2D nodes
    {
        vector<vector<int>> req(nranks);
        vector<vector<std::pair<int,size_t>>> back(nranks);
        for(int id : ids)
        {
            r3raw &R = *P[id]->raw;
            for(int a=0; a<=R.n[0]; ++a)
            for(int b=0; b<=R.n[1]; ++b)
            {
                const int I[3] = {R.lo[0]+a,R.lo[1]+b,0};
                const int C[3] = {MIN(I[0],gn0[0]-1),MIN(I[1],gn0[1]-1),0};
                const int r = owner0(C);
                req[r].insert(req[r].end(),I,I+2);
                back[r].push_back({id,(size_t)a*(R.n[1]+1)+b});
            }
        }
        vector<vector<double>> ans;
        exchange(req,2,2,[&](const int *I, double *v)
                 {
                     const int k = I[0]-org[0], q = I[1]-org[1];
                     if(p0->XN2D!=nullptr && k>=p0->imin && k<=p0->imin+p0->imax && q>=p0->jmin && q<=p0->jmin+p0->jmax)
                     {
                         v[0] = p0->XN2D[p0->nij(k,q)];
                         v[1] = p0->YN2D[p0->nij(k,q)];
                     }
                     else
                     {
                         v[0] = node(0,0,I[0]);
                         v[1] = node(0,1,I[1]);
                     }
                 },ans);
        for(int r=0; r<nranks; ++r)
        for(size_t q=0; q<back[r].size(); ++q)
        {
            r3raw &R = *P[back[r][q].first]->raw;
            R.xn2d[back[r][q].second] = ans[r][2*q];
            R.yn2d[back[r][q].second] = ans[r][2*q+1];
        }
    }
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

// refined tiles of every level: the seeds, the boxes G 15 on level G 1, every coarser level around
// the next finer one with a nesting band; no refinement in the wave generation and absorption zones
// (unless the model relaxes its patches there: zone_gen_refine, zone_beach_refine)
void reefamr3d::mark_tiles(const vector<vector<unsigned char>> *seed)
{
    tmask.assign(maxlev+1,vector<unsigned char>());

    for(int l=1; l<=maxlev; ++l)
    {
        tmask[l].assign((size_t)tn[l][0]*tn[l][1]*tn[l][2],0);
        if(seed!=nullptr && (int)seed->size()>l && (*seed)[l].size()==tmask[l].size())
        tmask[l] = (*seed)[l];
    }

    // the tiles of the cells lo..hi of level l (clipped to the domain; wrapped across a periodic
    // edge)
    auto mark_cells = [&](int l, const int *lo, const int *hi)
    {
        vector<int> t[3];
        for(int d=0; d<3; ++d)
        {
            if(lo[d]>hi[d])
            return;
            const int n = gn(l,d);
            if(!per[d])
            {
                if(hi[d]<0 || lo[d]>n-1)
                return;
                for(int u=MAX(lo[d],0)/tile[d]; u<=MIN(hi[d],n-1)/tile[d]; ++u)
                t[d].push_back(u);
                continue;
            }
            vector<char> on(tn[l][d],0);
            const int e = MIN(hi[d],lo[d]+n-1);
            for(int I=lo[d]; I<=e; ++I)
            on[(((I%n)+n)%n)/tile[d]] = 1;
            for(int u=0; u<tn[l][d]; ++u)
            if(on[u])
            t[d].push_back(u);
        }
        for(int ti : t[0])
        for(int tj : t[1])
        for(int tk : t[2])
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

    // wave generation and absorption zones (relaxation, B 98 2, B 99 1/2) stay on level 0, unless
    // the model relaxes its patches there
    double zx0 = 1.0e20, zx1 = -1.0e20;     // forbidden: x < zx0 and x > zx1
    const double xmin = gx[0][marge], xmax = gx[0][gn0[0]+marge];
    zx0 = -1.0e20; zx1 = 1.0e20;
    // (refined zones: the tiles at the in- and outflow edge stay on level 0, the patches have no
    // in- or outflow boundary; the relaxation there takes the target nearly fully)
    const double zeps = 1.0e-9*(xmax-xmin);
    if(p0->B98==2 && p0->B96_1>0.0)
    zx0 = zone_gen_refine ? xmin + zeps : xmin + p0->B96_1;
    if((p0->B99==1 || p0->B99==2) && p0->B96_2>0.0)
    zx1 = zone_beach_refine ? xmax - zeps : xmax - p0->B96_2;

    // in- and outflow boundaries (wave generation and absorption at the boundary, B 98 3/4,
    // B 99 3-5, the wavemakers, a current): the tiles whose patch with its margin of 3 fine cells
    // would reach the face stay on level 0, so the patches have no in- or outflow boundary.  The
    // faces of the domain with such boundary cells (gcb4 type 1, 6 inflow, 2, 7, 8 outflow; side
    // 1 x-, 4 x+, 3 y-, 2 y+, 5 z-, 6 z+)
    int ioface[6] = {0,0,0,0,0,0};      // x-, x+, y-, y+, z-, z+
    for(int n=0; n<p0->gcb4_count; ++n)
    {
        const int t = p0->gcb4[n][4], s = p0->gcb4[n][3];
        if(!(t==1 || t==2 || t==6 || t==7 || t==8))
        continue;
        const int f = (s==1) ? 0 : (s==4) ? 1 : (s==3) ? 2 : (s==2) ? 3 : (s==5) ? 4 : (s==6) ? 5 : -1;
        if(f>=0)
        ioface[f] = 1;
    }
    MPI_Allreduce(MPI_IN_PLACE,ioface,6,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    auto io_near = [&](int l, int ti, int tj, int tk) -> bool
    {
        const int t[3] = {ti,tj,tk};
        for(int d=0; d<3; ++d)
        {
            if(rr[d]!=2)
            continue;
            if(ioface[2*d] && t[d]*tile[d]-3<0)
            return true;
            if(ioface[2*d+1] && MIN((t[d]+1)*tile[d],gn(l,d))+3>gn(l,d))
            return true;
        }
        return false;
    };
    // the faces of the patch boundary conditions (B 440-442, the codes of patchBC_ini on the wall faces
    // of level 0, anywhere in the domain): the tiles whose patch with its margin of 3 fine cells would
    // reach the fluid cell of such a face stay on level 0 as well
    vector<vector<char>> pbc(maxlev+1);
    {
        int any = 0;
        for(int n=0; n<p0->gcb4_count; ++n)
        if(patch_bc(p0->gcb4[n][4]))
        any = 1;
        MPI_Allreduce(MPI_IN_PLACE,&any,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
        for(int l=1; l<=maxlev && any; ++l)
        {
            pbc[l].assign((size_t)tn[l][0]*tn[l][1]*tn[l][2],0);
            for(int n=0; n<p0->gcb4_count; ++n)
            {
                if(!patch_bc(p0->gcb4[n][4]))
                continue;
                int a[3], b[3];
                for(int d=0; d<3; ++d)
                {
                    const int I = p0->gcb4[n][d]+org[d], R = ratio(l,d);
                    const int m = (rr[d]==2) ? 3 : 0;
                    a[d] = MAX(I*R-m,0)/tile[d];
                    b[d] = MIN((I+1)*R-1+m,gn(l,d)-1)/tile[d];
                }
                for(int ti=a[0]; ti<=b[0]; ++ti)
                for(int tj=a[1]; tj<=b[1]; ++tj)
                for(int tk=a[2]; tk<=b[2]; ++tk)
                pbc[l][((size_t)ti*tn[l][1] + tj)*tn[l][2] + tk] = 1;
            }
            MPI_Allreduce(MPI_IN_PLACE,&pbc[l][0],(int)pbc[l].size(),MPI_CHAR,MPI_MAX,MPI_COMM_WORLD);
        }
    }

    // forbidden tiles, coarse to fine: at the face, or with their nesting band (nest cells of the
    // coarser level) on a forbidden tile of the coarser level - so the nesting of the tiles left
    // never asks for a forbidden one
    vector<vector<char>> fbd(maxlev+1);
    for(int l=1; l<=maxlev; ++l)
    {
        fbd[l].assign((size_t)tn[l][0]*tn[l][1]*tn[l][2],0);
        for(int ti=0; ti<tn[l][0]; ++ti)
        for(int tj=0; tj<tn[l][1]; ++tj)
        for(int tk=0; tk<tn[l][2]; ++tk)
        {
            bool f = io_near(l,ti,tj,tk) || (!pbc[l].empty() && pbc[l][((size_t)ti*tn[l][1] + tj)*tn[l][2] + tk]);
            if(!f && l>=2)
            {
                const int t[3] = {ti,tj,tk};
                int a[3], b[3];
                for(int d=0; d<3; ++d)
                {
                    const int lo = MAX(fdiv(t[d]*tile[d],rr[d]) - (rr[d]==2 ? nest : 0),0);
                    const int hi = MIN(fdiv((t[d]+1)*tile[d]-1,rr[d]) + (rr[d]==2 ? nest : 0),gn(l-1,d)-1);
                    a[d] = lo/tile[d];
                    b[d] = hi/tile[d];
                }
                for(int ui=a[0]; ui<=b[0] && !f; ++ui)
                for(int uj=a[1]; uj<=b[1] && !f; ++uj)
                for(int uk=a[2]; uk<=b[2] && !f; ++uk)
                if(fbd[l-1][((size_t)ui*tn[l-1][1] + uj)*tn[l-1][2] + uk])
                f = true;
            }
            fbd[l][((size_t)ti*tn[l][1] + tj)*tn[l][2] + tk] = f ? 1 : 0;
        }
    }
    auto io_tile = [&](int l, int ti, int tj, int tk) -> bool
    {
        return fbd[l][((size_t)ti*tn[l][1] + tj)*tn[l][2] + tk]!=0;
    };

    int cleared = 0, cleared_io = 0;
    for(int l=1; l<=maxlev; ++l)
    for(int ti=0; ti<tn[l][0]; ++ti)
    for(int tj=0; tj<tn[l][1]; ++tj)
    for(int tk=0; tk<tn[l][2]; ++tk)
    {
        unsigned char &m = tmask[l][((size_t)ti*tn[l][1] + tj)*tn[l][2] + tk];
        if(m && io_tile(l,ti,tj,tk))
        {
            m = 0;
            ++cleared_io;
        }
    }

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
            lo[d] = fdiv(t[d]*tile[d],rr[d]) - (rr[d]==2 ? 1 : 0);
            hi[d] = fdiv((t[d]+1)*tile[d]-1,rr[d]) + (rr[d]==2 ? 1 : 0);
            if(!per[d])
            {
                lo[d] = MAX(lo[d],0);
                hi[d] = MIN(hi[d],gn(l-1,d)-1);
            }
        }
        bool ok = true;
        int I[3], J[3];
        for(I[0]=lo[0]; I[0]<=hi[0] && ok; ++I[0])
        for(I[1]=lo[1]; I[1]<=hi[1] && ok; ++I[1])
        for(I[2]=lo[2]; I[2]<=hi[2] && ok; ++I[2])
        if(!wrap(l-1,I,J) || !refined(l-1,J))
        ok = false;

        if(!ok)
        {
            tmask[l][((size_t)ti*tn[l][1] + tj)*tn[l][2] + tk] = 0;
            ++dropped;
        }
    }

    if(myrank==0 && verbose && (cleared>0 || dropped>0))
    cout<<"CFD AMR: "<<cleared<<(zone_gen_refine || zone_beach_refine ? " tiles at the in- or outflow edge of the wave zones and " : " tiles in the wave zones and ")
        <<dropped<<" tiles without nesting not refined"<<endl;
    if(myrank==0 && verbose && cleared_io>0)
    cout<<"CFD AMR: "<<cleared_io<<" tiles at an in- or outflow boundary or a patch boundary condition (B 440-442) not refined (the patches keep off it)"<<endl;

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

    // a periodic edge: a level is refined on both sides of it alike (the cells next to it on one side
    // refined iff those on the other side are), so no patch edge lies on it - the margin across it is
    // filled by the patches at the other end (siblings); returns the tiles to add
    auto seam = [&](int l, bool fix) -> int
    {
        vector<size_t> add;
        for(int d=0; d<3; ++d)
        {
            if(!per[d])
            continue;
            const int n = gn(l,d);
            int I[3], J[3];
            int lo[3] = {0,0,0}, hi[3];
            for(int e=0; e<3; ++e)
            hi[e] = gn(l,e)-1;
            for(int side=0; side<2; ++side)
            {
                lo[d] = hi[d] = (side==0) ? 0 : n-1;
                for(I[0]=lo[0]; I[0]<=hi[0]; ++I[0])
                for(I[1]=lo[1]; I[1]<=hi[1]; ++I[1])
                for(I[2]=lo[2]; I[2]<=hi[2]; ++I[2])
                {
                    if(!refined(l,I))
                    continue;
                    for(int e=0; e<3; ++e)
                    J[e] = I[e];
                    J[d] = (side==0) ? n-1 : 0;
                    if(!refined(l,J))
                    add.push_back(((size_t)(J[0]/tile[0])*tn[l][1] + J[1]/tile[1])*tn[l][2] + J[2]/tile[2]);
                }
                lo[d] = 0; hi[d] = n-1;
            }
        }
        std::sort(add.begin(),add.end());
        add.erase(std::unique(add.begin(),add.end()),add.end());
        if(fix)
        for(size_t t : add)
        tmask[l][t] = 1;
        return (int)add.size();
    };

    // (G 40: the shadows give every level-1 patch its coarse cells, the finer levels are cut at
    // their units: no extension at the partition edges)
    int extended = 0, extended_seam = 0;
    for(int l=maxlev; l>=1 && (!place || per[0] || per[1]); --l)
    {
        for(int pass=0; pass<64; ++pass)
        {
            const int n = place ? 0 : edges(l,true);
            const int ns = seam(l,true);
            extended += n;
            extended_seam += ns;
            if(n+ns==0)
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
        for(int ti=0; ti<tn[l][0]; ++ti)
        for(int tj=0; tj<tn[l][1]; ++tj)
        for(int tk=0; tk<tn[l][2]; ++tk)
        if(tile_marked(l,ti,tj,tk) && io_tile(l,ti,tj,tk))
        ++bad_zone;
        if(!place)
        bad_edge += edges(l,false);
        bad_edge += seam(l,false);
    }

    if(myrank==0 && verbose && extended>0)
    cout<<"CFD AMR: "<<extended<<" tiles added: the edges of the refined region kept off the partition edges"<<endl;
    if(myrank==0 && verbose && extended_seam>0)
    cout<<"CFD AMR: "<<extended_seam<<" tiles added: the refined region alike on both sides of a periodic edge"<<endl;

    if(bad_zone>0 || bad_edge>0)
    {
        if(myrank==0)
        cout<<"CFD AMR: the refined region reaches into a wave zone or an in- or outflow boundary, or keeps an edge on a partition or periodic edge; "
              "move or resize the boxes G 15 (or the zones, the partition)"<<endl;
        MPI_Abort(MPI_COMM_WORLD,-3102);
    }
}

// tiles to boxes (greedy: runs along x, extended along y and z), cut at the rank boxes
void reefamr3d::make_boxes()
{
    GP.clear();
    glev.assign(maxlev+1,vector<int>());
    vector<vector<int>> mb(maxlev+1);       // G 40: the boxes of every level, not cut (6 per box)

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

            if(place)
            {
                for(int d=0; d<3; ++d)
                mb[l].push_back(blo[d]);
                for(int d=0; d<3; ++d)
                mb[l].push_back(bhi[d]);
                continue;
            }

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

    if(place)
    place_units(mb);

    // local ids on the owning rank, in the order of the global list
    vector<int> cnt(nranks,0);
    for(auto &q : GP)
    q.lid = cnt[q.rank]++;
    // restart (I 45) on the same decomposition: the ranks of the checkpoint (the same boxes)
    if(GPforce!=nullptr && GPforce->size()==GP.size())
    {
        bool same = true;
        for(size_t g=0; same && g<GP.size(); ++g)
        {
            if(GP[g].lev!=(*GPforce)[g].lev)
            same = false;
            for(int d=0; d<3; ++d)
            if(GP[g].lo[d]!=(*GPforce)[g].lo[d] || GP[g].hi[d]!=(*GPforce)[g].hi[d])
            same = false;
        }
        if(same)
        GP = *GPforce;
    }
}

// G 40: the work of the unit of the level-1 box lo-hi: its cells of all levels (the work per cell of
// a level, load_w); the boxes are on tile edges of level 1, the finer levels counted by their tiles
double reefamr3d::unit_work(const int *lo, const int *hi) const
{
    double w = 0.0;
    for(int l=1; l<=maxlev; ++l)
    {
        int a[3], b[3], ta[3], tb[3];
        for(int d=0; d<3; ++d)
        {
            const int R = ratio(l-1,d);
            a[d] = lo[d]*R;
            b[d] = MIN((hi[d]+1)*R,gn(l,d))-1;
            ta[d] = a[d]/tile[d];
            tb[d] = b[d]/tile[d];
        }
        double c = 0.0;
        for(int ti=ta[0]; ti<=tb[0]; ++ti)
        for(int tj=ta[1]; tj<=tb[1]; ++tj)
        for(int tk=ta[2]; tk<=tb[2]; ++tk)
        {
            if(!tile_marked(l,ti,tj,tk))
            continue;
            const int t[3] = {ti,tj,tk};
            double v = 1.0;
            for(int d=0; d<3; ++d)
            v *= double(MIN((t[d]+1)*tile[d]-1,b[d]) - MAX(t[d]*tile[d],a[d]) + 1);
            c += v;
        }
        w += load_w(l)*c;
    }
    return w;
}

// G 40: a cut of the level-1 box lo-hi at the level-1 cell index X along d (cells < X and >= X):
// no edge of a finer level on it (each finer level refined on both sides or on neither), so the
// finer patches of a unit keep their coarse cells across their edges in the unit
bool reefamr3d::conforming(int d, int X, const int *lo, const int *hi) const
{
    for(int l=2; l<=maxlev; ++l)
    {
        int ta[3], tb[3];
        for(int e=0; e<3; ++e)
        {
            const int R = ratio(l-1,e);
            ta[e] = (lo[e]*R)/tile[e];
            tb[e] = (MIN((hi[e]+1)*R,gn(l,e))-1)/tile[e];
        }
        const int XT = (X*ratio(l-1,d))/tile[d];
        ta[d] = tb[d] = XT;
        for(int ti=ta[0]; ti<=tb[0]; ++ti)
        for(int tj=ta[1]; tj<=tb[1]; ++tj)
        for(int tk=ta[2]; tk<=tb[2]; ++tk)
        {
            int u[3] = {ti,tj,tk};
            const bool right = tile_marked(l,u[0],u[1],u[2]);
            u[d] -= 1;
            const bool left = tile_marked(l,u[0],u[1],u[2]);
            if(left!=right)
            return false;
        }
    }
    return true;
}

// G 40: the units.  The work of a rank is its level-0 cells plus the work of its units (load_w).
// G 40 1 on several ranks: the boxes of level 1 (over the whole domain, not cut) in the order of
// their position fill the ranks in turn up to the mean work, a box cut at the conforming tile edge
// whose first part fills the rank best (few cuts).  G 40 2 (the test of the shadows): the boxes
// split, along the direction with the most tiles at the conforming tile edge that halves the work
// best, until a piece has at most an eighth of the mean work, placed largest first on the rank with
// the least work.  One rank: the boxes whole.  After a rebuild (G 40 1) the pieces of the same box
// as before keep their rank unless the new placement lowers the largest work by more than
// rebalance.  The finer levels are cut at the units, on their rank.  Every rank computes the same
// list.
void reefamr3d::place_units(vector<vector<int>> &mb)
{
    struct piece { int lo[3], hi[3]; double w; int rank; };
    vector<piece> pcs;

    vector<double> base(nranks,0.0);
    double wtot = 0.0;
    for(int r=0; r<nranks; ++r)
    {
        base[r] = load_w(0)*double(rbox[6*r+3])*double(rbox[6*r+4])*double(rbox[6*r+5]);
        wtot += base[r];
    }
    if(maxlev>=1)
    for(size_t b=0; b<mb[1].size()/6; ++b)
    wtot += unit_work(&mb[1][6*b],&mb[1][6*b+3]);
    const double wmean = wtot/nranks;
    // (one rank: the boxes whole, G 40 2 the pieces of the test as well)
    double wmax = (place==2 ? 0.125*wmean : 2.0*wtot);

    // REEF3D_AMR_PLACE_TEST=1 in the environment (a test of the shadows): the boxes whole, each on
    // the rank after the one of the level-0 cell at its centre
    // REEF3D_AMR_PLACE_TEST=2: the boxes cut at the rank boxes of level 0, on the rank below (the
    // patches of G 40 0, through the shadows)
    const int ptv = (getenv("REEF3D_AMR_PLACE_TEST")!=nullptr) ? atoi(getenv("REEF3D_AMR_PLACE_TEST")) : 0;
    const bool ptest = (ptv==1 || ptv==2);
    if(ptest)
    wmax = 2.0*wtot;
    // REEF3D_AMR_PLACE_TEST=3: the pieces of G 40 2 on two ranks, whatever the ranks
    if(ptv==3)
    wmax = 0.125*0.5*wtot;

    std::function<void(piece)> split = [&](piece c)
    {
        c.w = unit_work(c.lo,c.hi);
        if(c.w>wmax)
        {
            int nt[3], ord[3] = {0,1,2};
            for(int d=0; d<3; ++d)
            nt[d] = (rr[d]==2) ? (c.hi[d]-c.lo[d])/tile[d]+1 : 1;
            std::stable_sort(ord,ord+3,[&](int a, int b) { return (a==2)!=(b==2) ? b==2 : nt[a]>nt[b]; });   // vertical cuts first
            for(int oi=0; oi<3; ++oi)
            {
                const int d = ord[oi];
                if(nt[d]<2)
                continue;
                // the work of the tile slabs along d
                vector<double> cum(nt[d]+1,0.0);
                for(int m=0; m<nt[d]; ++m)
                {
                    int a[3] = {c.lo[0],c.lo[1],c.lo[2]}, b[3] = {c.hi[0],c.hi[1],c.hi[2]};
                    a[d] = c.lo[d]+m*tile[d];
                    b[d] = MIN(a[d]+tile[d]-1,c.hi[d]);
                    cum[m+1] = cum[m] + unit_work(a,b);
                }
                vector<int> cand;
                for(int m=1; m<nt[d]; ++m)
                cand.push_back(m);
                std::stable_sort(cand.begin(),cand.end(),[&](int a, int b) { return fabs(cum[a]-0.5*cum[nt[d]]) < fabs(cum[b]-0.5*cum[nt[d]]); });
                for(int m : cand)
                {
                    const int X = c.lo[d]+m*tile[d];
                    if(!conforming(d,X,c.lo,c.hi))
                    continue;
                    piece A = c, B = c;
                    A.hi[d] = X-1;
                    B.lo[d] = X;
                    split(A);
                    split(B);
                    return;
                }
            }
        }
        c.rank = -1;
        pcs.push_back(c);
    };

    // G 40 1 on several ranks: the ranks filled in turn to the mean work, the boxes in the order of
    // their position (x, y, z), a box cut at the conforming tile edge whose first part fills the rank
    // best: few cuts (the level V-cycles of the preconditioner are per patch across the cuts between
    // the ranks; cut into pieces of a quarter of the mean work, 10 pieces of a two-level dam break
    // on 2 ranks took 70-95 BiCGStab iterations instead of 8)
    const bool fill_ranks = (place==1 && !ptest && ptv!=3 && nranks>1);
    vector<int> rcap;           // the ranks of the pieces of fill_ranks
    if(maxlev>=1 && fill_ranks)
    {
        vector<piece> bx;
        for(size_t b=0; b<mb[1].size()/6; ++b)
        {
            piece c;
            for(int d=0; d<3; ++d)
            {
                c.lo[d] = mb[1][6*b+d];
                c.hi[d] = mb[1][6*b+3+d];
            }
            c.rank = -1;
            bx.push_back(c);
        }
        std::stable_sort(bx.begin(),bx.end(),[](const piece &a, const piece &b)
                         { return a.lo[0]!=b.lo[0] ? a.lo[0]<b.lo[0] : (a.lo[1]!=b.lo[1] ? a.lo[1]<b.lo[1] : a.lo[2]<b.lo[2]); });

        const double tol = 0.02*wmean;
        vector<double> cap(nranks);
        for(int q=0; q<nranks; ++q)
        cap[q] = MAX(0.0,wmean-base[q]);
        int r = 0;
        while(r<nranks-1 && cap[r]<=tol)
        ++r;
        double rem = cap[r];
        auto next = [&]()
        {
            if(r>=nranks-1)
            return;
            ++r;
            while(r<nranks-1 && cap[r]<=tol)
            ++r;
            rem = cap[r];
        };
        auto put = [&](piece c)
        {
            c.w = unit_work(c.lo,c.hi);
            pcs.push_back(c);
            rcap.push_back(r);
            rem -= c.w;
        };

        for(piece c : bx)
        for(;;)
        {
            const double w = unit_work(c.lo,c.hi);
            if(w<=rem+tol || r==nranks-1)
            {
                put(c);
                if(rem<=tol)
                next();
                break;
            }
            // the cut whose first part comes closest to the rest of the rank: vertical (x or y), a
            // horizontal one only without (a cut across the columns weakens the patch V-cycles of
            // the preconditioner: the dam break cut at the free surface took 21 BiCGStab
            // iterations instead of 7)
            int bd = -1, bX = 0;
            double berr = 1.0e300;
            for(int d=0; d<3; ++d)
            {
                if(rr[d]!=2 || (d==2 && bd>=0))
                continue;
                const int nt = (c.hi[d]-c.lo[d])/tile[d]+1;
                double cum = 0.0;
                for(int m=1; m<nt; ++m)
                {
                    int a[3] = {c.lo[0],c.lo[1],c.lo[2]}, b[3] = {c.hi[0],c.hi[1],c.hi[2]};
                    a[d] = c.lo[d]+(m-1)*tile[d];
                    b[d] = c.lo[d]+m*tile[d]-1;
                    cum += unit_work(a,b);
                    const int X = c.lo[d]+m*tile[d];
                    if(fabs(cum-rem)<berr && conforming(d,X,c.lo,c.hi))
                    {
                        berr = fabs(cum-rem);
                        bd = d;
                        bX = X;
                    }
                }
            }
            if(bd<0)
            {
                put(c);
                if(rem<=tol)
                next();
                break;
            }
            piece A = c, B = c;
            A.hi[bd] = bX-1;
            B.lo[bd] = bX;
            put(A);
            next();
            c = B;
        }
    }
    else if(maxlev>=1 && ptv==2)
    for(size_t b=0; b<mb[1].size()/6; ++b)
    for(int r=0; r<nranks; ++r)
    {
        piece c;
        bool empty = false;
        for(int d=0; d<3; ++d)
        {
            const int R = ratio(1,d);
            c.lo[d] = MAX(mb[1][6*b+d],rbox[6*r+d]*R);
            c.hi[d] = MIN(mb[1][6*b+3+d],(rbox[6*r+d]+rbox[6*r+3+d])*R-1);
            if(c.lo[d]>c.hi[d])
            empty = true;
        }
        if(empty)
        continue;
        c.w = unit_work(c.lo,c.hi);
        c.rank = r;
        pcs.push_back(c);
    }
    else if(maxlev>=1)
    for(size_t b=0; b<mb[1].size()/6; ++b)
    {
        piece c;
        for(int d=0; d<3; ++d)
        {
            c.lo[d] = mb[1][6*b+d];
            c.hi[d] = mb[1][6*b+3+d];
        }
        split(c);
    }

    // placement: largest first
    vector<int> order(pcs.size());
    for(size_t i=0; i<pcs.size(); ++i)
    order[i] = (int)i;
    std::stable_sort(order.begin(),order.end(),[&](int a, int b) { return pcs[a].w>pcs[b].w; });

    auto greedy = [&](vector<int> &rk, vector<double> &load)
    {
        for(int i : order)
        {
            if(rk[i]>=0)
            continue;
            int r = -1;
            if(place==1)
            {
                int C[3];
                for(int d=0; d<3; ++d)
                C[d] = fdiv((pcs[i].lo[d]+pcs[i].hi[d])/2,rr[d]);
                const int pr = owner0(C);
                if(load[pr]+pcs[i].w<=1.05*wmean)
                r = pr;
            }
            if(r<0)
            {
                r = 0;
                for(int q=1; q<nranks; ++q)
                if(load[q]<load[r])
                r = q;
            }
            rk[i] = r;
            load[r] += pcs[i].w;
        }
    };

    vector<int> rf(pcs.size(),-1), rkp(pcs.size(),-1);
    vector<double> lf = base, lk = base;
    if(fill_ranks)
    for(size_t i=0; i<pcs.size(); ++i)
    {
        rf[i] = rcap[i];
        lf[rf[i]] += pcs[i].w;
    }
    if(ptest)
    for(size_t i=0; i<pcs.size(); ++i)
    {
        int C[3];
        for(int d=0; d<3; ++d)
        C[d] = fdiv((pcs[i].lo[d]+pcs[i].hi[d])/2,rr[d]);
        rf[i] = (ptv==2) ? pcs[i].rank : (owner0(C)+1)%nranks;
        lf[rf[i]] += pcs[i].w;
    }
    greedy(rf,lf);

    bool keep = false;
    if(place==1 && !GPold.empty() && !ptest)
    {
        for(size_t i=0; i<pcs.size(); ++i)
        for(const r3gpatch &o : GPold)
        if(o.lev==1 && o.lo[0]==pcs[i].lo[0] && o.lo[1]==pcs[i].lo[1] && o.lo[2]==pcs[i].lo[2]
                    && o.hi[0]==pcs[i].hi[0] && o.hi[1]==pcs[i].hi[1] && o.hi[2]==pcs[i].hi[2])
        {
            rkp[i] = o.rank;
            lk[o.rank] += pcs[i].w;
            break;
        }
        greedy(rkp,lk);
        const double mf = *std::max_element(lf.begin(),lf.end()), mk = *std::max_element(lk.begin(),lk.end());
        keep = !(mf<(1.0-rebalance)*mk);
    }
    vector<int> rk = keep ? rkp : rf;

    // pieces that share a face on which a finer level has an edge (the boundary between two boxes of
    // the merge, not conforming) stay on one rank: the one with the most of their work (the coarse
    // cells across the edges of the finer patches local)
    {
        const int np = (int)pcs.size();
        vector<int> up(np);
        for(int i=0; i<np; ++i)
        up[i] = i;
        std::function<int(int)> root = [&](int g) { while(up[g]!=g) { up[g] = up[up[g]]; g = up[g]; } return g; };
        for(int i=0; i<np; ++i)
        for(int j=i+1; j<np; ++j)
        for(int d=0; d<3; ++d)
        {
            const piece &A = pcs[i], &B = pcs[j];
            int X = -1;
            if(A.hi[d]+1==B.lo[d]) X = B.lo[d];
            if(B.hi[d]+1==A.lo[d]) X = A.lo[d];
            if(X<0)
            continue;
            int lo[3], hi[3];
            bool ov = true;
            for(int e=0; e<3; ++e)
            {
                lo[e] = MAX(A.lo[e],B.lo[e]);
                hi[e] = MIN(A.hi[e],B.hi[e]);
                if(e!=d && lo[e]>hi[e])
                ov = false;
            }
            lo[d] = hi[d] = X;
            if(ov && !conforming(d,X,lo,hi))
            up[root(i)] = root(j);
        }
        std::map<int,vector<int>> groups;
        for(int i=0; i<np; ++i)
        groups[root(i)].push_back(i);
        for(auto &G : groups)
        {
            if(G.second.size()<2)
            continue;
            vector<double> wr(nranks,0.0);
            for(int i : G.second)
            wr[rk[i]] += pcs[i].w;
            int r = 0;
            for(int q=1; q<nranks; ++q)
            if(wr[q]>wr[r])
            r = q;
            for(int i : G.second)
            rk[i] = r;
        }
        vector<double> &L = keep ? lk : lf;
        L = base;
        for(int i=0; i<np; ++i)
        L[rk[i]] += pcs[i].w;
    }

    for(size_t i=0; i<pcs.size(); ++i)
    {
        r3gpatch q;
        q.lev = 1;
        q.rank = rk[i];
        q.lid = -1;
        for(int d=0; d<3; ++d)
        {
            q.lo[d] = pcs[i].lo[d];
            q.hi[d] = pcs[i].hi[d];
        }
        q.unit = (int)GP.size();
        glev[1].push_back((int)GP.size());
        GP.push_back(q);
    }

    // the finer levels cut at the units
    const vector<int> units = glev[1];
    for(int l=2; l<=maxlev; ++l)
    for(size_t b=0; b<mb[l].size()/6; ++b)
    for(int u : units)
    {
        r3gpatch q;
        q.lev = l;
        q.rank = GP[u].rank;
        q.lid = -1;
        q.unit = u;
        bool empty = false;
        for(int d=0; d<3; ++d)
        {
            const int R = ratio(l-1,d);
            q.lo[d] = MAX(mb[l][6*b+d],GP[u].lo[d]*R);
            q.hi[d] = MIN(mb[l][6*b+3+d],(GP[u].hi[d]+1)*R-1);
            if(q.lo[d]>q.hi[d])
            empty = true;
        }
        if(empty)
        continue;
        glev[l].push_back((int)GP.size());
        GP.push_back(q);
    }

    // (G 41 0: at every regrid, with the measured weight of a patch cell)
    if(myrank==0 && (verbose || (p0->G41<=0.0 && nranks>1)))
    {
        vector<double> &L = keep ? lk : lf;
        const double mx = *std::max_element(L.begin(),L.end());
        cout<<"CFD AMR: "<<pcs.size()<<" units placed (G 40 "<<place<<(keep ? ", kept" : "")<<"), predicted load max/mean "<<mx/wmean<<endl;
    }
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
    pp->bcside1 = p->bcside1; pp->bcside2 = p->bcside2; pp->bcside3 = p->bcside3;   // side types of the domain edges (y in 2D),
    pp->bcside4 = p->bcside4; pp->bcside5 = p->bcside5; pp->bcside6 = p->bcside6;   // as on the base grid (the lexer copy leaves them unset)
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
    // the runtime state of level 0 (time, wave parameters, water levels, extent): the lexer copy
    // (lexer_amr.cpp)

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

    // G 40: the level-0 data from the raw copy of a shadow (the shadow's own, or the one of the
    // unit of the patch), else from level 0 of this rank
    const r3raw *src = c.shadow ? c.raw : ((place && c.sh>=0) ? P[c.sh]->raw : nullptr);

    // level-0 cell below a patch cell (local index of level 0; with src: global)
    // (past a periodic edge on one rank: the cell at the other end, whose flags are those of the
    // fluid there; the halo of level 0 has the boundary flags)
    const bool wrap1 = (per[0]==1 || per[1]==1) && src==nullptr;
    auto below = [&](int ii, int jj, int kk, int *L)
    {
        const int a[3] = {ii,jj,kk};
        for(int d=0; d<3; ++d)
        L[d] = fdiv(c.lo[d]+a[d],R[d]) - (src ? 0 : org[d]);
        if(wrap1)
        {
            int Gb[3], Lw[3];
            for(int d=0; d<3; ++d)
            Gb[d] = L[d]+org[d];
            if(!in_domain(0,Gb) && wrap_local0(Gb,Lw))
            for(int d=0; d<3; ++d)
            L[d] = Lw[d];
        }
    };

    for(int ii=-m; ii<pp->knox+m; ++ii)
    for(int jj=-m; jj<pp->knoy+m; ++jj)
    for(int kk=-m; kk<pp->knoz+m; ++kk)
    {
        int L[3];
        below(ii,jj,kk,L);
        const int q = idx(pp,ii,jj,kk);
        pp->flag4[q] = src ? src->flag4[src->c3(L[0],L[1],L[2])] : p->flag4[idx(p,L[0],L[1],L[2])];
        pp->flag_solid[q] = 1.0e8;
        pp->flag_topo[q] = 1.0e8;
    }

    for(int ii=-m; ii<pp->knox+m; ++ii)
    for(int jj=-m; jj<pp->knoy+m; ++jj)
    {
        int L[3];
        below(ii,jj,0,L);
        const int q = lij(pp,ii,jj);
        if(src)
        {
            const size_t q0 = src->c2(L[0],L[1]);
            pp->flagslice1[q] = src->fs1[q0];
            pp->flagslice2[q] = src->fs2[q0];
            pp->flagslice4[q] = src->fs4[q0];
            pp->solidbed[q] = src->sbed[q0];
            pp->topobed[q] = src->tbed[q0];
            pp->geobed[q] = src->gbed[q0];
            pp->bed[q] = src->bed[q0];
            pp->wet[q] = src->wet[q0];
            pp->wet_n[q] = src->wetn[q0];
            pp->deep[q] = src->deep[q0];
            pp->depth[q] = src->depth[q0];
            pp->WL[q] = src->WL[q0];
            pp->data[q] = 0.0;
            continue;
        }
        const int q0 = lij(p,L[0],L[1]);
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
    if(src==nullptr)
    for(int q=0; q<p->gcb4_count; ++q)
    grp[key(p->gcb4[q][0],p->gcb4[q][1],p->gcb4[q][2],p->gcb4[q][3])] = p->gcb4[q][4];
    const unordered_map<long long,int> &gmap = src ? src->grp : grp;

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
            auto it = gmap.find(key(L[0],L[1],L[2],cs));
            int group = (it!=gmap.end()) ? it->second : 21;
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
    pp->Darray(pp->XN,pp->knox+1+4*marge);
    pp->Darray(pp->YN,pp->knoy+1+4*marge);
    pp->Darray(pp->ZN,pp->knoz+1+4*marge);
    for(int ii=-marge; ii<pp->knox+1+marge; ++ii)
    pp->XN[ii+marge] = node(l,0,c.lo[0]+ii);
    for(int jj=-marge; jj<pp->knoy+1+marge; ++jj)
    pp->YN[jj+marge] = node(l,1,c.lo[1]+jj);
    for(int kk=-marge; kk<pp->knoz+1+marge; ++kk)
    pp->ZN[kk+marge] = node(l,2,c.lo[2]+kk);

    // ---- gridini (without the solids of the grid file)
    pp->gridspacing(pgc0);

    // horizontal geometry layer: the patch nodes subdivide the level-0 2D nodes bilinearly (on a
    // Cartesian grid bit for bit the XN/YN above), then the metrics of the patch cells
    pp->geo_curv = p->geo_curv;
    pp->geometry_alloc();
    {
        auto node0 = [&](int K, int L, double &x, double &y)   // level-0 2D node, global index
        {
            if(src)
            {
                if(K>=src->lo[0] && K<=src->lo[0]+src->n[0] && L>=src->lo[1] && L<=src->lo[1]+src->n[1])
                {
                    x = src->xn2d[src->n2(K,L)];
                    y = src->yn2d[src->n2(K,L)];
                }
                else
                {
                    x = node(0,0,K);
                    y = node(0,1,L);
                }
                return;
            }
            const int k = K-org[0], q = L-org[1];
            if(p->XN2D!=nullptr && k>=p->imin && k<=p->imin+p->imax && q>=p->jmin && q<=p->jmin+p->jmax)
            {
                x = p->XN2D[p->nij(k,q)];
                y = p->YN2D[p->nij(k,q)];
            }
            else
            {
                x = node(0,0,K);
                y = node(0,1,L);
            }
        };

        for(int ii=pp->imin; ii<=pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<=pp->jmin+pp->jmax; ++jj)
        {
            const int NI = c.lo[0]+ii, NJ = c.lo[1]+jj;      // global level-l node
            const int K = fdiv(NI,R[0]), L = fdiv(NJ,R[1]);
            const double fr = double(NI-K*R[0])/double(R[0]), fs = double(NJ-L*R[1])/double(R[1]);
            double x00,y00,x10,y10,x01,y01,x11,y11;
            node0(K,L,x00,y00); node0(K+1,L,x10,y10); node0(K,L+1,x01,y01); node0(K+1,L+1,x11,y11);

            pp->XN2D[pp->nij(ii,jj)] = x00 + (x10-x00)*fr + (x01-x00)*fs + (x11-x10-x01+x00)*fr*fs;
            pp->YN2D[pp->nij(ii,jj)] = y00 + (y10-y00)*fr + (y01-y00)*fs + (y11-y10-y01+y00)*fr*fs;
        }
    }
    pp->geometry_metrics();
    {
        double err,gcl;
        const int mismatch = pp->geometry_check(err,gcl);
        if(mismatch!=0)
        cout<<"!!! AMR patch (level "<<l<<"): "<<mismatch<<" 2D nodes differ from XN/YN !!!"<<endl;
    }
    pp->DXM = p->DXM/double(R[0]);
    // the cells of the grid (the residual norm of the solvers of the patch, e.g. bicgstab_ijk)
    pp->cellnum = pp->cellnumtot = pp->knox*pp->knoy*pp->knoz;
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

    pgc0->flagx(pp,pp->flagsf1);
    pgc0->flagx(pp,pp->flagsf2);
    pgc0->flagx(pp,pp->flagsf3);
    pgc0->flagx(pp,pp->flagsf4);
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
                if(l==1 && place)
                add(c.sh,P[c.sh]->lo,P[c.sh]->hi);
                else if(l==1)
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
                f.depth = 0;
                {
                    const int a[3] = {ii,jj,kk}, n[3] = {pp->knox,pp->knoy,pp->knoz};
                    for(int d=0; d<3; ++d)
                    f.depth = MAX(f.depth, a[d]<0 ? -a[d] : (a[d]>=n[d] ? a[d]-n[d]+1 : 0));
                }
                f.kind = 3;
                f.g = -1;
                f.slot = -1;
                for(int d=0; d<3; ++d)
                f.s[d] = f.o[d] = 0;

                const int G[3] = {c.lo[0]+ii, c.lo[1]+jj, c.lo[2]+kk};
                int Gw[3];                  // (across a periodic edge: the cell at the other end)
                if(!wrap(l,G,Gw))
                {
                    c.fill.push_back(f);
                    continue;
                }

                const int gp = gpatch_at(l,Gw);
                if(gp>=0)
                {
                    const r3gpatch &q = GP[gp];
                    for(int d=0; d<3; ++d)
                    f.s[d] = Gw[d]-q.lo[d];
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

                if(l==1 && place)
                {
                    f.g = c.sh;
                    for(int d=0; d<3; ++d)
                    f.s[d] = Gp[d]-P[c.sh]->lo[d];
                }
                else if(l==1)
                {
                    // (past a periodic edge on one rank: the cell at the other end, else the halo)
                    int Lw[3];
                    const bool moved = (Gw[0]!=G[0] || Gw[1]!=G[1] || Gw[2]!=G[2]);
                    const bool w = (per[0]==1 || per[1]==1) && moved && wrap_local0(Gp,Lw);
                    f.g = -1;
                    for(int d=0; d<3; ++d)
                    f.s[d] = w ? Lw[d] : Gp[d]-org[d];
                }
                else
                {
                    // the parent's cell, also past a periodic edge (the patch at the other end)
                    int Gpw[3];
                    wrap(l-1,Gp,Gpw);
                    int best = -1, bestscore = -1;
                    for(int pass=0; pass<2 && best<0; ++pass)
                    {
                        const int *Q = (pass==0) ? Gp : Gpw;
                        if(pass==1 && Q[0]==Gp[0] && Q[1]==Gp[1] && Q[2]==Gp[2])
                        break;
                        for(int g : lev[l-1])
                        {
                            r3patch &q = *P[g];
                            int score = 2;
                            for(int d=0; d<3; ++d)
                            {
                                if(Q[d]<q.lo[d]-q.pp->margin+1 || Q[d]>q.hi[d]+q.pp->margin-1)
                                score = -1;
                                else if(score>=0 && (Q[d]<q.lo[d] || Q[d]>q.hi[d]))
                                score = 1;
                            }
                            if(score>bestscore)
                            {
                                bestscore = score;
                                best = g;
                            }
                        }
                        if(best>=0)
                        for(int d=0; d<3; ++d)
                        Gp[d] = Q[d];
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

    // G 40: every cell of the arrays of the shadows from the rank that holds it (collective)
    if(place)
    {
        vector<vector<int>> it(nranks);
        for(int id=npat; id<(int)P.size(); ++id)
        {
            r3patch *c = P[id];
            lexer *q = c->pp;
            for(int i=q->imin; i<q->imin+q->imax; ++i)
            for(int j=q->jmin; j<q->jmin+q->jmax; ++j)
            for(int k=q->kmin; k<q->kmin+q->kmax; ++k)
            {
                int I[3] = {c->lo[0]+i,c->lo[1]+j,c->lo[2]+k};
                wrap0(I);
                const int r = owner0(I);
                const int qq = (i-q->imin)*q->jmax*q->kmax + (j-q->jmin)*q->kmax + k-q->kmin;
                const int e[6] = {id,qq,I[0],I[1],I[2],0};
                it[r].insert(it[r].end(),e,e+6);
            }
        }
        cplan_build(sh_down,true,it);
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
