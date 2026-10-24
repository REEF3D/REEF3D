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

#include"cfd_amr.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"heaviside.h"
#include"interface_width.h"
#include"momentum_rk.h"
#include"pjm_corr.h"
#include"turbulence.h"
#include"cfd_amr_6dof.h"
#include<cmath>
#include<iostream>
#include<iomanip>
#include<algorithm>

//  Adaptive refinement (15c).  Every G 2 steps the cells of every level are flagged, each criterion
//  on its own or together (the union):
//    G 24 n   the interface band: |phi| < n psi (psi of the flagging level)
//    G 25 n   within n cells of a solid, the topography or a floating body (|solid|, |topo|, |fb| < n dx)
//    G 26 w   a vorticity |curl u| > w [1/s]
//    G 15     boxes, always refined to the finest level
//  A flagged cell of level l (with a buffer of G 3 cells) marks the tiles of its children on level
//  l+1; a tile stays refined G 5 regrids after its last flag.  reefamr3d nests the levels, keeps
//  the refinement out of the wave zones and the region edges off the partition edges, merges the
//  tiles into boxes.  The patches whose box is unchanged are kept; the new ones get their data from
//  the old patches of the same level where they overlap, else from their parents (prolongation,
//  level by level); then the coarser levels get the restriction of the finer ones.
//  The first hierarchy (setup) takes all levels from the flags of the initial level 0.

namespace
{
inline int fdiv(int a, int r) { return (a>=0) ? a/r : -((-a+r-1)/r); }

inline int cix(const lexer *q, int i, int j, int k)
{
    return (i-q->imin)*q->jmax*q->kmax + (j-q->jmin)*q->kmax + k-q->kmin;
}

struct pscope4
{
    ghostcell *g;
    fdm *back;
    bool oldc, oldl;
    pscope4(ghostcell *gg, fdm *a, fdm *a0) : g(gg), back(a0)
    {
        oldc = g->set_comms(false);
        oldl = g->set_local(true);
        g->fdm_update(a);
    }
    ~pscope4()
    {
        g->fdm_update(back);
        g->set_comms(oldc);
        g->set_local(oldl);
    }
};
}

// the criteria at cell (i,j,k) of grid g for the tiles of level lt (init: level 0 for all levels,
// with the interface thickness and the spacing of level lt-1)
bool cfd_amr::flagged(int g, int i, int j, int k, int lt, bool init)
{
    lexer *q = glex(g);
    fdm *a = gfd(g);
    const int lc = init ? lt-1 : glevel(g);
    const double scale = init ? 1.0/double(1<<lc) : 1.0;

    if(p0->G24>0.0)
    {
        const double psi = (g<0 ? p0->psi : q->psi)*scale;
        if(fabs(a->phi(i,j,k))<p0->G24*psi)
        return true;
    }

    if(p0->G25>0)
    {
        const double dx = MIN3(q->DXN[i+marge],q->DZN[k+marge],(q->j_dir==1 ? q->DYN[j+marge] : 1.0e20))*scale;
        const double w = double(p0->G25)*dx;
        if((p0->solidread!=0 && fabs(a->solid(i,j,k))<w) || (p0->toporead!=0 && fabs(a->topo(i,j,k))<w) || (bodies && fabs(a->fb(i,j,k))<w))
        return true;
    }

    if(p0->G26>0.0)
    {
        // cell-centred velocities, central differences
        auto uc = [&](int ii, int jj, int kk) { return 0.5*(a->u(ii,jj,kk)+a->u(ii-1,jj,kk)); };
        auto vc = [&](int ii, int jj, int kk) { return 0.5*(a->v(ii,jj,kk)+a->v(ii,jj-1,kk)); };
        auto wc = [&](int ii, int jj, int kk) { return 0.5*(a->w(ii,jj,kk)+a->w(ii,jj,kk-1)); };
        const double dx2 = q->XP[i+1+marge]-q->XP[i-1+marge];
        const double dz2 = q->ZP[k+1+marge]-q->ZP[k-1+marge];
        double ox = 0.0, oy, oz = 0.0;
        oy = (uc(i,j,k+1)-uc(i,j,k-1))/dz2 - (wc(i+1,j,k)-wc(i-1,j,k))/dx2;
        if(q->j_dir==1)
        {
            const double dy2 = q->YP[j+1+marge]-q->YP[j-1+marge];
            ox = (wc(i,j+1,k)-wc(i,j-1,k))/dy2 - (vc(i,j,k+1)-vc(i,j,k-1))/dz2;
            oz = (vc(i+1,j,k)-vc(i-1,j,k))/dx2 - (uc(i,j+1,k)-uc(i,j-1,k))/dy2;
        }
        if(sqrt(ox*ox+oy*oy+oz*oz)>p0->G26)
        return true;
    }

    return false;
}

// the seed tiles of every level from the flags
void cfd_amr::flags(vector<vector<unsigned char>> &seed, bool init)
{
    seed.assign(maxlev+1,vector<unsigned char>());
    for(int l=1; l<=maxlev; ++l)
    seed[l].assign((size_t)tn[l][0]*tn[l][1]*tn[l][2],0);

    if(!adapt)
    return;

    const int buf = MAX(p0->G3,0);
    int bb[3];
    for(int d=0; d<3; ++d)
    bb[d] = (rr[d]==2) ? buf : 0;

    for(int lt=1; lt<=maxlev; ++lt)
    {
        const int lc = init ? 0 : lt-1;
        vector<int> gs;
        if(lc==0)
        gs.push_back(-1);
        else
        gs = lev[lc];

        // level-(lt-1) cells per flagging cell and direction
        int R[3];
        for(int d=0; d<3; ++d)
        R[d] = ratio(lt-1-lc,d);

        for(int g : gs)
        {
            lexer *q = glex(g);
            int og[3];
            goff(g,og);
            for(int i=0; i<q->knox; ++i)
            for(int j=0; j<q->knoy; ++j)
            for(int k=0; k<q->knoz; ++k)
            {
                if(q->flag4[cix(q,i,j,k)]<=0)
                continue;
                if(!flagged(g,i,j,k,lt,init))
                continue;

                // the level-(lt-1) cells of this cell with the buffer, their children's tiles
                int lo[3], hi[3];
                const int I[3] = {i+og[0],j+og[1],k+og[2]};
                for(int d=0; d<3; ++d)
                {
                    lo[d] = MAX(I[d]*R[d]-bb[d],0);
                    hi[d] = MIN((I[d]+1)*R[d]-1+bb[d],gn(lt-1,d)-1);
                }
                int J[3];
                for(J[0]=lo[0]; J[0]<=hi[0]; ++J[0])
                for(J[1]=lo[1]; J[1]<=hi[1]; ++J[1])
                for(J[2]=lo[2]; J[2]<=hi[2]; ++J[2])
                seed[lt][tid(lt,J[0]*rr[0]/tile[0],J[1]*rr[1]/tile[1],J[2]*rr[2]/tile[2])] = 1;
            }
        }
    }

    for(int l=1; l<=maxlev; ++l)
    if(!seed[l].empty())
    MPI_Allreduce(MPI_IN_PLACE,seed[l].data(),(int)seed[l].size(),MPI_UNSIGNED_CHAR,MPI_MAX,MPI_COMM_WORLD);

    // hysteresis: a tile stays refined G 5 regrids after its last flag
    if(tage.empty())
    {
        tage.assign(maxlev+1,vector<int>());
        for(int l=1; l<=maxlev; ++l)
        tage[l].assign(seed[l].size(),1<<20);
    }
    for(int l=1; l<=maxlev; ++l)
    for(size_t t=0; t<seed[l].size(); ++t)
    {
        if(seed[l][t])
        tage[l][t] = 0;
        else if(tage[l][t]<(1<<20))
        ++tage[l][t];
        seed[l][t] = (tage[l][t]<=MAX(p0->G5,0)) ? 1 : 0;
    }
}

// G 40: the state of the cells of the new patches that lie in an old patch of their level, from that
// patch on its rank (the old placement, GPold); got: per new patch and cell (collective)
void cfd_amr::regrid_old(const vector<int> &fresh, vector<vector<char>> &got)
{
    vector<const cfd_amr_field*> st;
    for(const cfd_amr_field &e : reg)
    if(e.roles & AF_STATE)
    st.push_back(&e);
    const int nf = (int)st.size();

    vector<vector<int>> glo(maxlev+1);
    for(size_t g=0; g<GPold.size(); ++g)
    glo[GPold[g].lev].push_back((int)g);

    vector<vector<int>> req(nranks), back(nranks);
    for(int id : fresh)
    {
        cfd_amr_patch *c = CP(id);
        got[id].assign((size_t)c->n[0]*c->n[1]*c->n[2],0);
        int I[3];
        for(I[0]=c->lo[0]; I[0]<=c->hi[0]; ++I[0])
        for(I[1]=c->lo[1]; I[1]<=c->hi[1]; ++I[1])
        for(I[2]=c->lo[2]; I[2]<=c->hi[2]; ++I[2])
        {
            int go = -1;
            for(int g : glo[c->lev])
            {
                const r3gpatch &o = GPold[g];
                if(I[0]>=o.lo[0] && I[0]<=o.hi[0] && I[1]>=o.lo[1] && I[1]<=o.hi[1] && I[2]>=o.lo[2] && I[2]<=o.hi[2])
                {
                    go = g;
                    break;
                }
            }
            if(go<0)
            continue;
            const r3gpatch &o = GPold[go];
            const int li = I[0]-c->lo[0], lj = I[1]-c->lo[1], lk = I[2]-c->lo[2];
            const int oi = I[0]-o.lo[0], oj = I[1]-o.lo[1], ok = I[2]-o.lo[2];
            got[id][((size_t)li*c->n[1]+lj)*c->n[2]+lk] = 1;
            if(o.rank==myrank)
            {
                cfd_amr_patch *q = static_cast<cfd_amr_patch*>(Pold[o.lid]);
                for(const cfd_amr_field *e : st)
                rfp(*e,c)(li,lj,lk) = rfp(*e,q)(oi,oj,ok);
                continue;
            }
            const int a[4] = {o.lid,oi,oj,ok}, b[4] = {id,li,lj,lk};
            req[o.rank].insert(req[o.rank].end(),a,a+4);
            back[o.rank].insert(back[o.rank].end(),b,b+4);
        }
    }

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

    vector<double> ab((size_t)(nr/4)*nf+1);
    for(int m=0; m<nr/4; ++m)
    {
        cfd_amr_patch *q = static_cast<cfd_amr_patch*>(Pold[rb[4*m]]);
        for(int f=0; f<nf; ++f)
        ab[(size_t)m*nf+f] = rfp(*st[f],q)(rb[4*m+1],rb[4*m+2],rb[4*m+3]);
    }
    vector<int> ac(nranks), ad(nranks), bc(nranks), bd(nranks);
    for(int r=0; r<nranks; ++r)
    {
        ac[r] = rc[r]/4*nf; ad[r] = rd[r]/4*nf;
        bc[r] = sc[r]/4*nf; bd[r] = sd[r]/4*nf;
    }
    vector<double> bb((size_t)(ns/4)*nf+1);
    MPI_Alltoallv(ab.data(),ac.data(),ad.data(),MPI_DOUBLE,bb.data(),bc.data(),bd.data(),MPI_DOUBLE,MPI_COMM_WORLD);

    for(int r=0; r<nranks; ++r)
    for(size_t m=0; m<back[r].size()/4; ++m)
    {
        const int *b = &back[r][4*m];
        cfd_amr_patch *c = CP(b[0]);
        const double *v = &bb[(size_t)bd[r]+m*nf];
        for(int f=0; f<nf; ++f)
        rfp(*st[f],c)(b[1],b[2],b[3]) = v[f];
    }
}

// a new hierarchy from the flags, the data transferred
// G 41 0 on several ranks: the work of the ranks since the last regrid, from the timers (level 0: its
// level set, momentum, end of the stages (G 7 1: its step) and the level-0 V-cycles of the composite
// solves; the patches: the loops over the patches of the rank and their V-cycles, without the fills
// and the waits).  The work of a patch cell against a level-0 cell (both per step of their level) is
// the weight wpatch of the patch cells in the placement of G 40 1, 2 (the mean with the last one).
// The report: the measured work max/mean over the ranks, and the weight.
bool cfd_amr::lmeas() const
{
    return nranks>1 && p0->G41<=0.0;
}

void cfd_amr::lbar()
{
    if(lmeas())
    MPI_Barrier(MPI_COMM_WORLD);
}

void cfd_amr::measure_load()
{
    if(!lmeas())
    return;

    const double T0 = t0l[0]+t0l[1]+t0l[2]+tm[4]+tp[1];
    const double Tp = wtp+tp[0];
    double n0 = double(p0->knox)*p0->knoy*p0->knoz, np = 0.0;
    for(int qi=0; qi<npat; ++qi)
    np += ((p0->G7==1) ? double(1<<P[qi]->lev) : 1.0)*double(P[qi]->n[0])*P[qi]->n[1]*P[qi]->n[2];

    double d[5] = {T0-ml_t0, Tp-ml_tp, n0, np, 0.0};
    d[4] = d[0]+d[1];
    ml_t0 = T0;
    ml_tp = Tp;
    ++ml_n;
    double s[5], m[5];
    MPI_Allreduce(d,s,5,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(d,m,5,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    // (the first regrid: the timers since the start, with the initialisation)
    if(ml_n<2 || s[0]<=0.0 || s[2]<=0.0)
    return;

    double w = -1.0;
    if(s[1]>0.0 && s[3]>0.0)
    {
        w = (s[1]/s[3])/(s[0]/s[2]);
        w = MIN(MAX(w,0.25),16.0);
        wpatch = (ml_n==2) ? w : 0.5*(wpatch+w);
    }
    if(p0->mpirank==0)
    {
        const std::streamsize pr = cout.precision();
        cout<<"CFD AMR: measured since the last regrid: work max/mean over "<<nranks<<" ranks "<<setprecision(3)<<m[4]/(s[4]/nranks)
            <<" (level 0 "<<m[0]/(s[0]/nranks)<<", patches "<<(s[1]>0.0 ? m[1]/(s[1]/nranks) : 0.0)<<")";
        if(w>0.0)
        cout<<"; a patch cell "<<w<<" level-0 cells, the weight of the placement "<<wpatch;
        cout<<endl;
        cout.precision(pr);
    }
}

// G 40 1 with G 41 0 and a static refinement (no regrid): the work measured over steps 5-15, then
// the units placed once more with the measured weight of a patch cell
void cfd_amr::replace_static(lexer *p, ghostcell *pgc)
{
    if(adapt || place!=1 || !lmeas() || (p->count!=5 && p->count!=15))
    return;
    measure_load();
    if(p->count==15)
    regrid(p,pgc,true);
}

void cfd_amr::regrid(lexer *p, ghostcell *pgc, bool replace)
{
    amr_regrid_clock rclock;     // load statistics (G 43)
    const double t0 = MPI_Wtime();

    vector<r3patch*> gone;
    vector<int> fresh;
    verbose = false;
    if(replace)
    {
        if(!replace_units(gone,fresh))
        return;
    }
    else
    {
        vector<vector<unsigned char>> seed;

        // restart (I 45): the tile maps of the checkpoint (ckpt_read), no flags, no load sample
        if(!ck_seed.empty())
        seed = ck_seed;
        else
        {
        measure_load();
        flags(seed,false);
        }

        if(!rebuild(seed,gone,fresh))
        return;
    }
    ++nregrid;

    vector<char> isfresh(P.size(),0);
    for(int id : fresh)
    {
        isfresh[id] = 1;
        patch_fields_ini(id);
    }

    vector<fspec> fs = specs(AF_STATE);

    // G 40: the cells of the new patches in an old patch of their level (of any rank)
    vector<vector<char>> got(P.size());
    if(place)
    regrid_old(fresh,got);

    // G 40: the shadows of the level-1 patches (new ones: their porosity as well)
    if(place)
    {
        sh_fill_down(fs);
        vector<std::function<field&(int)>> por;
        por.push_back([this](int g) -> field& { return gfd(g)->porosity; });
        sh_down_f(por);
    }

    for(int l=1; l<=maxlev; ++l)
    {
        for(int id : lev[l])
        {
            if(!isfresh[id])
            continue;
            cfd_amr_patch *c = CP(id);

            int I[3];
            for(I[0]=c->lo[0]; I[0]<=c->hi[0]; ++I[0])
            for(I[1]=c->lo[1]; I[1]<=c->hi[1]; ++I[1])
            for(I[2]=c->lo[2]; I[2]<=c->hi[2]; ++I[2])
            {
                const int li = I[0]-c->lo[0], lj = I[1]-c->lo[1], lk = I[2]-c->lo[2];

                // G 40: from an old patch of the level on any rank (regrid_old)
                if(place && got[id][((size_t)li*c->n[1]+lj)*c->n[2]+lk])
                continue;

                // an old patch of the level
                cfd_amr_patch *o = nullptr;
                if(!place)
                for(r3patch *q : gone)
                if(q->lev==l && I[0]>=q->lo[0] && I[0]<=q->hi[0] && I[1]>=q->lo[1] && I[1]<=q->hi[1] && I[2]>=q->lo[2] && I[2]<=q->hi[2])
                {
                    o = static_cast<cfd_amr_patch*>(q);
                    break;
                }
                if(o!=nullptr)
                {
                    const int oi = I[0]-o->lo[0], oj = I[1]-o->lo[1], ok = I[2]-o->lo[2];
                    for(const cfd_amr_field &e : reg)
                    if(e.roles & AF_STATE)
                    rfp(e,c)(li,lj,lk) = rfp(e,o)(oi,oj,ok);
                    continue;
                }

                // the parent
                for(const r3patch::pblock &B : c->par)
                {
                    int og[3], Pg[3], s[3], oo[3];
                    goff(B.g,og);
                    bool in = true;
                    for(int d=0; d<3; ++d)
                    {
                        Pg[d] = fdiv(I[d],rr[d]);
                        oo[d] = I[d]-Pg[d]*rr[d];
                        s[d] = Pg[d]-og[d];
                        if(Pg[d]<B.lo[d] || Pg[d]>B.hi[d])
                        in = false;
                    }
                    if(!in)
                    continue;
                    for(auto &F : fs)
                    F.f(id)(li,lj,lk) = prolong(F,B.g,s,oo);
                    break;
                }
            }
        }

        // the solids and the topography at the resolution of the new patches
        for(int id : lev[l])
        if(isfresh[id])
        patch_geometry(id);

        fill(l,fs);

        // the boundary values of the new patches (the parents of the next level)
        for(int id : lev[l])
        if(isfresh[id])
        bc_state(id);
    }

    for(r3patch *q : gone)
    free_patch(static_cast<cfd_amr_patch*>(q));

    hierarchy(false);

    // density and viscosity; boundary values and the pressure objects of the new patches
    rovisc_all();
    for(int id : fresh)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        fdm *fa = c->a;
        pscope4 ps(pgc,fa,a0);
        pgc->start1(pp,fa->u,10);
        pgc->start2(pp,fa->v,11);
        pgc->start3(pp,fa->w,12);
        pgc->start4(pp,fa->phi,gcval_phi);
        pgc->start4(pp,fa->press,40);
        pgc->start4(pp,fa->ro,1);
        pgc->start4(pp,fa->visc,1);
        c->ppress->ini(pp,fa,pgc);
        if(solids)
        pgc->solid_forcing_ini(pp,fa);
        if(bodies)
        c->p6dof->body(pp,fa,pgc);
    }

    // the coarser levels from the finer ones
    {
        vector<std::function<field&(int)>> cells, faces;
        cells.push_back([this](int g) -> field& { return gfd(g)->phi; });
        cells.push_back([this](int g) -> field& { return gfd(g)->ro; });
        cells.push_back([this](int g) -> field& { return gfd(g)->visc; });
        if(!sub)
        cells.push_back([this](int g) -> field& { return gfd(g)->press; });
        if(turb)
        {
            cells.push_back([this](int g) -> field& { return gfd(g)->eddyv; });
            cells.push_back([this](int g) -> field& { return kfield(g,0); });
            cells.push_back([this](int g) -> field& { return kfield(g,1); });
        }
        faces.push_back([this](int g) -> field& { return gfd(g)->u; });
        faces.push_back([this](int g) -> field& { return gfd(g)->v; });
        faces.push_back([this](int g) -> field& { return gfd(g)->w; });
        for(int l=maxlev; l>=1; --l)
        {
            restrict_level(l,cells,faces);
            if(l>=2)
            fill(l-1,fs);
        }
        pgc0->start1(p0,a0->u,10);
        pgc0->start2(p0,a0->v,11);
        pgc0->start3(p0,a0->w,12);
        pgc0->start4(p0,a0->phi,gcval_phi);
        pgc0->start4(p0,a0->press,40);
        pgc0->start4(p0,a0->ro,1);
        pgc0->start4(p0,a0->visc,1);
        if(turb)
        pturb0->gcupdate(p0,a0,pgc0);
        for(int l=1; l<=maxlev; ++l)
        fill(l,fs);
        for(int id=0; id<npat; ++id)
        margin_rovisc(id);
    }

    // the old time level of kin, eps of the new patches (IM1: the time term; the kept patches
    // have theirs from the end of the step)
    if(turb)
    for(int id : fresh)
    {
        CP(id)->pturb->ktimesave(CP(id)->pp,CP(id)->a,pgc);
        CP(id)->pturb->etimesave(CP(id)->pp,CP(id)->a,pgc);
    }

    if(p->mpirank==0)
    {
        const std::streamsize pr = cout.precision();
        cout<<"CFD AMR: regrid "<<nregrid<<": "<<fresh.size()<<" new patches on rank 0, "<<setprecision(3)<<MPI_Wtime()-t0<<" s"<<endl;
        cout.precision(pr);
    }
}
