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

#include"sflow_amr.h"
#include"sflow_boussinesq.h"
#include"sflow_momentum_RK3.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"slice4.h"
#include"reefamr_krylov.h"
#include<cmath>
#include<limits>
#include<algorithm>

//  Boussinesq u_a (A 220 4) on the leaf cells of all levels.
//
//  Every grid builds the rows of its two line systems (sflow_boussinesq::invert_prepare: x-x and
//  y-y derivatives implicit, cross derivatives from the u_a at hand).  The unknowns are the leaf
//  cells, level 0 and the patches, as for the non-hydrostatic pressure (sflow_amr_nh.cpp): a row
//  takes its neighbours outside the grid from the finer or coarser level (covered cells: average
//  of the children; cells around a patch: sibling, owner rank across a partition edge, linear
//  interpolation of the coarser level).  Coarse and fine u_a are solved together in the stage,
//  instead of the patch taking the coarse u_a of the previous stage as boundary value, and the
//  result does not depend on how a level is split into patches and ranks.
//
//  BiCGStab (reefamr_bicgstab), right preconditioned by exact line solves on every grid (the
//  leaf cells of a line, zero outside), once for u_a and once for v_a.

namespace
{
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

typedef reefamr_comms_off comms_off;

// Krylov vectors per grid (in nh_vec storage; the non-hydrostatic solve is not used with A 220 4)
enum { BX=0, BB, BR, BRH, BPV, BVV, BS, BT, BPH, BSH, BNV };

// rows of a grid: dir 0: u_a, p + n (i+1) + s (i-1);  dir 1: v_a, p + w (j+1) + e (j-1)
struct bqrow { slice *p, *up, *dn, *r; };

bqrow rows(sflow_boussinesq *q, int dir)
{
    bqrow R;
    if(dir==0) { R.p=&q->Xp; R.up=&q->Xn; R.dn=&q->Xs; R.r=&q->Xr; }
    else       { R.p=&q->Yp; R.up=&q->Yw; R.dn=&q->Ye; R.r=&q->Yr; }
    return R;
}
}

// vectors, leaf cells (act: 1 leaf, -1 covered by a finer patch, -2 no cell), line segments of the
// preconditioner and the cells around the patches that the rows reach; once per hierarchy and
// level window (G 7 1: the levels wlo..wtop, the other grids without leaf cells)
void sflow_amr::bq_setup()
{
    while((int)nh0_v.size()<BNV)
    nh0_v.push_back(new slice4(p0));

    for(auto q : P)
    {
    sflow_amr_patch *c = SP(q);
    while((int)c->nv.size()<BNV)
    c->nv.push_back(new slice4(c->pp));
    }

    for(int n=0; n<(int)P.size(); ++n)
    SP(n)->bqid = n;

    if(bq_layout!=regrids)
    bqw.clear();
    bq_layout = regrids;

    const int wkey = wlo*1000 + wtop();
    auto it = bqw.find(wkey);
    if(it!=bqw.end())
    {
        bqc = &it->second;
        return;
    }

    vector<bqgrid> &bqg = bqw[wkey];
    bqc = &bqg;
    bqg.assign(P.size()+1,bqgrid());

    auto setup = [&](int g, int i0, int i1, int j0, int j1, int l)
    {
        nhg G = nh_grid(g);
        lexer *q = G.q;
        G.act->assign(q->imax*q->jmax,-2);

        if(l<wlo || l>wtop())
        return;

        for(int ii=i0; ii<=i1; ++ii)
        for(int jj=j0; jj<=j1; ++jj)
        {
            int c = lij(q,ii,jj);
            if(q->flagslice4[c]<0)
            continue;

            int I = (g<0) ? ii+O0i : ii-EXT+P[g]->I0;
            int J = (g<0) ? jj+O0j : jj-EXT+P[g]->J0;

            (*G.act)[c] = (l<wtop() && covered(l+1,2*I,2*J)) ? -1 : 1;
        }

        bqgrid &B = (*bqc)[g+1];
        for(int ii=i0; ii<=i1; ++ii)
        for(int jj=j0; jj<=j1; ++jj)
        if((*G.act)[lij(q,ii,jj)]>=0)
        B.leaf.push_back(lij(q,ii,jj));

        // runs of leaf cells along i (dir 0) and j (dir 1)
        for(int dir=0; dir<2; ++dir)
        {
            const int na = (dir==0) ? j1-j0+1 : i1-i0+1;
            const int nl = (dir==0) ? i1-i0+1 : j1-j0+1;
            for(int t=0; t<na; ++t)
            {
                int k=0;
                while(k<nl)
                {
                    int c = (dir==0) ? lij(q,i0+k,j0+t) : lij(q,i0+t,j0+k);
                    if((*G.act)[c]<0) { ++k; continue; }
                    int k0=k;
                    while(k<nl && (*G.act)[(dir==0) ? lij(q,i0+k,j0+t) : lij(q,i0+t,j0+k)]>=0) ++k;
                    B.seg[dir].push_back(c);
                    B.seg[dir].push_back(k-k0);
                }
            }
        }
    };

    setup(-1,0,NX0-1,0,NY0-1,0);
    for(int n=0; n<(int)P.size(); ++n)
    setup(n,EXT,EXT+P[n]->nx-1,EXT,EXT+P[n]->ny-1,P[n]->lev);

    // cells around a patch that the rows of its leaf cells reach
    for(int n=0; n<(int)P.size(); ++n)
    {
        sflow_amr_patch *c = SP(n);
        lexer *pp = c->pp;
        auto leaf = [&](int ii, int jj)
        {
            return ii>=EXT && ii<EXT+c->nx && jj>=EXT && jj<EXT+c->ny && c->act[lij(pp,ii,jj)]>=0;
        };

        for(int dir=0; dir<2; ++dir)
        {
            vector<int> &E = bqg[n+1].need[dir];
            const int di = (dir==0) ? 1 : 0, dj = 1-di;
            for(int k=0; k<(int)c->fill.size(); ++k)
            {
                const reefamr_fill &f = c->fill[k];
                if(leaf(f.di+di,f.dj+dj) || leaf(f.di-di,f.dj-dj))
                E.push_back(k);
            }
        }
    }
}

// the cells that the rows reach: covered cells (restriction), partition halo of level 0,
// cells around the patches next to a leaf cell along the lines
// (G 7 1, the lowest level of a window above level 0: the cells the parent fills keep their values
// in the unknowns, 0 in the Krylov vectors)
void sflow_amr::bq_sync(int k, int dir)
{
    nh_restrict_vec(k);

    if(p0->mpi_size>1 && wlo==0)
    {
        pgc0->gcslparax(p0,nh_vec(-1,k),4);
        pgc0->gcslparacox(p0,nh_vec(-1,k),10);
    }

    vector<bqgrid> &bqg = *bqc;

    for(int l=MAX(wlo,1); l<=wtop(); ++l)
    {
    const bool edge = (nh_edge && l==wlo);
    fill_run_sub(l,1,7300+l,
                 [&](const reefamr_fill &f, double *v)
                 {
                     if(edge && f.kind==1)
                     v[0] = std::numeric_limits<double>::quiet_NaN();
                     else
                     v[0] = nh_eval(f,k);
                 },
                 [&](reefamr_patch *c, int n, const reefamr_fill &f, const double *w)
                 {
                     if(edge && std::isnan(w[0]))
                     {
                         if(k!=nh_keep)
                         nh_vec(n,k)(f.di,f.dj) = 0.0;
                         return;
                     }
                     nh_vec(n,k)(f.di,f.dj) = w[0];
                 },
                 [&](reefamr_patch *c) { return &bqg[SP(c)->bqid+1].need[dir]; });
    }
}

double sflow_amr::bq_dot(int ka, int kb)
{
    double s=0.0;
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *a = nh_vec(g,ka).data(), *b = nh_vec(g,kb).data();
        for(int n : (*bqc)[g+1].leaf)
        s += a[n]*b[n];
    }
    return pgc0->globalsum(s);
}

// y = A x on the leaf cells
void sflow_amr::bq_apply(int kx, int ky, int dir)
{
    bq_sync(kx,dir);

    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = (g<0) ? p0 : P[g]->pp;
        bqrow R = rows((g<0 ? pmom0 : SP(g)->pmom)->boussinesq(),dir);
        const double *x = nh_vec(g,kx).data(), *rp = R.p->data(), *ru = R.up->data(), *rd = R.dn->data();
        double *y = nh_vec(g,ky).data();
        const int off = (dir==0) ? q->jmax : 1;

        for(int n : (*bqc)[g+1].leaf)
        y[n] = rp[n]*x[n] + ru[n]*x[n+off] + rd[n]*x[n-off];
    }
}

// z = M^-1 r: exact solves along the lines of every grid over its leaf cells, zero outside
void sflow_amr::bq_prec(int kr, int kz, int dir)
{
    vector<double> bb,d;

    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = (g<0) ? p0 : P[g]->pp;
        bqrow R = rows((g<0 ? pmom0 : SP(g)->pmom)->boussinesq(),dir);
        const double *r = nh_vec(g,kr).data(), *rp = R.p->data(), *ru = R.up->data(), *rd = R.dn->data();
        double *z = nh_vec(g,kz).data();
        const int off = (dir==0) ? q->jmax : 1;
        const vector<int> &S = (*bqc)[g+1].seg[dir];

        for(size_t k=0; k<S.size(); k+=2)
        {
            const int n0 = S[k], m = S[k+1];
            bb.resize(m); d.resize(m);

            bb[0] = rp[n0];
            d[0] = r[n0];
            for(int s=1; s<m; ++s)
            {
                const int n = n0+s*off;
                double w = rd[n]/bb[s-1];
                bb[s] = rp[n] - w*ru[n-off];
                d[s]  = r[n] - w*d[s-1];
            }
            d[m-1] /= bb[m-1];
            for(int s=m-2; s>=0; --s)
            d[s] = (d[s] - ru[n0+s*off]*d[s+1])/bb[s];

            for(int s=0; s<m; ++s)
            z[n0+s*off] = d[s];
        }
    }
}

// vector space of reefamr_bicgstab: element-wise operations on the leaf cells
void sflow_amr::bq_start(int dir)
{
    for(int g=0; g<(int)bqc->size(); ++g)
    {
        const double *vb = nh_vec(g-1,BB).data();
        double *vr = nh_vec(g-1,BR).data(), *vrh = nh_vec(g-1,BRH).data(), *vp = nh_vec(g-1,BPV).data(), *vv = nh_vec(g-1,BVV).data();
        for(int n : (*bqc)[g].leaf)
        {
            double r = vb[n] - vv[n];
            vr[n] = r; vrh[n] = r; vp[n] = 0.0; vv[n] = 0.0;
        }
    }
}

void sflow_amr::bq_p(double beta, double om)
{
    for(int g=0; g<(int)bqc->size(); ++g)
    {
        const double *vr = nh_vec(g-1,BR).data(), *vv = nh_vec(g-1,BVV).data();
        double *vp = nh_vec(g-1,BPV).data();
        for(int n : (*bqc)[g].leaf)
        vp[n] = vr[n] + beta*(vp[n] - om*vv[n]);
    }
}

void sflow_amr::bq_s(double alp)
{
    for(int g=0; g<(int)bqc->size(); ++g)
    {
        const double *vr = nh_vec(g-1,BR).data(), *vv = nh_vec(g-1,BVV).data();
        double *vs = nh_vec(g-1,BS).data();
        for(int n : (*bqc)[g].leaf)
        vs[n] = vr[n] - alp*vv[n];
    }
}

void sflow_amr::bq_x(double alp, double om)
{
    for(int g=0; g<(int)bqc->size(); ++g)
    {
        const double *vph = nh_vec(g-1,BPH).data(), *vsh = nh_vec(g-1,BSH).data(), *vs = nh_vec(g-1,BS).data(), *vt = nh_vec(g-1,BT).data();
        double *vx = nh_vec(g-1,BX).data(), *vr = nh_vec(g-1,BR).data();
        for(int n : (*bqc)[g].leaf)
        {
            vx[n] += alp*vph[n] + om*vsh[n];
            vr[n] = vs[n] - om*vt[n];
        }
    }
}

namespace
{
// the u_a (dir 0) or v_a (dir 1) line systems as the vector space of reefamr_bicgstab
struct bq_space
{
    static const int B=BB, R=BR, RH=BRH, PV=BPV, VV=BVV, S=BS, T=BT, PH=BPH, SH=BSH;
    sflow_amr *a;
    int dir;
    void apply(int x, int y) { a->bq_apply(x,y,dir); }
    void prec(int r, int z) { a->bq_prec(r,z,dir); }
    double dot(int x, int y) { return a->bq_dot(x,y); }
    void op_start() { a->bq_start(dir); }
    void op_p(double beta, double om) { a->bq_p(beta,om); }
    void op_s(double alp) { a->bq_s(alp); }
    void op_x(double alp, double om) { a->bq_x(alp,om); }
};
}

void sflow_amr::bous_solve(ghostcell *pgc)
{
    bous_window(pgc,-1,0.0);

    // the patch stages continue: velocities, flux M, relaxation zones
    const double t1 = MPI_Wtime();
    {
    comms_off guard(pgc);
    for(auto q : P)
    {
    sflow_amr_patch *c = SP(q);
    c->pmom->bous_resume(c->pp,c->b,pgc);
    }
    }
    tm[1] += MPI_Wtime()-t1;
}

// u_a of all levels (l < 0), or (G 7 1) of the level-l patches with the cells around them that the
// parent fills fixed at the parent u_a linear in time between the start and the end of its step (th)
void sflow_amr::bous_window(ghostcell *pgc, int l, double th)
{
    const double t0 = MPI_Wtime();

    if(l>=0)
    {
        wlo = whi = l;
        nh_edge = true;
        nh_keep = BX;
    }

    bq_setup();

    const int ng = (int)P.size()+1;
    vector<double*> vx(ng),vb(ng);
    for(int g=-1; g<ng-1; ++g)
    {
        vx[g+1]=nh_vec(g,BX).data(); vb[g+1]=nh_vec(g,BB).data();
    }

    auto ua = [&](int g, int dir) -> slice& { fdm2D *b = (g<0) ? b0 : SP(g)->b; return dir==0 ? b->UA : b->VA; };

    const int ndir = (p0->j_dir==1) ? 2 : 1;
    int itsum=0;

    for(int dir=0; dir<ndir; ++dir)
    {
        // initial guess: u_a of the stage; right-hand side
        for(int g=-1; g<ng-1; ++g)
        {
            std::ranges::copy(ua(g,dir),vx[g+1]);

            if((*bqc)[g+1].leaf.empty())
            continue;
            const double *rr = rows((g<0 ? pmom0 : SP(g)->pmom)->boussinesq(),dir).r->data();
            for(int n : (*bqc)[g+1].leaf)
            vb[g+1][n] = rr[n];
        }

        // G 7 1: the parent u_a around the level-l patches
        if(l>=0)
        fill_run(l,1,7360+l,
                 [&](const reefamr_fill &f, double *v) { v[0] = nh_eval_t(f,th,dir==0 ? 5 : 6); },
                 [&](reefamr_patch *c, int n, const reefamr_fill &f, const double *w) { nh_vec(n,BX)(f.di,f.dj) = w[0]; });

        bq_apply(BX,BVV,dir);

        double bn, rn;
        bq_space sp{this,dir};
        int it = reefamr_bicgstab(sp,p0->N44,p0->N46,bn,rn);
        itsum += it;

        // u_a everywhere: leaf cells, covered cells, around the patches, partition halo
        nh_sync(BX);

        for(int g=-1; g<ng-1; ++g)
        {
            const int lg = (g<0) ? 0 : P[g]->lev;
            if(lg<wlo || lg>wtop())
            continue;
            lexer *q = (g<0) ? p0 : P[g]->pp;
            slice &x = nh_vec(g,BX);
            slice &u = ua(g,dir);
            const int i1 = (g<0) ? NX0 : q->knox;
            const int j1 = (g<0) ? NY0 : q->knoy;

            // level 0: its cells (ghost cells follow in velcalc); patches: the computed cells
            for(int ii=0; ii<i1; ++ii)
            for(int jj=0; jj<j1; ++jj)
            if(q->flagslice4[lij(q,ii,jj)]>0)
            u(ii,jj) = (q->wet[lij(q,ii,jj)]==1) ? x(ii,jj) : 0.0;
        }
    }

    if(l>=0)
    {
        wlo = 0;
        whi = -1;
        nh_edge = false;
        nh_keep = -1;
        sub_lv_it += itsum;
        ++sub_lv_n;
    }

    bq_it_total += itsum;
    ++bq_solves;
    tm[9] += MPI_Wtime()-t0;
}

// G 7 1: u_a of the level-l stage s of step k (0, 1), then the rest of the stage
void sflow_amr::bous_level(lexer *p, ghostcell *pgc, int l, int k, int s)
{
    // times of the RK3 stage outputs in the step
    static const double rko[3] = {1.0, 0.5, 1.0};
    bous_window(pgc,l,0.5*(double(k)+rko[s]));

    const double t1 = MPI_Wtime();
    comms_off guard(pgc);
    for(int n : lev[l])
    {
        sflow_amr_patch *c = SP(n);
        c->pmom->bous_resume(c->pp,c->b,pgc);
    }
    tm[1] += MPI_Wtime()-t1;
}

// G 7 1: after refluxing and restriction into level l (in its step k): its u_a from the
// synchronised V, then M and U (level 0: its own inversion)
void sflow_amr::bous_sync(lexer *p, ghostcell *pgc, int l, int k)
{
    if(l==0)
    {
        const double t0 = MPI_Wtime();
        pmom0->velcalc(p0,b0,pgc,b0->UH,b0->VH,b0->WH,b0->WL,1);
        pmom0->velcalc(p0,b0,pgc,b0->UH,b0->VH,b0->WH,b0->WL,2);
        tm[9] += MPI_Wtime()-t0;
        return;
    }

    // the rows of the line systems from the synchronised V
    {
    comms_off guard(pgc);
    for(int n : lev[l])
    {
        sflow_amr_patch *c = SP(n);
        c->pmom->bous_prepare(c->pp,c->b,pgc,c->b->UH,c->b->VH,c->b->WL);
    }
    }

    bous_window(pgc,l,0.5*double(k+1));

    comms_off guard(pgc);
    for(int n : lev[l])
    {
        sflow_amr_patch *c = SP(n);
        c->pmom->velcalc(c->pp,c->b,pgc,c->b->UH,c->b->VH,c->b->WH,c->b->WL,2);
    }
}
