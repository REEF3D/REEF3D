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
// preconditioner; once per hierarchy and
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
}

// the cells that the rows reach: covered cells (restriction), partition halo of level 0,
// all cells around the patches (the prolongation into a finer patch reads the cells around its
// parent patch in both directions: with the cells along the lines only, the others kept the values
// of an earlier vector, an inconsistent operator; BiCGStab broke down or diverged)
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

    for(int l=MAX(wlo,1); l<=wtop(); ++l)
    {
    const bool edge = (nh_edge && l==wlo);
    fill_run(l,1,7300+l,
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
             });
    }
}

double sflow_amr::bq_dot(int ka, int kb)
{
    double s=0.0;
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *a = nh_vec(g,ka).V, *b = nh_vec(g,kb).V;
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
        const double *x = nh_vec(g,kx).V, *rp = R.p->V, *ru = R.up->V, *rd = R.dn->V;
        double *y = nh_vec(g,ky).V;
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
        const double *r = nh_vec(g,kr).V, *rp = R.p->V, *ru = R.up->V, *rd = R.dn->V;
        double *z = nh_vec(g,kz).V;
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
        const double *vb = nh_vec(g-1,BB).V;
        double *vr = nh_vec(g-1,BR).V, *vrh = nh_vec(g-1,BRH).V, *vp = nh_vec(g-1,BPV).V, *vv = nh_vec(g-1,BVV).V;
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
        const double *vr = nh_vec(g-1,BR).V, *vv = nh_vec(g-1,BVV).V;
        double *vp = nh_vec(g-1,BPV).V;
        for(int n : (*bqc)[g].leaf)
        vp[n] = vr[n] + beta*(vp[n] - om*vv[n]);
    }
}

void sflow_amr::bq_s(double alp)
{
    for(int g=0; g<(int)bqc->size(); ++g)
    {
        const double *vr = nh_vec(g-1,BR).V, *vv = nh_vec(g-1,BVV).V;
        double *vs = nh_vec(g-1,BS).V;
        for(int n : (*bqc)[g].leaf)
        vs[n] = vr[n] - alp*vv[n];
    }
}

void sflow_amr::bq_x(double alp, double om)
{
    for(int g=0; g<(int)bqc->size(); ++g)
    {
        const double *vph = nh_vec(g-1,BPH).V, *vsh = nh_vec(g-1,BSH).V, *vs = nh_vec(g-1,BS).V, *vt = nh_vec(g-1,BT).V;
        double *vx = nh_vec(g-1,BX).V, *vr = nh_vec(g-1,BR).V;
        for(int n : (*bqc)[g].leaf)
        {
            vx[n] += alp*vph[n] + om*vsh[n];
            vr[n] = vs[n] - om*vt[n];
        }
    }
}

// new shadow residual (reefamr_bicgstab: stagnation, breakdown)
void sflow_amr::bq_restart()
{
    for(int g=0; g<(int)bqc->size(); ++g)
    {
        const double *vr = nh_vec(g-1,BR).V;
        double *vrh = nh_vec(g-1,BRH).V, *vp = nh_vec(g-1,BPV).V, *vv = nh_vec(g-1,BVV).V;
        for(int n : (*bqc)[g].leaf)
        {
            vrh[n] = vr[n];
            vp[n] = 0.0;
            vv[n] = 0.0;
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
    void op_restart() { a->bq_restart(); }
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
    if(p0->geo_curv==1)
    {
        bous_window_curv(pgc,l,th);
        return;
    }

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
        vx[g+1]=nh_vec(g,BX).V; vb[g+1]=nh_vec(g,BB).V;
    }

    auto ua = [&](int g, int dir) -> slice& { fdm2D *b = (g<0) ? b0 : SP(g)->b; return dir==0 ? b->UA : b->VA; };

    const int ndir = (p0->j_dir==1) ? 2 : 1;
    int itsum=0;

    for(int dir=0; dir<ndir; ++dir)
    {
        // initial guess: u_a of the stage; right-hand side
        for(int g=-1; g<ng-1; ++g)
        {
            lexer *q = (g<0) ? p0 : P[g]->pp;
            const double *u = ua(g,dir).V;
            double *x = vx[g+1];
            for(int n=0; n<q->imax*q->jmax; ++n)
            x[n] = u[n];

            if((*bqc)[g+1].leaf.empty())
            continue;
            const double *rr = rows((g<0 ? pmom0 : SP(g)->pmom)->boussinesq(),dir).r->V;
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

// ---------------------------------------------------------------------------------------------
// curvilinear grid (G 50): the inversion V = H (u_a + V1'(u_a)) for the grid-aligned components
// a, b of all leaf cells together, per component with the other one lagged (as the Cartesian
// line systems), the own component with all its terms: BiCGStab with the operator of every grid
// (sflow_boussinesq::apply_curv, matrix free) and the line solves of the rows of prepare_curv as
// preconditioner (sflow_boussinesq::solve_curv on one grid).  The vectors are pairs (a: index k,
// b: k+CQN), as the covered cells and the cells around the patches take the components of the
// averaged (interpolated) vector a e1 + b e2 along their own basis: the bases of a parent and its
// children, or of the cells of a stencil, differ.  The operator is affine (in- and outflow ghost
// values, the lagged component, G 7 1 the parent values around the patches): A x0 with them, the
// Krylov vectors (other component 0) with the constant part C0 = A(0) taken off.
// ---------------------------------------------------------------------------------------------
namespace
{
enum { CQX=0, CQB, CQR, CQRH, CQPV, CQVV, CQS, CQT, CQPH, CQSH, CQC0, CQO, CQN };
}

slice& sflow_amr::cq_vec(int g, int k, int c)
{
    return nh_vec(g,k+c*CQN);
}

// covered cells, partition halo of level 0, cells around the patches of pair k (wlo..wtop)
void sflow_amr::cq_sync(int k)
{
    auto bous = [&](int g) { return (g<0 ? pmom0 : SP(g)->pmom)->boussinesq(); };

    auto vec = [&](int g, int ii, int jj, double &vx, double &vy)
    {
        double ex,ey;
        bous(g)->basis(glex(g),ii,jj,ex,ey);
        const double a = cq_vec(g,k,0)(ii,jj), c = cq_vec(g,k,1)(ii,jj);
        vx = a*ex - c*ey;
        vy = a*ey + c*ex;
    };
    auto put = [&](int g, int ii, int jj, double vx, double vy)
    {
        double ex,ey;
        bous(g)->basis(glex(g),ii,jj,ex,ey);
        cq_vec(g,k,0)(ii,jj) =  vx*ex + vy*ey;
        cq_vec(g,k,1)(ii,jj) = -vx*ey + vy*ex;
    };

    for(int l=wtop(); l>=wlo+1; --l)
    block_up(l,2,7040+l,
             [&](reefamr_patch *q, int n, int kk, double *v)
             {
                 sflow_amr_patch *c = SP(q);
                 const int nby = c->ny/2;
                 const int i0 = EXT+2*(kk/nby), j0 = EXT+2*(kk%nby);
                 v[0] = v[1] = 0.0;
                 for(int a=0; a<2; ++a)
                 for(int d=0; d<2; ++d)
                 {
                     double vx,vy;
                     vec(n,i0+a,j0+d,vx,vy);
                     v[0] += 0.25*vx;
                     v[1] += 0.25*vy;
                 }
             },
             [&](const reefamr_block &B, int key, const double *v) { put(B.g,B.ic,B.jc,v[0],v[1]); });

    if(p0->mpi_size>1 && wlo==0)
    for(int c=0; c<2; ++c)
    {
        pgc0->gcslparax(p0,cq_vec(-1,k,c),4);
        pgc0->gcslparacox(p0,cq_vec(-1,k,c),10);
    }

    for(int l=MAX(wlo,1); l<=wtop(); ++l)
    {
    const bool edge = (nh_edge && l==wlo);
    fill_run(l,2,7300+l,
             [&](const reefamr_fill &f, double *v)
             {
                 v[0] = v[1] = 0.0;
                 if(edge && f.kind==1)
                 {
                     v[0] = std::numeric_limits<double>::quiet_NaN();
                     return;
                 }
                 if(f.kind==0)
                 vec(f.g,f.si,f.sj,v[0],v[1]);
                 if(f.kind==1)
                 {
                     lexer *q = glex(f.g);
                     const int ic=f.si, jc=f.sj;
                     double x0,y0;
                     vec(f.g,ic,jc,x0,y0);
                     auto val = [&](int a, int bb, double &vx, double &vy)
                     {
                         if(q->flagslice4[lij(q,a,bb)]<0) { vx=x0; vy=y0; }
                         else vec(f.g,a,bb,vx,vy);
                     };
                     double xe,ye,xw,yw,xn,yn,xs,ys;
                     val(ic+1,jc,xe,ye); val(ic-1,jc,xw,yw);
                     val(ic,jc+1,xn,yn); val(ic,jc-1,xs,ys);
                     v[0] = x0 + 0.125*(xe-xw)*f.ox + 0.125*(xn-xs)*f.oy;
                     v[1] = y0 + 0.125*(ye-yw)*f.ox + 0.125*(yn-ys)*f.oy;
                 }
             },
             [&](reefamr_patch *c, int n, const reefamr_fill &f, const double *w)
             {
                 if(edge && std::isnan(w[0]))
                 {
                     // G 7 1: the parent values around the window patches stay in the start
                     if(k!=nh_keep)
                     cq_vec(n,k,0)(f.di,f.dj) = cq_vec(n,k,1)(f.di,f.dj) = 0.0;
                     return;
                 }
                 put(n,f.di,f.dj,w[0],w[1]);
             });
    }
}

// y = own component (dir) of A x on the computed cells of the window grids; lin: the other
// component of x 0, the constant part C0 taken off
void sflow_amr::cq_apply(int kx, int ky, bool lin, int dir)
{
    if(lin)
    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = (g<0) ? p0 : P[g]->pp;
        double *o = cq_vec(g,kx,1-dir).V;
        for(int n=0; n<q->imax*q->jmax; ++n)
        o[n] = 0.0;
    }

    cq_sync(kx);

    // A x0: the lagged component outside the leaf cells as in the stage (the state around the
    // patches, the restricted state in the covered cells), as the Cartesian line systems lag it
    if(kx==CQT)
    for(int g=-1; g<(int)P.size(); ++g)
    {
        if((*bqc)[g+1].leaf.empty())
        continue;
        lexer *q = (g<0) ? p0 : P[g]->pp;
        const vector<signed char> &act = *nh_grid(g).act;
        double *o = cq_vec(g,kx,1-dir).V;
        const double *os = cq_vec(g,CQO,1-dir).V;
        for(int n=0; n<q->imax*q->jmax; ++n)
        if(act[n]!=1)
        o[n] = os[n];
    }

    for(int g=-1; g<(int)P.size(); ++g)
    {
        // (level 0 on every rank while in the window: its operator exchanges)
        if(g<0 ? wlo>0 : (*bqc)[g+1].leaf.empty())
        continue;
        lexer *q = (g<0) ? p0 : P[g]->pp;
        fdm2D *b = (g<0) ? b0 : SP(g)->b;
        sflow_boussinesq *bq = (g<0 ? pmom0 : SP(g)->pmom)->boussinesq();
        slice &ya = (dir==0) ? cq_vec(g,ky,0) : cq_vec(g,CQS,0);
        slice &yb = (dir==1) ? cq_vec(g,ky,1) : cq_vec(g,CQS,1);
        slice &y = (dir==0) ? ya : yb;

        // patches: no partition exchange (level 0 exchanges, on all ranks in the same order)
        if(g<0)
        bq->apply_curv(q,b,pgc0,cq_vec(g,kx,0),cq_vec(g,kx,1),ya,yb,false);
        else
        {
        comms_off guard(pgc0);
        bq->apply_curv(q,b,pgc0,cq_vec(g,kx,0),cq_vec(g,kx,1),ya,yb,true);
        }

        // own component of ky
        if(dir==1)
        for(int n : (*bqc)[g+1].leaf)
        cq_vec(g,ky,1).V[n] = y.V[n];

        if(lin)
        {
            const double *c0 = cq_vec(g,CQC0,dir).V;
            double *yy = cq_vec(g,ky,dir).V;
            for(int n : (*bqc)[g+1].leaf)
            yy[n] -= c0[n];
        }
    }
}

double sflow_amr::cq_dot(int ka, int kb, int dir)
{
    double s=0.0;
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *a = cq_vec(g,ka,dir).V, *b = cq_vec(g,kb,dir).V;
        for(int n : (*bqc)[g+1].leaf)
        s += a[n]*b[n];
    }
    return pgc0->globalsum(s);
}

// z = M^-1 r: the line solves of the rows of a (dir 0) or b (dir 1)
void sflow_amr::cq_prec(int kr, int kz, int dir)
{
    bq_prec(kr+dir*CQN,kz+dir*CQN,dir);
}

// element-wise operations of reefamr_bicgstab on the leaf cells, component dir
void sflow_amr::cq_ops(int op, double s1, double s2, int dir)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        auto V = [&](int k) { return cq_vec(g,k,dir).V; };
        double *x=V(CQX), *bb=V(CQB), *r=V(CQR), *rh=V(CQRH), *pv=V(CQPV), *vv=V(CQVV), *sv=V(CQS), *t=V(CQT), *ph=V(CQPH), *sh=V(CQSH);
        for(int n : (*bqc)[g+1].leaf)
        {
            if(op==0)      { double rr = bb[n]-vv[n]; r[n]=rr; rh[n]=rr; pv[n]=0.0; vv[n]=0.0; }
            else if(op==1) pv[n] = r[n] + s1*(pv[n] - s2*vv[n]);
            else if(op==2) sv[n] = r[n] - s1*vv[n];
            else if(op==3) { x[n] += s1*ph[n] + s2*sh[n]; r[n] = sv[n] - s2*t[n]; }
            else           { rh[n] = r[n]; pv[n] = 0.0; vv[n] = 0.0; }
        }
    }
}

namespace
{
struct cq_space
{
    static const int B=CQB, R=CQR, RH=CQRH, PV=CQPV, VV=CQVV, S=CQS, T=CQT, PH=CQPH, SH=CQSH;
    sflow_amr *a;
    int dir;
    void apply(int x, int y) { a->cq_apply(x,y,true,dir); }
    void prec(int r, int z) { a->cq_prec(r,z,dir); }
    double dot(int x, int y) { return a->cq_dot(x,y,dir); }
    void op_start() { a->cq_ops(0,0.0,0.0,dir); }
    void op_p(double beta, double om) { a->cq_ops(1,beta,om,dir); }
    void op_s(double alp) { a->cq_ops(2,alp,0.0,dir); }
    void op_x(double alp, double om) { a->cq_ops(3,alp,om,dir); }
    void op_restart() { a->cq_ops(4,0.0,0.0,dir); }
};
}

void sflow_amr::bous_window_curv(ghostcell *pgc, int l, double th)
{
    const double t0 = MPI_Wtime();

    if(l>=0)
    {
        wlo = whi = l;
        nh_edge = true;
    }
    nh_keep = CQT;

    bq_setup();

    // the pairs of vectors
    while((int)nh0_v.size()<2*CQN)
    nh0_v.push_back(new slice4(p0));
    for(auto q : P)
    {
        sflow_amr_patch *c = SP(q);
        while((int)c->nv.size()<2*CQN)
        c->nv.push_back(new slice4(c->pp));
    }

    const int ng = (int)P.size()+1;
    auto bous = [&](int g) { return (g<0 ? pmom0 : SP(g)->pmom)->boussinesq(); };
    auto inwin = [&](int g) { const int lg = (g<0) ? 0 : P[g]->lev; return lg>=wlo && lg<=wtop(); };
    auto cells = [&](int g) { lexer *q = (g<0) ? p0 : P[g]->pp; return q->imax*q->jmax; };

    // start: a, b of prepare_curv (UA, VA hold them); G 7 1: the parent u_a around the level-l
    // patches (the components along the basis of the cell) for the own component
    for(int g=-1; g<ng-1; ++g)
    {
        fdm2D *b = (g<0) ? b0 : SP(g)->b;
        for(int n=0; n<cells(g); ++n)
        {
            cq_vec(g,CQX,0).V[n] = b->UA.V[n];
            cq_vec(g,CQX,1).V[n] = b->VA.V[n];
        }
    }

    // the lagged component: the start values (the stage state also around the patches)
    for(int g=-1; g<ng-1; ++g)
    for(int c=0; c<2; ++c)
    for(int n=0; n<cells(g); ++n)
    cq_vec(g,CQO,c).V[n] = cq_vec(g,CQX,c).V[n];

    if(l>=0)
    fill_run(l,2,7360+l,
             [&](const reefamr_fill &f, double *v) { v[0] = nh_eval_t(f,th,5); v[1] = nh_eval_t(f,th,6); },
             [&](reefamr_patch *c, int n, const reefamr_fill &f, const double *w)
             {
                 double ex,ey;
                 bous(n)->basis(SP(n)->pp,f.di,f.dj,ex,ey);
                 cq_vec(n,CQX,0)(f.di,f.dj) =  w[0]*ex + w[1]*ey;
                 cq_vec(n,CQX,1)(f.di,f.dj) = -w[0]*ey + w[1]*ex;
             });

    int itsum = 0;

    for(int dir=0; dir<2; ++dir)
    {
        // right-hand side V.e/H of the leaf cells
        for(int g=-1; g<ng-1; ++g)
        {
            const double *rh = (dir==0) ? bous(g)->Ra.V : bous(g)->Rb.V;
            for(int n : (*bqc)[g+1].leaf)
            cq_vec(g,CQB,dir).V[n] = rh[n];
        }

        // constant part C0 = A(0), A x0 with the lagged other component
        for(int g=-1; g<ng-1; ++g)
        for(int c=0; c<2; ++c)
        for(int n=0; n<cells(g); ++n)
        {
            cq_vec(g,CQPH,c).V[n] = 0.0;
            cq_vec(g,CQT,c).V[n] = (c==dir) ? cq_vec(g,CQX,c).V[n] : cq_vec(g,CQO,c).V[n];
        }
        cq_apply(CQPH,CQC0,false,dir);
        cq_apply(CQT,CQVV,false,dir);

        double bn, rn;
        cq_space sp{this,dir};
        int it = reefamr_bicgstab(sp,p0->N44,p0->N46,bn,rn);
        itsum += it;

        // at rest: 0
        if(bn<=0.0)
        for(int g=-1; g<ng-1; ++g)
        for(int n : (*bqc)[g+1].leaf)
        cq_vec(g,CQX,dir).V[n] = 0.0;
    }

    // u_a everywhere: leaf cells, covered cells, around the patches, partition halo; dry cells 0
    nh_keep = CQX;
    cq_sync(CQX);

    for(int g=-1; g<ng-1; ++g)
    {
        if(!inwin(g))
        continue;
        lexer *q = (g<0) ? p0 : P[g]->pp;
        fdm2D *b = (g<0) ? b0 : SP(g)->b;
        for(int n=0; n<cells(g); ++n)
        if(q->wet[n]==0)
        cq_vec(g,CQX,0).V[n] = cq_vec(g,CQX,1).V[n] = 0.0;
        bous(g)->set_curv(q,b,cq_vec(g,CQX,0),cq_vec(g,CQX,1),g>=0);
    }

    nh_keep = -1;
    if(l>=0)
    {
        wlo = 0;
        whi = -1;
        nh_edge = false;
        sub_lv_it += itsum;
        ++sub_lv_n;
    }

    bq_it_total += itsum;
    ++bq_solves;
    tm[9] += MPI_Wtime()-t0;
}
