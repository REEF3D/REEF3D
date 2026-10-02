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

#include"nhflow_amr.h"
#include"nhflow_amr_fill.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_pressure.h"
#include"reefmg_core.h"
#include"reefamr_krylov.h"
#include<cmath>
#include<mpi.h>
#include<iomanip>
#include<algorithm>

//  Composite pressure projection of REEF3D::NHFLOW on level 0 and the patches (A 520 1, 2).
//
//  Every grid assembles its own Poisson rows with its pressure object (nhflow_pjm or
//  nhflow_pjm_corr: amr_prepare, the rows of nhflow_poisson / nhflow_poisson_pcorr with the sigma
//  metrics of the grid).  The unknowns are the nodes of the leaf columns (level-0 and patch
//  columns not covered by a finer patch).  A leaf row is the grid's own matrix row; its neighbour
//  columns outside the grid come from the finer level (covered columns: restricted from the
//  children) or from the coarser level (the columns around a patch: bicubic in x,y, node by
//  node; a sibling patch or the owner rank across a partition edge).  Identity rows (dry or
//  shallow cells) are no unknowns.
//
//  BiCGStab (reefamr_bicgstab), right preconditioned by one FAC sweep: a REEFMG V-cycle on level 0
//  (the whole rank grid, covered columns included), then level by level the coarse correction
//  interpolated into the patches and patch-local REEFMG V-cycles (MPI_COMM_SELF, correction 0 at
//  the patch edge) on the residual it leaves.  The same structure as the composite Laplace of
//  FNPF AMR (fnpf_amr_lap.cpp).

using namespace nhflow_amr_detail;

namespace
{
typedef reefamr_comms_off comms_off;

enum { NR=REEFAMR_NR, NRH=REEFAMR_NRH, NPV=REEFAMR_NPV, NVV=REEFAMR_NVV, NS=REEFAMR_NS, NT=REEFAMR_NT,
       NPH=REEFAMR_NPH, NSH=REEFAMR_NSH, NVEC=REEFAMR_NVEC };
}

double* nhflow_amr::pvec(int g, int k)
{
    if(k<0)
    return (g<0) ? ptgt0 : NP(g)->ptgt;

    return (g<0) ? kv0[k] : NP(g)->kv[k];
}

// row numbering as the assembly in nhflow_poisson (LOOP: ILOOP JLOOP KLOOP, flag4>0); leaf and
// covered unknowns of every grid
void nhflow_amr::pr_rows()
{
    auto number = [&](lexer *q, vector<int> &row)
    {
        row.assign(q->imax*q->jmax*(q->kmax+2),-1);
        int n=0;
        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        for(int kk=0; kk<q->knoz; ++kk)
        if(q->flag4[cidx(q,ii,jj,kk)]>0)
        row[fidx(q,ii,jj,kk)] = n++;
    };

    number(p0,rowmap0);
    for(auto q : P)
    number(q->pp,NP(q)->row);

    pg.assign(P.size()+1,pgrid());

    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        vector<int> &row = (g<0) ? rowmap0 : NP(g)->row;
        pgrid &L = pg[g+1];
        const int l = (g<0) ? 0 : P[g]->lev;
        int i0,i1,j0,j1;
        if(g<0) { i0=0; i1=NX0-1; j0=0; j1=NY0-1; }
        else    { i0=EXT; i1=EXT+P[g]->nx-1; j0=EXT; j1=EXT+P[g]->ny-1; }

        for(int ii=i0; ii<=i1; ++ii)
        for(int jj=j0; jj<=j1; ++jj)
        {
            int I = (g<0) ? ii+O0i : ii-EXT+P[g]->I0;
            int J = (g<0) ? jj+O0j : jj-EXT+P[g]->J0;
            const bool covered = (l<maxlev && patch_at(l+1,2*I,2*J)>=0);

            for(int kk=0; kk<q->knoz; ++kk)
            {
                const int qq = fidx(q,ii,jj,kk);
                const int r = row[qq];
                if(r<0)
                continue;

                if(covered)
                L.cq.push_back(qq);
                else
                {
                    L.aq.push_back(qq);
                    L.ar.push_back(r);
                }
            }
        }
    }
}

// vectors, multigrids and their coefficients for this solve
void nhflow_amr::pr_prepare(ghostcell *pgc)
{
    if(pr_layout!=layout_id)
    {
        pr_rows();

        const int n0 = p0->imax*p0->jmax*(p0->kmax+2);
        if(kv0.empty())
        for(int k=0; k<NVEC; ++k)
        kv0.push_back(new double[n0]());

        for(auto q : P)
        {
            nhflow_amr_patch *c = NP(q);
            const int n = c->pp->imax*c->pp->jmax*(c->pp->kmax+2);
            if(c->kv.empty())
            for(int k=0; k<NVEC; ++k)
            c->kv.push_back(new double[n]());

            if(c->mg==nullptr)
            {
                lexer *pp = c->pp;
                c->mg = new reefmg_core;
                c->mg->set_precision(64);
                if(!c->mg->setup(MPI_COMM_SELF,c->nx,c->ny,pp->knoz,c->nx,c->ny,0,pp->DXN+marge+EXT-1,pp->DYN+marge+EXT-1))
                cout<<"NHFLOW AMR: patch multigrid "<<c->mg->err()<<endl;
            }
        }

        if(mg0==nullptr)
        {
            mg0 = new reefmg_core;
            mg0->set_precision(64);
            if(!mg0->setup(pgc->cart(),p0->knox,p0->knoy,p0->knoz,p0->gknox,p0->gknoy,p0->N13,p0->DXN+marge-1,p0->DYN+marge-1))
            {
                if(p0->mpirank==0)
                cout<<"NHFLOW AMR: level-0 multigrid "<<mg0->err()<<endl;
                MPI_Abort(MPI_COMM_WORLD,-2760);
            }
        }

        pr_layout = layout_id;
    }

    // identity rows (dry or shallow cells: P = 0) are no unknowns
    auto fixed = [&](const matrix_diag &M, int r)
    {
        return M.p[r]==1.0 && M.n[r]==0.0 && M.s[r]==0.0 && M.w[r]==0.0 && M.e[r]==0.0 && M.t[r]==0.0 && M.b[r]==0.0;
    };

    for(int g=-1; g<(int)P.size(); ++g)
    {
        pgrid &L = pg[g+1];
        const matrix_diag &M = gfd(g)->M;
        L.lq.clear();
        L.lr.clear();
        for(size_t n=0; n<L.aq.size(); ++n)
        if(!fixed(M,L.ar[n]))
        {
            L.lq.push_back(L.aq[n]);
            L.lr.push_back(L.ar[n]);
        }
    }

    auto clear = [](sc_level &L)
    {
        std::fill(L.p.begin(),L.p.end(),0.0);
        std::fill(L.n.begin(),L.n.end(),0.0); std::fill(L.s.begin(),L.s.end(),0.0);
        std::fill(L.w.begin(),L.w.end(),0.0); std::fill(L.e.begin(),L.e.end(),0.0);
        std::fill(L.t.begin(),L.t.end(),0.0); std::fill(L.b.begin(),L.b.end(),0.0);
        std::fill(L.u.begin(),L.u.end(),0.0); std::fill(L.f.begin(),L.f.end(),0.0);
        std::fill(L.act.begin(),L.act.end(),0);
    };

    // level 0: all rows of the rank grid
    {
        sc_level &L = mg0->fine();
        clear(L);
        const matrix_diag &M = d0->M;

        for(int ii=0; ii<p0->knox; ++ii)
        for(int jj=0; jj<p0->knoy; ++jj)
        for(int kk=0; kk<p0->knoz; ++kk)
        {
            const long lq = L.idx(ii,jj,kk);
            const int r = rowmap0[fidx(p0,ii,jj,kk)];
            if(r<0)
            {
                L.p[lq] = 1.0;
                continue;
            }
            L.p[lq]=M.p[r]; L.n[lq]=M.n[r]; L.s[lq]=M.s[r]; L.w[lq]=M.w[r]; L.e[lq]=M.e[r]; L.t[lq]=M.t[r]; L.b[lq]=M.b[r];
            L.act[lq] = fixed(M,r) ? 0 : 1;
        }
        mg0->coarsen();
    }

    // patches: interior rows, couplings to the cells around the patch dropped
    for(auto q : P)
    {
        nhflow_amr_patch *c = NP(q);
        lexer *pp = c->pp;
        sc_level &L = c->mg->fine();
        clear(L);
        const matrix_diag &M = c->d->M;

        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
        for(int kk=0; kk<pp->knoz; ++kk)
        {
            const long lq = L.idx(ii-EXT,jj-EXT,kk);
            const int r = c->row[fidx(pp,ii,jj,kk)];
            if(r<0)
            {
                L.p[lq] = 1.0;
                continue;
            }
            L.p[lq] = M.p[r];
            L.n[lq] = (ii+1<EXT+c->nx) ? M.n[r] : 0.0;
            L.s[lq] = (ii-1>=EXT) ? M.s[r] : 0.0;
            L.w[lq] = (jj+1<EXT+c->ny) ? M.w[r] : 0.0;
            L.e[lq] = (jj-1>=EXT) ? M.e[r] : 0.0;
            L.t[lq] = M.t[r];
            L.b[lq] = M.b[r];
            L.act[lq] = fixed(M,r) ? 0 : 1;
        }
        c->mg->coarsen();
    }
}

// covered columns, partition halo of level 0, columns around the patches
void nhflow_amr::pr_sync(int k)
{
    auto sel = [&](int g) -> double* { return pvec(g,k); };

    restrict_col(sel);

    if(p0->mpi_size>1)
    {
        double *x0 = pvec(-1,k);
        pgc0->gcparax7(p0,x0,7);
        pgc0->gcparax7co(p0,x0,7);
    }

    for(int l=1; l<=maxlev; ++l)
    fill_col(l,7600+l,sel);
}

// y = A x on the leaf unknowns
void nhflow_amr::pr_apply(int kx, int ky)
{
    pr_sync(kx);

    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        const matrix_diag &M = gfd(g)->M;
        const double *x = pvec(g,kx);
        double *y = pvec(g,ky);
        const int sI = q->jmax*q->kmaxF;
        const int sJ = q->kmaxF;
        const pgrid &L = pg[g+1];

        const double *Mp=M.p.data(), *Mn=M.n.data(), *Ms=M.s.data(), *Mw=M.w.data(), *Me=M.e.data(), *Mt=M.t.data(), *Mb=M.b.data();

        for(size_t n=0; n<L.lq.size(); ++n)
        {
            const int qq = L.lq[n], r = L.lr[n];
            y[qq] = Mp[r]*x[qq] + Mn[r]*x[qq+sI] + Ms[r]*x[qq-sI] + Mw[r]*x[qq+sJ] + Me[r]*x[qq-sJ] + Mt[r]*x[qq+1] + Mb[r]*x[qq-1];
        }

        for(int qq : L.cq)
        y[qq] = 0.0;
    }
}

double nhflow_amr::pr_dot(int ka, int kb)
{
    double s=0.0;
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *a = pvec(g,ka), *b = pvec(g,kb);
        for(int qq : pg[g+1].lq)
        s += a[qq]*b[qq];
    }
    return pgc0->globalsum(s);
}

void nhflow_amr::pr_start()
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vb = pvec(g,NS);
        double *vr = pvec(g,NR), *vrh = pvec(g,NRH), *vp = pvec(g,NPV), *vv = pvec(g,NVV);
        for(int qq : pg[g+1].lq)
        {
            double r = vb[qq] - vv[qq];
            vr[qq] = r; vrh[qq] = r; vp[qq] = 0.0; vv[qq] = 0.0;
        }
    }
}

void nhflow_amr::pr_p(double beta, double om)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vr = pvec(g,NR), *vv = pvec(g,NVV);
        double *vp = pvec(g,NPV);
        for(int qq : pg[g+1].lq)
        vp[qq] = vr[qq] + beta*(vp[qq] - om*vv[qq]);
    }
}

void nhflow_amr::pr_s(double alp)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vr = pvec(g,NR), *vv = pvec(g,NVV);
        double *vs = pvec(g,NS);
        for(int qq : pg[g+1].lq)
        vs[qq] = vr[qq] - alp*vv[qq];
    }
}

void nhflow_amr::pr_x(double alp, double om)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vph = pvec(g,NPH), *vsh = pvec(g,NSH), *vs = pvec(g,NS), *vt = pvec(g,NT);
        double *x = pvec(g,-1), *vr = pvec(g,NR);
        for(int qq : pg[g+1].lq)
        {
            x[qq] += alp*vph[qq] + om*vsh[qq];
            vr[qq] = vs[qq] - om*vt[qq];
        }
    }
}

// z = M^-1 r: one FAC sweep
void nhflow_amr::pr_prec(int kr, int kz)
{
    restrict_col([&](int g) -> double* { return pvec(g,kr); });

    // level 0: one V-cycle on the whole rank grid
    {
        sc_level &L = mg0->fine();
        const double *r = pvec(-1,kr);
        double *z = pvec(-1,kz);

        std::fill(L.u.begin(),L.u.end(),0.0);
        std::fill(L.f.begin(),L.f.end(),0.0);
        for(int ii=0; ii<p0->knox; ++ii)
        for(int jj=0; jj<p0->knoy; ++jj)
        for(int kk=0; kk<p0->knoz; ++kk)
        {
            const long lq = L.idx(ii,jj,kk);
            if(L.act[lq])
            L.f[lq] = r[fidx(p0,ii,jj,kk)];
        }

        mg0->vcycle(0,1,1);

        for(int ii=0; ii<p0->knox; ++ii)
        for(int jj=0; jj<p0->knoy; ++jj)
        for(int kk=0; kk<p0->knoz; ++kk)
        {
            const long lq = L.idx(ii,jj,kk);
            z[fidx(p0,ii,jj,kk)] = L.act[lq] ? L.u[lq] : 0.0;
        }

        if(p0->mpi_size>1)
        {
            pgc0->gcparax7(p0,z,7);
            pgc0->gcparax7co(p0,z,7);
        }
    }

    auto selz = [&](int g) -> double* { return pvec(g,kz); };

    for(int l=1; l<=maxlev; ++l)
    {
        // coarse correction interpolated into the patch interiors, then the columns around them
        for(int id : lev[l])
        prolong_interior_col(*NP(id),selz);

        fill_col(l,7700+l,selz);

        // patch-local corrections of the residual left by the coarse correction
        for(int id : lev[l])
        {
            nhflow_amr_patch *c = NP(id);
            lexer *pp = c->pp;
            const matrix_diag &M = c->d->M;
            const double *r = pvec(id,kr);
            double *z = pvec(id,kz);
            sc_level &L = c->mg->fine();
            const int sI = pp->jmax*pp->kmaxF;
            const int sJ = pp->kmaxF;

            std::fill(L.u.begin(),L.u.end(),0.0);
            std::fill(L.f.begin(),L.f.end(),0.0);
            for(int ii=EXT; ii<EXT+c->nx; ++ii)
            for(int jj=EXT; jj<EXT+c->ny; ++jj)
            for(int kk=0; kk<pp->knoz; ++kk)
            {
                const long lq = L.idx(ii-EXT,jj-EXT,kk);
                if(L.act[lq]==0)
                continue;

                const int qq = fidx(pp,ii,jj,kk);
                const int rw = c->row[qq];
                const double az = M.p[rw]*z[qq] + M.n[rw]*z[qq+sI] + M.s[rw]*z[qq-sI] + M.w[rw]*z[qq+sJ] + M.e[rw]*z[qq-sJ]
                                + M.t[rw]*z[qq+1] + M.b[rw]*z[qq-1];
                L.f[lq] = r[qq] - az;
            }

            c->mg->vcycle(0,1,1);

            for(int ii=EXT; ii<EXT+c->nx; ++ii)
            for(int jj=EXT; jj<EXT+c->ny; ++jj)
            for(int kk=0; kk<pp->knoz; ++kk)
            {
                const long lq = L.idx(ii-EXT,jj-EXT,kk);
                if(L.act[lq])
                z[fidx(pp,ii,jj,kk)] += L.u[lq];
            }
        }

        if(l<maxlev)
        fill_col(l,7700+l,selz);
    }
}

namespace
{
// the composite pressure equation as the vector space of reefamr_bicgstab
struct pr_space
{
    static const int B=NS, R=NR, RH=NRH, PV=NPV, VV=NVV, S=NS, T=NT, PH=NPH, SH=NSH;
    nhflow_amr *a;
    void apply(int x, int y) { a->pr_apply(x,y); }
    void prec(int r, int z) { a->pr_prec(r,z); }
    double dot(int x, int y) { return a->pr_dot(x,y); }
    void op_start() { a->pr_start(); }
    void op_p(double beta, double om) { a->pr_p(beta,om); }
    void op_s(double alp) { a->pr_s(alp); }
    void op_x(double alp, double om) { a->pr_x(alp,om); }
};
}

// assembled rows of all grids -> solution in the unknowns of the grids (initial guess: their
// values), then the covered columns, the columns around the patches and the halo
void nhflow_amr::pr_core(lexer *p, ghostcell *pgc)
{
    pr_prepare(pgc);

    // right-hand side into NS
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *R = gfd(g)->rhsvec.V.data();
        double *b = pvec(g,NS);
        const pgrid &L = pg[g+1];
        for(size_t n=0; n<L.lq.size(); ++n)
        b[L.lq[n]] = R[L.lr[n]];
        for(int qq : L.cq)
        b[qq] = 0.0;
    }

    pr_apply(-1,NVV);

    double bn, rn;
    pr_space sp{this};
    int it = reefamr_bicgstab(sp,p->N44,p->N46,bn,rn,&tm[5],&tm[6]);

    pr_it_last = it;
    pr_it_total += it;
    ++pr_solves;
    pr_res_last = bn>0.0 ? rn/bn : 0.0;
    p->solveriter = it;
    p->poissoniter = it;
    p->final_res = pr_res_last;

    pr_sync(-1);
}

// the projection of stage s on all grids: rows and right-hand sides of every grid, the composite
// solve, then the pressure and velocity update of every grid
void nhflow_amr::press_solve(lexer *p, ghostcell *pgc, int s)
{
    if(p->A520==0)
    return;

    const double t0 = MPI_Wtime();
    const double alpha = mom0->stage_alpha(s);

    {
    comms_off guard(pgc);
    for(auto q : P)
    {
        nhflow_amr_patch *c = NP(q);
        pscope ps(pgc,c->d,d0);
        c->ptgt = c->ppress->amr_prepare(c->pp,c->d,pgc,alpha);
    }
    }

    ptgt0 = S0p->ppress->amr_prepare(p,d0,pgc,alpha);

    pr_core(p,pgc);

    {
    comms_off guard(pgc);
    for(int n=0; n<(int)P.size(); ++n)
    {
        nhflow_amr_patch *c = NP(n);
        stg O = stage_out(n,s);
        pscope ps(pgc,c->d,d0);
        c->ppress->amr_finish(c->pp,c->d,pgc,*O.WL,O.UH,O.VH,O.WH,alpha);
    }
    }

    stg O = stage_out(-1,s);
    S0p->ppress->amr_finish(p,d0,pgc,*O.WL,O.UH,O.VH,O.WH,alpha);

    tm[4] += MPI_Wtime()-t0;

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"NHFLOW AMR pressure: iterations "<<pr_it_last<<"  res "<<setprecision(3)<<pr_res_last<<endl;
}
