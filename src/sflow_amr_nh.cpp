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
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"slice.h"
#include"slice4.h"
#include"vec2D.h"
#include"matrix2D.h"
#include"sflow_pressure_nh.h"
#include"sflow_momentum_RK3.h"
#include"reefmg_core.h"
#include"reefmg2D.h"
#include<cmath>
#include<mpi.h>
#include<iomanip>

//  Composite non-hydrostatic pressure (sflow_pjm_lin A 220 1, sflow_pjm_quad A 220 2/3) on level 0 and the patches.
//
//  Every grid assembles its own rows with sflow_pjm_lin (level 0 before the call, the patches
//  in nh_prepare).  The unknowns are the leaf cells: level-0 and patch cells not covered by a
//  finer patch.  A leaf row is the grid's own matrix row; its neighbours outside the grid come
//  from the finer or coarser level as for the flow variables (sibling copy, linear interpolation
//  from the parent, the owner rank across a partition edge).  A coarse face next to a patch
//  takes the mean of the two fine face gradients of q, sent across partition edges like the
//  fluxes.  Covered cells carry the average of their children.
//
//  BiCGStab, right preconditioned by one FAC sweep: a REEFMG V-cycle on level 0 (the whole rank
//  grid, covered cells included), then level by level the prolonged coarse correction and
//  patch-local REEFMG V-cycles (MPI_COMM_SELF, correction 0 at the patch edge) on the residual
//  it leaves, repeated NHPASS times with the siblings' corrections filled in between.  The
//  solve ends with the velocity correction on the patches and the rest of their stage; level 0
//  corrects itself in sflow_pjm_lin.

namespace
{
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

struct comms_off
{
    ghostcell *g;
    bool old;
    comms_off(ghostcell *gg) : g(gg) { old = g->set_comms(false); }
    ~comms_off() { g->set_comms(old); }
};

// Krylov vectors per grid
enum { NR=0, NRH, NPV, NVV, NS, NT, NPH, NSH, NTMP, NVEC };

// patch-local passes per level in the preconditioner: with one pass each patch corrects
// against the coarse interpolation only, and the iteration count grows with the number of
// sibling interfaces (dam break: 21 iterations on average; 2 passes: 9; 4 passes: 6.5,
// for the same run time as 2).  A level with a single patch needs one pass.
const int NHPASS = 4;
}

sflow_amr::nhg sflow_amr::nh_grid(int g)
{
    nhg r;
    if(g<0)
    {
        r.q = p0; r.b = b0; r.act = &nh0_act; r.row = &nh0_row; r.v = &nh0_v;
        return r;
    }
    sflow_amr_patch *c = P[g];
    r.q = c->pp; r.b = c->b; r.act = &c->act; r.row = &c->row; r.v = &c->nv;
    return r;
}

slice& sflow_amr::nh_vec(int g, int k)
{
    if(k<0)
    return (g<0) ? b0->press : P[g]->b->press;

    return (g<0) ? *nh0_v[k] : *P[g]->nv[k];
}

// rows, activity and the patch multigrids of this solve
void sflow_amr::nh_prepare(ghostcell *pgc)
{
    // vectors
    if(nh0_v.empty())
    for(int k=0; k<NVEC; ++k)
    nh0_v.push_back(new slice4(p0));

    if(nhmg0==nullptr)
    {
        nhmg0 = new reefmg2D(p0,pgc0);
        nhr0 = new vec2D(p0,p0->imax*p0->jmax);
    }
    nh_rebuild0 = true;

    {
    comms_off guard(pgc);

    for(auto c : P)
    {
        if(c->nv.empty())
        for(int k=0; k<NVEC; ++k)
        c->nv.push_back(new slice4(c->pp));

        // the patch rows (right-hand side from the state after the momentum update)
        c->pnh->assemble(c->pp,c->b,pgc,*c->pmom->nhUH,*c->pmom->nhVH,*c->pmom->nhWL,*c->pmom->nhUn,*c->pmom->nhVn,c->pmom->nh_alpha);
    }
    }

    // row numbers and activity
    auto rows = [&](int g, int i0, int i1, int j0, int j1, int l)
    {
        nhg G = nh_grid(g);
        lexer *q = G.q;
        const int nsl = q->imax*q->jmax;
        G.act->assign(nsl,-2);
        G.row->assign(nsl,-1);

        int n=0;
        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        if(q->flagslice4[lij(q,ii,jj)]>0)
        {
            (*G.row)[lij(q,ii,jj)] = n;
            ++n;
        }

        matrix2D &M = G.b->M;
        for(int ii=i0; ii<=i1; ++ii)
        for(int jj=j0; jj<=j1; ++jj)
        {
            int c = lij(q,ii,jj);
            int r = (*G.row)[c];
            if(r<0)
            continue;

            int I = (g<0) ? ii+O0i : ii-EXT+P[g]->I0;
            int J = (g<0) ? jj+O0j : jj-EXT+P[g]->J0;

            if(l<maxlev && patch_at(l+1,2*I,2*J)>=0)
            {
                (*G.act)[c] = -1;
                continue;
            }

            (*G.act)[c] = (M.n[r]==0.0 && M.s[r]==0.0 && M.w[r]==0.0 && M.e[r]==0.0) ? 0 : 1;
        }
    };

    rows(-1,0,NX0-1,0,NY0-1,0);
    for(int n=0; n<(int)P.size(); ++n)
    rows(n,EXT,EXT+P[n]->nx-1,EXT,EXT+P[n]->ny-1,P[n]->lev);

    // patch multigrids: the patch rows with q = 0 at the patch edge
    for(int n=0; n<(int)P.size(); ++n)
    {
        sflow_amr_patch *c = P[n];
        lexer *pp = c->pp;
        matrix2D &M = c->b->M;

        if(c->mg==nullptr)
        {
            c->mg = new reefmg_core;
            c->mg->set_precision(64);
            if(!c->mg->setup(MPI_COMM_SELF,c->nx,c->ny,1,c->nx,c->ny,0,pp->DXN+marge+EXT-1,pp->DYN+marge+EXT-1))
            cout<<"SFLOW AMR: patch multigrid "<<c->mg->err()<<endl;
        }

        sc_level &L = c->mg->fine();
        std::fill(L.p.begin(),L.p.end(),0.0);
        std::fill(L.n.begin(),L.n.end(),0.0); std::fill(L.s.begin(),L.s.end(),0.0);
        std::fill(L.w.begin(),L.w.end(),0.0); std::fill(L.e.begin(),L.e.end(),0.0);
        std::fill(L.t.begin(),L.t.end(),0.0); std::fill(L.b.begin(),L.b.end(),0.0);
        std::fill(L.u.begin(),L.u.end(),0.0); std::fill(L.f.begin(),L.f.end(),0.0);
        std::fill(L.act.begin(),L.act.end(),0);

        double ldx = 0.0;
        auto inner = [&](int ii, int jj) { return ii>=EXT && ii<EXT+c->nx && jj>=EXT && jj<EXT+c->ny; };

        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
        {
            long qd = L.idx(ii-EXT,jj-EXT,0);
            int r = c->row[lij(pp,ii,jj)];

            if(r<0 || c->act[lij(pp,ii,jj)]==0 || (M.n[r]==0.0 && M.s[r]==0.0 && M.w[r]==0.0 && M.e[r]==0.0))
            {
                L.p[qd] = 1.0;
                continue;
            }

            double cf[4] = {M.n[r],M.s[r],M.w[r],M.e[r]};
            const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};
            for(int d=0; d<4; ++d)
            {
                int a = ii+di[d], bb = jj+dj[d];
                int rn = inner(a,bb) ? c->row[lij(pp,a,bb)] : -1;
                if(rn<0 || (M.n[rn]==0.0 && M.s[rn]==0.0 && M.w[rn]==0.0 && M.e[rn]==0.0))
                cf[d] = 0.0;
            }

            const double mp = M.p[r];
            double react = mp+M.n[r]+M.s[r]+M.w[r]+M.e[r];
            double diff = fabs(cf[0])+fabs(cf[1])+fabs(cf[2])+fabs(cf[3]);
            if(react>1.0e-12*fabs(mp))
            ldx = MAX(ldx,sqrt(diff/(4.0*react)));

            L.p[qd]=mp; L.n[qd]=cf[0]; L.s[qd]=cf[1]; L.w[qd]=cf[2]; L.e[qd]=cf[3];
            L.act[qd]=1;
        }

        int need = 2 + (int)ceil(log2(ldx>1.0 ? ldx : 1.0));
        c->mg->set_active_levels(need);
        c->mg->coarsen();
    }
}

// covered cells: average of their children, finest first
void sflow_amr::nh_restrict_vec(int k)
{
    for(int l=maxlev; l>=1; --l)
    for(int n : lev[l])
    {
        sflow_amr_patch *c = P[n];
        slice &x = nh_vec(n,k);
        const int nby = c->ny/2;

        for(int bi=0; bi<c->nx/2; ++bi)
        for(int bj=0; bj<nby; ++bj)
        {
            int kk = bi*nby+bj;
            int g = c->rgrid[kk];
            if(g<-1)
            continue;

            int i0 = EXT+2*bi, j0 = EXT+2*bj;
            nh_vec(g,k)(c->ric[kk],c->rjc[kk]) = 0.25*(x(i0,j0)+x(i0+1,j0)+x(i0,j0+1)+x(i0+1,j0+1));
        }
    }
}

// value of a filled cell: sibling copy or linear interpolation of the coarser level
double sflow_amr::nh_eval(const sflow_amr_fill &f, int k)
{
    if(f.kind==0)
    return nh_vec(f.g,k)(f.si,f.sj);

    if(f.kind==1)
    {
        slice &x = nh_vec(f.g,k);
        lexer *q = (f.g<0) ? p0 : P[f.g]->pp;
        const int ic=f.si, jc=f.sj;
        const double x0 = x(ic,jc);

        auto val = [&](int a, int bb) { return q->flagslice4[lij(q,a,bb)]<0 ? x0 : x(a,bb); };

        return x0 + 0.125*(val(ic+1,jc)-val(ic-1,jc))*f.ox + 0.125*(val(ic,jc+1)-val(ic,jc-1))*f.oy;
    }

    return 0.0;
}

// cells of the level-l patches outside their interior
void sflow_amr::nh_qfill(int l, int k)
{
    sflow_amr_xplan &X = gplan[l];

    for(size_t kk=0; kk<X.speer.size(); ++kk)
    {
        vector<double> &sb = X.sbuf[kk];
        sb.resize(X.sitem[kk].size());
        for(size_t m=0; m<X.sitem[kk].size(); ++m)
        sb[m] = nh_eval(gserve[l][X.sitem[kk][m]],k);
    }

    xrun(X,1,7300+l);

    for(int n : lev[l])
    {
        sflow_amr_patch *c = P[n];
        slice &x = nh_vec(n,k);
        for(auto &f : c->fill)
        {
            if(f.kind==2)
            x(f.di,f.dj) = X.rbuf[f.si][f.slot];
            else
            x(f.di,f.dj) = nh_eval(f,k);
        }
    }
}

void sflow_amr::nh_sync(int k)
{
    nh_restrict_vec(k);

    if(p0->mpi_size>1)
    {
        pgc0->gcslparax(p0,nh_vec(-1,k),4);
        pgc0->gcslparacox(p0,nh_vec(-1,k),10);
    }

    for(int l=1; l<=maxlev; ++l)
    nh_qfill(l,k);
}

// y = A x on the leaf cells
void sflow_amr::nh_apply(int kx, int ky)
{
    nh_sync(kx);

    // fine face gradients on the patch edges
    for(int n=0; n<(int)P.size(); ++n)
    {
        sflow_amr_patch *c = P[n];
        lexer *pp = c->pp;
        slice &x = nh_vec(n,kx);
        const int il=EXT-1, ih=EXT+c->nx-1, jl=EXT-1, jh=EXT+c->ny-1;

        c->qrec[0].resize(c->ny); c->qrec[1].resize(c->ny);
        c->qrec[2].resize(c->nx); c->qrec[3].resize(c->nx);

        for(int r=0; r<c->ny; ++r)
        {
            int jj=EXT+r;
            c->qrec[0][r] = (x(il+1,jj)-x(il,jj))/pp->DXP[il+marge];
            c->qrec[1][r] = (x(ih+1,jj)-x(ih,jj))/pp->DXP[ih+marge];
        }
        for(int r=0; r<c->nx; ++r)
        {
            int ii=EXT+r;
            c->qrec[2][r] = (x(ii,jl+1)-x(ii,jl))/pp->DYP[jl+marge];
            c->qrec[3][r] = (x(ii,jh+1)-x(ii,jh))/pp->DYP[jh+marge];
        }
    }

    // across partition edges
    for(int l=1; l<=maxlev; ++l)
    {
        sflow_amr_xplan &F = fplan[l];
        for(size_t k=0; k<F.speer.size(); ++k)
        {
            vector<double> &sb = F.sbuf[k];
            sb.resize(F.sitem[k].size());
            for(size_t m=0; m<F.sitem[k].size(); ++m)
            {
                int it = F.sitem[k][m];
                sflow_amr_patch *c = P[fsend[l][3*it]];
                int side = fsend[l][3*it+1], r = fsend[l][3*it+2];
                sb[m] = 0.5*(c->qrec[side][r]+c->qrec[side][r+1]);
            }
        }

        xrun(F,1,7400+l);

        size_t it=0;
        for(size_t k=0; k<F.rpeer.size(); ++k)
        for(int m=0; m<F.rcount[k]; ++m)
        {
            int tgt = frecv[l][2*it], idx = frecv[l][2*it+1];
            if(tgt>=0)
            rmatch[tgt][idx].val[5] = F.rbuf[k][m];
            ++it;
        }
    }

    // rows
    for(int g=-1; g<(int)P.size(); ++g)
    {
        nhg G = nh_grid(g);
        lexer *q = G.q;
        matrix2D &M = G.b->M;
        slice &x = nh_vec(g,kx);
        slice &y = nh_vec(g,ky);
        int i0,i1,j0,j1;
        if(g<0) { i0=0; i1=NX0-1; j0=0; j1=NY0-1; }
        else    { i0=EXT; i1=EXT+P[g]->nx-1; j0=EXT; j1=EXT+P[g]->ny-1; }

        for(int ii=i0; ii<=i1; ++ii)
        for(int jj=j0; jj<=j1; ++jj)
        {
            int c = lij(q,ii,jj);
            int a = (*G.act)[c];
            if(a<-1)
            continue;
            if(a==-1)
            {
                y(ii,jj) = 0.0;
                continue;
            }
            int r = (*G.row)[c];
            y(ii,jj) = M.p[r]*x(ii,jj) + M.n[r]*x(ii+1,jj) + M.s[r]*x(ii-1,jj) + M.w[r]*x(ii,jj+1) + M.e[r]*x(ii,jj-1);
        }

        // coarse faces next to finer patches: the fine face gradients replace the coarse one
        auto over = [&](const sflow_amr_match &m, double gf)
        {
            int ii = m.fi + (m.side==1 ? 1 : 0);
            int jj = m.fj + (m.side==3 ? 1 : 0);
            int c = lij(q,ii,jj);
            if((*G.act)[c]<0)
            return;
            int r = (*G.row)[c];
            double x0 = x(ii,jj);
            if(m.side==0) y(ii,jj) += M.n[r]*(x0 + q->DXP[ii+marge]*gf - x(ii+1,jj));
            if(m.side==1) y(ii,jj) += M.s[r]*(x0 - q->DXP[ii-1+marge]*gf - x(ii-1,jj));
            if(m.side==2) y(ii,jj) += M.w[r]*(x0 + q->DYP[jj+marge]*gf - x(ii,jj+1));
            if(m.side==3) y(ii,jj) += M.e[r]*(x0 - q->DYP[jj-1+marge]*gf - x(ii,jj-1));
        };

        for(auto &m : match[g+1])
        {
            sflow_amr_patch *c = P[m.child];
            over(m,0.5*(c->qrec[m.side][m.r]+c->qrec[m.side][m.r+1]));
        }
        for(auto &m : rmatch[g+1])
        over(m,m.val[5]);
    }
}

double sflow_amr::nh_dot(int ka, int kb)
{
    double s=0.0;
    for(int g=-1; g<(int)P.size(); ++g)
    {
        nhg G = nh_grid(g);
        lexer *q = G.q;
        slice &a = nh_vec(g,ka);
        slice &bv = nh_vec(g,kb);
        int i0,i1,j0,j1;
        if(g<0) { i0=0; i1=NX0-1; j0=0; j1=NY0-1; }
        else    { i0=EXT; i1=EXT+P[g]->nx-1; j0=EXT; j1=EXT+P[g]->ny-1; }

        for(int ii=i0; ii<=i1; ++ii)
        for(int jj=j0; jj<=j1; ++jj)
        if((*G.act)[lij(q,ii,jj)]>=0)
        s += a(ii,jj)*bv(ii,jj);
    }
    return pgc0->globalsum(s);
}

// z = M^-1 r: one FAC sweep
void sflow_amr::nh_prec(int kr, int kz)
{
    nh_restrict_vec(kr);

    // level 0: one V-cycle on the whole rank grid
    {
        slice &r = nh_vec(-1,kr);
        for(int ii=0; ii<NX0; ++ii)
        for(int jj=0; jj<NY0; ++jj)
        {
            int rw = nh0_row[lij(p0,ii,jj)];
            if(rw>=0)
            nhr0->V[rw] = r(ii,jj);
        }

        nhmg0->vcycle(p0,pgc0,nh_vec(-1,kz),b0->M,*nhr0,nh_rebuild0);
        nh_rebuild0 = false;

        if(p0->mpi_size>1)
        {
            pgc0->gcslparax(p0,nh_vec(-1,kz),4);
            pgc0->gcslparacox(p0,nh_vec(-1,kz),10);
        }
    }

    for(int l=1; l<=maxlev; ++l)
    {
        // coarse correction interpolated into the patch interior, then the cells around it
        for(int n : lev[l])
        {
            sflow_amr_patch *c = P[n];
            slice &z = nh_vec(n,kz);
            const int nby = c->ny/2;

            for(int bi=0; bi<c->nx/2; ++bi)
            for(int bj=0; bj<nby; ++bj)
            {
                int kk = bi*nby+bj;
                int g = c->rgrid[kk];
                for(int a=0;a<2;++a)
                for(int d=0;d<2;++d)
                {
                    int ii=EXT+2*bi+a, jj=EXT+2*bj+d;
                    if(g<-1)
                    {
                        z(ii,jj) = 0.0;
                        continue;
                    }
                    sflow_amr_fill f;
                    f.kind=1; f.g=g; f.si=c->ric[kk]; f.sj=c->rjc[kk];
                    f.ox = a==0 ? -1 : 1; f.oy = d==0 ? -1 : 1;
                    z(ii,jj) = nh_eval(f,kz);
                }
            }
        }

        nh_qfill(l,kz);

        // patch-local corrections of the residual left by the coarse correction; repeated
        // passes let siblings see each other's corrections through the filled cells
        const int npass = (nlevg[l]>1) ? NHPASS : 1;
        for(int pass=0; pass<npass; ++pass)
        {
        if(pass>0)
        nh_qfill(l,kz);

        for(int n : lev[l])
        {
            sflow_amr_patch *c = P[n];
            lexer *pp = c->pp;
            matrix2D &M = c->b->M;
            slice &z = nh_vec(n,kz);
            slice &r = nh_vec(n,kr);
            sc_level &L = c->mg->fine();

            for(int ii=EXT; ii<EXT+c->nx; ++ii)
            for(int jj=EXT; jj<EXT+c->ny; ++jj)
            {
                long qd = L.idx(ii-EXT,jj-EXT,0);
                L.u[qd] = 0.0;
                L.f[qd] = 0.0;
                if(L.act[qd]==0)
                continue;
                int rw = c->row[lij(pp,ii,jj)];
                double az = M.p[rw]*z(ii,jj) + M.n[rw]*z(ii+1,jj) + M.s[rw]*z(ii-1,jj) + M.w[rw]*z(ii,jj+1) + M.e[rw]*z(ii,jj-1);
                L.f[qd] = r(ii,jj) - az;
            }

            c->mg->vcycle(0,1,1);

            for(int ii=EXT; ii<EXT+c->nx; ++ii)
            for(int jj=EXT; jj<EXT+c->ny; ++jj)
            {
                long qd = L.idx(ii-EXT,jj-EXT,0);
                if(L.act[qd])
                z(ii,jj) += L.u[qd];
                else
                z(ii,jj) = r(ii,jj);
            }
        }
        }

        if(l<maxlev)
        nh_qfill(l,kz);
    }

    // q = 0 rows
    for(int g=-1; g<(int)P.size(); ++g)
    {
        nhg G = nh_grid(g);
        lexer *q = G.q;
        slice &z = nh_vec(g,kz);
        slice &r = nh_vec(g,kr);
        int i0,i1,j0,j1;
        if(g<0) { i0=0; i1=NX0-1; j0=0; j1=NY0-1; }
        else    { i0=EXT; i1=EXT+P[g]->nx-1; j0=EXT; j1=EXT+P[g]->ny-1; }

        for(int ii=i0; ii<=i1; ++ii)
        for(int jj=j0; jj<=j1; ++jj)
        if((*G.act)[lij(q,ii,jj)]==0)
        z(ii,jj) = r(ii,jj);
    }
}

void sflow_amr::nh_solve(lexer *p, fdm2D *b, ghostcell *pgc, slice &UH, slice &VH, slice &WH, slice &WL, double alpha)
{
    const double t0 = MPI_Wtime();
    nh_prepare(pgc);

    auto each = [&](auto fn)
    {
        for(int g=-1; g<(int)P.size(); ++g)
        {
            nhg G = nh_grid(g);
            lexer *q = G.q;
            int i0,i1,j0,j1;
            if(g<0) { i0=0; i1=NX0-1; j0=0; j1=NY0-1; }
            else    { i0=EXT; i1=EXT+P[g]->nx-1; j0=EXT; j1=EXT+P[g]->ny-1; }

            for(int ii=i0; ii<=i1; ++ii)
            for(int jj=j0; jj<=j1; ++jj)
            {
                int c = lij(q,ii,jj);
                if((*G.act)[c]>=-1)
                fn(g,G,ii,jj,c);
            }
        }
    };

    // right-hand side into NS (temporarily), q = 0 rows start at 0
    each([&](int g, nhg &G, int ii, int jj, int c)
    {
        int a = (*G.act)[c];
        nh_vec(g,NS)(ii,jj) = (a>=0) ? G.b->rhsvec.V[(*G.row)[c]] : 0.0;
        if(a==0)
        nh_vec(g,-1)(ii,jj) = 0.0;
    });

    // BiCGStab, right preconditioned; initial guess: the previous pressure
    nh_apply(-1,NVV);
    each([&](int g, nhg &G, int ii, int jj, int c)
    {
        double r = ((*G.act)[c]>=0) ? nh_vec(g,NS)(ii,jj) - nh_vec(g,NVV)(ii,jj) : 0.0;
        nh_vec(g,NR)(ii,jj) = r;
        nh_vec(g,NRH)(ii,jj) = r;
        nh_vec(g,NPV)(ii,jj) = 0.0;
        nh_vec(g,NVV)(ii,jj) = 0.0;
    });

    double bn = sqrt(nh_dot(NS,NS));
    double rn = sqrt(nh_dot(NR,NR));
    double rho=1.0, alp=1.0, om=1.0;
    int it=0;

    if(bn>0.0)
    while(rn/bn>p->N44 && it<p->N46)
    {
        double rhon = nh_dot(NRH,NR);
        double beta = (rhon/rho)*(alp/om);
        rho = rhon;

        each([&](int g, nhg &G, int ii, int jj, int c)
        {
            nh_vec(g,NPV)(ii,jj) = nh_vec(g,NR)(ii,jj) + beta*(nh_vec(g,NPV)(ii,jj) - om*nh_vec(g,NVV)(ii,jj));
        });

        double ta = MPI_Wtime();
        nh_prec(NPV,NPH);
        double tb = MPI_Wtime();
        nh_apply(NPH,NVV);
        tm[6] += tb-ta;
        tm[7] += MPI_Wtime()-tb;
        alp = rho/nh_dot(NRH,NVV);

        each([&](int g, nhg &G, int ii, int jj, int c)
        {
            nh_vec(g,NS)(ii,jj) = nh_vec(g,NR)(ii,jj) - alp*nh_vec(g,NVV)(ii,jj);
        });

        ta = MPI_Wtime();
        nh_prec(NS,NSH);
        tb = MPI_Wtime();
        nh_apply(NSH,NT);
        tm[6] += tb-ta;
        tm[7] += MPI_Wtime()-tb;
        double tt = nh_dot(NT,NT);
        om = tt>0.0 ? nh_dot(NT,NS)/tt : 0.0;

        each([&](int g, nhg &G, int ii, int jj, int c)
        {
            if((*G.act)[c]<0)
            return;
            nh_vec(g,-1)(ii,jj) += alp*nh_vec(g,NPH)(ii,jj) + om*nh_vec(g,NSH)(ii,jj);
            nh_vec(g,NR)(ii,jj) = nh_vec(g,NS)(ii,jj) - om*nh_vec(g,NT)(ii,jj);
        });

        rn = sqrt(nh_dot(NR,NR));
        ++it;
    }

    nh_it_last = it;
    nh_it_total += it;
    ++nh_solves;
    p->solveriter = it;
    p->final_res = bn>0.0 ? rn/bn : 0.0;

    // q = 0 rows exactly, then the values around the patches
    each([&](int g, nhg &G, int ii, int jj, int c)
    {
        if((*G.act)[c]==0)
        nh_vec(g,-1)(ii,jj) = 0.0;
    });
    nh_sync(-1);

    // velocity correction and the rest of the stage on the patches
    comms_off guard(pgc);
    for(auto c : P)
    {
        sflow_momentum_RK3 *m = c->pmom;
        pgc->gcsl_start4(c->pp,c->b->press,c->pnh->gcval());
        c->pnh->correct(c->pp,c->b,*m->nhUH,*m->nhVH,*m->nhWH,*m->nhWL,m->nh_alpha);
        m->stage_finish(c->pp,c->b,pgc,*m->nhUH,*m->nhVH,*m->nhWH,*m->nhWL);
    }

    tm[5] += MPI_Wtime()-t0;
}
