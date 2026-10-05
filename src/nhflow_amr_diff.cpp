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

#include"nhflow_amr.h"
#include"nhflow_amr_fill.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"solver.h"
#include"nhflow_diffusion.h"
#include"reefamr_krylov.h"
#include<cmath>
#include<mpi.h>
#include<iostream>

//  Composite implicit diffusion of REEF3D::NHFLOW on level 0 and the patches (G 31 1, A 512 2).
//
//  Every grid assembles its own rows with its diffusion object (nhflow_idiff, called through
//  nhflow_momentum_func::phase_D): (1/dt - D) X = rhs over its computed cells, D with the
//  viscosity, the eddy viscosity and the breaking viscosity of the grid.  The unknowns are the
//  leaf cells (level-0 cells and patch box cells not covered by a finer patch), solved for the
//  increment d = X - X0 over the stage input X0 (UHDIFF as nhflow_idiff leaves it before the
//  solve): A d = rhs - A X0, A X0 with every grid's own neighbours.  A leaf row is the grid's own
//  row; the increment of its neighbours outside the leaf cells of the grid is
//    - covered cells: the 2x2 mean of the children per layer (with G 6 also of both fine layers),
//      as restrict_momentum
//    - the cells around a patch: a sibling patch (copy, any rank), or the coarser grid as the
//      stage fill interpolates UH: the velocity d/WL bicubic in x,y from wet cells, times the WL
//      of the patch (with G 6 split into the two fine layers with the central slope: no limiter,
//      the operator of the Krylov solve is linear)
//    - level 0 across a partition edge: the halo (gcparaxijk_single, as bicgstab_ijk)
//  Other neighbours (solids, beyond the domain) are 0, as in bicgstab_ijk.  Every cell of a grid
//  then gets X0 + d (the cells around the patches and the covered cells their interpolated or
//  restricted increment): where nothing diffuses a grid keeps its stage values.
//
//  BiCGStab (reefamr_bicgstab) with the Jacobi preconditioner of bicgstab_ijk, residual N 43
//  relative to rhs, at most N 46 iterations.  A patch split in pieces (partition edge, placed
//  patches) gives the same solution as the whole patch, up to that tolerance.  The coupling is not
//  conservative across a coarse-fine interface: the coarse cell next to a patch sees the
//  restricted fine cells (as the composite pressure).

using namespace nhflow_amr_detail;

namespace
{
enum { DR=REEFAMR_NR, DRH=REEFAMR_NRH, DPV=REEFAMR_NPV, DVV=REEFAMR_NVV, DS=REEFAMR_NS, DT=REEFAMR_NT,
       DPH=REEFAMR_NPH, DSH=REEFAMR_NSH, DX=REEFAMR_NVEC, DNV=REEFAMR_NVEC+1 };

typedef reefamr_comms_off comms_off;

// the solver of phase_D: records the system of the grid instead of solving it
class df_capture_solver : public solver
{
public:
    explicit df_capture_solver(nhflow_amr *a) : amr(a) {}

    void startV(lexer*, ghostcell*, double *f, vec &rhs, matrix_diag &M, int var) override
    {
        amr->df_capture(f,rhs,M,var);
    }

    void start(lexer*, fdm*, ghostcell*, field&, vec&, int) override { fail(); }
    void startf(lexer*, ghostcell*, field&, vec&, matrix_diag&, int) override { fail(); }
    void startF(lexer*, ghostcell*, double*, vec&, matrix_diag&, int) override { fail(); }
    void startM(lexer*, ghostcell*, double*, double*, double*, int) override { fail(); }

private:
    void fail()
    {
        std::cout<<"NHFLOW AMR: composite diffusion, unexpected solver call"<<std::endl;
        MPI_Abort(MPI_COMM_WORLD,-2871);
    }
    nhflow_amr *amr;
};

// the diffusion of phase_M after the composite solve: UHDIFF, VHDIFF, WHDIFF stay as they are
class df_keep : public nhflow_diffusion
{
public:
    void diff_u(lexer*, fdm_nhf*, ghostcell*, ioflow*, solver*, double*, double*, double*, double*, double*, slice&, double) override {}
    void diff_v(lexer*, fdm_nhf*, ghostcell*, ioflow*, solver*, double*, double*, double*, double*, double*, slice&, double) override {}
    void diff_w(lexer*, fdm_nhf*, ghostcell*, ioflow*, solver*, double*, double*, double*, double*, double*, slice&, double) override {}
    void diff_scalar(lexer*, fdm_nhf*, ghostcell*, solver*, double*, double, double) override {}
};

// G 6: a coarse layer into its two fine layers, linear in sigma with the central slope (one-sided
// at the bed and the surface); the mean of the halves is the coarse value
void vsplit_lin(lexer *qc, const double *c, int nc, double *f, int marge)
{
    const double *ZP = qc->ZP + marge;
    const double *DZ = qc->DZN + marge;
    for(int k=0; k<nc; ++k)
    {
        double s = 0.0;
        if(nc>1)
        {
            if(k==0)
            s = (c[1]-c[0])/(ZP[1]-ZP[0]);
            else if(k==nc-1)
            s = (c[k]-c[k-1])/(ZP[k]-ZP[k-1]);
            else
            s = (c[k+1]-c[k-1])/(ZP[k+1]-ZP[k-1]);
        }
        const double h = 0.25*DZ[k]*s;
        f[2*k] = c[k] - h;
        f[2*k+1] = c[k] + h;
    }
}

// the composite diffusion equation as the vector space of reefamr_bicgstab
struct df_space
{
    static const int B=DS, R=DR, RH=DRH, PV=DPV, VV=DVV, S=DS, T=DT, PH=DPH, SH=DSH;
    nhflow_amr *a;
    void apply(int x, int y) { a->df_apply(x,y); }
    void prec(int r, int z) { a->df_prec(r,z); }
    double dot(int x, int y) { return a->df_dot(x,y); }
    void op_start() { a->df_start(); }
    void op_p(double beta, double om) { a->df_p(beta,om); }
    void op_s(double alp) { a->df_s(alp); }
    void op_x(double alp, double om) { a->df_x(alp,om); }
};
}

double* nhflow_amr::dvec(int g, int k)
{
    if(k<0)
    k = DX;
    return (g<0) ? dkv0[k].data() : NP(g)->dkv[k].data();
}

void nhflow_amr::df_capture(double *f, vec &rhs, matrix_diag &M, int var)
{
    dgrid &D = dg[dcur+1];
    D.x = f;
    D.rhs = &rhs;
    D.M = &M;
    D.var = var;
}

// the implicit diffusion of stage s on all grids: per component the rows of every grid (phase_D),
// one composite solve, the result into UHDIFF/VHDIFF/WHDIFF of every grid with its ghost cells
void nhflow_amr::diff_solve(lexer *p, ghostcell *pgc, int s)
{
    if(dcap==nullptr)
    {
        dcap = new df_capture_solver(this);
        dkeep = new df_keep();
    }

    dg.assign(P.size()+1,dgrid());
    df_fs.clear();
    df_stage = s;

    // vectors and leaf lists (the layout and the solid flags do not change during the step)
    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        const size_t n = (size_t)q->imax*q->jmax*(q->kmax+2);
        vector<vector<double>> &V = (g<0) ? dkv0 : NP(g)->dkv;
        if(V.size()!=(size_t)DNV || V[0].size()!=n)
        V.assign(DNV,vector<double>(n,0.0));

        dgrid &D = dg[g+1];
        const int l = (g<0) ? 0 : P[g]->lev;
        int i0,i1,j0,j1;
        if(g<0) { i0=0; i1=NX0-1; j0=0; j1=NY0-1; }
        else    { i0=EXT; i1=EXT+P[g]->nx-1; j0=EXT; j1=EXT+P[g]->ny-1; }

        int r=0;
        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        {
            const bool box = (ii>=i0 && ii<=i1 && jj>=j0 && jj<=j1);
            bool cov = false;
            if(box && l<maxlev)
            {
                const int I = (g<0) ? ii+O0i : ii-EXT+P[g]->I0;
                const int J = (g<0) ? jj+O0j : jj-EXT+P[g]->J0;
                cov = covered(l+1,2*I,2*J);
            }
            for(int kk=0; kk<q->knoz; ++kk)
            {
                const int qq = cidx(q,ii,jj,kk);
                if(q->flag4[qq]<=0)
                continue;
                D.aq.push_back(qq);
                if(box && !cov)
                {
                    D.lq.push_back(qq);
                    D.lr.push_back(r);
                }
                ++r;
            }
        }
    }

    for(int m=0; m<3; ++m)
    {
        // the rows of every grid
        {
        comms_off guard(pgc);
        for(int n=0; n<(int)P.size(); ++n)
        {
            nhflow_amr_patch *c = NP(n);
            pscope ps(pgc,c->d,d0);
            nhflow_stage_obj S = c->S;
            S.psolv = dcap;
            dcur = n;
            c->pmom->phase_D(c->pp,c->d,pgc,S,s,m);
        }
        }
        {
        nhflow_stage_obj S = *S0p;
        S.psolv = dcap;
        dcur = -1;
        mom0->phase_D(p,d0,pgc,S,s,m);
        }

        df_core(p);

        // the increment with the covered cells and the cells around the patches onto the arrays
        // of the grids, then their ghost cells
        df_sync(DX);
        for(int g=-1; g<(int)P.size(); ++g)
        {
            const dgrid &D = dg[g+1];
            const double *d = dvec(g,DX);
            for(int qq : D.aq)
            D.x[qq] += d[qq];
        }

        {
        comms_off guard(pgc);
        for(int n=0; n<(int)P.size(); ++n)
        {
            nhflow_amr_patch *c = NP(n);
            pscope ps(pgc,c->d,d0);
            pgc->start4V(c->pp,dg[n+1].x,14+m);
        }
        }
        pgc->start4V(p,dg[0].x,14+m);
    }
}

// one component: right-hand side, initial guess, BiCGStab
void nhflow_amr::df_core(lexer *p)
{
    // the stage input of every grid with its own cells around the patches and covered cells,
    // level 0 with its halo; outside the computed cells 0, as in bicgstab_ijk
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const dgrid &D = dg[g+1];
        double *x = dvec(g,DX);
        for(int qq : D.aq)
        x[qq] = D.x[qq];
    }
    if(p0->mpi_size>1)
    pgc0->gcparaxijk_single(p0,dvec(-1,DX),dg[0].var);
    pgc0->gc_periodic_ijk(p0,dvec(-1,DX));

    // the increment d = X - X0 solves A d = b - A X0 (A X0 with the grid's own neighbours): the
    // cells around the patches and the covered cells keep their stage values plus the
    // interpolated or restricted increment, so without diffusion nothing changes
    double sb=0.0, sr=0.0;
    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        const dgrid &D = dg[g+1];
        const matrix_diag &M = *D.M;
        const double *R = D.rhs->V.data();
        const double *x = dvec(g,DX);
        double *b = dvec(g,DS);
        const int sI = q->jmax*q->kmax;
        const int sJ = q->kmax;
        for(size_t n=0; n<D.lq.size(); ++n)
        {
            const int qq = D.lq[n], r = D.lr[n];
            const double ax = M.p[r]*x[qq] + M.n[r]*x[qq+sI] + M.s[r]*x[qq-sI] + M.w[r]*x[qq+sJ] + M.e[r]*x[qq-sJ]
                            + M.t[r]*x[qq+1] + M.b[r]*x[qq-1];
            b[qq] = R[r] - ax;
            sb += R[r]*R[r];
            sr += b[qq]*b[qq];
        }
    }
    const double bfull = sqrt(pgc0->globalsum(sb));
    const double r0 = sqrt(pgc0->globalsum(sr));

    for(int g=-1; g<(int)P.size(); ++g)
    {
        double *d = dvec(g,DX);
        double *vv = dvec(g,DVV);
        for(int qq : dg[g+1].lq)
        d[qq] = vv[qq] = 0.0;
    }

    // relative residual N 43 of the whole system (bicgstab_ijk divides the residual norm by the
    // cell count, which over all grids would loosen the patch solves)
    int it = 0;
    if(r0>0.0)
    {
        double bn, rn;
        df_space sp{this};
        it = reefamr_bicgstab(sp,p->N43*bfull/r0,p->N46,bn,rn);
    }

    df_it_last = it;
    df_it_total += it;
    ++df_solves;
    p->solveriter = it;
}

// covered cells, level-0 halo, cells around the patches of vector k
void nhflow_amr::df_sync(int k)
{
    // covered cells: 2x2 mean per layer (with G 6 of both fine layers), finest first
    for(int l=maxlev; l>=1; --l)
    {
        const int Kc = klev(l-1);
        block_up(l,Kc,7800+l,
                 [&](reefamr_patch *q, int id, int kb, double *v)
                 {
                     lexer *pp = q->pp;
                     const double *f = dvec(id,k);
                     const int nby = q->ny/2;
                     const int i0 = EXT+2*(kb/nby), j0 = EXT+2*(kb%nby);
                     const int fz = pp->knoz/Kc;
                     for(int kk=0; kk<Kc; ++kk)
                     {
                         double r = 0.0;
                         for(int kf=fz*kk; kf<fz*kk+fz; ++kf)
                         r += f[cidx(pp,i0,j0,kf)]+f[cidx(pp,i0+1,j0,kf)]+f[cidx(pp,i0,j0+1,kf)]+f[cidx(pp,i0+1,j0+1,kf)];
                         v[kk] = 0.25*r/double(fz);
                     }
                 },
                 [&](const reefamr_block &B, int key, const double *v)
                 {
                     lexer *q = glex(B.g);
                     double *dst = dvec(B.g,k);
                     for(int kk=0; kk<Kc; ++kk)
                     dst[cidx(q,B.ic,B.jc,kk)] = v[kk];
                 });
    }

    if(p0->mpi_size>1)
    {
        double *x0 = dvec(-1,k);
        pgc0->gcparaxijk_single(p0,x0,dg[0].var);
    }
    pgc0->gc_periodic_ijk(p0,dvec(-1,k));

    for(int l=1; l<=maxlev; ++l)
    {
        const int K = klev(l);
        const int Kc = klev(l-1);
        vector<double> uc(Kc);

        fill_run(l,K,7820+l,
                 [&](const reefamr_fill &f, double *v)
                 {
                     lexer *q = glex(f.g);
                     const double *src = dvec(f.g,k);
                     if(f.kind==0)
                     {
                         for(int kk=0; kk<K; ++kk)
                         v[kk] = src[cidx(q,f.si,f.sj,kk)];
                         return;
                     }
                     if(f.kind!=1)
                     {
                         for(int kk=0; kk<K; ++kk)
                         v[kk] = 0.0;
                         return;
                     }

                     // the stencil of the stage fill (wet cells), the weights divided by the WL of
                     // their column: velocities of the coarser grid
                     auto it = df_fs.find(&f);
                     if(it==df_fs.end())
                     {
                         double w1[25];
                         pweights(q,f.si,f.sj,f.ox,f.oy,w1,1);
                         stg C = stage_in(f.g,df_stage);
                         dstencil st;
                         st.nw = 0;
                         const int n0 = cidx(q,f.si,f.sj,0);
                         const int sI = q->jmax*q->kmax, sJ = q->kmax;
                         for(int di=-2; di<=2; ++di)
                         for(int dj=-2; dj<=2; ++dj)
                         {
                             const double c = w1[(di+2)*5+(dj+2)];
                             if(c==0.0)
                             continue;
                             const double wl = (*C.WL)(f.si+di,f.sj+dj);
                             const double wlvl = fabs(wl)>p0->A544 ? wl : 1.0e20;
                             st.off[st.nw] = n0 + di*sI + dj*sJ;
                             st.w[st.nw] = c/wlvl;
                             ++st.nw;
                         }
                         it = df_fs.emplace(&f,st).first;
                     }
                     const dstencil &st = it->second;
                     for(int kk=0; kk<Kc; ++kk)
                     {
                         double r = 0.0;
                         for(int m=0; m<st.nw; ++m)
                         r += st.w[m]*src[st.off[m]+kk];
                         uc[kk] = r;
                     }
                     if(K==Kc)
                     for(int kk=0; kk<K; ++kk)
                     v[kk] = uc[kk];
                     else
                     vsplit_lin(q,uc.data(),Kc,v,marge);
                 },
                 [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
                 {
                     lexer *pp = q->pp;
                     double *dst = dvec(id,k);
                     double wl = 1.0;
                     if(f.kind==1)
                     {
                         stg S = stage_in(id,df_stage);
                         wl = (*S.WL)(f.di,f.dj);
                         // G 30: a dry cell around the patch has no momentum (dry_cell)
                         if(shore && pp->wet[lij(pp,f.di,f.dj)]==0)
                         wl = 0.0;
                     }
                     for(int kk=0; kk<K; ++kk)
                     dst[cidx(pp,f.di,f.dj,kk)] = wl*w[kk];
                 });
    }
}

// y = A x on the leaf cells
void nhflow_amr::df_apply(int kx, int ky)
{
    df_sync(kx);

    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        const dgrid &D = dg[g+1];
        const matrix_diag &M = *D.M;
        const double *x = dvec(g,kx);
        double *y = dvec(g,ky);
        const int sI = q->jmax*q->kmax;
        const int sJ = q->kmax;
        const double *Mp=M.p.data(), *Mn=M.n.data(), *Ms=M.s.data(), *Mw=M.w.data(), *Me=M.e.data(), *Mt=M.t.data(), *Mb=M.b.data();

        for(size_t n=0; n<D.lq.size(); ++n)
        {
            const int qq = D.lq[n], r = D.lr[n];
            y[qq] = Mp[r]*x[qq] + Mn[r]*x[qq+sI] + Ms[r]*x[qq-sI] + Mw[r]*x[qq+sJ] + Me[r]*x[qq-sJ] + Mt[r]*x[qq+1] + Mb[r]*x[qq-1];
        }
    }
}

// z = r / diagonal (the preconditioner of bicgstab_ijk)
void nhflow_amr::df_prec(int kr, int kz)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const dgrid &D = dg[g+1];
        const double *Mp = D.M->p.data();
        const double *r = dvec(g,kr);
        double *z = dvec(g,kz);
        for(size_t n=0; n<D.lq.size(); ++n)
        {
            const int qq = D.lq[n];
            z[qq] = r[qq]/(Mp[D.lr[n]]+1.0e-20);
        }
    }
}

double nhflow_amr::df_dot(int ka, int kb)
{
    double s=0.0;
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *a = dvec(g,ka), *b = dvec(g,kb);
        for(int qq : dg[g+1].lq)
        s += a[qq]*b[qq];
    }
    return pgc0->globalsum(s);
}

void nhflow_amr::df_start()
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vb = dvec(g,DS);
        double *vr = dvec(g,DR), *vrh = dvec(g,DRH), *vp = dvec(g,DPV), *vv = dvec(g,DVV);
        for(int qq : dg[g+1].lq)
        {
            double r = vb[qq] - vv[qq];
            vr[qq] = r; vrh[qq] = r; vp[qq] = 0.0; vv[qq] = 0.0;
        }
    }
}

void nhflow_amr::df_p(double beta, double om)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vr = dvec(g,DR), *vv = dvec(g,DVV);
        double *vp = dvec(g,DPV);
        for(int qq : dg[g+1].lq)
        vp[qq] = vr[qq] + beta*(vp[qq] - om*vv[qq]);
    }
}

void nhflow_amr::df_s(double alp)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vr = dvec(g,DR), *vv = dvec(g,DVV);
        double *vs = dvec(g,DS);
        for(int qq : dg[g+1].lq)
        vs[qq] = vr[qq] - alp*vv[qq];
    }
}

void nhflow_amr::df_x(double alp, double om)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vph = dvec(g,DPH), *vsh = dvec(g,DSH), *vs = dvec(g,DS), *vt = dvec(g,DT);
        double *x = dvec(g,DX), *vr = dvec(g,DR);
        for(int qq : dg[g+1].lq)
        {
            x[qq] += alp*vph[qq] + om*vsh[qq];
            vr[qq] = vs[qq] - om*vt[qq];
        }
    }
}

// the objects of the composite diffusion (destructor)
void nhflow_amr::df_free()
{
    delete static_cast<df_capture_solver*>(dcap);
    delete static_cast<df_keep*>(dkeep);
    dcap = nullptr;
    dkeep = nullptr;
}
