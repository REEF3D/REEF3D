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
#include"nhflow_breaking.h"
#include"reefamr_krylov.h"
#include<cmath>
#include<mpi.h>
#include<iostream>
#include<iomanip>

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
//  patches) gives the same solution as the whole patch, up to that tolerance.  The coupling is
//  conservative across a coarse-fine interface: the coarse face of a leaf cell next to a patch takes
//  the sum of the viscous fluxes of the fine faces there (df_fix, below; single-block grids).

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
    void op_restart() { a->df_restart(); }
};

// safety net of the BiCGStab (reefamr_krylov.h, as the composite pressure): restart after a stagnation,
// stop after a second one below STAGTOL times the tolerance
const double STAGTOL = 100.0;
}

double* nhflow_amr::dvec(int g, int k)
{
    if(k<0)
    k = DX;
    // a grid without vectors (the patches of levels 1.. in the turbulence of level 0, turb_composite0): no cells
    vector<vector<double>> &V = (g<0) ? dkv0 : NP(g)->dkv;
    return ((size_t)k<V.size()) ? V[k].data() : nullptr;
}

void nhflow_amr::df_capture(double *f, vec &rhs, matrix_diag &M, int var)
{
    dgrid &D = dg[dcur+1];
    D.x = f;
    D.rhs = &rhs;
    D.M = &M;
    D.var = var;
}

// the capture solver of the systems of the grids (created with the first use)
solver* nhflow_amr::df_capsolver()
{
    if(dcap==nullptr)
    {
        dcap = new df_capture_solver(this);
        dkeep = new df_keep();
    }
    return dcap;
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
        else    { i0=EXT; i1=EXT+P[g]->nx-1; j0=EY; j1=EY+P[g]->ny-1; }

        int r=0;
        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        {
            // multi-block: the ring of the native grid stands for cells of the other blocks
            const bool box = (ii>=i0 && ii<=i1 && jj>=j0 && jj<=j1) && !(g<0 && ring0(ii,jj));
            bool cov = false;
            if(box && l<maxlev)
            {
                const int I = (g<0) ? ii+O0i : ii-EXT+P[g]->I0;
                const int J = (g<0) ? jj+O0j : jj-EY+P[g]->J0;
                cov = covered(l+1,2*I,ychild(J));
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

        D.row.assign((size_t)q->imax*q->jmax*q->kmax,-1);
        for(size_t n=0; n<D.aq.size(); ++n)
        D.row[D.aq[n]] = (int)n;
    }

    dfid.clear();
    for(int n=0; n<(int)P.size(); ++n)
    dfid[P[n]] = n;

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
        if(D.M==nullptr)       // a grid outside the system (the patches of the turbulence of level 0)
        continue;
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
    // the coarse faces next to the finer patches take the sum of the fine faces (df_fix), for X0 with the grids'
    // own cells around the patches and covered cells as for the rows above
    df_fix(DX,DS,-1.0);
    sr = 0.0;
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *b = dvec(g,DS);
        for(int qq : dg[g+1].lq)
        sr += b[qq]*b[qq];
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
        // the system of level 0 alone (the turbulence of the blocks, df_l0): the residual per cell as bicgstab_ijk
        // (the turbulence of one grid), so that the blocks solve as tightly as one grid
        double tol = p->N43*bfull/r0;
        if(df_l0)
        {
            double nl = 0.0;
            for(const dgrid &D : dg)
            nl += double(D.lq.size());
            tol = p->N43*pgc0->globalsum(nl)/r0;
            // at least one iteration, as bicgstab_ijk
            tol = MIN(tol,1.0-1.0e-9);
        }
        it = reefamr_bicgstab(sp,tol,p->N46,bn,rn,nullptr,nullptr,&df_restarts,STAGTOL*tol,&df_stalls);

        if(it>=p->N46 && rn>tol*bn)
        {
            ++df_capped;
            if(p->mpirank==0)
            cout<<"NHFLOW AMR WARNING - diffusion iteration limit N 46 = "<<p->N46<<" reached, res "<<setprecision(3)<<(bn>0.0 ? rn/bn : 0.0)
                <<" of the initial residual (tolerance "<<tol<<")"<<endl;
        }
    }

    df_it_last = it;
    df_it_total += it;
    ++df_solves;
    p->solveriter = it;
}

// covered cells, level-0 halo, cells around the patches of vector k
void nhflow_amr::df_sync(int k)
{
    // covered cells: 2x2 mean per layer (with G 6 of both fine layers), finest first (not in a system of level 0
    // alone: df_l0)
    for(int l=(df_l0 ? 0 : maxlev); l>=1; --l)
    {
        const int Kc = klev(l-1);
        block_up(l,Kc,7800+l,
                 [&](reefamr_patch *q, int id, int kb, double *v)
                 {
                     lexer *pp = q->pp;
                     const double *f = dvec(id,k);
                     const int nby = q->ny>>sy;
                     const int i0 = EXT+2*(kb/nby), j0 = EY+ychild(kb%nby);
                     const int fz = pp->knoz/Kc;
                     double w[4];
                     cw4(pp,i0,j0,w);    // curvilinear grid: per area of the children
                     for(int kk=0; kk<Kc; ++kk)
                     {
                         double r = 0.0;
                         if(pp->geo_curv==1)
                         {
                         for(int kf=fz*kk; kf<fz*kk+fz; ++kf)
                         r += (sy==0 ? w[0]*f[cidx(pp,i0,j0,kf)]+w[1]*f[cidx(pp,i0+1,j0,kf)] : w[0]*f[cidx(pp,i0,j0,kf)]+w[1]*f[cidx(pp,i0+1,j0,kf)]+w[2]*f[cidx(pp,i0,j0+1,kf)]+w[3]*f[cidx(pp,i0+1,j0+1,kf)]);
                         v[kk] = r/double(fz);
                         continue;
                         }
                         if(sy==0)      // 2D grid: the 2x1 children
                         {
                         for(int kf=fz*kk; kf<fz*kk+fz; ++kf)
                         r += f[cidx(pp,i0,j0,kf)]+f[cidx(pp,i0+1,j0,kf)];
                         v[kk] = 0.5*r/double(fz);
                         continue;
                         }
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

    // multi-block: the cells around the blocks and the ring from the cells they stand for
    if(nblk>0)
    {
        const int K0 = p0->knoz;
        fill_run0(K0,7819,
                  [&](const reefamr_fill &f, double *v)
                  {
                      if(f.kind!=0)
                      {
                          for(int kk=0; kk<K0; ++kk)
                          v[kk] = 0.0;
                          return;
                      }
                      const double *src = dvec(f.g,k);
                      lexer *q = glex(f.g);
                      for(int kk=0; kk<K0; ++kk)
                      v[kk] = src[cidx(q,f.si,f.sj,kk)];
                  },
                  [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
                  {
                      double *dst = dvec(id,k);
                      for(int kk=0; kk<K0; ++kk)
                      dst[cidx(q->pp,f.di,f.dj,kk)] = w[kk];
                  },
                  [&](const reefamr_fill &f, const double *w)
                  {
                      double *dst = dvec(-1,k);
                      for(int kk=0; kk<K0; ++kk)
                      dst[cidx(p0,f.di,f.dj,kk)] = w[kk];
                  });
    }

    if(p0->mpi_size>1)
    {
        double *x0 = dvec(-1,k);
        pgc0->gcparaxijk_single(p0,x0,dg[0].var);
    }
    pgc0->gc_periodic_ijk(p0,dvec(-1,k));

    for(int l=1; l<=(df_l0 ? 0 : maxlev); ++l)
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
        if(D.M==nullptr)
        continue;
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
    // the coarse faces next to the finer patches: the sum of the fine faces
    df_fix(kx,ky,1.0);
}

// z = r / diagonal (the preconditioner of bicgstab_ijk)
void nhflow_amr::df_prec(int kr, int kz)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const dgrid &D = dg[g+1];
        if(D.M==nullptr)
        continue;
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

// restart of the BiCGStab recurrence from the current residual (reefamr_bicgstab)
void nhflow_amr::df_restart()
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vr = dvec(g,DR);
        double *vrh = dvec(g,DRH), *vp = dvec(g,DPV), *vv = dvec(g,DVV);
        for(int qq : dg[g+1].lq)
        {
            vrh[qq] = vr[qq];
            vp[qq] = 0.0;
            vv[qq] = 0.0;
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

// --------------------------------------------------------------------- breaking on split patches
// G 31 1 with A 550 1: the parts of nhflow_breaking on every patch, coarse to fine.  A cell around
// a patch filled from a sibling patch (fill kind 0) takes the sibling's detection (after part 0)
// and its flags (after part 1) of this stage, and its breaking viscosity at the end, so that the
// cells of a patch next to the cut see what the whole patch would; the cells filled from the
// coarser grid keep the breaking of their source (lexer::amrvb, the stage fill).
void nhflow_amr::breaking_patches(ghostcell *pgc)
{
    for(int part=0; part<3; ++part)
    {
        // multi-block: level 0 first (block 0, the blocks 1.., their exchange)
        if(brk0)
        breaking_level0(pgc,part);

        {
        comms_off guard(pgc);
        for(int n=0; n<(int)P.size(); ++n)
        {
            nhflow_amr_patch *c = NP(n);
            if(brk0 && c->lev==0)
            continue;
            nhflow_breaking *b = dynamic_cast<nhflow_breaking*>(c->pmom);
            if(b==nullptr)
            continue;
            pscope ps(pgc,c->d,d0);
            b->breaking_part(c->pp,c->d,pgc,part);
        }
        }

        for(int l=1; l<=maxlev; ++l)
        brk_exchange(l,part);
    }
}

// --------------------------------------------------------------------- breaking on a multi-block grid
// A 550 1 with G 50 1: part 'part' (0 detection, 1 flags, 2 breaking viscosity; -1 all three in turn) of
// nhflow_breaking on block 0 and the blocks 1.., then the ring of block 0 and the cells around the blocks
// take the values of the cells they stand for (brk_exchange0).  The detection looks at most two cells
// beyond a cell, the flags spread two cells and the viscosity is read one cell beyond by the diffusion:
// the EXT = 3 cells around the blocks and the ring hold what one grid would have there.  (The filter of
// A 552 runs on every block alone.)
void nhflow_amr::breaking_level0(ghostcell *pgc, int part)
{
    // A 550 3: the eddy viscosity of Kennedy et al. (2000) in its parts
    if(p0->A550==3)
    {
        if(part<=0)
        kennedy_level0(pgc);
        return;
    }

    const int p0s = (part<0) ? 0 : part;
    const int p1s = (part<0) ? 2 : part;

    for(int q=p0s; q<=p1s; ++q)
    {
        if(nhflow_breaking *b = dynamic_cast<nhflow_breaking*>(mom0))
        b->breaking_part(p0,d0,pgc,q);

        {
        comms_off guard(pgc);
        for(int id : lev[0])
        {
            nhflow_amr_patch *c = NP(id);
            nhflow_breaking *b = dynamic_cast<nhflow_breaking*>(c->pmom);
            if(b==nullptr)
            continue;
            pscope ps(pgc,c->d,d0);
            b->breaking_part(c->pp,c->d,pgc,q);
        }
        }

        brk_exchange0(pgc,q);
    }
}

// A 550 3 with G 50 1: breaking_kennedy on block 0 and the blocks 1.. as one grid.  The ages of the last step, then
// every pass of the band of the viscosity read the cells around the blocks and the ring; these take the values of the
// cells they stand for after part 0 (ages), before every pass (viscosity, wavenumber) and at the end (with the
// flags), so that the band spreads across the interfaces as on one grid
void nhflow_amr::kennedy_level0(ghostcell *pgc)
{
    auto run = [&](int part, int qn)
    {
        int band = 0, cnt = 0;
        if(nhflow_breaking *b = dynamic_cast<nhflow_breaking*>(mom0))
        {
            b->kennedy_run(p0,d0,pgc,part,qn);
            band = MAX(band,b->kn_band);
            cnt += b->kn_count;
        }
        {
        comms_off guard(pgc);
        for(int id : lev[0])
        {
            nhflow_amr_patch *c = NP(id);
            nhflow_breaking *b = dynamic_cast<nhflow_breaking*>(c->pmom);
            if(b==nullptr)
            continue;
            pscope ps(pgc,c->d,d0);
            b->kennedy_run(c->pp,c->d,pgc,part,qn);
            band = MAX(band,b->kn_band);
            cnt += b->kn_count;
        }
        }
        return std::make_pair(band,cnt);
    };

    run(0,0);
    ken_exchange0(pgc,0);

    const int nband = pgc->globalimax(run(1,0).first);

    for(int qn=0; qn<nband; ++qn)
    {
        run(2,qn);
        ken_exchange0(pgc,1);
        run(3,qn);
    }

    run(4,0);
    ken_exchange0(pgc,2);

    // the breaking columns of the cells of the grids (not the ring, not the cells around the blocks)
    int cnt = 0;
    for(int ii=0; ii<p0->knox; ++ii)
    for(int jj=0; jj<p0->knoy; ++jj)
    if(!ring0(ii,jj) && d0->breaking(ii,jj)>0)
    ++cnt;
    for(int id : lev[0])
    {
        const reefamr_patch *q = P[id];
        fdm_nhf *d = NP(id)->d;
        for(int ii=EXT; ii<EXT+q->nx; ++ii)
        for(int jj=EY; jj<EY+q->ny; ++jj)
        if(d->breaking(ii,jj)>0)
        ++cnt;
    }
    const int count = pgc->globalisum(cnt);
    if(p0->mpirank==0 && (p0->count%p0->P12==0))
    cout<<"breaking: "<<count<<endl;
}

// A 550 3: the ring of block 0 and the cells around the blocks 1.. from the cells they stand for: 0 the ages, the
// wavenumbers and the thresholds of the last step, 1 the viscosity and its wavenumber, 2 these and the flags
void nhflow_amr::ken_exchange0(ghostcell *pgc, int what)
{
    const int nv = (what==1) ? 2 : 3;
    auto brk = [&](int g) { return dynamic_cast<nhflow_breaking*>(gmom(g)); };

    auto put = [&](int g, const reefamr_fill &f, const double *w)
    {
        if(w[0]<-0.5)
        return;
        nhflow_breaking *b = brk(g);
        if(b==nullptr)
        return;
        fdm_nhf *d = gfd(g);
        if(what==0)
        {
            b->brk_age_n()(f.di,f.dj) = w[1];
            b->brk_k_n()(f.di,f.dj) = w[2];
            b->brk_r_n()(f.di,f.dj) = w[3];
            return;
        }
        d->vb(f.di,f.dj) = w[1];
        d->vbk(f.di,f.dj) = w[2];
        if(what==2)
        d->breaking(f.di,f.dj) = (int)w[3];
    };

    // the first value marks a cell with a source (1) or without (-1)
    fill_run0(nv+1,7880+what,
              [&](const reefamr_fill &f, double *v)
              {
                  nhflow_breaking *b = (f.kind==0) ? brk(f.g) : nullptr;
                  for(int m=0; m<=nv; ++m)
                  v[m] = 0.0;
                  if(b==nullptr)
                  {
                      v[0] = -1.0;
                      return;
                  }
                  v[0] = 1.0;
                  fdm_nhf *d = gfd(f.g);
                  if(what==0)
                  {
                      v[1] = b->brk_age_n()(f.si,f.sj);
                      v[2] = b->brk_k_n()(f.si,f.sj);
                      v[3] = b->brk_r_n()(f.si,f.sj);
                      return;
                  }
                  v[1] = d->vb(f.si,f.sj);
                  v[2] = d->vbk(f.si,f.sj);
                  if(what==2)
                  v[3] = d->breaking(f.si,f.sj);
              },
              [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w) { put(id,f,w); },
              [&](const reefamr_fill &f, const double *w) { put(-1,f,w); });

    // the ring cells in the partition halos of block 0
    if(p0->mpi_size>1)
    if(nhflow_breaking *b = brk(-1))
    {
        if(what==0)
        {
            pgc->gcsl_start4(p0,b->brk_age_n(),1);
            pgc->gcsl_start4(p0,b->brk_k_n(),1);
            pgc->gcsl_start4(p0,b->brk_r_n(),1);
            return;
        }
        pgc->gcsl_start4(p0,d0->vb,1);
        pgc->gcsl_start4(p0,d0->vbk,1);
        if(what==2)
        pgc->gcsl_start4int(p0,d0->breaking,50);
    }
}

// values of the ring of block 0 and of the cells around the blocks 1..: 0 detection (bx, by, bd), 1 flags,
// 2 breaking viscosity; then the partition halos of block 0
void nhflow_amr::brk_exchange0(ghostcell *pgc, int what)
{
    // A 552: the water level that the filter of part 2 changed, with the viscosity
    const bool wl = (what==2 && p0->A552>0);
    const int nv = (what==0) ? 3 : (wl ? 2 : 1);
    auto brk = [&](int g) { return dynamic_cast<nhflow_breaking*>(gmom(g)); };

    auto put = [&](int g, const reefamr_fill &f, const double *w)
    {
        if(w[0]<0.0)
        return;
        nhflow_breaking *b = brk(g);
        if(b==nullptr)
        return;
        if(what==0)
        {
            b->brk_bx()(f.di,f.dj) = (int)w[0];
            b->brk_by()(f.di,f.dj) = (int)w[1];
            b->brk_bd()(f.di,f.dj) = (int)w[2];
        }
        else if(what==1)
        b->brk_flag()(f.di,f.dj) = (int)w[0];
        else
        {
        gfd(g)->vb(f.di,f.dj) = w[0];
        if(wl && b->brk_wl()!=nullptr)
        (*b->brk_wl())(f.di,f.dj) = w[1];
        }
    };

    fill_run0(nv,7870+what,
              [&](const reefamr_fill &f, double *v)
              {
                  nhflow_breaking *b = (f.kind==0) ? brk(f.g) : nullptr;
                  if(b==nullptr)
                  {
                      for(int m=0; m<nv; ++m)
                      v[m] = -1.0;
                      return;
                  }
                  if(what==0)
                  {
                      v[0] = b->brk_bx()(f.si,f.sj);
                      v[1] = b->brk_by()(f.si,f.sj);
                      v[2] = b->brk_bd()(f.si,f.sj);
                  }
                  else if(what==1)
                  v[0] = b->brk_flag()(f.si,f.sj);
                  else
                  {
                  v[0] = gfd(f.g)->vb(f.si,f.sj);
                  if(wl)
                  v[1] = (b->brk_wl()!=nullptr) ? (*b->brk_wl())(f.si,f.sj) : 0.0;
                  }
              },
              [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w) { put(id,f,w); },
              [&](const reefamr_fill &f, const double *w) { put(-1,f,w); });

    // the ring cells in the partition halos of block 0
    if(p0->mpi_size>1)
    if(nhflow_breaking *b = brk(-1))
    {
        if(what==0)
        {
            pgc->gcsl_start4int(p0,b->brk_bx(),50);
            pgc->gcsl_start4int(p0,b->brk_by(),50);
            pgc->gcsl_start4int(p0,b->brk_bd(),50);
        }
        else if(what==1)
        pgc->gcsl_start4int(p0,b->brk_flag(),50);
        else
        {
        pgc->gcsl_start4(p0,d0->vb,1);
        if(wl && b->brk_wl()!=nullptr)
        pgc->gcsl_start4(p0,*b->brk_wl(),1);
        }
    }
}

// values of the sibling-filled cells around the level-l patches: 0 detection (bx, by, bd), 1
// flags, 2 breaking viscosity
void nhflow_amr::brk_exchange(int l, int what)
{
    const int nv = (what==0) ? 3 : 1;

    fill_run(l,nv,7840+10*what+l,
             [&](const reefamr_fill &f, double *v)
             {
                 if(f.kind!=0 || f.g<0)
                 {
                     for(int m=0; m<nv; ++m)
                     v[m] = -1.0;
                     return;
                 }
                 nhflow_amr_patch *c = NP(f.g);
                 nhflow_breaking *b = dynamic_cast<nhflow_breaking*>(c->pmom);
                 if(what==0)
                 {
                     v[0] = b->brk_bx()(f.si,f.sj);
                     v[1] = b->brk_by()(f.si,f.sj);
                     v[2] = b->brk_bd()(f.si,f.sj);
                 }
                 else if(what==1)
                 v[0] = b->brk_flag()(f.si,f.sj);
                 else
                 v[0] = c->d->vb(f.si,f.sj);
             },
             [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
             {
                 if(w[0]<0.0)
                 return;
                 nhflow_amr_patch *c = NP(q);
                 nhflow_breaking *b = dynamic_cast<nhflow_breaking*>(c->pmom);
                 if(what==0)
                 {
                     b->brk_bx()(f.di,f.dj) = (int)w[0];
                     b->brk_by()(f.di,f.dj) = (int)w[1];
                     b->brk_bd()(f.di,f.dj) = (int)w[2];
                 }
                 else if(what==1)
                 b->brk_flag()(f.di,f.dj) = (int)w[0];
                 else
                 c->d->vb(f.di,f.dj) = w[0];
             });
}

// ---- the viscous flux across the coarse-fine interfaces
//
// A leaf cell c of grid g next to a level-(l+1) patch (a coarse face of the flux registry, match / rmatch) has the
// row term T_c = a_c (x_cov - x_c) for that face, x_cov the covered cell (the mean of the children), a_c the
// coefficient of the face (-M.n, -M.s, -M.w or -M.e of its row).  The fine cells f of the patch at that face have
// T_f = a_f (x_out - x_f), x_out their cell around the patch (interpolated from the coarser grid).  The transfer
// over the interface is V_c T_c from the coarse side and -sum V_f T_f from the fine side (V the cell volume, the
// rows are per unit volume); they differ, and the diffusion neither conserves momentum there nor sees the fine
// gradient.  df_fix replaces T_c by -sum V_f T_f / V_c: the coarse row gets the fine flux, as the convective
// fluxes of the continuity and momentum equations (flux_hook).  The fine rows are unchanged.  Single-block grids
// (Cartesian and curvilinear); on a multi-block grid not.  Applied in the
// operator (df_apply, on the increment) and to the right-hand side (df_core, on X0 with the grids' own cells around
// the patches and covered cells), so the composite system is that of the matched fluxes for X = X0 + d.
//
// df_fine: sum over the fine cells of coarse layer k at the coarse face (side, first fine index r) of V_f T_f
double nhflow_amr::df_fine(int id, int side, int r, int K, int k, int kx)
{
    nhflow_amr_patch *c = NP(id);
    lexer *pp = c->pp;
    const dgrid &D = dg[id+1];
    const matrix_diag &M = *D.M;
    const double *x = dvec(id,kx);
    const int fz = pp->knoz/K;
    const int sI = pp->jmax*pp->kmax, sJ = pp->kmax;

    double sum = 0.0;
    const int nr = (fr1(side)==0) ? 1 : 2;
    for(int a=0; a<nr; ++a)
    {
        const int rr = r + a;
        int ii, jj;
        if(side==0)      { ii = EXT;           jj = EY+rr; }
        else if(side==1) { ii = EXT+c->nx-1;   jj = EY+rr; }
        else if(side==2) { ii = EXT+rr;        jj = EY; }
        else             { ii = EXT+rr;        jj = EY+c->ny-1; }

        const int n2 = lij(pp,ii,jj);
        if(pp->wet[n2]!=1)
        continue;
        const double area = (pp->geo_curv==1) ? pp->AREA2D[n2] : pp->DXN[ii+marge]*pp->DYN[jj+marge];

        for(int kf=fz*k; kf<fz*k+fz; ++kf)
        {
            const int qq = cidx(pp,ii,jj,kf);
            const int rw = D.row[qq];
            if(rw<0)
            continue;
            // the coefficient and the cell across the patch box
            double a_f; int qo;
            if(side==0)      { a_f = -M.s[rw]; qo = qq-sI; }
            else if(side==1) { a_f = -M.n[rw]; qo = qq+sI; }
            else if(side==2) { a_f = -M.e[rw]; qo = qq-sJ; }
            else             { a_f = -M.w[rw]; qo = qq+sJ; }
            sum += area*pp->DZN[kf+marge]*a_f*(x[qo]-x[qq]);
        }
    }
    return sum;
}

// y += sg * (T_c + sum V_f T_f / V_c) for every coarse face next to a finer patch (x of vector kx with its cells
// around the patches and covered cells set)
void nhflow_amr::df_fix(int kx, int ky, double sg)
{
    // (multi-block grids: the coarse faces of a patch at a block face may lie in a turned block, where the
    // velocity components of both sides differ: not matched there, as before)
    if(maxlev<1 || nblk>0)
    return;

    // the fine sums of the coarse faces whose patch is on another rank
    dfr.assign(rmatch.size(),vector<double>());
    for(size_t t=0; t<rmatch.size(); ++t)
    if(!rmatch[t].empty())
    dfr[t].assign(rmatch[t].size()*(size_t)glex((int)t-1)->knoz,0.0);

    for(int l=0; l<maxlev; ++l)
    {
        const int K = klev(l);
        face_run_at(l+1,K,7860+l+1,
                    [&](reefamr_patch *q, int side, int r, double *sb)
                    {
                        const int id = dfid[q];
                        for(int k=0; k<K; ++k)
                        sb[k] = df_fine(id,side,r,K,k,kx);
                    },
                    [&](int tgt, int idx, const double *v)
                    {
                        for(int k=0; k<K; ++k)
                        dfr[tgt][(size_t)idx*K+k] = v[k];
                    });
    }

    for(int g=-1; g<(int)P.size(); ++g)
    {
        const int t = g+1;
        if(t>=(int)match.size())
        continue;
        lexer *q = glex(g);
        const dgrid &D = dg[t];
        const matrix_diag &M = *D.M;
        const double *x = dvec(g,kx);
        double *y = dvec(g,ky);
        const int K = q->knoz;
        const int sI = q->jmax*q->kmax, sJ = q->kmax;

        auto fix = [&](const reefamr_match &m, auto fine)
        {
            // the leaf cell outside the patch and the direction to the covered cell
            const int ii = m.fi + (m.side==1 ? 1 : 0);
            const int jj = m.fj + (m.side==3 ? 1 : 0);
            const int n2 = lij(q,ii,jj);
            if(q->wet[n2]!=1)
            return;
            const double area = (q->geo_curv==1) ? q->AREA2D[n2] : q->DXN[ii+marge]*q->DYN[jj+marge];
            for(int k=0; k<K; ++k)
            {
                const int qq = cidx(q,ii,jj,k);
                const int rw = D.row[qq];
                if(rw<0)
                continue;
                double a_c; int qn;
                if(m.side==0)      { a_c = -M.n[rw]; qn = qq+sI; }
                else if(m.side==1) { a_c = -M.s[rw]; qn = qq-sI; }
                else if(m.side==2) { a_c = -M.w[rw]; qn = qq+sJ; }
                else               { a_c = -M.e[rw]; qn = qq-sJ; }
                const double Tc = a_c*(x[qn]-x[qq]);
                const double Vc = area*q->DZN[k+marge];
                y[qq] += sg*(Tc + fine(k)/Vc);
            }
        };

        for(size_t n=0; n<match[t].size(); ++n)
        {
            const reefamr_match &m = match[t][n];
            fix(m,[&](int k) { return df_fine(m.child,m.side,m.r,K,k,kx); });
        }
        if(t<(int)rmatch.size())
        for(size_t n=0; n<rmatch[t].size(); ++n)
        {
            const reefamr_match &m = rmatch[t][n];
            fix(m,[&](int k) { return dfr[t][n*(size_t)K+k]; });
        }
    }
}
