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

#include"fnpf_amr.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"slice4.h"
#include"fnpf_fsfbc.h"
#include"fnpf_sigma.h"
#include"fnpf_fsf_update.h"
#include"fnpf_bed_update.h"
#include"fnpf_laplace_cds2.h"
#include"reefmg_core.h"
#include"reefamr_krylov.h"
#include"fnpf_amr_fill.h"
#include"fnpf_body.h"
#include<cmath>
#include<mpi.h>
#include<iomanip>
#include<algorithm>

//  Composite Laplace equation of REEF3D::FNPF on level 0 and the patches.
//
//  Every grid assembles its own rows with fnpf_laplace_cds2 (level 0 in fnpf_RK3, the patches
//  here), with its own sigma metrics, bed condition and free-surface Dirichlet value.  The
//  unknowns are the nodes of the leaf columns (level-0 and patch columns not covered by a finer
//  patch).  A leaf row is the grid's own matrix row; its neighbour columns outside the grid come
//  from the finer level (covered columns: average of the 2x2 children, node by node) or from
//  the coarser level (the columns around a patch: biquadratic, vertically linear in sigma; a
//  sibling patch or the owner rank across a partition edge).
//
//  BiCGStab (reefamr_bicgstab), right preconditioned by one multiplicative FAC sweep (lap_prec):
//  patch-local REEFMG V-cycles (MPI_COMM_SELF, correction 0 at the patch edge) as pre-smoothing,
//  a REEFMG V-cycle on level 0 (the whole rank grid, covered columns included) on the composite
//  residual they leave, then level by level the coarse increment interpolated into the patches
//  and patch-local V-cycles as post-smoothing (LPASS passes, the siblings' corrections filled in
//  between).

namespace
{
typedef reefamr_comms_off comms_off;

enum { NR=REEFAMR_NR, NRH=REEFAMR_NRH, NPV=REEFAMR_NPV, NVV=REEFAMR_NVV, NS=REEFAMR_NS, NT=REEFAMR_NT,
       NPH=REEFAMR_NPH, NSH=REEFAMR_NSH,
       NRP=REEFAMR_NTMP,        // preconditioner: residual after the pre-smoothing, prolonged increment
       NPRE=REEFAMR_NVEC,       // preconditioner: z after the pre-smoothing (patches)
       NVEC=REEFAMR_NVEC+1 };

// patch-local passes per level: more passes let sibling patches see each other's corrections
// before the next Krylov step; for patches split at partition edges the iteration count did
// not change with 2 or 4 passes (the level-0 correction couples them), only the cost
const int LPASS = 1;

// smoothing sweeps (pre and post) of the patch-local V-cycles of the FAC sweep.  With the pre-
// smoothing one sweep is enough (moored box, 700 steps: max 4 iterations, 8 % less Laplace time
// than two sweeps)
const int PSWEEP = 1;

// safety net: a solve that stagnates again after a BiCGStab restart stops if its relative
// residual is below STAGTOL*N 44 (0: run to N 46).  With the coarse-first FAC sweep of the
// first version about 2 % of the psi solves of a moored box stagnated at 1e-8..3e-7 relative;
// with the pre-smoothing neither the restart nor the stop has triggered
const double STAGTOL = 100.0;

inline long ijk4(const lexer *q, int ii, int jj, int kk)
{
    return (long)(ii-q->imin)*q->jmax*q->kmax + (long)(jj-q->jmin)*q->kmax + kk-q->kmin;
}
}

double* fnpf_amr::lvec(int g, int k)
{
    if(k<0)
    return ltgt[g+1];

    return (g<0) ? kv0[k] : FP(g)->kv[k];
}

// row numbering as the assembly in fnpf_laplace_cds2 (ILOOP JLOOP KLOOP, flag4>0); leaf and
// covered unknowns of every grid
void fnpf_amr::lap_rows()
{
    auto number = [&](lexer *q, vector<int> &row)
    {
        row.assign(q->imax*q->jmax*(q->kmax+2),-1);
        int n=0;
        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        for(int kk=0; kk<q->knoz; ++kk)
        if(q->flag4[ijk4(q,ii,jj,kk)]>0)
        row[fidx(q,ii,jj,kk)] = n++;
    };

    number(p0,rowmap0);
    for(auto q : P)
    number(q->pp,FP(q)->row);

    lg.assign(P.size()+1,lgrid());

    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        vector<int> &row = (g<0) ? rowmap0 : FP(g)->row;
        lgrid &L = lg[g+1];
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
void fnpf_amr::lap_prepare(ghostcell *pgc)
{
    if(lap_layout!=layout_id)
    {
        const double ts = MPI_Wtime();

        lap_rows();

        const int n0 = p0->imax*p0->jmax*(p0->kmax+2);
        if(kv0.empty())
        for(int k=0; k<NVEC; ++k)
        kv0.push_back(new double[n0]());

        for(auto q : P)
        {
            fnpf_amr_patch *c = FP(q);
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
                cout<<"FNPF AMR: patch multigrid "<<c->mg->err()<<endl;
            }
        }

        if(mg0==nullptr)
        {
            mg0 = new reefmg_core;
            mg0->set_precision(64);
            if(!mg0->setup(pgc->cart(),p0->knox,p0->knoy,p0->knoz,p0->gknox,p0->gknoy,p0->N13,p0->DXN+marge-1,p0->DYN+marge-1))
            {
                if(p0->mpirank==0)
                cout<<"FNPF AMR: level-0 multigrid "<<mg0->err()<<endl;
                MPI_Abort(MPI_COMM_WORLD,-2760);
            }
        }

        lap_layout = layout_id;
        tm[5] += MPI_Wtime()-ts;
    }

    // fixed rows (identity rows of the assembly: nodes inside a resolved body) are no unknowns:
    // they keep their value, and no fluid row couples to them (Neumann faces).  Left in the
    // composite system they took up the coarse correction across the hull and the solve
    // needed 200 iterations instead of 5
    const bool hasbody = (body!=nullptr);
    auto fixed = [&](const matrix_diag &M, int r)
    {
        return hasbody && M.p[r]==1.0 && M.n[r]==0.0 && M.s[r]==0.0 && M.w[r]==0.0 && M.e[r]==0.0 && M.t[r]==0.0 && M.b[r]==0.0;
    };

    for(int g=-1; g<(int)P.size(); ++g)
    {
        lgrid &L = lg[g+1];
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

    // level 0: all rows of the rank grid, 7-point part
    {
        sc_level &L = mg0->fine();
        clear(L);
        const matrix_diag &M = c0->M;

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
        fnpf_amr_patch *c = FP(q);
        lexer *pp = c->pp;
        sc_level &L = c->mg->fine();
        clear(L);
        const matrix_diag &M = c->c->M;

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
void fnpf_amr::lap_sync(int k)
{
    auto sel = [&](int g) -> double* { return lvec(g,k); };

    restrict_col(sel);

    if(p0->mpi_size>1)
    {
        double *x0 = lvec(-1,k);
        pgc0->gcparax7(p0,x0,7);
        pgc0->gcparax7co(p0,x0,7);
    }

    for(int l=1; l<=maxlev; ++l)
    fill_col(l,7600+l,sel);
}

// y = A x on the leaf unknowns
void fnpf_amr::lap_apply(int kx, int ky)
{
    lap_sync(kx);

    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        const matrix_diag &M = gfd(g)->M;
        const double *x = lvec(g,kx);
        double *y = lvec(g,ky);
        const int sI = q->jmax*q->kmaxF;
        const int sJ = q->kmaxF;
        const lgrid &L = lg[g+1];

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

double fnpf_amr::lap_dot(int ka, int kb)
{
    double s=0.0;
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *a = lvec(g,ka), *b = lvec(g,kb);
        for(int qq : lg[g+1].lq)
        s += a[qq]*b[qq];
    }
    return pgc0->globalsum(s);
}

void fnpf_amr::lap_start()
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vb = lvec(g,NS);
        double *vr = lvec(g,NR), *vrh = lvec(g,NRH), *vp = lvec(g,NPV), *vv = lvec(g,NVV);
        for(int qq : lg[g+1].lq)
        {
            double r = vb[qq] - vv[qq];
            vr[qq] = r; vrh[qq] = r; vp[qq] = 0.0; vv[qq] = 0.0;
        }
    }
}

// restart of the BiCGStab recurrence from the current residual (reefamr_bicgstab)
void fnpf_amr::lap_restart()
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vr = lvec(g,NR);
        double *vrh = lvec(g,NRH), *vp = lvec(g,NPV), *vv = lvec(g,NVV);
        for(int qq : lg[g+1].lq)
        {
            vrh[qq] = vr[qq];
            vp[qq] = 0.0;
            vv[qq] = 0.0;
        }
    }
}

void fnpf_amr::lap_p(double beta, double om)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vr = lvec(g,NR), *vv = lvec(g,NVV);
        double *vp = lvec(g,NPV);
        for(int qq : lg[g+1].lq)
        vp[qq] = vr[qq] + beta*(vp[qq] - om*vv[qq]);
    }
}

void fnpf_amr::lap_s(double alp)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vr = lvec(g,NR), *vv = lvec(g,NVV);
        double *vs = lvec(g,NS);
        for(int qq : lg[g+1].lq)
        vs[qq] = vr[qq] - alp*vv[qq];
    }
}

void fnpf_amr::lap_x(double alp, double om)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vph = lvec(g,NPH), *vsh = lvec(g,NSH), *vs = lvec(g,NS), *vt = lvec(g,NT);
        double *x = lvec(g,-1), *vr = lvec(g,NR);
        for(int qq : lg[g+1].lq)
        {
            x[qq] += alp*vph[qq] + om*vsh[qq];
            vr[qq] = vs[qq] - om*vt[qq];
        }
    }
}

// patch-local correction: z += V-cycle on (r - A z) over the interior of patch id, with the
// current z around the patch (correction 0 at the patch edge)
void fnpf_amr::lap_local(int id, int kr, int kz)
{
    fnpf_amr_patch *c = FP(id);
    lexer *pp = c->pp;
    const matrix_diag &M = c->c->M;
    const double *r = lvec(id,kr);
    double *z = lvec(id,kz);
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

    c->mg->vcycle(0,PSWEEP,PSWEEP);

    for(int ii=EXT; ii<EXT+c->nx; ++ii)
    for(int jj=EXT; jj<EXT+c->ny; ++jj)
    for(int kk=0; kk<pp->knoz; ++kk)
    {
        const long lq = L.idx(ii-EXT,jj-EXT,kk);
        if(L.act[lq])
        z[fidx(pp,ii,jj,kk)] += L.u[lq];
    }
}

// z = M^-1 r: one multiplicative FAC sweep
//  1. pre-smoothing: z = 0, patch-local V-cycles on r (all levels)
//  2. composite residual r' = r - A z, restricted into the covered columns
//  3. level 0: one V-cycle on r' over the whole rank grid (covered columns included)
//  4. coarse to fine: the increment of the coarser grid since its pre-smoothing interpolated
//     into the patch interiors and added, the columns around the patches filled, then
//     post-smoothing with patch-local V-cycles on r - A z
//  Without step 1 (coarse correction first) the residual restricted next to a resolved body
//  carried the unsmoothed fine error at the hull into the level-0 correction, and the roll and
//  pitch unit-mode psi solves stagnated (1e-8..3e-7 relative, up to N 46 iterations)
void fnpf_amr::lap_prec(int kr, int kz)
{
    auto selz = [&](int g) -> double* { return lvec(g,kz); };

    // 1. pre-smoothing
    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        double *z = lvec(g,kz);
        std::fill(z,z+q->imax*q->jmax*(q->kmax+2),0.0);
    }

    for(int id=0; id<(int)P.size(); ++id)
    lap_local(id,kr,kz);

    // 2. composite residual (lap_apply also fills the covered columns and the columns around
    //    the patches of z), z after the pre-smoothing kept in NPRE
    lap_apply(kz,NRP);
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *r = lvec(g,kr);
        double *t = lvec(g,NRP);
        for(int qq : lg[g+1].lq)
        t[qq] = r[qq] - t[qq];
    }

    for(int id=0; id<(int)P.size(); ++id)
    {
        lexer *q = glex(id);
        std::copy(lvec(id,kz),lvec(id,kz)+q->imax*q->jmax*(q->kmax+2),lvec(id,NPRE));
    }

    restrict_col([&](int g) -> double* { return lvec(g,NRP); });

    // 3. level 0: one V-cycle on the whole rank grid
    {
        sc_level &L = mg0->fine();
        const double *r = lvec(-1,NRP);
        double *z = lvec(-1,kz);

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

    // 4. coarse to fine
    for(int l=1; l<=maxlev; ++l)
    {
        // increment of the level l-1 patches since their pre-smoothing (level 0: z itself)
        if(l>1)
        for(int id : lev[l-1])
        {
            lexer *q = glex(id);
            const double *z = lvec(id,kz);
            double *d = lvec(id,NPRE);
            const int n = q->imax*q->jmax*(q->kmax+2);
            for(int m=0; m<n; ++m)
            d[m] = z[m] - d[m];
        }

        // the increment interpolated into the patch interiors (NRP) and added
        for(int id : lev[l])
        {
            fnpf_amr_patch *c = FP(id);
            lexer *pp = c->pp;
            prolong_interior_col(*c,[&](int g) -> double*
            {
                if(g<0)
                return lvec(-1,kz);
                return (P[g]->lev==l) ? lvec(g,NRP) : lvec(g,NPRE);
            });

            double *z = lvec(id,kz);
            const double *t = lvec(id,NRP);
            for(int ii=EXT; ii<EXT+c->nx; ++ii)
            for(int jj=EXT; jj<EXT+c->ny; ++jj)
            for(int kk=0; kk<=pp->knoz; ++kk)
            {
                const int qq = fidx(pp,ii,jj,kk);
                z[qq] += t[qq];
            }
        }

        fill_col(l,7700+l,selz);

        // post-smoothing; repeated passes let siblings see each other's corrections through
        // the filled columns
        const int npass = (nlevg[l]>1) ? LPASS : 1;
        for(int pass=0; pass<npass; ++pass)
        {
            if(pass>0)
            fill_col(l,7700+l,selz);

            for(int id : lev[l])
            lap_local(id,kr,kz);
        }

        if(l<maxlev)
        fill_col(l,7700+l,selz);
    }
}

namespace
{
// the composite Laplace as the vector space of reefamr_bicgstab
struct lap_space
{
    static const int B=NS, R=NR, RH=NRH, PV=NPV, VV=NVV, S=NS, T=NT, PH=NPH, SH=NSH;
    fnpf_amr *a;
    void apply(int x, int y) { a->lap_apply(x,y); }
    void prec(int r, int z) { a->lap_prec(r,z); }
    double dot(int x, int y) { return a->lap_dot(x,y); }
    void op_start() { a->lap_start(); }
    void op_p(double beta, double om) { a->lap_p(beta,om); }
    void op_s(double alp) { a->lap_s(alp); }
    void op_x(double alp, double om) { a->lap_x(alp,om); }
    void op_restart() { a->lap_restart(); }
};
}

// assembled rows of all grids -> solution in ltgt: right-hand side, BiCGStab (initial guess:
// the values in ltgt), then the covered columns, the columns around the patches and the halo
void fnpf_amr::lap_core(lexer *p, ghostcell *pgc)
{
    lap_prepare(pgc);

    // right-hand side into NS
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *R = gfd(g)->rhsvec.V.data();
        double *b = lvec(g,NS);
        const lgrid &L = lg[g+1];
        for(size_t n=0; n<L.lq.size(); ++n)
        b[L.lq[n]] = R[L.lr[n]];
        for(int qq : L.cq)
        b[qq] = 0.0;
    }

    lap_apply(-1,NVV);

    double bn, rn;
    lap_space sp{this};
    int it = reefamr_bicgstab(sp,p->N44,p->N46,bn,rn,&tm[2],&tm[3],&lap_restarts,STAGTOL*p->N44,&lap_stalls);

    const double res = bn>0.0 ? rn/bn : 0.0;
    lap_it_total += it;
    ++lap_solves;
    lap_it_step += it;
    ++lap_solves_step;
    if(lap_kind==0)
    {
        lap_it_last = it;
        lap_res_last = res;
        lap_it_phi_max = MAX(lap_it_phi_max,it);
    }
    else
    lap_it_psi_max = MAX(lap_it_psi_max,it);
    p->solveriter = it;
    p->final_res = res;

    // the uniform REEFMG path warns at the iteration limit; a stalled composite solve was silent
    if(it>=p->N46 && res>p->N44)
    {
        ++lap_capped;
        if(p->mpirank==0)
        cout<<"FNPF AMR WARNING - "<<(lap_kind==0 ? "phi" : "psi")<<" Laplace iteration limit N 46 = "<<p->N46
            <<" reached, res "<<setprecision(3)<<res<<" (N 44 = "<<p->N44<<")"<<endl;
    }

    lap_sync(-1);
}

// the Laplace equation of the RK stage on all grids (fnpf_RK3, through fnpf_laplace_amr)
void fnpf_amr::lap_solve(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf, double *f, slice &Fifsf)
{
    const double t0 = MPI_Wtime();

    // the patches: free-surface derivatives, sigma grid, boundary conditions
    {
    comms_off guard(pgc);
    for(int id=0; id<(int)P.size(); ++id)
    {
        fnpf_amr_patch *pc = FP(id);
        lexer *pp = pc->pp;
        fdm_fnpf *cc = pc->c;
        slice &Se = sval(id,stg,0);
        slice &Sf = sval(id,stg,1);

        pc->pf->fsfdisc(pp,cc,pgc,Se,Sf);
        pc->psig->sigma_update(pp,cc,pgc,pc->pf,Se);
        pc->pfu->fsfbc_sig(pp,cc,pgc,Sf,cc->Fi);
        pc->pbu->bedbc_sig(pp,cc,pgc,cc->Fi,pc->pf);
        walls_fi(*pc,cc->Fi);
    }
    }

    // resolved body on the new sigma grids
    if(body!=nullptr)
    body->amr_geometry(p,c,pgc);

    // rows: level 0 with its own solver, the patches
    plap0->assemble_only = true;
    plap0->start(p,c,pgc,psolv,pf,f,Fifsf);
    plap0->assemble_only = false;

    {
    comms_off guard(pgc);
    for(int id=0; id<(int)P.size(); ++id)
    {
        fnpf_amr_patch *pc = FP(id);
        pc->plap->start(pc->pp,pc->c,pgc,nullptr,pc->pf,pc->c->Fi,sval(id,stg,1));
    }
    }

    // initial guess: Fi of the previous stage
    ltgt.assign(P.size()+1,nullptr);
    ltgt[0] = f;
    for(int id=0; id<(int)P.size(); ++id)
    ltgt[id+1] = FP(id)->c->Fi;

    lap_kind = 0;
    lap_core(p,pgc);

    // body band
    if(body!=nullptr)
    body->amr_post_solve(p,c,pgc,f);

    // the patches: wall nodes, vertical velocity at the free surface
    {
    comms_off guard(pgc);
    for(int id=0; id<(int)P.size(); ++id)
    {
        fnpf_amr_patch *pc = FP(id);
        walls_fi(*pc,pc->c->Fi);
        pc->pf->fsfwvel(pc->pp,pc->c,pgc,sval(id,stg,0),sval(id,stg,1));
    }
    }

    tm[1] += MPI_Wtime()-t0;

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"FNPF AMR Laplace: iterations "<<lap_it_last<<"  res "<<setprecision(3)<<lap_res_last<<endl;
}

// a psi solve of the body loads on all grids, with the operator of the current sigma grids
void fnpf_amr::lap_solve_psi(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf, double **f, slice **D)
{
    const double t0 = MPI_Wtime();

    plap0->assemble_only = true;
    plap0->start(p,c,pgc,psolv,pf,f[0],*D[0]);
    plap0->assemble_only = false;

    {
    comms_off guard(pgc);
    for(int id=0; id<(int)P.size(); ++id)
    {
        fnpf_amr_patch *pc = FP(id);
        pc->plap->start(pc->pp,pc->c,pgc,nullptr,pc->pf,f[id+1],*D[id+1]);
    }
    }

    ltgt.assign(f,f+P.size()+1);

    lap_kind = 1;
    lap_core(p,pgc);
    lap_kind = 0;

    tm[1] += MPI_Wtime()-t0;
}
