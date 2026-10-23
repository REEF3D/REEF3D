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
#include<cstring>

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

// float coefficients of the V-cycles of the preconditioner: reefmg_core wants the exact fine
// operator for its own Krylov solver and residual(); the FAC sweep only calls vcycle(), so this
// one is never used
struct lap_fineop : public sc_operator
{
    void fine_apply(const sc_level&, const double*, double*) override
    {
        cout<<"FNPF AMR: fine operator of the preconditioner called"<<endl;
        MPI_Abort(MPI_COMM_WORLD,-2762);
    }
};
lap_fineop lap_fineop_none;

// coefficients of row lq of the fine level, in the storage precision of the multigrid
inline void lap_coef(reefmg_core *mg, sc_level &L, long lq, double p, double n, double s, double w, double e, double t, double b)
{
    if(mg->precision()==32)
    {
        L.pf[lq]=(float)p; L.nf[lq]=(float)n; L.sf[lq]=(float)s; L.wf[lq]=(float)w; L.ef[lq]=(float)e; L.tf[lq]=(float)t; L.bf[lq]=(float)b;
    }
    else
    {
        L.p[lq]=p; L.n[lq]=n; L.s[lq]=s; L.w[lq]=w; L.e[lq]=e; L.t[lq]=t; L.b[lq]=b;
    }
}

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

    // (X 36: a window solve, its patches; the patches outside are not touched)
    if(g>=0 && lap_nol0 && !FP(g)->kvw.empty())
    return FP(g)->kvw[k];

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

    // (X 36, a window without level 0: level 0 keeps its numbering of the layout)
    if(!(lap_nol0 && rowmap0_layout==layout_id))
    {
        number(p0,rowmap0);
        rowmap0_layout = layout_id;
    }
    // (X 36: not the patches outside the window of the solve; a solve on all grids numbers them again)
    for(int id=0; id<(int)P.size(); ++id)
    if(!skipped(id))
    number(P[id]->pp,FP(id)->row);

    // the fill entries next to the patch: the columns the rows of the patch reach
    for(int id=0; id<(int)P.size(); ++id)
    {
        if(skipped(id))
        continue;
        fnpf_amr_patch *c = FP(id);
        c->fill1.clear();
        c->fill2.clear();
        for(size_t e=0; e<c->fill.size(); ++e)
        {
            const reefamr_fill &f = c->fill[e];
            if(f.di>=EXT-1 && f.di<=EXT+c->nx && f.dj>=EY-1 && f.dj<=EY+c->ny)
            c->fill1.push_back((int)e);
            if(f.di>=EXT-2 && f.di<=EXT+c->nx+1 && f.dj>=EY-2 && f.dj<=EY+c->ny+1)
            c->fill2.push_back((int)e);
        }
    }

    lg.assign(P.size()+1,lgrid());

    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        vector<int> &row = (g<0) ? rowmap0 : FP(g)->row;
        lgrid &L = lg[g+1];
        const int l = (g<0) ? 0 : P[g]->lev;

        // G 7 1: the grids of the level window only (wlo..wtop: all levels without subcycling); X 36:
        // not the patches outside the window of the unit-mode solve
        if(l<wlo || l>wtop() || (g>=0 && skipped(g)))
        continue;

        int i0,i1,j0,j1;
        if(g<0) { i0=0; i1=NX0-1; j0=0; j1=NY0-1; }
        else    { i0=EXT; i1=EXT+P[g]->nx-1; j0=EY; j1=EY+P[g]->ny-1; }

        for(int ii=i0; ii<=i1; ++ii)
        for(int jj=j0; jj<=j1; ++jj)
        {
            int I = (g<0) ? ii+O0i : ii-EXT+P[g]->I0;
            int J = (g<0) ? jj+O0j : jj-EY+P[g]->J0;
            // (X 36: covered by a patch of the window only; under a patch outside it the columns are
            // leaf columns of the window, the hulls there on this grid's resolution)
            const int gc = (l<wtop()) ? gpatch_at(l+1,2*I,ychild(J)) : -1;
            const bool covered = (gc>=0 && (skip==nullptr || lap_wgp==nullptr || (*lap_wgp)[gc]));

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
    const int wkey = wlo*1000 + wtop() + lap_wextra;
    if(lap_layout==layout_id && lap_wkey!=wkey)
    {
        lap_rows();
        lap_wkey = wkey;
    }

    if(lap_layout!=layout_id)
    {
        const double ts = MPI_Wtime();

        lap_rows();
        lap_wkey = wkey;

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
                c->mg->set_precision(pp->N14==32 ? 32 : 64);     // (N 14, as the uniform solver)
                c->mg->set_fine_operator(&lap_fineop_none);
                if(!c->mg->setup(MPI_COMM_SELF,c->nx,c->ny,pp->knoz,c->nx,c->ny,0,pp->DXN+marge+EXT-1,pp->DYN+marge+EY-1))
                cout<<"FNPF AMR: patch multigrid "<<c->mg->err()<<endl;
            }
        }

        if(mg0==nullptr)
        {
            mg0 = new reefmg_core;
            mg0->set_precision(p0->N14==32 ? 32 : 64);
            mg0->set_fine_operator(&lap_fineop_none);
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

    // the rows of the last preparation: with a body (psi solves) or Picard passes (A 324) a stage
    // has several solves with the same rows
    if((body!=nullptr || p0->A324>0) && lap_same(pgc,wkey))
    {
        ++lap_kept;
        return;
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
        std::fill(L.pf.begin(),L.pf.end(),0.0f);
        std::fill(L.nf.begin(),L.nf.end(),0.0f); std::fill(L.sf.begin(),L.sf.end(),0.0f);
        std::fill(L.wf.begin(),L.wf.end(),0.0f); std::fill(L.ef.begin(),L.ef.end(),0.0f);
        std::fill(L.tf.begin(),L.tf.end(),0.0f); std::fill(L.bf.begin(),L.bf.end(),0.0f);
        std::fill(L.u.begin(),L.u.end(),0.0); std::fill(L.f.begin(),L.f.end(),0.0);
        std::fill(L.act.begin(),L.act.end(),0);
    };

    // level 0: all rows of the rank grid, 7-point part (X 36, a window without level 0: none)
    if(!lap_nol0 && !(sy==0 && xz_keep(xz0key,xz0on,wkey)))
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
                lap_coef(mg0,L,lq,1.0,0.0,0.0,0.0,0.0,0.0,0.0);
                continue;
            }
            lap_coef(mg0,L,lq,M.p[r],M.n[r],M.s[r],M.w[r],M.e[r],M.t[r],M.b[r]);
            L.act[lq] = fixed(M,r) ? 0 : 1;
        }
        // a 2D grid: level 0 is one x-z problem, solved exactly (reefamr_xzlu; several ranks: the
        // problem of all ranks on every rank)
        xz0on = false;
        if(sy==0 && p0->mpi_size==1)
        {
            xz0.factor(L,mg0->precision());
            xz0on = true;
        }
        else
        if(sy==0)
        xz0on = xz0.factor_mpi(L,mg0->precision(),MPI_COMM_WORLD,p0->origin_i,p0->gknox);
        if(!xz0on)
        mg0->coarsen();
        xz0key[0] = layout_id; xz0key[1] = wkey; xz0key[2] = p0->count;
    }

    // patches: interior rows, couplings to the cells around the patch dropped
    for(int id=0; id<(int)P.size(); ++id)
    {
        reefamr_patch *q = P[id];
        if(q->lev<wlo || q->lev>wtop() || skipped(id))
        continue;

        fnpf_amr_patch *c = FP(q);
        if(sy==0 && xz_keep(c->xzkey,true,wkey))
        continue;
        lexer *pp = c->pp;
        sc_level &L = c->mg->fine();
        clear(L);
        const matrix_diag &M = c->c->M;

        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EY; jj<EY+c->ny; ++jj)
        for(int kk=0; kk<pp->knoz; ++kk)
        {
            const long lq = L.idx(ii-EXT,jj-EY,kk);
            const int r = c->row[fidx(pp,ii,jj,kk)];
            if(r<0)
            {
                lap_coef(c->mg,L,lq,1.0,0.0,0.0,0.0,0.0,0.0,0.0);
                continue;
            }
            lap_coef(c->mg,L,lq,M.p[r],(ii+1<EXT+c->nx) ? M.n[r] : 0.0,(ii-1>=EXT) ? M.s[r] : 0.0,
                     (jj+1<EY+c->ny) ? M.w[r] : 0.0,(jj-1>=EY) ? M.e[r] : 0.0,M.t[r],M.b[r]);
            L.act[lq] = fixed(M,r) ? 0 : 1;
        }
        // a 2D grid: the patch-local problem is an x-z problem, solved exactly (reefamr_xzlu) instead
        // of the V-cycle
        if(sy==0)
        {
            c->xz.factor(L,c->mg->precision());
            c->xzkey[0] = layout_id; c->xzkey[1] = wkey; c->xzkey[2] = p0->count;
        }
        else
        c->mg->coarsen();
    }
}

// the solves of a stage share the rows: phi, then psi_0 and the unit modes of the body loads
// (the Neumann data of the body faces enter the right-hand side only), and the Picard passes
// of A 324.  The 7-point rows of level 0 and of the patches of the window are compared bit
// for bit with those of the last preparation; unchanged on every rank, the preparation is
// skipped and the multigrid coefficients, their coarse levels and the x-z factorisations
// are those of the last solve, the same preconditioner.  With a body on a 2D grid this
// replaces the factorisation of every solve (xz_keep).  Collective: mg0 and the level-0
// factorisation are set up by all ranks together, so the decision is the same on all ranks
bool fnpf_amr::lap_same(ghostcell *pgc, int wkey)
{
    if(lmr.size()!=P.size()+1)
    lmr.assign(P.size()+1,lmrows());

    int same = 1;

    for(int g=-1; g<(int)P.size(); ++g)
    {
        if(g>=0 && (P[g]->lev<wlo || P[g]->lev>wtop() || skipped(g)))
        continue;

        // (X 36, a window without level 0: its rows are not part of the solve)
        if(g<0 && lap_nol0)
        continue;

        const matrix_diag &M = gfd(g)->M;
        const size_t n = M.p.size();
        const double *a[7] = {M.p.data(),M.n.data(),M.s.data(),M.w.data(),M.e.data(),M.t.data(),M.b.data()};
        lmrows &R = lmr[g+1];

        if(R.layout!=layout_id || R.wkey!=wkey || R.m.size()!=7*n)
        {
            R.layout = layout_id;
            R.wkey = wkey;
            R.m.resize(7*n);
            same = 0;
        }

        for(int k=0; k<7; ++k)
        {
            double *s = R.m.data() + k*n;

            if(same && memcmp(s,a[k],n*sizeof(double))!=0)
            same = 0;

            if(!same)
            memcpy(s,a[k],n*sizeof(double));
        }
    }

    return pgc->globalimin(same)==1;
}

// a 2D grid: the exact factorisation of a grid is kept for the other solves of the step (key: the
// layout, the window and the step of the factorisation).  The coefficients follow the free
// surface: the factorisation of the first solve of the step is the preconditioner for its other
// solves (2D flume, two levels: 0.3 % more iterations, 14 % less time; kept for 4 steps: 1.7 %
// more iterations and no further gain).  With a body (fixed rows) refactorised when the rows
// change (lap_same keeps it for the solves of a stage)
bool fnpf_amr::xz_keep(long *key, bool on, int wkey)
{
    if(!on || body!=nullptr)
    return false;
    return key[0]==layout_id && key[1]==wkey && key[2]==p0->count;
}

// covered columns, partition halo of level 0, columns around the patches
void fnpf_amr::lap_sync(int k, bool d1)
{
    auto sel = [&](int g) -> double* { return lvec(g,k); };

    lap_dir = (k>=0);

    restrict_col(sel);

    if(p0->mpi_size>1 && wlo==0)
    {
        double *x0 = lvec(-1,k);
        pgc0->gcparax7(p0,x0,7);
        pgc0->gcparax7co(p0,x0,7);
    }

    for(int l=MAX(wlo,1); l<=wtop(); ++l)
    fill_col(l,7600+l,sel,d1);
}

// y = A x on the leaf unknowns
void fnpf_amr::lap_apply(int kx, int ky, bool d1)
{
    lap_sync(kx,d1);

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

        // A 328 1 (N 10 1): the sigma cross-derivatives in the matrix (fnpf_laplace_cds2), the couplings
        // to the nodes above and below the four neighbour nodes; those in the columns around a patch
        // and the covered columns are filled (lap_sync) as for the 7-point part
        if(lap_mix && M.nt.size()==M.p.size())
        {
            const double *Ant=M.nt.data(), *Anb=M.nb.data(), *Ast=M.st.data(), *Asb=M.sb.data();
            const double *Awt=M.wt.data(), *Awb=M.wb.data(), *Aet=M.et.data(), *Aeb=M.eb.data();
            for(size_t n=0; n<L.lq.size(); ++n)
            {
                const int qq = L.lq[n], r = L.lr[n];
                y[qq] += Ant[r]*x[qq+sI+1] + Anb[r]*x[qq+sI-1] + Ast[r]*x[qq-sI+1] + Asb[r]*x[qq-sI-1]
                       + Awt[r]*x[qq+sJ+1] + Awb[r]*x[qq+sJ-1] + Aet[r]*x[qq-sJ+1] + Aeb[r]*x[qq-sJ-1];
            }
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
        
        // the composite residual of the preconditioner: the restriction to level 0 reads entries
        // the residual of this solve does not write (they held those of an earlier solve, e.g. one
        // of another level window of the subcycling): the solve did not depend on its own data only,
        // and a restart from a checkpoint (P 48, I 45) did not continue bitwise.  (X 36: a window
        // without level 0 restricts and reads neither level 0 nor the patches outside it)
        if((g<0 && lap_nol0) || (g>=0 && skipped(g)))
        continue;
        lexer *q = glex(g);
        double *vt = lvec(g,NRP);
        std::fill(vt,vt+q->imax*q->jmax*(q->kmax+2),0.0);
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

    // (the exact solve of a 2D grid reads f of the active rows only and writes u)
    if(sy==1)
    {
    std::fill(L.u.begin(),L.u.end(),0.0);
    std::fill(L.f.begin(),L.f.end(),0.0);
    }
    for(int ii=EXT; ii<EXT+c->nx; ++ii)
    for(int jj=EY; jj<EY+c->ny; ++jj)
    for(int kk=0; kk<pp->knoz; ++kk)
    {
        const long lq = L.idx(ii-EXT,jj-EY,kk);
        if(L.act[lq]==0)
        continue;

        const int qq = fidx(pp,ii,jj,kk);
        const int rw = c->row[qq];
        const double az = M.p[rw]*z[qq] + M.n[rw]*z[qq+sI] + M.s[rw]*z[qq-sI] + M.w[rw]*z[qq+sJ] + M.e[rw]*z[qq-sJ]
                        + M.t[rw]*z[qq+1] + M.b[rw]*z[qq-1];
        L.f[lq] = r[qq] - az;
    }

    if(sy==0)
    c->xz.solve(L,L.f,L.u);
    else
    c->mg->vcycle(0,PSWEEP,PSWEEP);

    for(int ii=EXT; ii<EXT+c->nx; ++ii)
    for(int jj=EY; jj<EY+c->ny; ++jj)
    for(int kk=0; kk<pp->knoz; ++kk)
    {
        const long lq = L.idx(ii-EXT,jj-EY,kk);
        if(L.act[lq])
        z[fidx(pp,ii,jj,kk)] += L.u[lq];
    }
}

// the level-0 cycle of the preconditioner on mg0->fine(): a V-cycle, on a 2D grid the
// exact solve
void fnpf_amr::lap_l0solve()
{
    sc_level &L = mg0->fine();
    if(xz0on)
    xz0.solve(L,L.f,L.u);
    else
    mg0->vcycle(0,1,1);
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
    if(wlo>0)
    {
        lap_prec_win(kr,kz);
        return;
    }

    auto selz = [&](int g) -> double* { return lvec(g,kz); };

    // 1. pre-smoothing (z 0: level 0; the patches and the columns next to them, the cells the
    //    patch-local cycles read; lap_apply below fills the columns around the patches)
    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        double *z = lvec(g,kz);
        if(g<0)
        {
            std::fill(z,z+q->imax*q->jmax*(q->kmax+2),0.0);
            continue;
        }
        // (the columns jj of one ii are contiguous: one fill per ii)
        for(int ii=EXT-1; ii<=EXT+P[g]->nx; ++ii)
        {
            const int a = fidx(q,ii,EY-1,q->kmin);
            std::fill(z+a,z+a+(P[g]->ny+2)*q->kmaxF,0.0);
        }

        // NRP of the patches: the composite residual (step 2) writes the leaf and covered rows
        // only, but restrict_col reads every child node.  The other nodes (resolved-body and dry
        // rows) held the increment interpolated in step 4 of the previous sweep, which entered the
        // covered columns of the coarser grid: the preconditioner changed from one Krylov step to
        // the next by an error of the size of z/L^2 against r, negligible at model scale (dx 0.1 m)
        // and dominant at full scale (dx 2.5 m and more), where the psi solves diverged to NaN
        std::fill(lvec(g,NRP),lvec(g,NRP)+q->imax*q->jmax*(q->kmax+2),0.0);
    }

    // (only with a resolved body: without one, the sweep without the pre-smoothing takes as many
    // or fewer iterations and half the patch-local cycles and no composite residual here)
    const bool pre = (body!=nullptr);
    if(pre)
    {
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
    }
    else
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *r = lvec(g,kr);
        double *t = lvec(g,NRP);
        for(int qq : lg[g+1].lq)
        t[qq] = r[qq];
    }

    // (a 2D grid: the rows next to the strip only, the prolongation reads the strip)
    for(int id=0; id<(int)P.size(); ++id)
    {
        lexer *q = glex(id);
        if(sy==0)
        for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
        {
            const int a = fidx(q,ii,EY-1,q->kmin);
            std::copy(lvec(id,kz)+a,lvec(id,kz)+a+(P[id]->ny+2)*q->kmaxF,lvec(id,NPRE)+a);
        }
        else
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

        lap_l0solve();

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
            if(sy==0)
            for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
            {
                const int a = fidx(q,ii,EY-1,q->kmin), n = (P[id]->ny+2)*q->kmaxF;
                for(int m=a; m<a+n; ++m)
                d[m] = z[m] - d[m];
            }
            else
            {
            const int n = q->imax*q->jmax*(q->kmax+2);
            for(int m=0; m<n; ++m)
            d[m] = z[m] - d[m];
            }
        }

        // the increment interpolated into the patch interiors (NRP, from the parent on its rank)
        // and added
        // (biquadratic: the bicubic stencil gave the same iterations in 3D and more on a 2D grid,
        // 4.6 against 3.7 with two levels, at 16 instead of 9 coarse columns per fine column)
        const int pord_keep = pord;
        pord = 3;
        prolong_interior_col(l,[](reefamr_patch*) { return true; },
                             [&](int g) -> double* { return (g<0) ? lvec(-1,kz) : lvec(g,NPRE); },
                             [&](int id) -> double* { return lvec(id,NRP); });
        pord = pord_keep;

        for(int id : lev[l])
        {
            fnpf_amr_patch *c = FP(id);
            lexer *pp = c->pp;

            double *z = lvec(id,kz);
            const double *t = lvec(id,NRP);
            for(int ii=EXT; ii<EXT+c->nx; ++ii)
            for(int jj=EY; jj<EY+c->ny; ++jj)
            for(int kk=0; kk<=pp->knoz; ++kk)
            {
                const int qq = fidx(pp,ii,jj,kk);
                z[qq] += t[qq];
            }
        }

        // (the columns next to the patches, as the patch-local cycles read them; the full fill
        // below for the next level)
        fill_col(l,7700+l,selz,true);

        // post-smoothing; repeated passes let siblings see each other's corrections through
        // the filled columns
        const int npass = (nlevg[l]>1) ? LPASS : 1;
        for(int pass=0; pass<npass; ++pass)
        {
            if(pass>0)
            fill_col(l,7700+l,selz,true);

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
    void apply(int x, int y) { a->lap_apply(x,y,true); }
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
    const double ts = MPI_Wtime();
    lap_prepare(pgc);
    lap_tsetup(MPI_Wtime()-ts);

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
    // the psi solves of the floating bodies (fnpf_6DOF) stop at the half step once below N 44.  The
    // phi solves keep the full steps: there an earlier stop changes the iterates within N 44, which
    // the breaking detection of the AMR breaking cases turns into differences of up to 1e-3
    const bool half = (lap_kind==1);
    int it = reefamr_bicgstab(sp,p->N44,p->N46,bn,rn,&tm[2],&tm[3],&lap_restarts,STAGTOL*p->N44,&lap_stalls,half);

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
    {
        lap_it_psi_max = MAX(lap_it_psi_max,it);
        lap_it_psi += it;
        ++lap_solves_psi;
    }
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

// A 324 > 0: the composite solve repeated with the rows assembled from the new iterate (the explicit
// sigma cross-derivatives), until the change of the leaf unknowns is below A 325 times that of the
// first solve, at most A 324 times (as fnpf_laplace_cds2 on a uniform grid); between the solves the
// bed ghost nodes and the halo of level 0 and the bed ghost and wall nodes of the patches from the
// iterate (the columns around the patches and the covered ones: lap_sync at the end of lap_core)
template<class ASM>
void fnpf_amr::lap_picard(lexer *p, fdm_fnpf *c, ghostcell *pgc, fnpf_fsf *pf, ASM assemble)
{
    const int nouter = 1 + MAX(p->A324,0);
    double dfirst = 0.0;

    for(int qo=0; qo<nouter; ++qo)
    {
        if(qo>0)
        {
            plap0->picard_bc(p,c,pgc,ltgt[0],pf);
            comms_off guard(pgc);
            for(int id=0; id<(int)P.size(); ++id)
            {
                fnpf_amr_patch *pc = FP(id);
                pc->pbu->bedbc_sig(pc->pp,pc->c,pgc,ltgt[id+1],pc->pf);
                walls_fi(*pc,ltgt[id+1]);
            }
        }

        const double ts = MPI_Wtime();
        assemble();
        lap_tsetup(MPI_Wtime()-ts);

        if(nouter>1)
        {
            pold.resize(P.size()+1);
            for(int g=-1; g<(int)P.size(); ++g)
            {
                lexer *q = glex(g);
                pold[g+1].assign(ltgt[g+1],ltgt[g+1]+(size_t)q->imax*q->jmax*(q->kmax+2));
            }
        }

        lap_core(p,pgc);

        if(nouter>1)
        {
            double dmax = 0.0;
            for(int g=-1; g<(int)P.size(); ++g)
            {
                const vector<int> &lq = lg[g+1].lq;
                const vector<double> &o = pold[g+1];
                for(int qq : lq)
                dmax = MAX(dmax,fabs(ltgt[g+1][qq]-o[qq]));
            }
            dmax = pgc->globalmax(dmax);

            if(qo==0)
            dfirst = dmax;
            else if(dmax<=p->A325*dfirst)
            break;
        }
    }
}

// the Laplace equation of the RK stage on all grids (fnpf_RK3, through fnpf_laplace_amr)
void fnpf_amr::lap_solve(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf, double *f, slice &Fifsf)
{
    const double t0 = MPI_Wtime();

    // breaking (A 350 1): detection on the patch stage values, the viscosity of the cells around
    // the patches from their source grid (level 0 has detected in fnpf_RK3 before this solve,
    // after the stage values: alpha 1, 1/2, 1)
    if(brk)
    patch_breaking(-1,pgc,stg,(stg==1) ? 0.5 : 1.0);

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

    // initial guess: Fi of the previous stage
    ltgt.assign(P.size()+1,nullptr);
    ltgt[0] = f;
    for(int id=0; id<(int)P.size(); ++id)
    ltgt[id+1] = FP(id)->c->Fi;

    // rows: level 0 with its own solver, the patches; A 324 > 0: assembled again with the new
    // iterate (Picard on the explicit sigma cross-derivatives, as fnpf_laplace_cds2)
    lap_kind = 0;
    lap_picard(p,c,pgc,pf,[&]()
    {
        plap0->assemble_only = true;
        plap0->start(p,c,pgc,psolv,pf,f,Fifsf);
        plap0->assemble_only = false;

        comms_off guard(pgc);
        for(int id=0; id<(int)P.size(); ++id)
        {
            fnpf_amr_patch *pc = FP(id);
            pc->plap->start(pc->pp,pc->c,pgc,nullptr,pc->pf,pc->c->Fi,sval(id,stg,1));
        }
    });

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

    ltgt.assign(f,f+P.size()+1);

    lap_kind = 1;
    lap_picard(p,c,pgc,pf,[&]()
    {
        plap0->assemble_only = true;
        plap0->start(p,c,pgc,psolv,pf,f[0],*D[0]);
        plap0->assemble_only = false;

        comms_off guard(pgc);
        for(int id=0; id<(int)P.size(); ++id)
        {
            fnpf_amr_patch *pc = FP(id);
            pc->plap->start(pc->pp,pc->c,pgc,nullptr,pc->pf,f[id+1],*D[id+1]);
        }
    });
    lap_kind = 0;

    tm[1] += MPI_Wtime()-t0;
    tm[6] += MPI_Wtime()-t0;
}

// X 36: a unit-mode psi solve of one body on the patches of its window (out[id] 1: outside).  Level 0 and
// the patches outside hold 0: the columns around the level-1 patches of the window are filled from them
// before the solve (0, a sibling of the window its value) and kept during it (the lowest level of the
// window, lap_edge), level 0 has no rows, and the preconditioner has no level-0 cycle (lap_prec_win with
// lap_nol0: the patch-local V-cycles of the window's level-1 patches carry the coarse correction).  The
// right-hand side and the rows of the window's patches only; the rows and the preconditioner are kept
// between the solves of the window with the same rows (lap_same, wkey with key)
void fnpf_amr::lap_solve_psi_window(lexer *p, ghostcell *pgc, double **f, slice **D, const vector<char> &out, const vector<char> &wgp, int key)
{
    const double t0 = MPI_Wtime();

    lap_skip = out;
    skip = &lap_skip;
    lap_wgp = &wgp;

    // the Krylov vectors of the window's patches (kvw)
    for(int id=0; id<(int)P.size(); ++id)
    if(!out[id] && FP(id)->kvw.empty())
    {
        lexer *q = FP(id)->pp;
        const size_t n = (size_t)q->imax*q->jmax*(q->kmax+2);
        for(int k=0; k<NVEC; ++k)
        FP(id)->kvw.push_back(new double[n]());
    }

    lap_kind = 1;
    {
    comms_off guard(pgc);
    for(int id=0; id<(int)P.size(); ++id)
    if(!lap_skip[id])
    {
        fnpf_amr_patch *pc = FP(id);
        pc->plap->start(pc->pp,pc->c,pgc,nullptr,pc->pf,f[id+1],*D[id+1]);
    }
    }
    lap_tsetup(MPI_Wtime()-t0);

    ltgt.assign(f,f+P.size()+1);

    // psi 0 at the edge of the window
    if(maxlev>=1)
    fill_col(1,7995,[&](int g) -> double* { return ltgt[g+1]; });

    wlo = 1;
    whi = -1;
    lap_edge = true;
    lap_nol0 = true;
    lap_wextra = 100000*(key+1);
    lap_core(p,pgc);
    lap_wextra = 0;
    lap_nol0 = false;
    lap_edge = false;
    lap_kind = 0;
    wlo = 0;
    whi = -1;

    skip = nullptr;
    lap_wgp = nullptr;

    tm[1] += MPI_Wtime()-t0;
    tm[6] += MPI_Wtime()-t0;
}

void fnpf_amr::window_gp(const vector<char> &win, const vector<int> &first, vector<char> &wgp)
{
    wgp.assign(GP.size(),0);
    for(size_t k=0; k<GP.size(); ++k)
    {
        const size_t e = (size_t)first[GP[k].rank] + GP[k].lid;
        wgp[k] = (e<win.size()) ? win[e] : 0;
    }
}

void fnpf_amr::patch_box(int n, double &x0, double &x1, double &y0, double &y1)
{
    lexer *pp = P[n]->pp;
    const int m = marge;
    x0 = pp->XN[EXT+m];
    x1 = pp->XN[EXT+P[n]->nx+m];
    y0 = pp->YN[EY+m];
    y1 = pp->YN[EY+P[n]->ny+m];
}

bool fnpf_amr::zone_box(int nb, double &x0, double &x1, double &y0, double &y1)
{
    if(nb<0 || nb>=(int)zones.size() || zones[nb].bx0>zones[nb].bx1)
    return false;
    x0 = zones[nb].bx0 - par.zr;
    x1 = zones[nb].bx1 + par.zr;
    y0 = zones[nb].by0 - par.zr;
    y1 = zones[nb].by1 + par.zr;
    return true;
}
