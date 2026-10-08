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
#include<limits>

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
//  BiCGStab (reefamr_bicgstab), right preconditioned by one FAC sweep: a REEFMG V-cycle (coefficients
//  in float with N 14 32, the default, as the single-grid REEFMG) on level 0
//  (the whole rank grid, covered columns included), then level by level the coarse correction
//  interpolated into the patches and patch-local REEFMG V-cycles (MPI_COMM_SELF, correction 0 at
//  the patch edge) on the residual it leaves.  The same structure as the composite Laplace of
//  FNPF AMR (fnpf_amr_lap.cpp).

using namespace nhflow_amr_detail;

namespace
{
typedef reefamr_comms_off comms_off;

// N 14 32: the V-cycles of the preconditioner keep their coefficients in float.  reefmg_core then
// wants the exact fine operator for its own Krylov solver and residual(); the FAC sweep only calls
// vcycle(), which smooths with the float coefficients, so this one is never used.
struct pr_fineop : public sc_operator
{
    void fine_apply(const sc_level&, const double*, double*) override
    {
        cout<<"NHFLOW AMR: fine operator of the preconditioner called"<<endl;
        MPI_Abort(MPI_COMM_WORLD,-2761);
    }
};
pr_fineop pr_fineop_none;

// coefficients of row lq of the fine level, in the storage precision of the multigrid
inline void pr_coef(reefmg_core *mg, sc_level &L, long lq, double p, double n, double s, double w, double e, double t, double b)
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

        // G 7 1: the grids of the level window only (wlo..wtop: all levels without subcycling)
        if(l<wlo || l>wtop())
        continue;

        int i0,i1,j0,j1;
        if(g<0) { i0=0; i1=NX0-1; j0=0; j1=NY0-1; }
        else    { i0=EXT; i1=EXT+P[g]->nx-1; j0=EXT; j1=EXT+P[g]->ny-1; }

        for(int ii=i0; ii<=i1; ++ii)
        for(int jj=j0; jj<=j1; ++jj)
        {
            int I = (g<0) ? ii+O0i : ii-EXT+P[g]->I0;
            int J = (g<0) ? jj+O0j : jj-EXT+P[g]->J0;
            const bool cov = (l<wtop() && covered(l+1,2*I,2*J));

            for(int kk=0; kk<q->knoz; ++kk)
            {
                const int qq = fidx(q,ii,jj,kk);
                const int r = row[qq];
                if(r<0)
                continue;

                if(cov)
                L.cq.push_back(qq);
                else
                {
                    L.aq.push_back(qq);
                    L.ar.push_back(r);
                    L.aw.push_back((ii-q->imin)*q->jmax + (jj-q->jmin));
                }
            }
        }
    }
}

// vectors, multigrids and their coefficients for this solve
void nhflow_amr::pr_prepare(ghostcell *pgc)
{
    // row lists of the layout and of the level window (G 7 1: level solves and synchronisation
    // projections over part of the levels)
    const int wkey = wlo*1000 + wtop();
    if(pr_layout!=layout_id || pr_wkey!=wkey)
    {
        pr_rows();
        pr_wkey = wkey;
    }

    if(pr_layout!=layout_id)
    {

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
                c->mg->set_precision(pp->N14==32 ? 32 : 64);
                c->mg->set_fine_operator(&pr_fineop_none);
                if(!c->mg->setup(MPI_COMM_SELF,c->nx,c->ny,pp->knoz,c->nx,c->ny,0,pp->DXN+marge+EXT-1,pp->DYN+marge+EXT-1))
                cout<<"NHFLOW AMR: patch multigrid "<<c->mg->err()<<endl;
            }
        }

        if(mg0==nullptr)
        {
            mg0 = new reefmg_core;
            mg0->set_precision(p0->N14==32 ? 32 : 64);
            mg0->set_fine_operator(&pr_fineop_none);
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
        lexer *q = glex(g);
        L.lq.clear();
        L.lr.clear();
        L.fq.clear();
        for(size_t n=0; n<L.aq.size(); ++n)
        if(!fixed(M,L.ar[n]))
        {
            L.lq.push_back(L.aq[n]);
            L.lr.push_back(L.ar[n]);
        }
        else if(q->wet[L.aw[n]]==1)
        L.fq.push_back(L.aq[n]);
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

    // level 0: all rows of the rank grid (G 7 1, a level solve: the operator of the coarse
    // correction, from the matrix of the last level-0 solve)
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
                pr_coef(mg0,L,lq,1.0,0.0,0.0,0.0,0.0,0.0,0.0);
                continue;
            }
            pr_coef(mg0,L,lq,M.p[r],M.n[r],M.s[r],M.w[r],M.e[r],M.t[r],M.b[r]);
            L.act[lq] = fixed(M,r) ? 0 : 1;
        }
        mg0->coarsen();
    }

    // patches: interior rows, couplings to the cells around the patch dropped
    for(auto q : P)
    {
        if(q->lev<wlo || q->lev>wtop())
        continue;

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
                pr_coef(c->mg,L,lq,1.0,0.0,0.0,0.0,0.0,0.0,0.0);
                continue;
            }
            pr_coef(c->mg,L,lq,M.p[r],
                    (ii+1<EXT+c->nx) ? M.n[r] : 0.0,
                    (ii-1>=EXT) ? M.s[r] : 0.0,
                    (jj+1<EXT+c->ny) ? M.w[r] : 0.0,
                    (jj-1>=EXT) ? M.e[r] : 0.0,
                    M.t[r],M.b[r]);
            L.act[lq] = fixed(M,r) ? 0 : 1;
        }
        c->mg->coarsen();
    }

    pr_stencils();
}

// pcol with the weights of a solve kept: the source offsets and weights of the 5x5 stencil
void nhflow_amr::pst_make(int g, int ic, int jc, int ox, int oy, pstencil &st)
{
    lexer *q = glex(g);
    const int sI = q->jmax*q->kmaxF;
    const int sJ = q->kmaxF;

    double w[25];
    pweights(q,ic,jc,ox,oy,w,2);

    st.nw = 0;
    const int n0 = fidx(q,ic,jc,0);
    for(int di=-2; di<=2; ++di)
    for(int dj=-2; dj<=2; ++dj)
    {
        const double c = w[(di+2)*5+(dj+2)];
        if(c==0.0)
        continue;
        st.off[st.nw] = n0 + di*sI + dj*sJ;
        st.w[st.nw] = c;
        ++st.nw;
    }
}

void nhflow_amr::pst_col(const pstencil &st, const double *src, int kc, int knf, double *v) const
{
    const int fz = knf/kc;

    for(int kk=0; kk<=kc; ++kk)
    {
        double r = 0.0;
        for(int m=0; m<st.nw; ++m)
        r += st.w[m]*src[st.off[m]+kk];
        v[fz*kk] = r;
    }

    if(fz==2)
    for(int kk=0; kk<kc; ++kk)
    v[2*kk+1] = 0.5*(v[2*kk]+v[2*kk+2]);
}

// the restriction blocks, interior prolongation stencils and active index lists of this solve
void nhflow_amr::pr_stencils()
{
    // restriction: the blocks of the local patches (their parent on any rank)
    pr_rbk.assign(P.size(),vector<prblock>());
    for(int id=0; id<(int)P.size(); ++id)
    {
        reefamr_patch *c = P[id];
        lexer *pp = c->pp;
        const int nby = c->ny/2;
        pr_rbk[id].resize(c->rgrid.size());

        for(int bi=0; bi<c->nx/2; ++bi)
        for(int bj=0; bj<nby; ++bj)
        {
            const int k = bi*nby+bj;
            if(c->rgrid[k]==-2)
            continue;

            const int i0 = EXT+2*bi, j0 = EXT+2*bj;

            bool hi = (bi>0 && bi<c->nx/2-1 && bj>0 && bj<nby-1) && rcubic(pp,i0,j0) && sub==0;

            int na = 4;
            double wa[4] = {0.25,0.25,0.25,0.25};
            if(shore)
            {
                for(int a=-1; a<=2 && hi; ++a)
                for(int b=-1; b<=2 && hi; ++b)
                if(!wet_at(pp,i0+a,j0+b,2))
                hi = false;

                na = 0;
                for(int a=0; a<2; ++a)
                for(int b=0; b<2; ++b)
                {
                    wa[2*a+b] = wet_at(pp,i0+a,j0+b,2) ? 1.0 : 0.0;
                    na += (int)wa[2*a+b];
                }
                for(int m=0; m<4; ++m)
                wa[m] = (na>0) ? wa[m]/double(na) : 0.0;
            }

            prblock &B = pr_rbk[id][k];
            B.mode = hi ? 1 : (na==4 ? 2 : 3);
            B.src = fidx(pp,i0,j0,0);
            for(int m=0; m<4; ++m)
            B.wa[m] = wa[m];
        }
    }

    // interior prolongation: the stencils of the parents this rank holds, made at first use
    pr_pik.assign(maxlev+1,vector<pstencil>());
    pr_pik_ok.assign(maxlev+1,vector<char>());
    for(int l=1; l<=maxlev; ++l)
    {
        pr_pik[l].resize(4*(size_t)block_keys(l));
        pr_pik_ok[l].assign(4*(size_t)block_keys(l),0);
    }

    pr_fs.clear();

    // the level-0 rows of the V-cycle (G 7 1: also for the coarse correction of a level solve)
    {
        const sc_level &L = mg0->fine();
        pr_l0a.clear();
        pr_l0f.clear();
        pr_l0z.clear();
        for(int ii=0; ii<p0->knox; ++ii)
        for(int jj=0; jj<p0->knoy; ++jj)
        for(int kk=0; kk<p0->knoz; ++kk)
        {
            const long lq = L.idx(ii,jj,kk);
            if(L.act[lq])
            {
                pr_l0a.push_back(lq);
                pr_l0f.push_back(fidx(p0,ii,jj,kk));
            }
            else
            pr_l0z.push_back(fidx(p0,ii,jj,kk));
        }
    }

    pr_pa.assign(P.size(),vector<pact>());
    for(int id=0; id<(int)P.size(); ++id)
    {
        nhflow_amr_patch *c = NP(id);
        lexer *pp = c->pp;
        const sc_level &L = c->mg->fine();
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
        for(int kk=0; kk<pp->knoz; ++kk)
        {
            const long lq = L.idx(ii-EXT,jj-EXT,kk);
            if(L.act[lq]==0)
            continue;
            const int qq = fidx(pp,ii,jj,kk);
            pr_pa[id].push_back({lq,qq,c->row[qq]});
        }
    }
}

// restrict_col with the blocks of pr_stencils
template<class SEL>
void nhflow_amr::pr_restrict(SEL sel)
{
    using namespace nhflow_amr_detail;

    for(int l=wtop(); l>=wlo+1; --l)
    {
        const int Kc = klev(l-1);
        block_up(l,Kc+1,7500+l,
                 [&](reefamr_patch *c, int id, int k, double *v)
                 {
                     const prblock &B = pr_rbk[id][k];
                     lexer *pp = c->pp;
                     const double *src = sel(id);
                     const int sI = pp->jmax*pp->kmaxF, sJ = pp->kmaxF;
                     const int fz = pp->knoz/Kc;

                     for(int K=0; K<=Kc; ++K)
                     {
                         const int n0 = B.src + fz*K;
                         if(B.mode==1)
                         v[K] = rc4([&](int a, int b) { return src[n0+a*sI+b*sJ]; });
                         else if(B.mode==2)
                         v[K] = 0.25*(src[n0] + src[n0+sI] + src[n0+sJ] + src[n0+sI+sJ]);
                         else
                         v[K] = B.wa[0]*src[n0] + B.wa[2]*src[n0+sI] + B.wa[1]*src[n0+sJ] + B.wa[3]*src[n0+sI+sJ];
                     }
                 },
                 [&](const reefamr_block &B, int key, const double *v)
                 {
                     double *dst = sel(B.g);
                     const int d0 = fidx(glex(B.g),B.ic,B.jc,0);
                     for(int K=0; K<=Kc; ++K)
                     dst[d0+K] = v[K];
                 });
    }
}

// fill_col with the stencils of the fills from the coarser grid kept for the solve
template<class SEL>
void nhflow_amr::pr_fill(int l, int tag, SEL sel)
{
    const int knf = klev(l);
    const int nv = knf+1;

    // G 7 1, the lowest level of the window above level 0: its parent columns are fixed (the
    // parent pressure of the level solve, 0 for a correction); they keep their values in the
    // unknowns and are 0 in the Krylov vectors
    const bool edge = (pr_edge && l==wlo);

    fill_run(l,nv,tag,
             [&](const reefamr_fill &f, double *v)
             {
                 if(f.kind==0)
                 {
                     lexer *q = glex(f.g);
                     const double *src = sel(f.g);
                     for(int kk=0; kk<=knf; ++kk)
                     v[kk] = src[fidx(q,f.si,f.sj,kk)];
                 }
                 else if(f.kind==1 && edge)
                 {
                     for(int kk=0; kk<=knf; ++kk)
                     v[kk] = std::numeric_limits<double>::quiet_NaN();
                 }
                 else if(f.kind==1)
                 {
                     auto it = pr_fs.find(&f);
                     if(it==pr_fs.end())
                     {
                         pstencil st;
                         pst_make(f.g,f.si,f.sj,f.ox,f.oy,st);
                         it = pr_fs.emplace(&f,st).first;
                     }
                     pst_col(it->second,sel(f.g),glex(f.g)->knoz,knf,v);
                 }
                 else
                 for(int kk=0; kk<=knf; ++kk)
                 v[kk] = 0.0;
             },
             [&](reefamr_patch *c, int id, const reefamr_fill &f, const double *w)
             {
                 lexer *pp = c->pp;
                 double *dst = sel(id);
                 if(edge && std::isnan(w[0]))
                 {
                     if(pr_dir)
                     for(int kk=0; kk<=knf; ++kk)
                     dst[fidx(pp,f.di,f.dj,kk)] = 0.0;
                     return;
                 }
                 for(int kk=0; kk<=knf; ++kk)
                 dst[fidx(pp,f.di,f.dj,kk)] = w[kk];
             });
}

// the coarse correction of vector k into the interior of the level-l patches: from the parent
// of every 2x2 block (on its rank) with the stencils of the solve
void nhflow_amr::pr_prolong(int l, int k)
{
    const int knf = klev(l);
    const int nc = knf+1;

    block_down(l,4*nc,7520+l,
               [&](const reefamr_block &B, int key, double *v)
               {
                   lexer *q = glex(B.g);
                   for(int a=0; a<2; ++a)
                   for(int d=0; d<2; ++d)
                   {
                       const size_t m = 4*(size_t)key + 2*a+d;
                       if(!pr_pik_ok[l][m])
                       {
                           pst_make(B.g,B.ic,B.jc,a==0?-1:1,d==0?-1:1,pr_pik[l][m]);
                           pr_pik_ok[l][m] = 1;
                       }
                       pst_col(pr_pik[l][m],pvec(B.g,k),q->knoz,knf,&v[(2*a+d)*nc]);
                   }
               },
               [&](reefamr_patch *c, int id, int kb, const double *v)
               {
                   lexer *pp = c->pp;
                   const int nby = c->ny/2;
                   const int bi = kb/nby, bj = kb%nby;
                   double *dst = pvec(id,k);
                   for(int a=0; a<2; ++a)
                   for(int d=0; d<2; ++d)
                   {
                       const int d0 = fidx(pp,EXT+2*bi+a,EXT+2*bj+d,0);
                       for(int kk=0; kk<=knf; ++kk)
                       dst[d0+kk] = v[(2*a+d)*nc+kk];
                   }
               });
}

// covered columns, partition halo of level 0, columns around the patches
void nhflow_amr::pr_sync(int k)
{
    auto sel = [&](int g) -> double* { return pvec(g,k); };

    pr_dir = (k>=0);

    pr_restrict(sel);

    if(p0->mpi_size>1 && wlo==0)
    {
        double *x0 = pvec(-1,k);
        pgc0->gcparax7(p0,x0,7);
        pgc0->gcparax7co(p0,x0,7);
    }

    for(int l=MAX(wlo,1); l<=wtop(); ++l)
    pr_fill(l,7600+l,sel);
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
    pr_restrict([&](int g) -> double* { return pvec(g,kr); });
    pr_dir = true;

    // level 0: one V-cycle on the whole rank grid
    if(wlo==0)
    {
        sc_level &L = mg0->fine();
        const double *r = pvec(-1,kr);
        double *z = pvec(-1,kz);

        std::fill(L.u.begin(),L.u.end(),0.0);
        std::fill(L.f.begin(),L.f.end(),0.0);
        for(size_t n=0; n<pr_l0a.size(); ++n)
        L.f[pr_l0a[n]] = r[pr_l0f[n]];

        mg0->vcycle(0,1,1);

        for(size_t n=0; n<pr_l0a.size(); ++n)
        z[pr_l0f[n]] = L.u[pr_l0a[n]];
        for(long qq : pr_l0z)
        z[qq] = 0.0;

        if(p0->mpi_size>1)
        {
            pgc0->gcparax7(p0,z,7);
            pgc0->gcparax7co(p0,z,7);
        }
    }

    auto selz = [&](int g) -> double* { return pvec(g,kz); };

    // G 7 1, the level solve of level wlo > 0 with a coarse correction: its residual restricted
    // through the levels below down to level 0 (in the vectors of the grids below the window,
    // which have no rows in this solve: scratch, 0 outside the covered cells), one level-0
    // V-cycle, the correction interpolated back up into the patch interiors; then the patch-local
    // passes below.  The patch-local V-cycles alone do not reduce the error that spans several
    // patches or a large patch (the ring wave of 0019: 11 BiCGStab iterations against 5 of the
    // composite solve with its level-0 cycle)
    if(wlo>0)
    {
        const int wl = wlo;

        for(int g=-1; g<(int)P.size(); ++g)
        {
            if(g>=0 && P[g]->lev>=wl)
            continue;
            lexer *q = glex(g);
            const size_t n = (size_t)q->imax*q->jmax*(q->kmax+2);
            std::fill(pvec(g,kr),pvec(g,kr)+n,0.0);
            std::fill(pvec(g,kz),pvec(g,kz)+n,0.0);
        }

        wlo = 0;
        pr_restrict([&](int g) -> double* { return pvec(g,kr); });

        sc_level &L = mg0->fine();
        const double *r = pvec(-1,kr);
        double *z = pvec(-1,kz);
        std::fill(L.u.begin(),L.u.end(),0.0);
        std::fill(L.f.begin(),L.f.end(),0.0);
        for(size_t n=0; n<pr_l0a.size(); ++n)
        L.f[pr_l0a[n]] = r[pr_l0f[n]];

        mg0->vcycle(0,1,1);

        for(size_t n=0; n<pr_l0a.size(); ++n)
        z[pr_l0f[n]] = L.u[pr_l0a[n]];
        for(long qq : pr_l0z)
        z[qq] = 0.0;

        if(p0->mpi_size>1)
        {
            pgc0->gcparax7(p0,z,7);
            pgc0->gcparax7co(p0,z,7);
        }

        for(int l=1; l<=wl; ++l)
        {
            pr_prolong(l,kz);
            if(l<wl)
            pr_fill(l,7800+l,selz);
        }

        wlo = wl;
    }

    for(int l=MAX(wlo,1); l<=wtop(); ++l)
    {
        // coarse correction interpolated into the patch interiors, then the columns around them
        if(l>wlo)
        pr_prolong(l,kz);

        // G 7 1, the level solve of a level with several patches: the patch-local corrections
        // are repeated once with the siblings' corrections filled in between (with one pass every
        // patch corrects against zero at its sibling edges: the ring wave with 3-7 patches needs
        // 21 BiCGStab iterations with one pass, 11 with two, 7.7 with four, for the same time)
        const int npass = (l==wlo && wlo>0 && nlevg[l]>1) ? 2 : 1;

        for(int pass=0; pass<npass; ++pass)
        {
        pr_fill(l,7700+l,selz);

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
            for(const pact &A : pr_pa[id])
            {
                const int qq = A.qq;
                const int rw = A.rw;
                const double az = M.p[rw]*z[qq] + M.n[rw]*z[qq+sI] + M.s[rw]*z[qq-sI] + M.w[rw]*z[qq+sJ] + M.e[rw]*z[qq-sJ]
                                + M.t[rw]*z[qq+1] + M.b[rw]*z[qq-1];
                L.f[A.lq] = r[qq] - az;
            }

            c->mg->vcycle(0,1,1);

            for(const pact &A : pr_pa[id])
            z[A.qq] += L.u[A.lq];
        }
        }

        if(l<wtop())
        pr_fill(l,7700+l,selz);
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

    // identity rows of wet but shallow columns: P = 0 as the single-grid solvers leave them
    // (PCORR is 0 there already); dry columns keep their pressure
    for(int g=-1; g<(int)P.size(); ++g)
    {
        double *x = pvec(g,-1);
        for(int qq : pg[g+1].fq)
        x[qq] = 0.0;
    }

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
