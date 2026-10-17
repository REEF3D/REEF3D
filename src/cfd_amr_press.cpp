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
#include<map>
#include<tuple>
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"heaviside.h"
#include"interface_width.h"
#include"momentum_rk.h"
#include"pjm_corr.h"
#include"ioflow.h"
#include"reefmg_core.h"
#include"reefamr_krylov.h"
#include<cmath>
#include<iostream>
#include<iomanip>
#include<algorithm>
#include<limits>
#include"definitions.h"

//  The composite pressure projection of REEF3D::CFD with mesh refinement.
//
//  Every grid assembles its Poisson rows with its own objects (pjm_corr::amr_prepare: boundary
//  velocities, right-hand side -div u / (alpha dt), poisson_pcorr).  The unknowns are the leaf
//  cells: the level-0 cells and patch cells not under a finer patch.  A face between a patch cell f
//  and a cell C of the coarser level (the patch edge) is one face of the composite grid:
//
//    - before the solve, its velocity is the coarse face velocity of the stage, prolonged linearly
//      in the transverse directions (the mean of the fine faces is the coarse value), so the
//      right-hand sides of both sides see the same flux
//    - the flux of the pressure correction through the fine face is (x_C - x_f) / (rho d) with d
//      the distance of the cell centres: the fine row takes it with its own cell width, the coarse
//      row C the sum over the fine faces of its face (area weighted) instead of its coupling to the
//      covered cell; volume weighted the matrix is symmetric
//    - after the solve the fine faces are corrected with that flux, the coarse face is the mean of
//      the fine faces (restriction): the velocities of all leaf cells are divergence free together
//
//  BiCGStab (reefamr_bicgstab), right preconditioned by one FAC sweep: the residual restricted
//  onto level 0, a REEFMG V-cycle on level 0 (the whole rank grid with its own Poisson matrix),
//  then level by level the correction prolonged into the patches and a patch-local REEFMG V-cycle
//  on the residual it leaves (zero correction around the patch), as the composite solves of
//  NHFLOW and FNPF with mesh refinement.

namespace
{
inline int fdiv(int a, int r) { return (a>=0) ? a/r : -((-a+r-1)/r); }

inline int cix(const lexer *q, int i, int j, int k)
{
    return (i-q->imin)*q->jmax*q->kmax + (j-q->jmin)*q->kmax + k-q->kmin;
}

struct pr_fineop : public sc_operator
{
    void fine_apply(const sc_level&, const double*, double*) override
    {
        cout<<"CFD AMR: fine operator of the preconditioner called"<<endl;
        MPI_Abort(MPI_COMM_WORLD,-3121);
    }
};
pr_fineop pr_fineop_none;

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

void pr_clear(sc_level &L)
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
}

// the coefficient of row n towards the neighbour at side dir (+1, -1) of direction d
inline double& mcoef(matrix_diag &M, int n, int d, int dir)
{
    if(d==0) return dir>0 ? M.n[n] : M.s[n];
    if(d==1) return dir>0 ? M.w[n] : M.e[n];
    return dir>0 ? M.t[n] : M.b[n];
}

// the coefficient of multigrid row q towards its neighbour in direction d (dir > 0: up)
template<class T>
inline T& mcoef_mg(std::vector<T> &n, std::vector<T> &s, std::vector<T> &w, std::vector<T> &e, std::vector<T> &t, std::vector<T> &b, long q, int d, int dir)
{
    if(d==0) return dir>0 ? n[q] : s[q];
    if(d==1) return dir>0 ? w[q] : e[q];
    return dir>0 ? t[q] : b[q];
}

// G 40: one row of ks_apply on its rank (the other cell on another rank): h 1 the coarse cell c (its
// coupling towards -dir), h 2 the covered cell c (towards dir)
void ks_half(reefmg_core *mg, const int *c, int d, int dir, int h, const int *n, double K)
{
    for(int e=0; e<3; ++e)
    if(c[e]<0 || c[e]>=n[e])
    return;
    sc_level &L = mg->fine();
    const long lq = L.idx(c[0],c[1],c[2]);
    if(!L.act[lq])
    return;
    const int s = (h==1) ? -dir : dir;
    if(mg->precision()==32)
    {
        float &a = mcoef_mg(L.nf,L.sf,L.wf,L.ef,L.tf,L.bf,lq,d,s);
        if(a==0.0f)
        return;
        L.pf[lq] += a + (float)K;
        a = -(float)K;
    }
    else
    {
        double &a = mcoef_mg(L.n,L.s,L.w,L.e,L.t,L.b,lq,d,s);
        if(a==0.0)
        return;
        L.p[lq] += a + K;
        a = -K;
    }
}

// the couplings of a coarse cell ci of multigrid mg next to a patch to the covered cell across the
// face (side dir of direction d) from the fine faces: the sum K of their kc (n: the cells of the grid)
void ks_apply(reefmg_core *mg, const int *ci, int d, int dir, const int *n, double K)
{
    int pi[3] = {ci[0],ci[1],ci[2]};
    pi[d] -= dir;
    if(pi[d]<0 || pi[d]>=n[d] || ci[d]<0 || ci[d]>=n[d])
    return;

    sc_level &L = mg->fine();
    const long lc = L.idx(ci[0],ci[1],ci[2]), lp = L.idx(pi[0],pi[1],pi[2]);
    if(!L.act[lc] || !L.act[lp])
    return;
    if(mg->precision()==32)
    {
        float &c = mcoef_mg(L.nf,L.sf,L.wf,L.ef,L.tf,L.bf,lc,d,-dir);
        float &b = mcoef_mg(L.nf,L.sf,L.wf,L.ef,L.tf,L.bf,lp,d,dir);
        if(c==0.0f || b==0.0f)
        return;
        L.pf[lc] += c + (float)K;
        L.pf[lp] += b + (float)K;
        c = -(float)K;
        b = -(float)K;
    }
    else
    {
        double &c = mcoef_mg(L.n,L.s,L.w,L.e,L.t,L.b,lc,d,-dir);
        double &b = mcoef_mg(L.n,L.s,L.w,L.e,L.t,L.b,lp,d,dir);
        if(c==0.0 || b==0.0)
        return;
        L.p[lc] += c + K;
        L.p[lp] += b + K;
        c = -K;
        b = -K;
    }
}

enum { NR=REEFAMR_NR, NRH=REEFAMR_NRH, NPV=REEFAMR_NPV, NVV=REEFAMR_NVV, NS=REEFAMR_NS, NT=REEFAMR_NT,
       NPH=REEFAMR_NPH, NSH=REEFAMR_NSH, NVEC=REEFAMR_NVEC };
}

field& cfd_amr::pfield(int g, int k)
{
    if(k<0)
    return (g<0) ? press0->pcorr : CP(g)->ppress->pcorr;

    return (g<0) ? *kv0[k] : *CP(g)->kv[k];
}

double* cfd_amr::pvec(int g, int k)
{
    return pfield(g,k).V;
}

// rows, leaf cells, the faces at the patch edges, vectors and multigrids
void cfd_amr::pr_setup()
{
    auto number = [&](lexer *q, vector<int> &row)
    {
        row.assign(q->imax*q->jmax*q->kmax,-1);
        int n=0;
        for(int i=0; i<q->knox; ++i)
        for(int j=0; j<q->knoy; ++j)
        for(int k=0; k<q->knoz; ++k)
        if(q->flag4[cix(q,i,j,k)]>0)
        row[cix(q,i,j,k)] = n++;
    };

    number(p0,row0);
    for(int id=0; id<(int)P.size(); ++id)
    number(CP(id)->pp,CP(id)->row);

    lq.assign(P.size()+1,vector<int>());
    lr.assign(P.size()+1,vector<int>());
    cq.assign(P.size()+1,vector<int>());
    aq.assign(P.size()+1,vector<int>());
    ar.assign(P.size()+1,vector<int>());
    prl = -1;
    prg.clear();
    for(int g=-1; g<npat; ++g)
    prg.push_back(g);

    for(int g=-1; g<npat; ++g)
    {
        lexer *q = glex(g);
        const vector<int> &row = (g<0) ? row0 : CP(g)->row;
        const int l = glevel(g);
        int og[3];
        goff(g,og);

        for(int i=0; i<q->knox; ++i)
        for(int j=0; j<q->knoy; ++j)
        for(int k=0; k<q->knoz; ++k)
        {
            const int qq = cix(q,i,j,k);
            if(row[qq]<0)
            continue;
            const int I[3] = {i+og[0],j+og[1],k+og[2]};
            aq[g+1].push_back(qq);
            ar[g+1].push_back(row[qq]);
            if(covered(l,I))
            cq[g+1].push_back(qq);
            else
            {
                lq[g+1].push_back(qq);
                lr[g+1].push_back(row[qq]);
            }
        }
    }

    // the local grid of level lv whose interior holds the global cell I (-1: level 0; G 40: the
    // shadow sh of the level-1 patch)
    int sh = -1;
    auto grid_of = [&](int lv, const int *I, int *loc) -> int
    {
        if(lv==0 && place)
        {
            const r3patch *c = P[sh];
            for(int d=0; d<3; ++d)
            {
                loc[d] = I[d]-c->lo[d];
                if(loc[d]<0 || loc[d]>=c->n[d])
                return -2;
            }
            return sh;
        }
        if(lv==0)
        {
            for(int d=0; d<3; ++d)
            loc[d] = I[d]-org[d];
            return -1;
        }
        for(int g : lev[lv])
        {
            r3patch *c = P[g];
            if(I[0]>=c->lo[0] && I[0]<=c->hi[0] && I[1]>=c->lo[1] && I[1]<=c->hi[1] && I[2]>=c->lo[2] && I[2]<=c->hi[2])
            {
                for(int d=0; d<3; ++d)
                loc[d] = I[d]-c->lo[d];
                return g;
            }
        }
        return -2;
    };

    const double *cc[3], *dn[3];
    cf.clear();
    int bad = 0;

    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        const int l = c->lev;
        sh = c->sh;

        for(int i=0; i<pp->knox; ++i)
        for(int j=0; j<pp->knoy; ++j)
        for(int k=0; k<pp->knoz; ++k)
        {
            const int qf = cix(pp,i,j,k);
            if(c->row[qf]<0)
            continue;

            const int fi[3] = {i,j,k};
            for(int d=0; d<3; ++d)
            {
                if(rr[d]!=2)
                continue;
                for(int dir=-1; dir<=1; dir+=2)
                {
                    int mi[3] = {i,j,k};
                    mi[d] += dir;
                    const int qm = cix(pp,mi[0],mi[1],mi[2]);
                    if(c->mkind[qm]!=1 || pp->flag4[qm]<=0)
                    continue;

                    cfd_amr_cf F;
                    F.pid = id;
                    F.d = d;
                    F.dir = dir;
                    F.qf = qf;
                    F.nf = c->row[qf];
                    F.qm = qm;

                    int I[3], Pg[3], Cg[3];
                    for(int e=0; e<3; ++e)
                    {
                        I[e] = c->lo[e]+fi[e];
                        Pg[e] = fdiv(I[e],rr[e]);
                        F.o[e] = I[e]-Pg[e]*rr[e];
                        Cg[e] = Pg[e];
                        F.ff[e] = fi[e];
                    }
                    Cg[d] += dir;
                    if(dir<0)
                    F.ff[d] -= 1;

                    int pl[3], cl[3];
                    F.gp = grid_of(l-1,Pg,pl);
                    F.gc = grid_of(l-1,Cg,cl);
                    if(F.gp<-1 || F.gc<-1)
                    {
                        if(bad<4)
                        {
                            const int gq = gpatch_at(l-1,Cg), gpp = gpatch_at(l-1,Pg);
                            cout<<"CFD AMR: rank "<<myrank<<": the face "<<d<<" "<<dir<<" of cell ("<<I[0]<<","<<I[1]<<","<<I[2]<<") of the level-"<<l
                                <<" patch "<<c->gid<<" (unit "<<GP[c->gid].unit<<"): its coarse cell ("<<Cg[0]<<","<<Cg[1]<<","<<Cg[2]<<") in patch "<<gq
                                <<(gq>=0 ? " of unit " : "")<<(gq>=0 ? GP[gq].unit : -1)<<" on rank "<<(gq>=0 ? GP[gq].rank : -1)<<", its parent in "<<gpp<<endl;
                        }
                        ++bad;
                        continue;
                    }
                    lexer *qp = glex(F.gp);
                    lexer *qc = glex(F.gc);
                    F.qp = cix(qp,pl[0],pl[1],pl[2]);
                    F.qc = cix(qc,cl[0],cl[1],cl[2]);
                    const vector<int> &rowc = (F.gc<0) ? row0 : CP(F.gc)->row;     // (a shadow: the cells with fluid)
                    F.nc = rowc[F.qc];
                    if(F.nc<0)
                    {
                        ++bad;
                        continue;
                    }

                    for(int e=0; e<3; ++e)
                    {
                        F.fi[e] = fi[e];
                        F.ci[e] = cl[e];
                        F.pi[e] = pl[e];
                    }

                    F.gface = (dir>0) ? F.gp : F.gc;
                    for(int e=0; e<3; ++e)
                    F.cf[e] = (dir>0) ? pl[e] : cl[e];

                    cc[0]=pp->XP; cc[1]=pp->YP; cc[2]=pp->ZP;
                    dn[0]=pp->DXN; dn[1]=pp->DYN; dn[2]=pp->DZN;
                    const double xf = cc[d][fi[d]+marge];
                    F.dxnf = dn[d][fi[d]+marge];
                    cc[0]=qc->XP; cc[1]=qc->YP; cc[2]=qc->ZP;
                    dn[0]=qc->DXN; dn[1]=qc->DYN; dn[2]=qc->DZN;
                    const double xc = cc[d][cl[d]+marge];
                    F.dxnc = dn[d][cl[d]+marge];
                    F.dist = fabs(xc-xf);

                    F.area = 1.0;
                    for(int e=0; e<3; ++e)
                    if(e!=d)
                    F.area /= double(rr[e]);

                    F.ustar = F.kf = F.kc = 0.0;
                    cf.push_back(F);
                }
            }
        }
    }

    pr_stencils();
    pr_components();

    bad = pgc0->globalimax(bad);
    if(bad>0)
    {
        if(myrank==0 && !place)
        cout<<"CFD AMR: a patch edge lies on a partition edge of level 0 (the coarse cell across it is on another rank); "
              "move the refinement box (G 15) or change the partition"<<endl;
        if(myrank==0 && place)
        cout<<"CFD AMR (G 40): the coarse cell across a patch edge is not on the rank of the patch (the units)"<<endl;
        MPI_Abort(MPI_COMM_WORLD,-3120);
    }

    // vectors (a regrid keeps those of the kept patches)
    if(kv0.empty())
    for(int k=0; k<NVEC+1; ++k)
    kv0.push_back(new field4(p0));
    for(int id=0; id<(int)P.size(); ++id)
    if(CP(id)->kv.empty())
    for(int k=0; k<NVEC+1; ++k)
    CP(id)->kv.push_back(new field4(CP(id)->pp));

    // multigrids
    if(mg0==nullptr)
    {
    mg0 = new reefmg_core;
    mg0->set_precision(p0->N14==32 ? 32 : 64);
    mg0->set_fine_operator(&pr_fineop_none);
    if(!mg0->setup(pgc0->cart(),p0->knox,p0->knoy,p0->knoz,p0->gknox,p0->gknoy,p0->N13,p0->DXN+marge-1,p0->DYN+marge-1,
                   p0->periodic1>0,p0->periodic2>0 && p0->j_dir==1))
    {
        if(p0->mpirank==0)
        cout<<"CFD AMR: level-0 multigrid "<<mg0->err()<<endl;
        MPI_Abort(MPI_COMM_WORLD,-3122);
    }
    }

    // the joint multigrids of the levels (the patches of a rank next to each other), then the overlap
    // of the patch multigrids (G 40; periodic: the patches not in a joint multigrid)
    pr_level_mg_setup();
    pr_overlap();

    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        bool same = true;
        for(int d=0; d<3; ++d)
        if(c->mlo[d]!=c->olo[d] || c->mhi[d]!=c->ohi[d])
        same = false;
        if(c->mg!=nullptr && same)
        continue;
        delete c->mg;
        for(int d=0; d<3; ++d)
        {
            c->mlo[d] = c->olo[d];
            c->mhi[d] = c->ohi[d];
        }
        const int nx = pp->knox+c->olo[0]+c->ohi[0], ny = pp->knoy+c->olo[1]+c->ohi[1], nz = pp->knoz+c->olo[2]+c->ohi[2];
        // the cell widths of the box with one cell beyond it (the array of the patch: its margin,
        // the last width repeated beyond)
        auto widths = [&](const double *DN, int n, int lo, int nb) -> vector<double>
        {
            vector<double> w(nb+2);
            for(int b=0; b<nb+2; ++b)
            {
                const int i = MIN(MAX(b-1-lo,-marge),n+marge-1);
                w[b] = DN[i+marge];
            }
            return w;
        };
        const vector<double> wx = widths(pp->DXN,pp->knox,c->olo[0],nx), wy = widths(pp->DYN,pp->knoy,c->olo[1],ny);
        c->mg = new reefmg_core;
        c->mg->set_precision(pp->N14==32 ? 32 : 64);
        c->mg->set_fine_operator(&pr_fineop_none);
        // (a patch over the whole length of a periodic direction: periodic there)
        if(!c->mg->setup(MPI_COMM_SELF,nx,ny,nz,nx,ny,0,wx.data(),wy.data(),pr_wrap(id,0),pr_wrap(id,1)))
        cout<<"CFD AMR: patch multigrid "<<c->mg->err()<<endl;
    }
    if(pr_ras())
    for(int id=0; id<npat; ++id)
    if(CP(id)->kc7.empty())
    for(int m=0; m<7; ++m)
    CP(id)->kc7.push_back(new field4(CP(id)->pp));

    l0a.clear();
    l0f.clear();
    {
        const sc_level &L = mg0->fine();
        for(int i=0; i<p0->knox; ++i)
        for(int j=0; j<p0->knoy; ++j)
        for(int k=0; k<p0->knoz; ++k)
        if(row0[cix(p0,i,j,k)]>=0)
        {
            l0a.push_back(L.idx(i,j,k));
            l0f.push_back(cix(p0,i,j,k));
        }
    }
}

// The joint multigrids of the levels: the patches of a level on this rank that touch each other
// (a component) get one multigrid over the bounding box of the component (the cells of the level
// index space; the cells outside the patches identity rows), so the couplings between the patches
// stay in the preconditioner.  With a patch-local V-cycle on each patch the preconditioner was a
// block Jacobi over the patches: many small patches of an adaptive level took 30+ iterations in the
// level projections of G 7 1 and 50-90 in the synchronisation.  A component whose box is more than
// four times its cells keeps the patch-local V-cycles.  The cell widths of the box from the patches
// (those of level 0 refined where no patch has the column).
void cfd_amr::pr_level_mg_setup()
{
    for(cfd_amr_lmg &M : lmg)
    delete M.mg;
    lmg.clear();
    lmg_of.assign(P.size(),-1);

    // G 40: the patch multigrids with their overlap instead
    if(place)
    return;

    for(int l=1; l<=maxlev; ++l)
    {
        // the patches of the level on this rank next to each other (sharing a face): together
        vector<int> up(P.size());
        for(int id : lev[l])
        up[id] = id;
        std::function<int(int)> root = [&](int g) { while(up[g]!=g) { up[g] = up[up[g]]; g = up[g]; } return g; };
        // (and across a periodic edge: the shift of b against a, +-1 length of the level)
        struct jn { int a, b, d, s; };
        vector<jn> join;
        for(int a : lev[l])
        for(int b : lev[l])
        {
            if(b<=a)
            continue;
            const r3patch *A = P[a], *B = P[b];
            for(int d=0; d<3; ++d)
            {
                int sft = 0;
                if(A->hi[d]+1==B->lo[d] || B->hi[d]+1==A->lo[d])
                sft = 0;
                else if(d<2 && per[d] && A->hi[d]==gn(l,d)-1 && B->lo[d]==0)
                sft = 1;
                else if(d<2 && per[d] && B->hi[d]==gn(l,d)-1 && A->lo[d]==0)
                sft = -1;
                else
                continue;
                bool ov = true;
                for(int e=0; e<3; ++e)
                if(e!=d && (A->lo[e]>B->hi[e] || B->lo[e]>A->hi[e]))
                ov = false;
                if(ov)
                {
                    up[root(a)] = root(b);
                    join.push_back({a,b,d,sft});
                }
            }
        }
        std::map<int,vector<int>> groups;
        for(int id : lev[l])
        groups[root(id)].push_back(id);

        for(auto &G : groups)
        {
            if(G.second.size()<2)
            continue;

            cfd_amr_lmg M;
            M.ids = G.second;

            // the positions: the patches joined across a periodic edge shifted by the length of
            // the level, so the box is contiguous; a group around the whole length (the shifts do not
            // agree) or with a patch over it: the box over the whole length, periodic
            std::map<int,vector<int>> shift;
            for(int id : M.ids)
            for(int d=0; d<2; ++d)
            if(per[d] && P[id]->lo[d]==0 && P[id]->hi[d]==gn(l,d)-1)
            M.wrap[d] = true;
            for(int pass=0; pass<2; ++pass)
            {
                shift.clear();
                shift[M.ids[0]] = vector<int>(3,0);
                bool grown = true, clash = false;
                while(grown)
                {
                    grown = false;
                    for(const jn &J : join)
                    {
                        const bool ha = shift.count(J.a)>0, hb = shift.count(J.b)>0;
                        if(ha==hb && !ha)
                        continue;
                        const int L = M.wrap[J.d] ? 0 : J.s*gn(l,J.d);
                        if(ha && hb)
                        {
                            if(shift[J.b][J.d]!=shift[J.a][J.d]+L)
                            {
                                clash = true;
                                M.wrap[J.d] = true;
                            }
                            continue;
                        }
                        if(ha)
                        {
                            shift[J.b] = shift[J.a];
                            shift[J.b][J.d] += L;
                        }
                        else
                        {
                            shift[J.a] = shift[J.b];
                            shift[J.a][J.d] -= L;
                        }
                        grown = true;
                    }
                }
                if(!clash)
                break;
            }

            int hi[3];
            size_t ncell = 0;
            for(int d=0; d<3; ++d)
            {
                M.lo[d] = std::numeric_limits<int>::max();
                hi[d] = std::numeric_limits<int>::min();
            }
            for(int id : M.ids)
            {
                ncell += (size_t)P[id]->n[0]*P[id]->n[1]*P[id]->n[2];
                for(int d=0; d<3; ++d)
                {
                    M.lo[d] = MIN(M.lo[d],P[id]->lo[d]+shift[id][d]);
                    hi[d] = MAX(hi[d],P[id]->hi[d]+shift[id][d]);
                }
            }
            for(int d=0; d<3; ++d)
            {
                if(M.wrap[d])
                {
                    M.lo[d] = 0;
                    hi[d] = gn(l,d)-1;
                }
                M.n[d] = hi[d]-M.lo[d]+1;
            }
            for(int id : M.ids)
            {
                M.off[id] = vector<int>(3);
                for(int d=0; d<3; ++d)
                M.off[id][d] = P[id]->lo[d]+shift[id][d]-M.lo[d];
            }

            // a box mostly outside the patches (a long thin component): the patch-local V-cycles
            if((size_t)M.n[0]*M.n[1]*M.n[2] > 4*ncell)
            continue;

            // periodic: the overlap cells (the margin of a member next to a patch of another rank, OV
            // cells deep across a side); the box grows to hold them
            struct ovc { int id, qq, b[3]; };
            vector<ovc> ovl;
            if(pr_ras())
            {
                const int OV = 3;
                for(int id : M.ids)
                {
                    cfd_amr_patch *c = CP(id);
                    lexer *pp = c->pp;
                    const int nn[3] = {pp->knox,pp->knoy,pp->knoz};
                    for(int i=-OV; i<nn[0]+OV; ++i)
                    for(int j=-OV; j<nn[1]+OV; ++j)
                    for(int k=-OV; k<nn[2]+OV; ++k)
                    {
                        const int q[3] = {i,j,k};
                        int out = 0;
                        for(int d=0; d<3; ++d)
                        if(q[d]<0 || q[d]>=nn[d])
                        ++out;
                        if(out!=1)
                        continue;
                        bool skip = false;
                        for(int d=0; d<3; ++d)
                        if((q[d]<0 || q[d]>=nn[d]) && (rr[d]!=2 || pr_wrap(id,d) || M.wrap[d]))
                        skip = true;
                        const int qq = cix(pp,i,j,k);
                        if(skip || c->mkind[qq]!=2 || pp->flag4[qq]<=0)
                        continue;
                        ovc v;
                        v.id = id;
                        v.qq = qq;
                        for(int d=0; d<3; ++d)
                        v.b[d] = M.o(id,d)+q[d];
                        ovl.push_back(v);
                    }
                }
                int glo[3] = {0,0,0}, ghi[3];
                for(int d=0; d<3; ++d)
                ghi[d] = M.n[d]-1;
                for(const ovc &v : ovl)
                for(int d=0; d<3; ++d)
                {
                    glo[d] = MIN(glo[d],v.b[d]);
                    ghi[d] = MAX(ghi[d],v.b[d]);
                }
                for(int d=0; d<3; ++d)
                {
                    M.lo[d] += glo[d];
                    M.n[d] = ghi[d]-glo[d]+1;
                    for(int id : M.ids)
                    M.off[id][d] -= glo[d];
                    for(ovc &v : ovl)
                    v.b[d] -= glo[d];
                }
            }

            // the multigrid is 2D in x, y with full columns in z: the box layout of sc_level
            M.own.assign((size_t)M.n[0]*M.n[1]*M.n[2],-1);
            auto bix = [&](int i, int j, int k) { return ((size_t)i*M.n[1] + j)*M.n[2] + k; };
            vector<double> wx(M.n[0]+2,-1.0), wy(M.n[1]+2,-1.0);
            for(int id : M.ids)
            {
                r3patch *c = P[id];
                lexer *pp = c->pp;
                for(int i=0; i<c->n[0]; ++i)
                for(int j=0; j<c->n[1]; ++j)
                for(int k=0; k<c->n[2]; ++k)
                M.own[bix(M.o(id,0)+i,M.o(id,1)+j,M.o(id,2)+k)] = id;
                for(int i=-1; i<=c->n[0]; ++i)
                {
                    const int b = M.o(id,0)+i;
                    if(b>=-1 && b<=M.n[0] && (wx[b+1]<0.0 || (i>=0 && i<c->n[0])))
                    wx[b+1] = pp->DXN[i+marge];
                }
                for(int j=-1; j<=c->n[1]; ++j)
                {
                    const int b = M.o(id,1)+j;
                    if(b>=-1 && b<=M.n[1] && (wy[b+1]<0.0 || (j>=0 && j<c->n[1])))
                    wy[b+1] = pp->DYN[j+marge];
                }
            }
            if(!ovl.empty())
            {
                M.ovm.assign(M.own.size(),-1);
                M.ovq.assign(M.own.size(),-1);
                for(const ovc &v : ovl)
                {
                    const size_t b = bix(v.b[0],v.b[1],v.b[2]);
                    if(M.own[b]!=-1)
                    continue;
                    M.own[b] = -2;
                    M.ovm[b] = v.id;
                    M.ovq[b] = v.qq;
                }
            }

            // the columns of no patch: the width of level 0 refined
            double f = 1.0, fy = 1.0;
            for(int m=0; m<l; ++m)
            {
                f /= rr[0];
                fy /= rr[1];
            }
            for(int b=-1; b<=M.n[0]; ++b)
            if(wx[b+1]<0.0)
            {
                const int I = (int)std::floor(double(M.lo[0]+b)*f) - org[0];
                wx[b+1] = p0->DXN[MIN(MAX(I,-1),p0->knox)+marge]*f;
            }
            for(int b=-1; b<=M.n[1]; ++b)
            if(wy[b+1]<0.0)
            {
                const int J = (int)std::floor(double(M.lo[1]+b)*fy) - org[1];
                wy[b+1] = p0->DYN[MIN(MAX(J,-1),p0->knoy)+marge]*fy;
            }

            M.mg = new reefmg_core;
            M.mg->set_precision(p0->N14==32 ? 32 : 64);
            M.mg->set_fine_operator(&pr_fineop_none);
            if(!M.mg->setup(MPI_COMM_SELF,M.n[0],M.n[1],M.n[2],M.n[0],M.n[1],0,wx.data(),wy.data(),M.wrap[0],M.wrap[1]))
            {
                cout<<"CFD AMR: level multigrid "<<M.mg->err()<<" (level "<<l<<"), the patch-local ones"<<endl;
                delete M.mg;
                continue;
            }
            for(int id : M.ids)
            lmg_of[id] = (int)lmg.size();
            lmg.push_back(std::move(M));
        }
    }
}

// the coefficients of the level multigrid from the rows of the patches (pr_level: the couplings
// across the patch edges to the coarser level on the diagonal): a coupling to a cell of another
// patch of the level on this rank stays, the others (another rank, the boundary) are dropped as in
// the patch-local multigrids
void cfd_amr::pr_level_mg_fill(int m)
{
    cfd_amr_lmg &M = lmg[m];
    sc_level &L = M.mg->fine();
    pr_clear(L);
    auto bix = [&](int i, int j, int k) { return ((size_t)i*M.n[1] + j)*M.n[2] + k; };

    for(int i=0; i<M.n[0]; ++i)
    for(int j=0; j<M.n[1]; ++j)
    for(int k=0; k<M.n[2]; ++k)
    pr_coef(M.mg,L,L.idx(i,j,k),1.0,0.0,0.0,0.0,0.0,0.0,0.0);

    for(int id : M.ids)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        const matrix_diag &A = c->a->M;
        const int o[3] = {M.o(id,0), M.o(id,1), M.o(id,2)};
        // the neighbour at side dir of direction d of patch cell (i,j,k): in the patch, or a cell
        // of another patch of the level with a row
        auto nb = [&](int i, int j, int k, int d, int dir) -> bool
        {
            int q[3] = {i,j,k};
            q[d] += dir;
            if(q[d]>=0 && q[d]<pp->knox*(d==0) + pp->knoy*(d==1) + pp->knoz*(d==2))
            return true;
            int b[3] = {q[0]+o[0],q[1]+o[1],q[2]+o[2]};
            if(M.wrap[d])
            b[d] = (b[d]+M.n[d])%M.n[d];
            if(b[d]<0 || b[d]>=M.n[d])
            return false;
            const int g = M.own[bix(b[0],b[1],b[2])];
            if(g==-2)
            return true;
            if(g<0 || g==id)
            return false;
            cfd_amr_patch *e = CP(g);
            const int qi = b[0]-M.o(g,0), qj = b[1]-M.o(g,1), qk = b[2]-M.o(g,2);
            return e->row[cix(e->pp,qi,qj,qk)]>=0;
        };

        for(int i=0; i<pp->knox; ++i)
        for(int j=0; j<pp->knoy; ++j)
        for(int k=0; k<pp->knoz; ++k)
        {
            const long q = L.idx(i+o[0],j+o[1],k+o[2]);
            const int r = c->row[cix(pp,i,j,k)];
            if(r<0)
            continue;
            pr_coef(M.mg,L,q,A.p[r]*(1.0+1.0e-6),
                    nb(i,j,k,0,1) ? A.n[r] : 0.0, nb(i,j,k,0,-1) ? A.s[r] : 0.0,
                    nb(i,j,k,1,1) ? A.w[r] : 0.0, nb(i,j,k,1,-1) ? A.e[r] : 0.0,
                    nb(i,j,k,2,1) ? A.t[r] : 0.0, nb(i,j,k,2,-1) ? A.b[r] : 0.0);
            L.act[q] = 1;
        }
    }
    pr_level_mg_ovrows(m,1.0+1.0e-6);
    M.mg->coarsen();
}

// periodic: the rows of the overlap cells of a joint multigrid (the rows of the patches of the other
// rank there, from the kc7 fields of the members), the couplings to cells of the box with a row
void cfd_amr::pr_level_mg_ovrows(int m, double dshift)
{
    cfd_amr_lmg &M = lmg[m];
    if(M.ovm.empty())
    return;
    sc_level &L = M.mg->fine();
    auto bix = [&](int i, int j, int k) { return ((size_t)i*M.n[1] + j)*M.n[2] + k; };
    auto has_row = [&](int i, int j, int k) -> bool
    {
        int b[3] = {i,j,k};
        for(int d=0; d<3; ++d)
        {
            if(M.wrap[d])
            b[d] = (b[d]+M.n[d])%M.n[d];
            if(b[d]<0 || b[d]>=M.n[d])
            return false;
        }
        const int g = M.own[bix(b[0],b[1],b[2])];
        if(g==-2)
        return true;
        if(g<0)
        return false;
        cfd_amr_patch *e = CP(g);
        return e->row[cix(e->pp,b[0]-M.o(g,0),b[1]-M.o(g,1),b[2]-M.o(g,2))]>=0;
    };
    for(int i=0; i<M.n[0]; ++i)
    for(int j=0; j<M.n[1]; ++j)
    for(int k=0; k<M.n[2]; ++k)
    {
        const size_t b = bix(i,j,k);
        if(M.own[b]!=-2)
        continue;
        cfd_amr_patch *c = CP(M.ovm[b]);
        const int qq = M.ovq[b];
        double cc[7];
        for(int n=0; n<7; ++n)
        cc[n] = c->kc7[n]->V[qq];
        if(cc[0]==0.0)
        continue;
        const long q = L.idx(i,j,k);
        pr_coef(M.mg,L,q,cc[0]*dshift,
                has_row(i+1,j,k) ? cc[1] : 0.0, has_row(i-1,j,k) ? cc[2] : 0.0,
                has_row(i,j+1,k) ? cc[3] : 0.0, has_row(i,j-1,k) ? cc[4] : 0.0,
                has_row(i,j,k+1) ? cc[5] : 0.0, has_row(i,j,k-1) ? cc[6] : 0.0);
        L.act[q] = 1;
    }
}

// the level multigrid of the composite projection (pr_matrix): the rows of the patch multigrids
// (the couplings across the patch edges on the diagonal, those of a coarse cell next to a finer
// patch as there), the couplings between the patches of the level on this rank added
void cfd_amr::pr_level_mg_fill_comp(int m)
{
    cfd_amr_lmg &M = lmg[m];
    reefmg_core *U = M.mg;
    sc_level &L = U->fine();
    pr_clear(L);
    const bool f32 = U->precision()==32;
    auto bix = [&](int i, int j, int k) { return ((size_t)i*M.n[1] + j)*M.n[2] + k; };

    for(int i=0; i<M.n[0]; ++i)
    for(int j=0; j<M.n[1]; ++j)
    for(int k=0; k<M.n[2]; ++k)
    pr_coef(U,L,L.idx(i,j,k),1.0,0.0,0.0,0.0,0.0,0.0,0.0);

    for(int id : M.ids)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        const matrix_diag &A = c->a->M;
        const sc_level &S = c->mg->fine();
        const int o[3] = {M.o(id,0), M.o(id,1), M.o(id,2)};
        const int nn[3] = {pp->knox,pp->knoy,pp->knoz};
        auto sib = [&](int i, int j, int k, int d, int dir) -> bool
        {
            int q[3] = {i,j,k};
            q[d] += dir;
            if(q[d]>=0 && q[d]<nn[d])
            return false;
            int b[3] = {q[0]+o[0],q[1]+o[1],q[2]+o[2]};
            if(M.wrap[d])
            b[d] = (b[d]+M.n[d])%M.n[d];
            if(b[d]<0 || b[d]>=M.n[d])
            return false;
            const int g = M.own[bix(b[0],b[1],b[2])];
            if(g==-2)
            return true;
            if(g<0 || g==id)
            return false;
            cfd_amr_patch *e = CP(g);
            return e->row[cix(e->pp,b[0]-M.o(g,0),b[1]-M.o(g,1),b[2]-M.o(g,2))]>=0;
        };

        for(int i=0; i<nn[0]; ++i)
        for(int j=0; j<nn[1]; ++j)
        for(int k=0; k<nn[2]; ++k)
        {
            const long q = L.idx(i+o[0],j+o[1],k+o[2]), s = S.idx(i,j,k);
            double cc[7];
            if(f32)
            {
                cc[0]=S.pf[s]; cc[1]=S.nf[s]; cc[2]=S.sf[s]; cc[3]=S.wf[s]; cc[4]=S.ef[s]; cc[5]=S.tf[s]; cc[6]=S.bf[s];
            }
            else
            {
                cc[0]=S.p[s]; cc[1]=S.n[s]; cc[2]=S.s[s]; cc[3]=S.w[s]; cc[4]=S.e[s]; cc[5]=S.t[s]; cc[6]=S.b[s];
            }
            const int r = c->row[cix(pp,i,j,k)];
            if(r>=0)
            {
                if(sib(i,j,k,0,1))  cc[1] = A.n[r];
                if(sib(i,j,k,0,-1)) cc[2] = A.s[r];
                if(sib(i,j,k,1,1))  cc[3] = A.w[r];
                if(sib(i,j,k,1,-1)) cc[4] = A.e[r];
                if(sib(i,j,k,2,1))  cc[5] = A.t[r];
                if(sib(i,j,k,2,-1)) cc[6] = A.b[r];
            }
            pr_coef(U,L,q,cc[0],cc[1],cc[2],cc[3],cc[4],cc[5],cc[6]);
            L.act[q] = S.act[s];
        }
    }
    pr_level_mg_ovrows(m,1.0);
    U->coarsen();
}

// one V-cycle of the level multigrid of level l on the right-hand side kr of its patches, added
// to z
void cfd_amr::pr_level_vcycle(int m, int kr, int kz)
{
    const double t0 = MPI_Wtime();
    cfd_amr_lmg &M = lmg[m];
    sc_level &L = M.mg->fine();
    std::fill(L.u.begin(),L.u.end(),0.0);
    std::fill(L.f.begin(),L.f.end(),0.0);
    for(int g : M.ids)
    {
        cfd_amr_patch *c = CP(g);
        lexer *pp = c->pp;
        const double *t = pvec(g,kr);
        const int o[3] = {M.o(g,0), M.o(g,1), M.o(g,2)};
        for(int i=0; i<pp->knox; ++i)
        for(int j=0; j<pp->knoy; ++j)
        for(int k=0; k<pp->knoz; ++k)
        {
            const int qq = cix(pp,i,j,k);
            if(c->row[qq]>=0)
            L.f[L.idx(i+o[0],j+o[1],k+o[2])] = t[qq];
        }
    }
    // (periodic: the right-hand side in the overlap, filled into the margins of the members)
    if(!M.ovm.empty())
    for(size_t b=0; b<M.own.size(); ++b)
    if(M.own[b]==-2)
    {
        const int i = (int)(b/((size_t)M.n[1]*M.n[2])), j = (int)((b/M.n[2])%M.n[1]), k = (int)(b%M.n[2]);
        L.f[L.idx(i,j,k)] = pvec(M.ovm[b],kr)[M.ovq[b]];
    }
    M.mg->vcycle(0,2,2);
    for(int g : M.ids)
    {
        cfd_amr_patch *c = CP(g);
        lexer *pp = c->pp;
        double *z = pvec(g,kz);
        const int o[3] = {M.o(g,0), M.o(g,1), M.o(g,2)};
        for(int i=0; i<pp->knox; ++i)
        for(int j=0; j<pp->knoy; ++j)
        for(int k=0; k<pp->knoz; ++k)
        {
            const int qq = cix(pp,i,j,k);
            if(c->row[qq]>=0)
            z[qq] += L.u[L.idx(i+o[0],j+o[1],k+o[2])];
        }
    }
    tp[0] += MPI_Wtime()-t0;
}

void cfd_amr::pr_level_solve(int l, int kr, int kz)
{
    // G 40, periodic: the right-hand side in the overlap of the patch multigrids (the patches next
    // to them)
    if(pr_ras())
    {
        vector<fspec> fs;
        fs.push_back({CELL,2,[this,kr](int g) -> field& { return pfield(g,kr); }});
        fill(l,fs,true);
    }
    for(int m=0; m<(int)lmg.size(); ++m)
    if(P[lmg[m].ids[0]]->lev==l)
    pr_level_vcycle(m,kr,kz);
    for(int id : lev[l])
    if(lmg_of[id]<0)
    pr_patch_solve(id,kr,kz,-1);
}

// the coarse value at the position of the fine cell of a face at a patch edge: from the coarse cell
// and its leaf neighbours in the transverse directions (pr_weights)
void cfd_amr::pr_stencils()
{
    for(cfd_amr_cf &F : cf)
    {
        const int lc = CP(F.pid)->lev-1;
        lexer *qc = glex(F.gc);
        int og[3];
        goff(F.gc,og);

        // (past a periodic edge: the cell at the other end; on level 0 of one rank its own index)
        auto leaf = [&](int *cl) -> bool
        {
            const int G[3] = {cl[0]+og[0],cl[1]+og[1],cl[2]+og[2]};
            int Gw[3];
            if(!wrap(lc,G,Gw))
            return false;
            if(lc>0 && !refined(lc,Gw))
            return false;
            if(covered(lc,Gw))
            return false;
            if(Gw[0]!=G[0] || Gw[1]!=G[1] || Gw[2]!=G[2])
            {
                int Lw[3];
                if(lc>0 || F.gc!=-1)
                return false;
                if(wrap_local0(G,Lw))
                for(int d=0; d<3; ++d)
                cl[d] = Lw[d];
            }
            if(cl[0]<-1 || cl[0]>qc->knox || cl[1]<-1 || cl[1]>qc->knoy || cl[2]<-1 || cl[2]>qc->knoz)
            return false;
            return qc->flag4[cix(qc,cl[0],cl[1],cl[2])]>0;
        };

        for(int t=0; t<3; ++t)
        {
            F.tq[t][0] = F.tq[t][1] = -1;
            F.tdel[t] = 0.0;
            if(t==F.d || rr[t]!=2)
            continue;
            F.tdel[t] = (F.o[t]==0) ? -0.25 : 0.25;
            int cm[3] = {F.ci[0],F.ci[1],F.ci[2]}, cp[3] = {F.ci[0],F.ci[1],F.ci[2]};
            cm[t] -= 1;
            cp[t] += 1;
            if(leaf(cm))
            F.tq[t][0] = cix(qc,cm[0],cm[1],cm[2]);
            if(leaf(cp))
            F.tq[t][1] = cix(qc,cp[0],cp[1],cp[2]);
        }
    }
}

namespace
{
// the porosity factor CPOR n of the pressure gradient (and of every momentum source) at a face of
// porosity n (B 200 1; looping.h, pjm_corr): 1 in open flow
inline double porfac(const lexer *p, double n)
{
    return (n<1.0 ? 1.0/(1.0 + p->B260*(p->B266==1 ? (1.0-n)/n : 1.0)) : 1.0)*n;
}

// R(phi) = int_0^phi rho ds with the density of the smoothed level set (half width e)
inline double rint(double phi, double e, double row, double roa)
{
    auto I = [&](double x) -> double
    {
        if(x>=e) return x;
        if(x<=-e) return 0.0;
        return 0.5*((x+e) + (x*x-e*e)/(2.0*e) - (e/(PI*PI))*(cos(PI*x/e)+1.0));
    };
    return roa*phi + (row-roa)*(I(phi)-I(0.0));
}

// int rho ds between two points at the distance ds with the level set phi1, phi2 (linear between)
inline double rmass(double phi1, double phi2, double ds, double e, double row, double roa)
{
    if(fabs(phi2-phi1)>1.0e-6*fabs(ds))
    return ds*(rint(phi2,e,row,roa)-rint(phi1,e,row,roa))/(phi2-phi1);
    const double H = heaviside(0.5*(phi1+phi2),e);
    return ds*(row*H + roa*(1.0-H));
}
}

// the weights of the coarse value at the position of the fine cell, for the stage: linear in the
// transverse directions in the density-weighted distance (the acceleration (dp/ds)/rho is taken
// constant between the cells, not dp/ds): at a free surface across a patch edge, the hydrostatic
// pressure and the pressure correction keep their kink, the faces in the air no spurious gradient
void cfd_amr::pr_weights()
{
    // face density of the fine face (the level set of the fine cell and of the cell across)
    for(cfd_amr_cf &F : cf)
    {
        cfd_amr_patch *c = CP(F.pid);
        int o[3] = {0,0,0};
        o[F.d] = F.dir;
        const double w = interface_width_face(c->pp,c->a->phi,F.fi[0],F.fi[1],F.fi[2],o[0],o[1],o[2]);
        const double H = heaviside(0.5*(c->a->phi.V[F.qf] + c->a->phi.V[F.qm]),w);
        F.rof = rmix(c->pp,H,0.5*(sval(F.pid,F.qf) + sval(F.pid,F.qm)));

        // porous structures: the face density divided by the porosity factor of the face (all
        // the fluxes of the face, kf, kc and the corrections, carry it as pjm_corr / poisson_pcorr)
        F.rof /= porfac(c->pp,0.5*(c->a->porosity.V[F.qm] + c->a->porosity.V[F.qf]));
    }

    for(cfd_amr_cf &F : cf)
    {
        lexer *qc = glex(F.gc);
        fdm *ac = gfd(F.gc);
        const double e = qc->psi, row = qc->W1, roa = qc->W3;
        const double *dn[3] = {qc->DXN,qc->DYN,qc->DZN};
        const double phC = ac->phi.V[F.qc];
        const double phf = CP(F.pid)->a->phi.V[F.qf];

        F.ns = 1;
        F.sg[0] = F.gc;
        F.sq[0] = F.qc;
        F.sw[0] = 1.0;
        for(int e=0; e<3; ++e)
        F.sc[0][e] = F.ci[e];

        auto put = [&](int q, double w, int t, int off)
        {
            F.sg[F.ns] = F.gc; F.sq[F.ns] = q; F.sw[F.ns] = w;
            for(int e=0; e<3; ++e)
            F.sc[F.ns][e] = F.ci[e];
            F.sc[F.ns][t] += off;
            ++F.ns;
        };

        for(int t=0; t<3; ++t)
        {
            if(F.tdel[t]==0.0)
            continue;
            const double h = dn[t][F.ci[t]+marge];
            const double mf = rmass(phC,phf,F.tdel[t]*h,e,row,roa);
            const int qm = F.tq[t][0], qp = F.tq[t][1];

            if(qm>=0 && qp>=0)
            {
                const double mm = rmass(ac->phi.V[qm],ac->phi.V[qp],2.0*h,e,row,roa);
                put(qp,mf/mm,t,1);
                put(qm,-mf/mm,t,-1);
            }
            else if(qp>=0)
            {
                const double mp = rmass(phC,ac->phi.V[qp],h,e,row,roa);
                put(qp,mf/mp,t,1);
                F.sw[0] -= mf/mp;
            }
            else if(qm>=0)
            {
                const double mm = rmass(phC,ac->phi.V[qm],-h,e,row,roa);
                put(qm,mf/mm,t,-1);
                F.sw[0] -= mf/mm;
            }
        }
    }
}

// The hydrostatic pressure of the fine cell and at its position on the coarse side, for the faces
// at a patch edge (gravity along z): the discrete hydrostatic pressure of each grid (the face
// densities of its own w equation), integrated down from the plane of the top of the patch.  The
// predictor takes the gradient across such a face without it (and a face normal to z without the
// gravity it balances): the hydrostatic pressures of the coarse and the fine grid differ inside the
// interface band (where the free surface crosses a patch edge), at rest the patch edge then has no
// gradient of the rest.
void cfd_amr::pr_hydro(int lsub)
{
    // the level set of grid gg: with subcycling the coarser grids at the time of the stage (as
    // their pressure, tsf), the grids of the level their own
    auto lsf = [&](int gg) -> field& { return (lsub>=0 && glevel(gg)<lsub) ? tsf(gg,tl_phi) : gfd(gg)->phi; };

    const bool zgrav = (p0->W20==0.0 && p0->W21==0.0 && p0->W22!=0.0);
    const double g = fabs(p0->W22);


    // grid gg, column (i,j): p_h of cell kc, from the plane on top of cell ktop-1 (kc <= ktop)
    auto column = [&](int gg, int i, int j, int ktop, int kc, double &val) -> bool
    {
        lexer *q = glex(gg);
        field &phi = lsf(gg);
        const int m = q->margin;
        const int k0 = MIN(kc,ktop-1);
        if(k0<-m || ktop>q->knoz+m-1 || i<-m || i>q->knox+m-1 || j<-m || j>q->knoy+m-1)
        return false;
        // densities as the grid has them: the local interface width (interface_width.h)
        auto rho = [&](double ph, double e) { const double H = heaviside(ph,e); return q->W1*H + q->W3*(1.0-H); };
        for(int k=k0; k<=MIN(kc+1,ktop); ++k)
        if(q->flag4[cix(q,i,j,k)]<=0 && k<ktop)
        return false;

        double h = g*rho(phi(i,j,ktop-1),interface_width(q,phi,i,j,ktop-1))*0.5*q->DZN[ktop-1+marge];   // cell ktop-1
        if(kc==ktop)
        {
            if(q->flag4[cix(q,i,j,ktop)]<=0)
            return false;
            val = h - g*rho(0.5*(phi(i,j,ktop-1)+phi(i,j,ktop)),interface_width_face(q,phi,i,j,ktop-1,0,0,1))*q->DZP[ktop-1+marge];
            return true;
        }
        for(int k=ktop-2; k>=kc; --k)
        {
            if(q->flag4[cix(q,i,j,k)]<=0)
            return false;
            h += g*rho(0.5*(phi(i,j,k)+phi(i,j,k+1)),interface_width_face(q,phi,i,j,k,0,0,1))*q->DZP[k+marge];
        }
        val = h;
        return true;
    };

    for(cfd_amr_cf &F : cf)
    {
        F.hf = F.hs = 0.0;
        F.hok = false;
        if(!zgrav || rr[2]!=2)
        continue;

        cfd_amr_patch *c = CP(F.pid);
        lexer *pp = c->pp;

        // fine column
        double hfine;
        if(!column(F.pid,F.fi[0],F.fi[1],pp->knoz,F.fi[2],hfine))
        continue;

        // the coarse hydrostatic pressure at the position of the fine cell, up to the same plane:
        // for a face normal to x or y in the column of the covered parent (its level set restricted
        // from the fine cells: the same profile, only the discretisation differs; the difference
        // of the hydrostatic pressure between the coarse and the fine column, the slope of the
        // free surface, stays in the gradient), for a face normal to z in the column of the
        // coarse cell (the same column)
        const int gh = (F.d==2) ? F.gc : F.gp;
        const int *ch = (F.d==2) ? F.ci : F.pi;
        int og[3];
        goff(gh,og);
        const int ktop = fdiv(c->hi[2]+1,rr[2]) - og[2];
        double hs = 0.0;
        bool ok = true;
        for(int m=0; m<F.ns && ok; ++m)
        {
            int sc[3];
            for(int e=0; e<3; ++e)
            sc[e] = ch[e] + F.sc[m][e] - F.ci[e];
            if(sc[2]>ktop)
            {
                ok = false;
                break;
            }
            double hc;
            if(!column(gh,sc[0],sc[1],ktop,sc[2],hc))
            ok = false;
            else
            hs += F.sw[m]*hc;
        }
        if(!ok)
        continue;

        F.hf = hfine;
        F.hs = hs;
        F.hok = true;
    }
}

double cfd_amr::pr_xs(const cfd_amr_cf &F, int k)
{
    double v = 0.0;
    for(int m=0; m<F.ns; ++m)
    v += F.sw[m]*pvec(F.sg[m],k)[F.sq[m]];
    return v;
}

double cfd_amr::pr_ps(const cfd_amr_cf &F)
{
    double v = 0.0;
    for(int m=0; m<F.ns; ++m)
    v += F.sw[m]*gfd(F.sg[m])->press.V[F.sq[m]];
    return v;
}

// the stage predictor at the faces of the patch edges: the coarse face velocity, linear in the
// transverse directions
void cfd_amr::pr_predict(int s)
{
    // G 40: the coarse data of the faces at the edges of the level-1 patches in the shadows
    if(place)
    {
        vector<std::function<field&(int)>> f;
        f.push_back([this](int g) -> field& { return gfd(g)->phi; });
        f.push_back([this](int g) -> field& { return gfd(g)->press; });
        for(int c=0; c<3; ++c)
        f.push_back([this,s,c](int g) -> field& { return gmom(g)->amr_velout(gfd(g),c,s); });
        if(scal)
        f.push_back([this](int g) -> field& { return sfield(g); });
        sh_down_f(f);
    }

    pr_weights();
    pr_hydro();

    for(cfd_amr_cf &F : cf)
    {
        field &fc = gmom(F.gface)->amr_velout(gfd(F.gface),F.d,s);
        double v = fc(F.cf[0],F.cf[1],F.cf[2]);
        for(int t=0; t<3; ++t)
        if(t!=F.d)
        v += cell_slope(fc,F.gface,F.cf,t,1)*(F.o[t]==0 ? -0.25 : 0.25);

        // the pressure gradient of the coarse face (in the coarse stage velocity) replaced by the
        // one across the fine face: the fine cell and the coarse cell, the coarse pressure taken
        // to the height (position) of the fine cell; the face density of the fine face.  At rest
        // the hydrostatic pressure has no gradient across the patch edge then.
        {
            cfd_amr_patch *c = CP(F.pid);
            lexer *pp = c->pp;
            fdm *af = c->a;
            fdm *ac = gfd(F.gc);
            fdm *ap = gfd(F.gp);
            lexer *qc = glex(F.gc);
            const double bs = mom0->amr_alpha(s);

            const double pC = ac->press(F.ci[0],F.ci[1],F.ci[2]);
            const double pP = ap->press(F.pi[0],F.pi[1],F.pi[2]);
            const double pCs = pr_ps(F) - F.hs;
            const double pf = af->press(F.fi[0],F.fi[1],F.fi[2]) - F.hf;

            const double phC = ac->phi(F.ci[0],F.ci[1],F.ci[2]);
            const double phP = ap->phi(F.pi[0],F.pi[1],F.pi[2]);
            lexer *qp = glex(F.gp);
            const double wc = 0.5*(interface_width(qc,ac->phi,F.ci[0],F.ci[1],F.ci[2]) + interface_width(qp,ap->phi,F.pi[0],F.pi[1],F.pi[2]));
            const double Hc = heaviside(0.5*(phC+phP),wc);
            // the coarse face with its porosity factor (the coarse stage velocity has it on its
            // pressure gradient and its gravity)
            const double fpc = porfac(qc,0.5*(ac->porosity(F.ci[0],F.ci[1],F.ci[2]) + ap->porosity(F.pi[0],F.pi[1],F.pi[2])));
            const double roc = rmix(qc,Hc,0.5*(sval(F.gc,F.qc) + sval(F.gp,F.qp)))/fpc;
            const double rof = F.rof;

            // the spacing of the centres across the coarse face, on the grid that holds it (the
            // parent and the coarse cell may lie on different patches of the coarser level)
            lexer *ql = (F.dir>0) ? qp : qc;
            const double *dp[3] = {ql->DXP,ql->DYP,ql->DZP};
            const int Lc = (F.dir>0) ? F.pi[F.d] : F.ci[F.d];
            const double Gc = double(F.dir)*(pC-pP)/(dp[F.d][Lc+marge]*roc);
            const double Gf = double(F.dir)*(pCs-pf)/(F.dist*rof);

            v += bs*p0->dt*(Gc-Gf);

            // a face normal to gravity: the hydrostatic part of the pressure gradient balances the
            // gravity of the coarse stage velocity, both out
            if(F.hok && F.d==2)
            v -= bs*p0->dt*p0->W22*fpc;
        }

        F.ustar = v;
        field &ff = CP(F.pid)->pmom->amr_velout(CP(F.pid)->a,F.d,s);
        ff(F.ff[0],F.ff[1],F.ff[2]) = v;
    }

    // the coarse faces: the mean of their fine faces (the right-hand sides of both sides see the
    // same flux).  At a high patch edge the face belongs to the grid of the parent; when the coarse
    // cell lies on another patch of the coarser level, its right-hand side takes the face from its
    // own copy (in the cells around that patch): set as well
    for(cfd_amr_cf &F : cf)
    {
        gmom(F.gface)->amr_velout(gfd(F.gface),F.d,s)(F.cf[0],F.cf[1],F.cf[2]) = 0.0;
        if(F.dir>0 && F.gc!=F.gface)
        {
            int c[3] = {F.ci[0],F.ci[1],F.ci[2]};
            c[F.d] -= 1;
            gmom(F.gc)->amr_velout(gfd(F.gc),F.d,s)(c[0],c[1],c[2]) = 0.0;
        }
    }
    for(cfd_amr_cf &F : cf)
    {
        gmom(F.gface)->amr_velout(gfd(F.gface),F.d,s)(F.cf[0],F.cf[1],F.cf[2]) += F.area*F.ustar;
        if(F.dir>0 && F.gc!=F.gface)
        {
            int c[3] = {F.ci[0],F.ci[1],F.ci[2]};
            c[F.d] -= 1;
            gmom(F.gc)->amr_velout(gfd(F.gc),F.d,s)(c[0],c[1],c[2]) += F.area*F.ustar;
        }
    }

    // G 40: the coarse faces of the shadows onto level 0
    if(place)
    for(int d=0; d<3; ++d)
    {
        vector<std::function<field&(int)>> f(1,[this,s,d](int g) -> field& { return gmom(g)->amr_velout(gfd(g),d,s); });
        sh_run(pu_cf[d],f,false);
    }
}

// multigrid coefficients from the assembled rows, then the couplings across the patch edges
void cfd_amr::pr_matrix(int s)
{
    // level 0: its own matrix, covered cells included
    {
        sc_level &L = mg0->fine();
        pr_clear(L);
        const matrix_diag &M = a0->M;
        for(int i=0; i<p0->knox; ++i)
        for(int j=0; j<p0->knoy; ++j)
        for(int k=0; k<p0->knoz; ++k)
        {
            const long q = L.idx(i,j,k);
            const int r = row0[cix(p0,i,j,k)];
            if(r<0)
            {
                pr_coef(mg0,L,q,1.0,0.0,0.0,0.0,0.0,0.0,0.0);
                continue;
            }
            pr_coef(mg0,L,q,M.p[r],M.n[r],M.s[r],M.w[r],M.e[r],M.t[r],M.b[r]);
            L.act[q] = 1;
        }
    }

    // patches: their own matrix (before the couplings across the patch edges below), the couplings
    // to the cells around the patch dropped (G 40: those into the overlap kept, the rows of the
    // patches next to it there)
    if(pr_ras())
    pr_coef_fill(-1);
    for(int id=0; id<npat; ++id)
    pr_patch_mg(id,1.0);

    // the faces at the patch edges
    for(cfd_amr_cf &F : cf)
    {
        cfd_amr_patch *c = CP(F.pid);
        lexer *pp = c->pp;
        fdm *a = c->a;

        const double rof = F.rof;
        F.kf = 1.0/(rof*F.dist*F.dxnf);
        F.kc = F.area/(rof*F.dist*F.dxnc);

        // the couplings across the face are in pr_apply (the fine row to the coarse value at its
        // position, the coarse row to the fluxes of its fine faces); the diagonal of the fine row
        // keeps its part
        matrix_diag &Mf = a->M;
        double &cfn = mcoef(Mf,F.nf,F.d,F.dir);
        Mf.p[F.nf] += cfn + F.kf;
        cfn = 0.0;

        // (G 40: the coarse cell in a shadow; its rank folds the row, sh_folds)
        if(is_shadow(F.gc))
        continue;
        matrix_diag &Mc = gfd(F.gc)->M;
        double &ccn = mcoef(Mc,F.nc,F.d,-F.dir);
        Mc.p[F.nc] += ccn;
        ccn = 0.0;
    }
    if(place)
    sh_folds();

    // the multigrid of a coarser grid (level 0 or a patch with finer patches): a coarse cell next
    // to a patch couples to the covered cell across the face with the couplings of the fine faces
    // (sum of kc, their densities and distances; the coarse face may lie on the other side of the
    // free surface, its density off by the density ratio), the covered cell back the same way.
    // The preconditioner then sees the operator of the composite rows (2 levels with the free
    // surface across the patches: 27 -> 11 BiCGStab iterations)
    {
        std::map<std::tuple<int,int,int,int,int,int>,double> ks;
        for(const cfd_amr_cf &F : cf)
        if(F.gc==F.gp)
        ks[std::make_tuple(F.gc,F.ci[0],F.ci[1],F.ci[2],F.d,F.dir)] += F.kc;

        // G 40: the keys of the shadows to the ranks of their coarse cells (sh_ks)
        if(place)
        {
            for(vector<double> &v : kval)
            std::fill(v.begin(),v.end(),0.0);
            for(const auto &kv : ks)
            {
                const int g = std::get<0>(kv.first);
                if(!is_shadow(g))
                continue;
                const int *og = P[g]->lo;
                const auto key = std::make_tuple(g,std::get<1>(kv.first)+og[0],std::get<2>(kv.first)+og[1],std::get<3>(kv.first)+og[2],
                                                 std::get<4>(kv.first),std::get<5>(kv.first));
                kval[g-npat][kidx[key]] = kv.second;
            }
            sh_ks();
        }

        for(const auto &kv : ks)
        {
            const int g = std::get<0>(kv.first);
            if(is_shadow(g))
            continue;
            int ci[3] = {std::get<1>(kv.first),std::get<2>(kv.first),std::get<3>(kv.first)};
            const int d = std::get<4>(kv.first), dir = std::get<5>(kv.first);
            lexer *q = glex(g);
            int n[3] = {q->knox,q->knoy,q->knoz};
            if(g>=0)        // (G 40: the overlap of the patch multigrid)
            for(int e=0; e<3; ++e)
            {
                ci[e] += CP(g)->olo[e];
                n[e] += CP(g)->olo[e]+CP(g)->ohi[e];
            }
            ks_apply((g<0) ? mg0 : CP(g)->mg,ci,d,dir,n,kv.second);
        }
    }

    pr_covered_coef();

    mg0->coarsen();
    for(int id=0; id<npat; ++id)
    if(lmg_of[id]<0)
    CP(id)->mg->coarsen();
    for(int m=0; m<(int)lmg.size(); ++m)
    pr_level_mg_fill_comp(m);
}

// the multigrid of a coarser grid under the finer patches: the coupling between two covered cells
// from the densities of the fine faces at the position of the coarse face (1/rho their mean, area
// weighted) instead of the density of the coarse face (its level set smeared over the interface
// band of the coarser level).  The coarse correction of the FAC sweep then sees the free surface as
// sharp as the finer level does.  The covered rows are not in the composite operator: the
// preconditioner alone changes
void cfd_amr::pr_covered_coef()
{
    if(p0->F30<=0 || place)
    return;

    // the face density of grid g at the + face of its local cell I along d
    auto face_rho = [&](int g, const int *I, int d) -> double
    {
        lexer *q = glex(g);
        field &phi = gfd(g)->phi;
        int o[3] = {0,0,0};
        o[d] = 1;
        const int qa = cix(q,I[0],I[1],I[2]), qb = cix(q,I[0]+o[0],I[1]+o[1],I[2]+o[2]);
        const double H = heaviside_face(q,phi,0.5*(phi.V[qa]+phi.V[qb]),I[0],I[1],I[2],o[0],o[1],o[2]);
        return rmix(q,H,0.5*(sval(g,qa)+sval(g,qb)));
    };
    auto getc = [](reefmg_core *U, sc_level &L, long lq, int m) -> double
    {
        if(U->precision()==32)
        {
            float *a[7] = {L.pf.data(),L.nf.data(),L.sf.data(),L.wf.data(),L.ef.data(),L.tf.data(),L.bf.data()};
            return a[m][lq];
        }
        double *a[7] = {L.p.data(),L.n.data(),L.s.data(),L.w.data(),L.e.data(),L.t.data(),L.b.data()};
        return a[m][lq];
    };
    auto setc = [](reefmg_core *U, sc_level &L, long lq, int m, double v)
    {
        if(U->precision()==32)
        {
            float *a[7] = {L.pf.data(),L.nf.data(),L.sf.data(),L.wf.data(),L.ef.data(),L.tf.data(),L.bf.data()};
            a[m][lq] = (float)v;
            return;
        }
        double *a[7] = {L.p.data(),L.n.data(),L.s.data(),L.w.data(),L.e.data(),L.t.data(),L.b.data()};
        a[m][lq] = v;
    };
    const int mplus[3] = {1,3,5}, mminus[3] = {2,4,6};      // n s, w e, t b

    for(int g=-1; g<npat; ++g)
    {
        const int lc = glevel(g);
        if(lc>=maxlev)
        continue;
        lexer *q = glex(g);
        reefmg_core *U = (g<0) ? mg0 : CP(g)->mg;
        sc_level &L = U->fine();
        const int *ol = (g<0) ? nullptr : CP(g)->olo;
        const vector<int> &row = (g<0) ? row0 : CP(g)->row;
        int og[3];
        goff(g,og);
        const int nn[3] = {q->knox,q->knoy,q->knoz};
        int last = -1;      // the fine patch of the last face (the next face mostly in the same one)

        for(int i=0; i<q->knox; ++i)
        for(int j=0; j<q->knoy; ++j)
        for(int k=0; k<q->knoz; ++k)
        {
            const int I[3] = {i,j,k};
            const int G[3] = {i+og[0],j+og[1],k+og[2]};
            if(row[cix(q,i,j,k)]<0 || !covered(lc,G))
            continue;
            for(int d=0; d<3; ++d)
            {
                if(rr[d]!=2 || I[d]+1>=nn[d])
                continue;
                int J[3] = {i,j,k}, GJ[3] = {G[0],G[1],G[2]};
                J[d] += 1;
                GJ[d] += 1;
                if(row[cix(q,J[0],J[1],J[2])]<0 || !covered(lc,GJ))
                continue;

                // the fine faces at the coarse face: the + faces of the fine cells of the last layer of I
                double inv = 0.0;
                int nf = 0;
                bool ok = true;
                int a[3], b[3];
                for(int e=0; e<3; ++e)
                {
                    a[e] = G[e]*rr[e];
                    b[e] = a[e]+rr[e]-1;
                }
                a[d] = b[d];
                int F[3];
                for(F[0]=a[0]; F[0]<=b[0] && ok; ++F[0])
                for(F[1]=a[1]; F[1]<=b[1] && ok; ++F[1])
                for(F[2]=a[2]; F[2]<=b[2] && ok; ++F[2])
                {
                    auto holds = [&](int c)
                    {
                        const r3patch *P2 = P[c];
                        return F[0]>=P2->lo[0] && F[0]<=P2->hi[0] && F[1]>=P2->lo[1] && F[1]<=P2->hi[1] && F[2]>=P2->lo[2] && F[2]<=P2->hi[2];
                    };
                    int id = (last>=0 && holds(last)) ? last : -1;
                    if(id<0)
                    for(int c : lev[lc+1])
                    if(holds(c))
                    {
                        id = c;
                        break;
                    }
                    last = id;
                    if(id<0)
                    {
                        ok = false;
                        break;
                    }
                    const int Fl[3] = {F[0]-P[id]->lo[0],F[1]-P[id]->lo[1],F[2]-P[id]->lo[2]};
                    inv += 1.0/face_rho(id,Fl,d);
                    ++nf;
                }
                if(!ok || nf==0)
                continue;
                inv /= double(nf);
                const double f = face_rho(g,I,d)*inv;      // the coupling scales with 1/rho

                const long qi = (g<0) ? L.idx(i,j,k) : L.idx(i+ol[0],j+ol[1],k+ol[2]);
                const long qj = (g<0) ? L.idx(J[0],J[1],J[2]) : L.idx(J[0]+ol[0],J[1]+ol[1],J[2]+ol[2]);
                const double ci = getc(U,L,qi,mplus[d]), cj = getc(U,L,qj,mminus[d]);
                setc(U,L,qi,mplus[d],ci*f);
                setc(U,L,qi,0,getc(U,L,qi,0)+ci-ci*f);
                setc(U,L,qj,mminus[d],cj*f);
                setc(U,L,qj,0,getc(U,L,qj,0)+cj-cj*f);
            }
        }
    }
}

// G 40: the level-0 rows of the coarse cells at the edges of the level-1 patches in the shadows:
// their coupling to the covered cell across the face folded into the diagonal (as pr_matrix for a
// local coarse cell; the composite operator adds the fluxes of the fine faces, pr_apply)
void cfd_amr::sh_folds()
{
    matrix_diag &M = a0->M;
    for(size_t r=0; r<pu_key.rq.size(); ++r)
    for(size_t q=0; q<pu_key.rq[r].size(); ++q)
    {
        const int n = row0[pu_key.rq[r][q]];
        if(n<0)
        continue;
        const int t = pu_key.rt[r][q];
        if(t%3==2)
        continue;
        const int d = t/6, dir = ((t/3)%2) ? 1 : -1;
        double &ccn = mcoef(M,n,d,-dir);
        M.p[n] += ccn;
        ccn = 0.0;
    }
}

// G 40: the couplings of the level-0 multigrid at the edges of the level-1 patches, from the sums of
// kc of the shadows (collective)
void cfd_amr::sh_ks()
{
    cplan_xfer(pu_key,1,[this](int, int g) -> double* { return kval[g-npat].data(); });
    const int n[3] = {p0->knox,p0->knoy,p0->knoz};
    const int jk = p0->jmax*p0->kmax;
    for(size_t r=0; r<pu_key.rq.size(); ++r)
    for(size_t q=0; q<pu_key.rq[r].size(); ++q)
    {
        const int qq = pu_key.rq[r][q];
        const int ci[3] = {qq/jk + p0->imin, (qq%jk)/p0->kmax + p0->jmin, qq%p0->kmax + p0->kmin};
        const int t = pu_key.rt[r][q];
        const int d = t/6, dir = ((t/3)%2) ? 1 : -1, h = t%3;
        if(h==0)
        ks_apply(mg0,ci,d,dir,n,pu_key.rb[r][q]);
        else
        ks_half(mg0,ci,d,dir,h,n,pu_key.rb[r][q]);
    }
}

// the cells around the patches: the vector of the patch of the same level, or the value of the
// coarser cell (the patch rows couple to the coarse cell across a patch edge); the level-0 halo
void cfd_amr::pr_sync(int k)
{
    vector<fspec> fs;
    fs.push_back({CELL,2,[this,k](int g) -> field& { return pfield(g,k); }});

    // level projection (G 7 1): the patches of the level next to each other
    if(prl>0)
    {
        fill(prl,fs);
        return;
    }

    pr_halo0(pfield(-1,k));

    // (the operator reaches one cell around a patch: the cells next to it)
    fill_all(fs,1);
}

// the halo of a level-0 vector: the neighbour ranks, the periodic edges
void cfd_amr::pr_halo0(field &f)
{
    if(p0->mpi_size>1)
    pgc0->gcparax(p0,f,4);

    if(p0->periodic1>0 || p0->periodic2>0)
    {
        if(p0->mpi_size>1)
        pgc0->gcperiodicx(p0,f,4);
        if(p0->periodic1==1)
        pgc0->gc_periodic(p0,f,4,1);
        if(p0->periodic2==1 && p0->j_dir==1)
        pgc0->gc_periodic(p0,f,4,2);
    }
}

// y = A x on the leaf cells
void cfd_amr::pr_apply(int kx, int ky)
{
    pr_sync(kx);

    for(int g : prg)
    {
        lexer *q = glex(g);
        const matrix_diag &M = gfd(g)->M;
        const double *x = pvec(g,kx);
        double *y = pvec(g,ky);
        const int sI = q->jmax*q->kmax;
        const int sJ = q->kmax;
        const vector<int> &LQ = UQ(g), &LR = UR(g);

        const double *Mp=M.p.data(), *Mn=M.n.data(), *Ms=M.s.data(), *Mw=M.w.data(), *Me=M.e.data(), *Mt=M.t.data(), *Mb=M.b.data();

        for(size_t n=0; n<LQ.size(); ++n)
        {
            const int qq = LQ[n], r = LR[n];
            y[qq] = Mp[r]*x[qq] + Mn[r]*x[qq+sI] + Ms[r]*x[qq-sI] + Mw[r]*x[qq+sJ] + Me[r]*x[qq-sJ] + Mt[r]*x[qq+1] + Mb[r]*x[qq-1];
        }

        if(prl<0)
        for(int qq : cq[g+1])
        y[qq] = 0.0;
    }

    // the level projection has no couplings across the patch edges (prescribed flux)
    if(prl>0)
    return;

    // G 40: the fluxes into the coarse cells of the shadows go to the ranks of those cells
    const std::function<field&(int)> fy = [this,ky](int g) -> field& { return pfield(g,ky); };
    if(place)
    sh_zero(pu_cc,fy);

    for(const cfd_amr_cf &F : cf)
    {
        const double xs = pr_xs(F,kx);
        const double xf = pvec(F.pid,kx)[F.qf];
        pvec(F.pid,ky)[F.qf] -= F.kf*xs;
        pvec(F.gc,ky)[F.qc] += F.kc*(xs-xf);
    }

    if(place)
    {
        vector<std::function<field&(int)>> f(1,fy);
        sh_run(pu_cc,f,true);
    }
}

double cfd_amr::pr_dot(int ka, int kb)
{
    double s=0.0;
    for(int g : prg)
    {
        const double *a = pvec(g,ka), *b = pvec(g,kb);
        for(int qq : UQ(g))
        s += a[qq]*b[qq];
    }
    MPI_Allreduce(MPI_IN_PLACE,&s,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return s;
}

void cfd_amr::pr_start()
{
    for(int g : prg)
    {
        const double *vb = pvec(g,NS);
        double *vr = pvec(g,NR), *vrh = pvec(g,NRH), *vp = pvec(g,NPV), *vv = pvec(g,NVV);
        for(int qq : UQ(g))
        {
            double r = vb[qq] - vv[qq];
            vr[qq] = r; vrh[qq] = r; vp[qq] = 0.0; vv[qq] = 0.0;
        }
    }
}

void cfd_amr::pr_p(double beta, double om)
{
    for(int g : prg)
    {
        const double *vr = pvec(g,NR), *vv = pvec(g,NVV);
        double *vp = pvec(g,NPV);
        for(int qq : UQ(g))
        vp[qq] = vr[qq] + beta*(vp[qq] - om*vv[qq]);
    }
}

void cfd_amr::pr_s(double alp)
{
    for(int g : prg)
    {
        const double *vr = pvec(g,NR), *vv = pvec(g,NVV);
        double *vs = pvec(g,NS);
        for(int qq : UQ(g))
        vs[qq] = vr[qq] - alp*vv[qq];
    }
}

void cfd_amr::pr_x(double alp, double om)
{
    for(int g : prg)
    {
        const double *vph = pvec(g,NPH), *vsh = pvec(g,NSH), *vs = pvec(g,NS), *vt = pvec(g,NT);
        double *x = pvec(g,-1), *vr = pvec(g,NR);
        for(int qq : UQ(g))
        {
            x[qq] += alp*vph[qq] + om*vsh[qq];
            vr[qq] = vs[qq] - om*vt[qq];
        }
    }
}

// the composite residual r - A z into kt: leaf cells, then the cells under the patches of the levels
// above lmin from their children (finest first)
void cfd_amr::pr_residual(int kr, int kz, int kt, int lmin)
{
    pr_apply(kz,kt);
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *r = pvec(g,kr);
        double *t = pvec(g,kt);
        for(int qq : lq[g+1])
        t[qq] = r[qq] - t[qq];
    }
    vector<std::function<field&(int)>> cells, none;
    cells.push_back([this,kt](int g) -> field& { return pfield(g,kt); });
    for(int l=maxlev; l>lmin; --l)
    restrict_level(l,cells,none);
}

// the half width of the interface band of grid g (psi of its level)
double cfd_amr::psi_band(int g)
{
    return (g<0) ? p0->psi : glex(g)->psi;
}

// z = M^-1 r: one FAC V-cycle over the levels.  Down: on every level, finest first, a patch-local
// V-cycle on the composite residual (the cells under a finer patch take the residual of the finer
// level); level 0: a V-cycle on the residual left; up: the correction of the coarser level prolonged
// into the patches (not in the interface band of the coarser level), a patch-local V-cycle on the
// residual left.
void cfd_amr::pr_prec(int kr, int kz)
{
    const int KT = REEFAMR_NTMP;        // residual
    const int KC = REEFAMR_NVEC;        // correction of a level, prolonged to the next finer one

    // level projection (G 7 1): one patch-local V-cycle on every patch of the level
    if(prl>0)
    {
        for(int g : prg)
        {
            lexer *q = glex(g);
            std::fill(pvec(g,kz),pvec(g,kz)+(size_t)q->imax*q->jmax*q->kmax,0.0);
        }
        pr_level_solve(prl,kr,kz);
        return;
    }

    // (a patch: its cells and the layer around it, the cells the sweep reads and writes)
    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        double *z = pvec(g,kz), *c = pvec(g,KC);
        if(g<0 || g>=npat)
        {
            const size_t n = (size_t)q->imax*q->jmax*q->kmax;
            std::fill(z,z+n,0.0);
            std::fill(c,c+n,0.0);
            continue;
        }
        for(int i=-1; i<=q->knox; ++i)
        for(int j=-1; j<=q->knoy; ++j)
        {
            const size_t a = cix(q,i,j,-1);
            std::fill(z+a,z+a+q->knoz+2,0.0);
            std::fill(c+a,c+a+q->knoz+2,0.0);
        }
    }

    // a level: its multigrid over the patches together (several patches on the rank), else a
    // patch-local V-cycle on each
    auto level_solve = [&](int l, int kt) { pr_level_solve(l,kt,kz); };

    vector<fspec> fc;
    fc.push_back({CELL,2,[this](int g) -> field& { return pfield(g,REEFAMR_NVEC); }});
    const fspec cl = {CELL,0,[this](int g) -> field& { return pfield(g,REEFAMR_NVEC); }};

    // down
    for(int l=maxlev; l>=1; --l)
    {
        if(l==maxlev)
        {
            // z = 0: the residual is r
            level_solve(l,kr);
        }
        else
        {
            pr_residual(kr,kz,KT,l);
            level_solve(l,KT);
        }
    }

    // level 0
    pr_residual(kr,kz,KT,0);
    {
        sc_level &L = mg0->fine();
        double *t = pvec(-1,KT);
        double *z = pvec(-1,kz);
        double *c = pvec(-1,KC);
        std::fill(L.u.begin(),L.u.end(),0.0);
        std::fill(L.f.begin(),L.f.end(),0.0);
        for(size_t n=0; n<l0a.size(); ++n)
        L.f[l0a[n]] = t[l0f[n]];

        const double t0 = MPI_Wtime();
        mg0->vcycle(0,1,1);
        tp[1] += MPI_Wtime()-t0;

        for(size_t n=0; n<l0a.size(); ++n)
        {
            z[l0f[n]] += L.u[l0a[n]];
            c[l0f[n]] = L.u[l0a[n]];
        }

        pr_halo0(pfield(-1,KC));
    }

    // G 40: the correction of level 0 in the shadows
    if(place)
    {
        vector<std::function<field&(int)>> f(1,[this](int g) -> field& { return pfield(g,REEFAMR_NVEC); });
        sh_down_f(f);
    }

    // up
    for(int l=1; l<=maxlev; ++l)
    {
        for(int id : lev[l])
        {
            cfd_amr_patch *c = CP(id);
            field &z = pfield(id,kz);
            field &t = pfield(id,KC);
            for(const r3patch::pblock &B : c->par)
            {
                int og[3];
                goff(B.g,og);
                int Pg[3];
                for(Pg[0]=B.lo[0]; Pg[0]<=B.hi[0]; ++Pg[0])
                for(Pg[1]=B.lo[1]; Pg[1]<=B.hi[1]; ++Pg[1])
                for(Pg[2]=B.lo[2]; Pg[2]<=B.hi[2]; ++Pg[2])
                {
                    const int sp[3] = {Pg[0]-og[0],Pg[1]-og[1],Pg[2]-og[2]};
                    for(int a=0; a<rr[0]; ++a)
                    for(int b=0; b<rr[1]; ++b)
                    for(int e=0; e<rr[2]; ++e)
                    {
                        const int o[3] = {a,b,e};
                        const int i = Pg[0]*rr[0]+a-c->lo[0], j = Pg[1]*rr[1]+b-c->lo[1], k = Pg[2]*rr[2]+e-c->lo[2];
                        if(c->row[cix(c->pp,i,j,k)]<0)
                        continue;
                        // not inside the interface band of the coarser level: there its density is
                        // smeared over its own band (psi of each level), the fine cell is water or air,
                        // and the coarse correction does not fit the fine rows (with 4 levels and a
                        // free surface in the patches the composite BiCGStab diverged, 250 iterations);
                        // the level solve below takes these cells from the residual alone
                        if(p0->F30>0 && fabs(c->a->phi(i,j,k))<psi_band(B.g))
                        continue;
                        const double v = prolong(cl,B.g,sp,o);
                        t(i,j,k) = v;
                        z(i,j,k) += v;
                    }
                }
            }
        }

        pr_residual(kr,kz,KT,l);
        level_solve(l,KT);

        // the correction of this level for the next finer one: all of z (zero before the cycle;
        // the part of the down sweep is the coarse correction of the finer level's residual)
        if(l<maxlev)
        {
            for(int id : lev[l])
            {
                lexer *q = CP(id)->pp;
                std::copy(pvec(id,kz),pvec(id,kz)+q->imax*q->jmax*q->kmax,pvec(id,KC));
            }
            fill(l,fc);
        }
    }
}

// G 40: the overlap of every patch multigrid: OV cells (the cells around the patch) across a side of
// the patch where it meets a patch of its level (copied from that patch, kind 0 or 2).  The
// preconditioner is then a restricted additive Schwarz across the cuts between the patches (the
// right-hand side and the rows in the overlap from the patch there, the correction kept inside
// the patch): the patch-local V-cycles with zero correction at the cut (block Jacobi) took 60
// BiCGStab iterations with four units side by side in the dam break, 11 with an overlap of 2
// cells, 9 with 3 (G 40 0 on four ranks: 10)
void cfd_amr::pr_overlap()
{
    const int OV = 3;
    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        for(int d=0; d<3; ++d)
        c->olo[d] = c->ohi[d] = 0;
        if(!place && !(pr_ras() && lmg_of[id]<0))
        continue;
        const int n[3] = {pp->knox,pp->knoy,pp->knoz};
        for(int d=0; d<3; ++d)
        {
            if(rr[d]!=2 || pr_wrap(id,d))
            continue;
            for(int side=0; side<2; ++side)
            {
                bool sib = false;
                int a[3] = {0,0,0}, b[3] = {n[0]-1,n[1]-1,n[2]-1};
                a[d] = b[d] = (side==0) ? -1 : n[d];
                int I[3];
                for(I[0]=a[0]; I[0]<=b[0] && !sib; ++I[0])
                for(I[1]=a[1]; I[1]<=b[1] && !sib; ++I[1])
                for(I[2]=a[2]; I[2]<=b[2] && !sib; ++I[2])
                {
                    const int k = c->mkind[cix(pp,I[0],I[1],I[2])];
                    if(k==0 || k==2)
                    sib = true;
                }
                if(sib)
                (side==0 ? c->olo[d] : c->ohi[d]) = OV;
            }
        }
    }
}

// G 40: the matrix rows of the patches (level l, all for l < 0) in their kc7 fields, filled into the
// cells around the patches next to them (collective)
void cfd_amr::pr_coef_fill(int l)
{
    for(int lv=1; lv<=maxlev; ++lv)
    {
        if(l>=0 && lv!=l)
        continue;
        for(int id : lev[lv])
        {
            cfd_amr_patch *c = CP(id);
            lexer *pp = c->pp;
            const matrix_diag &M = c->a->M;
            for(field4 *f : c->kc7)
            std::fill(f->V,f->V+(size_t)pp->imax*pp->jmax*pp->kmax,0.0);
            for(int i=0; i<pp->knox; ++i)
            for(int j=0; j<pp->knoy; ++j)
            for(int k=0; k<pp->knoz; ++k)
            {
                const int qq = cix(pp,i,j,k);
                const int r = c->row[qq];
                if(r<0)
                continue;
                c->kc7[0]->V[qq] = M.p[r]; c->kc7[1]->V[qq] = M.n[r]; c->kc7[2]->V[qq] = M.s[r];
                c->kc7[3]->V[qq] = M.w[r]; c->kc7[4]->V[qq] = M.e[r]; c->kc7[5]->V[qq] = M.t[r]; c->kc7[6]->V[qq] = M.b[r];
            }
        }
        vector<fspec> fs;
        for(int m=0; m<7; ++m)
        fs.push_back({CELL,2,[this,m](int g) -> field& { return *CP(g)->kc7[m]; }});
        fill(lv,fs,true);
    }
}

// a patch over the whole length of a periodic direction d (x, y): its multigrid is periodic there
// (no overlap there), its margin across the edge is the patch itself
bool cfd_amr::pr_wrap(int id, int d) const
{
    const r3patch *c = P[id];
    if(d>1 || !per[d])
    return false;
    return c->lo[d]==0 && c->hi[d]==gn(c->lev,d)-1;
}

void cfd_amr::pr_patch_mg(int id, double dshift)
{
    cfd_amr_patch *c = CP(id);
    lexer *pp = c->pp;
    sc_level &L = c->mg->fine();
    pr_clear(L);
    const matrix_diag &M = c->a->M;
    const int nx=pp->knox, ny=pp->knoy, nz=pp->knoz;
    const int *ol = c->olo, *oh = c->ohi;

    if(ol[0]+ol[1]+ol[2]+oh[0]+oh[1]+oh[2]==0)
    {
        const bool wx = pr_wrap(id,0), wy = pr_wrap(id,1);   // (periodic over the whole patch)
        for(int i=0; i<nx; ++i)
        for(int j=0; j<ny; ++j)
        for(int k=0; k<nz; ++k)
        {
            const long q = L.idx(i,j,k);
            const int r = c->row[cix(pp,i,j,k)];
            if(r<0)
            {
                pr_coef(c->mg,L,q,1.0,0.0,0.0,0.0,0.0,0.0,0.0);
                continue;
            }
            pr_coef(c->mg,L,q,(dshift==1.0) ? M.p[r] : M.p[r]*dshift,
                    (i+1<nx || wx) ? M.n[r] : 0.0, (i>0 || wx) ? M.s[r] : 0.0,
                    (j+1<ny || wy) ? M.w[r] : 0.0, (j>0 || wy) ? M.e[r] : 0.0,
                    (k+1<nz) ? M.t[r] : 0.0, (k>0) ? M.b[r] : 0.0);
            L.act[q] = 1;
        }
        return;
    }

    // with the overlap: the cells of the patch and those of the patches next to it in the box
    auto in = [&](int i, int j, int k) { return i>=0 && i<nx && j>=0 && j<ny && k>=0 && k<nz; };
    auto row_o = [&](int i, int j, int k) -> bool     // an overlap cell with a row of the patch there
    {
        if(i<-ol[0] || i>=nx+oh[0] || j<-ol[1] || j>=ny+oh[1] || k<-ol[2] || k>=nz+oh[2] || in(i,j,k))
        return false;
        const int qq = cix(pp,i,j,k);
        return (c->mkind[qq]==0 || c->mkind[qq]==2) && pp->flag4[qq]>0;
    };
    // (a direction over which the patch is periodic over itself: no overlap there, the couplings
    // across the edge kept)
    const bool wx = pr_wrap(id,0), wy = pr_wrap(id,1);
    auto keep = [&](int i, int j, int k)
    {
        if(wx) i = (i+nx)%nx;
        if(wy) j = (j+ny)%ny;
        return in(i,j,k) || row_o(i,j,k);
    };

    for(int i=-ol[0]; i<nx+oh[0]; ++i)
    for(int j=-ol[1]; j<ny+oh[1]; ++j)
    for(int k=-ol[2]; k<nz+oh[2]; ++k)
    {
        const long q = L.idx(i+ol[0],j+ol[1],k+ol[2]);
        const int qq = cix(pp,i,j,k);
        double cc[7];
        if(in(i,j,k))
        {
            const int r = c->row[qq];
            if(r<0)
            {
                pr_coef(c->mg,L,q,1.0,0.0,0.0,0.0,0.0,0.0,0.0);
                continue;
            }
            cc[0]=M.p[r]; cc[1]=M.n[r]; cc[2]=M.s[r]; cc[3]=M.w[r]; cc[4]=M.e[r]; cc[5]=M.t[r]; cc[6]=M.b[r];
        }
        else if(row_o(i,j,k))
        {
            for(int m=0; m<7; ++m)
            cc[m] = c->kc7[m]->V[qq];
        }
        else
        {
            pr_coef(c->mg,L,q,1.0,0.0,0.0,0.0,0.0,0.0,0.0);
            continue;
        }
        pr_coef(c->mg,L,q,cc[0]*dshift,
                keep(i+1,j,k) ? cc[1] : 0.0, keep(i-1,j,k) ? cc[2] : 0.0,
                keep(i,j+1,k) ? cc[3] : 0.0, keep(i,j-1,k) ? cc[4] : 0.0,
                keep(i,j,k+1) ? cc[5] : 0.0, keep(i,j,k-1) ? cc[6] : 0.0);
        L.act[q] = 1;
    }
}

// patch-local V-cycle with the right-hand side kt on the patch rows, added to z (and to kc)
void cfd_amr::pr_patch_solve(int id, int kt, int kz, int kc)
{
    const double t0 = MPI_Wtime();
    cfd_amr_patch *c = CP(id);
    lexer *pp = c->pp;
    const double *t = pvec(id,kt);
    double *z = pvec(id,kz);
    double *cc = (kc>=0) ? pvec(id,kc) : nullptr;
    sc_level &L = c->mg->fine();
    const int *ol = c->olo, *oh = c->ohi;

    std::fill(L.u.begin(),L.u.end(),0.0);
    std::fill(L.f.begin(),L.f.end(),0.0);
    for(int i=0; i<pp->knox; ++i)
    for(int j=0; j<pp->knoy; ++j)
    for(int k=0; k<pp->knoz; ++k)
    {
        const int qq = cix(pp,i,j,k);
        if(c->row[qq]>=0)
        L.f[L.idx(i+ol[0],j+ol[1],k+ol[2])] = t[qq];
    }
    // G 40: the overlap (the cells of the patches next to it, their right-hand side filled)
    if(ol[0]+ol[1]+ol[2]+oh[0]+oh[1]+oh[2]>0)
    for(int i=-ol[0]; i<pp->knox+oh[0]; ++i)
    for(int j=-ol[1]; j<pp->knoy+oh[1]; ++j)
    for(int k=-ol[2]; k<pp->knoz+oh[2]; ++k)
    {
        if(i>=0 && i<pp->knox && j>=0 && j<pp->knoy && k>=0 && k<pp->knoz)
        continue;
        const int qq = cix(pp,i,j,k);
        if((c->mkind[qq]==0 || c->mkind[qq]==2) && pp->flag4[qq]>0)
        L.f[L.idx(i+ol[0],j+ol[1],k+ol[2])] = t[qq];
    }

    c->mg->vcycle(0,2,2);

    for(int i=0; i<pp->knox; ++i)
    for(int j=0; j<pp->knoy; ++j)
    for(int k=0; k<pp->knoz; ++k)
    {
        const int qq = cix(pp,i,j,k);
        if(c->row[qq]<0)
        continue;
        const double u = L.u[L.idx(i+ol[0],j+ol[1],k+ol[2])];
        z[qq] += u;
        if(cc)
        cc[qq] += u;
    }
    tp[0] += MPI_Wtime()-t0;
}

namespace
{
struct pr_space
{
    static const int B=NS, R=NR, RH=NRH, PV=NPV, VV=NVV, S=NS, T=NT, PH=NPH, SH=NSH;
    cfd_amr *a;
    void apply(int x, int y) { a->pr_apply(x,y); }
    void prec(int r, int z) { a->pr_prec(r,z); }
    double dot(int x, int y) { return a->pr_dot(x,y); }
    void op_start() { a->pr_start(); }
    void op_p(double beta, double om) { a->pr_p(beta,om); }
    void op_s(double alp) { a->pr_s(alp); }
    void op_x(double alp, double om) { a->pr_x(alp,om); }
};

// patch work: no exchange, rank-local reductions, the fdm of the patch in the ghostcell object
struct pscope2
{
    ghostcell *g;
    fdm *back;
    bool oldc, oldl;
    pscope2(ghostcell *gg, fdm *a, fdm *a0) : g(gg), back(a0)
    {
        oldc = g->set_comms(false);
        oldl = g->set_local(true);
        g->fdm_update(a);
    }
    ~pscope2()
    {
        g->fdm_update(back);
        g->set_comms(oldc);
        g->set_local(oldl);
    }
};
}

// the faces at the patch edges after the solve: corrected with the flux of the composite rows
void cfd_amr::pr_correct(int s)
{
    const double alpha = mom0->amr_alpha(s);

    for(const cfd_amr_cf &F : cf)
    {
        const double xf = pvec(F.pid,-1)[F.qf];
        const double xc = pr_xs(F,-1);
        const double rof = 1.0/(F.kf*F.dist*F.dxnf);
        field &ff = CP(F.pid)->pmom->amr_velout(CP(F.pid)->a,F.d,s);
        ff(F.ff[0],F.ff[1],F.ff[2]) = F.ustar - double(F.dir)*alpha*p0->dt*(xc-xf)/(F.dist*rof);
    }
}

// the projection of stage s on all grids
void cfd_amr::pr_project(lexer *p, ghostcell *pgc, int s)
{
    const double alpha = mom0->amr_alpha(s);
    const double t0 = MPI_Wtime();

    // the stage velocities of the patches of the same level next to each other, then the faces at
    // the patch edges
    fill_vel_after(s);
    pr_predict(s);

    field &u0 = mom0->amr_velout(a0,0,s), &v0 = mom0->amr_velout(a0,1,s), &w0 = mom0->amr_velout(a0,2,s);

    pflow0->pressure_io(p,a0,pgc);
    press0->amr_prepare(p,a0,pois0,pgc,u0,v0,w0,alpha);

    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        pscope2 ps(pgc,c->a,a0);
        momentum_rk *m = c->pmom;
        c->pflow->pressure_io(c->pp,c->a,pgc);
        c->ppress->amr_prepare(c->pp,c->a,c->ppois,pgc,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s),alpha);
    }

    pr_matrix(s);

    // right-hand side into NS, A x0 into NVV (x0 = 0: pcorr cleared by pjm_corr::rhs)
    for(int g=-1; g<npat; ++g)
    {
        const double *R = gfd(g)->rhsvec.V.data();
        double *b = pvec(g,NS);
        for(size_t n=0; n<lq[g+1].size(); ++n)
        b[lq[g+1][n]] = R[lr[g+1][n]];
        for(int qq : cq[g+1])
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

    // pressure and velocities of every grid
    press0->amr_finish(p,a0,pgc,u0,v0,w0,alpha);

    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        pscope2 ps(pgc,c->a,a0);
        momentum_rk *m = c->pmom;
        c->ppress->amr_finish(c->pp,c->a,pgc,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s),alpha);
    }

    fill_vel_after(s);
    pr_correct(s);

    mom0->amr_project_after(p,a0,pgc,u0,v0,w0);
    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        pscope2 ps(pgc,c->a,a0);
        momentum_rk *m = c->pmom;
        patch_relax_vel(id,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s));
        m->amr_project_after(c->pp,c->a,pgc,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s));
    }

    restrict_vel(s);

    p->poissontime = MPI_Wtime()-t0;
}

// --------------------------------------------------------------------- level projection (G 7 1)
// the connected components of the patches of every level (patches of a level that share a face)
void cfd_amr::pr_components()
{
    const int ng = (int)GP.size();
    vector<int> up(ng);
    for(int g=0; g<ng; ++g)
    up[g] = g;
    std::function<int(int)> root = [&](int g) { while(up[g]!=g) { up[g] = up[up[g]]; g = up[g]; } return g; };

    for(int a=0; a<ng; ++a)
    for(int b=a+1; b<ng; ++b)
    {
        const r3gpatch &A = GP[a], &B = GP[b];
        if(A.lev!=B.lev)
        continue;
        for(int d=0; d<3; ++d)
        {
            if(!(A.hi[d]+1==B.lo[d] || B.hi[d]+1==A.lo[d]))
            continue;
            bool ov = true;
            for(int e=0; e<3; ++e)
            if(e!=d && (A.lo[e]>B.hi[e] || B.lo[e]>A.hi[e]))
            ov = false;
            if(ov)
            up[root(a)] = root(b);
        }
    }

    vector<int> num(ng,-1);
    ncomp = 0;
    for(int g=0; g<ng; ++g)
    {
        const int r = root(g);
        if(num[r]<0)
        num[r] = ncomp++;
    }
    for(int id=0; id<npat; ++id)
    CP(id)->comp = num[root(P[id]->gid)];
}

// the projection of stage s on the patches of level l alone: the faces at the patch edges keep the
// face velocity of the coarser level at the time of the stage output (in tsf, the transverse
// distribution linear as the prolongation), so the patches of a level take the flux of the coarser
// level and the composite field stays divergence free.  Neumann at the patch edges: a component of
// the patches without a boundary of fixed pressure has a singular problem, the right-hand side is
// made compatible (the round-off of the coarser solve) and the correction gets the mean 0.
void cfd_amr::pr_level(lexer *p, ghostcell *pgc, int l, int s)
{
    const double alpha = mom0->amr_alpha(s);
    const double t0 = MPI_Wtime();

    auto vel_after = [&]()
    {
        vector<fspec> fs;
        for(int c=0; c<3; ++c)
        fs.push_back({c,1,[this,s,c,l](int g) -> field& { return glevel(g)==l ? gmom(g)->amr_velout(gfd(g),c,s) : tsf(g,c); },
                          [this,s,c](int g) -> field& { return gmom(g)->amr_velout(gfd(g),c,s); }});
        fill(l,fs);
    };

    vel_after();

    // the faces at the patch edges: the coarse face velocity at the output time of the stage
    // (tsf), its pressure gradient replaced by the one across the fine face as in pr_predict: the
    // fine cell and the coarse pressure at its position, each without the discrete hydrostatic
    // pressure of its grid (the free surface across a patch edge; the coarse one with its level set
    // at the time of the stage), the face density of the stage.
    // The correction of the coarser level is final (Dirichlet 0): the fine pressure answers the
    // coarse one across the edge, the patch takes its own flux there (the coarser level takes it
    // in the synchronisation at the end of the level-0 step)
    // G 40: the coarse data of the level-1 edges in the shadows (level 0 at the time of the stage)
    if(place && l==1)
    {
        vector<std::function<field&(int)>> f;
        f.push_back([this](int g) -> field& { return gfd(g)->phi; });
        vector<int> ms = {0,1,2,tl_press,tl_phi};
        if(scal)
        ms.push_back(tl_s);
        for(int m : ms)
        f.push_back([this,m](int g) -> field& { return tsf(g,m); });
        sh_down_f(f);
    }

    pr_weights();
    pr_hydro(l);
    for(cfd_amr_cf &F : cf)
    {
        cfd_amr_patch *c = CP(F.pid);
        if(c->lev!=l)
        continue;
        field &fc = tsf(F.gface,F.d);
        double v = fc(F.cf[0],F.cf[1],F.cf[2]);
        for(int t=0; t<3; ++t)
        if(t!=F.d)
        v += cell_slope(fc,F.gface,F.cf,t,1)*(F.o[t]==0 ? -0.25 : 0.25);

        lexer *pp = c->pp;
        lexer *qc = glex(F.gc), *qp = glex(F.gp);
        const double pC = tsf(F.gc,tl_press).V[F.qc];
        const double pP = tsf(F.gp,tl_press).V[F.qp];
        double pCs = -F.hs;
        for(int m=0; m<F.ns; ++m)
        pCs += F.sw[m]*tsf(F.sg[m],tl_press).V[F.sq[m]];
        const double pf = c->a->press.V[F.qf] - F.hf;

        field &lc = tsf(F.gc,tl_phi), &lp = tsf(F.gp,tl_phi);
        const double phC = lc(F.ci[0],F.ci[1],F.ci[2]);
        const double phP = lp(F.pi[0],F.pi[1],F.pi[2]);
        const double wc = 0.5*(interface_width(qc,lc,F.ci[0],F.ci[1],F.ci[2]) + interface_width(qp,lp,F.pi[0],F.pi[1],F.pi[2]));
        const double Hc = heaviside(0.5*(phC+phP),wc);
        const double fpc = porfac(qc,0.5*(gfd(F.gc)->porosity.V[F.qc] + gfd(F.gp)->porosity.V[F.qp]));
        const double roc = rmix(qc,Hc,scal ? 0.5*(tsf(F.gc,tl_s).V[F.qc] + tsf(F.gp,tl_s).V[F.qp]) : 0.0)/fpc;

        lexer *ql = (F.dir>0) ? qp : qc;
        const double *dp[3] = {ql->DXP,ql->DYP,ql->DZP};
        const int Lc = (F.dir>0) ? F.pi[F.d] : F.ci[F.d];
        const double Gc = double(F.dir)*(pC-pP)/(dp[F.d][Lc+marge]*roc);
        const double Gf = double(F.dir)*(pCs-pf)/(F.dist*F.rof);

        v += alpha*pp->dt*(Gc-Gf);
        if(F.hok && F.d==2)
        v -= alpha*pp->dt*p0->W22*fpc;

        F.ustar = v;
        c->pmom->amr_velout(c->a,F.d,s)(F.ff[0],F.ff[1],F.ff[2]) = v;
    }

    for(int id : lev[l])
    {
        cfd_amr_patch *c = CP(id);
        pscope2 ps(pgc,c->a,a0);
        momentum_rk *m = c->pmom;
        c->pflow->pressure_io(c->pp,c->a,pgc);
        c->ppress->amr_prepare(c->pp,c->a,c->ppois,pgc,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s),alpha);
    }

    // the faces at the patch edges Dirichlet: the coupling of the fine cell to the coarse value at
    // its position stays on the diagonal, the coarse correction is 0; multigrids from the patch rows
    // (a small shift of the diagonal keeps a closed patch, one without such faces and without a
    // boundary of fixed pressure, regular)
    for(cfd_amr_cf &F : cf)
    {
        if(CP(F.pid)->lev!=l)
        continue;
        F.kf = 1.0/(F.rof*F.dist*F.dxnf);
        matrix_diag &Mf = CP(F.pid)->a->M;
        double &cfn = mcoef(Mf,F.nf,F.d,F.dir);
        Mf.p[F.nf] += cfn + F.kf;
        cfn = 0.0;
    }
    if(pr_ras())
    pr_coef_fill(l);
    for(int id : lev[l])
    {
        if(lmg_of[id]>=0)
        continue;
        pr_patch_mg(id,1.0+1.0e-6);
        CP(id)->mg->coarsen();
    }
    for(int m=0; m<(int)lmg.size(); ++m)
    if(P[lmg[m].ids[0]]->lev==l)
    pr_level_mg_fill(m);

    prl = l;
    prg = lev[l];

    // components with a boundary of fixed pressure (a coupling to a cell outside the domain)
    vector<double> cs(3*ncomp,0.0);     // sum b V, sum V, fixed-pressure faces
    for(int id : lev[l])
    {
        cfd_amr_patch *c = CP(id);
        lexer *q = c->pp;
        const matrix_diag &M = c->a->M;
        const int sI = q->jmax*q->kmax, sJ = q->kmax;
        for(size_t n=0; n<aq[id+1].size(); ++n)
        {
            const int qq = aq[id+1][n], r = ar[id+1][n];
            const double cc[6] = {M.n[r],M.s[r],M.w[r],M.e[r],M.t[r],M.b[r]};
            const int nb[6] = {qq+sI,qq-sI,qq+sJ,qq-sJ,qq+1,qq-1};
            for(int m=0; m<6; ++m)
            if(cc[m]!=0.0 && q->flag4[nb[m]]<0)
            cs[3*c->comp+2] += 1.0;
        }
    }
    for(const cfd_amr_cf &F : cf)
    if(CP(F.pid)->lev==l)
    cs[3*CP(F.pid)->comp+2] += 1.0;

    // right-hand side, compatible on the closed components
    auto cvol = [&](lexer *q, int qq) -> double
    {
        const int jk = q->jmax*q->kmax;
        const int i = qq/jk + q->imin, r = qq%jk;
        const int j = r/q->kmax + q->jmin, k = r%q->kmax + q->kmin;
        return q->DXN[i+marge]*q->DYN[j+marge]*q->DZN[k+marge];
    };

    for(int id : lev[l])
    {
        cfd_amr_patch *c = CP(id);
        const double *R = c->a->rhsvec.V.data();
        double *b = pvec(id,NS);
        for(size_t n=0; n<aq[id+1].size(); ++n)
        {
            const int qq = aq[id+1][n];
            b[qq] = R[ar[id+1][n]];
            const double V = cvol(c->pp,qq);
            cs[3*c->comp] += b[qq]*V;
            cs[3*c->comp+1] += V;
        }
    }
    if(ncomp>0)
    MPI_Allreduce(MPI_IN_PLACE,cs.data(),3*ncomp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    vector<char> closed(ncomp,0);
    for(int m=0; m<ncomp; ++m)
    closed[m] = (cs[3*m+2]==0.0 && cs[3*m+1]>0.0) ? 1 : 0;

    for(int id : lev[l])
    {
        cfd_amr_patch *c = CP(id);
        const int m = c->comp;
        if(!closed[m])
        continue;
        const double mean = cs[3*m]/cs[3*m+1];
        double *b = pvec(id,NS);
        for(int qq : aq[id+1])
        b[qq] -= mean;
    }

    pr_apply(-1,NVV);

    double bn, rn;
    pr_space sp{this};
    int it = reefamr_bicgstab(sp,p->N44,p->N46,bn,rn,&tm[5],&tm[6]);
    pr_it_last = it;
    pr_it_total += it;
    ++pr_solves;
    pr_res_last = bn>0.0 ? rn/bn : 0.0;
    sync_bref = MAX(sync_bref,bn);

    // the correction of a closed component: mean 0
    std::fill(cs.begin(),cs.end(),0.0);
    for(int id : lev[l])
    {
        cfd_amr_patch *c = CP(id);
        const double *x = pvec(id,-1);
        for(int qq : aq[id+1])
        {
            const double V = cvol(c->pp,qq);
            cs[3*c->comp] += x[qq]*V;
            cs[3*c->comp+1] += V;
        }
    }
    if(ncomp>0)
    MPI_Allreduce(MPI_IN_PLACE,cs.data(),3*ncomp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(int id : lev[l])
    {
        cfd_amr_patch *c = CP(id);
        const int m = c->comp;
        if(!closed[m])
        continue;
        const double mean = cs[3*m]/cs[3*m+1];
        double *x = pvec(id,-1);
        for(int qq : aq[id+1])
        x[qq] -= mean;
    }

    pr_sync(-1);

    for(int id : lev[l])
    {
        cfd_amr_patch *c = CP(id);
        pscope2 ps(pgc,c->a,a0);
        momentum_rk *m = c->pmom;
        c->ppress->amr_finish(c->pp,c->a,pgc,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s),alpha);
    }

    // the faces at the patch edges: the predictor corrected with the fine correction against the
    // coarse one (0), as the rows of the solve
    for(cfd_amr_cf &F : cf)
    {
        cfd_amr_patch *c = CP(F.pid);
        if(c->lev!=l)
        continue;
        const double xf = pvec(F.pid,-1)[F.qf];
        F.ufin = F.ustar + double(F.dir)*alpha*c->pp->dt*xf/(F.dist*F.rof);
    }
    vel_after();
    for(const cfd_amr_cf &F : cf)
    if(CP(F.pid)->lev==l)
    CP(F.pid)->pmom->amr_velout(CP(F.pid)->a,F.d,s)(F.ff[0],F.ff[1],F.ff[2]) = F.ufin;

    for(int id : lev[l])
    {
        cfd_amr_patch *c = CP(id);
        pscope2 ps(pgc,c->a,a0);
        momentum_rk *m = c->pmom;
        patch_relax_vel(id,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s));
        m->amr_project_after(c->pp,c->a,pgc,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s));
    }

    // back to the composite projection
    prl = -1;
    prg.clear();
    for(int g=-1; g<npat; ++g)
    prg.push_back(g);

    tm[2] += MPI_Wtime()-t0;
}

// G 7 1: the synchronisation projection at the end of the level-0 step.  The patches took their
// own flux at their edges (Dirichlet from the coarser level in the level projections), the
// coarser levels kept theirs; after the restriction the coarse faces at the patch edges are the
// mean of the fine ones and the coarse cells next to the patches have the divergence of the
// mismatch.  One composite correction solve over all levels (as G 7 0, the faces at the patch edges
// with their fine values as the predictor), the velocities of all grids corrected, the
// correction added to the pressure.  The step of level 0, alpha 1.
void cfd_amr::pr_sync_project(lexer *p, ghostcell *pgc)
{
    const int s = mom0->amr_stages()-1;     // the end state of the step: velout of the last stage
    const double alpha = 1.0;
    const double t0 = MPI_Wtime();

    vector<double> dts(npat);
    for(int id=0; id<npat; ++id)
    {
        dts[id] = CP(id)->pp->dt;
        CP(id)->pp->dt = p->dt;
    }

    fill_vel_after(s);
    if(place)
    {
        vector<std::function<field&(int)>> f(1,[this](int g) -> field& { return gfd(g)->phi; });
        sh_down_f(f);
    }
    pr_weights();

    // the predictor at the faces of the patch edges: their own value; the coarse faces their mean
    for(cfd_amr_cf &F : cf)
    {
        F.hok = false;
        F.ustar = CP(F.pid)->pmom->amr_velout(CP(F.pid)->a,F.d,s)(F.ff[0],F.ff[1],F.ff[2]);
    }
    for(cfd_amr_cf &F : cf)
    {
        gmom(F.gface)->amr_velout(gfd(F.gface),F.d,s)(F.cf[0],F.cf[1],F.cf[2]) = 0.0;
        if(F.dir>0 && F.gc!=F.gface)
        {
            int c[3] = {F.ci[0],F.ci[1],F.ci[2]};
            c[F.d] -= 1;
            gmom(F.gc)->amr_velout(gfd(F.gc),F.d,s)(c[0],c[1],c[2]) = 0.0;
        }
    }
    for(cfd_amr_cf &F : cf)
    {
        gmom(F.gface)->amr_velout(gfd(F.gface),F.d,s)(F.cf[0],F.cf[1],F.cf[2]) += F.area*F.ustar;
        if(F.dir>0 && F.gc!=F.gface)
        {
            int c[3] = {F.ci[0],F.ci[1],F.ci[2]};
            c[F.d] -= 1;
            gmom(F.gc)->amr_velout(gfd(F.gc),F.d,s)(c[0],c[1],c[2]) += F.area*F.ustar;
        }
    }
    if(place)
    for(int d=0; d<3; ++d)
    {
        vector<std::function<field&(int)>> f(1,[this,s,d](int g) -> field& { return gmom(g)->amr_velout(gfd(g),d,s); });
        sh_run(pu_cf[d],f,false);
    }

    field &u0 = mom0->amr_velout(a0,0,s), &v0 = mom0->amr_velout(a0,1,s), &w0 = mom0->amr_velout(a0,2,s);

    pflow0->pressure_io(p,a0,pgc);
    press0->amr_prepare(p,a0,pois0,pgc,u0,v0,w0,alpha);

    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        pscope2 ps(pgc,c->a,a0);
        momentum_rk *m = c->pmom;
        c->pflow->pressure_io(c->pp,c->a,pgc);
        c->ppress->amr_prepare(c->pp,c->a,c->ppois,pgc,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s),alpha);
    }

    pr_matrix(s);

    for(int g=-1; g<npat; ++g)
    {
        const double *R = gfd(g)->rhsvec.V.data();
        double *b = pvec(g,NS);
        for(size_t n=0; n<lq[g+1].size(); ++n)
        b[lq[g+1][n]] = R[lr[g+1][n]];
        for(int qq : cq[g+1])
        b[qq] = 0.0;
    }

    // the start: the correction of the last synchronisation (the mismatch at the patch edges
    // changes slowly from step to step), as long as the patches are the same
    if(sync_x_ok && sync_x.size()==npat+1)
    {
        for(int g=-1; g<npat; ++g)
        if(sync_x[g+1].size()==lq[g+1].size())
        {
            double *x = pvec(g,-1);
            for(size_t n=0; n<lq[g+1].size(); ++n)
            x[lq[g+1][n]] = sync_x[g+1][n];
        }
        pr_sync(-1);
    }

    pr_apply(-1,NVV);

    // the tolerance: N 44 relative to the right-hand side of the level projections of the step
    // (the divergence they leave), when that is larger than the one of the synchronisation
    double tol = p->N44;
    {
        const double b0 = sqrt(pr_dot(NS,NS));
        if(b0>0.0 && sync_bref>b0)
        tol *= sync_bref/b0;
    }
    sync_bref = 0.0;

    double bn, rn;
    pr_space sp{this};
    int it = reefamr_bicgstab(sp,tol,p->N46,bn,rn,&tm[5],&tm[6]);
    pr_it_sync = it;

    sync_x.assign(npat+1,vector<double>());
    for(int g=-1; g<npat; ++g)
    {
        const double *x = pvec(g,-1);
        sync_x[g+1].reserve(lq[g+1].size());
        for(int qq : lq[g+1])
        sync_x[g+1].push_back(x[qq]);
    }
    sync_x_ok = true;

    pr_sync(-1);

    press0->amr_finish(p,a0,pgc,u0,v0,w0,alpha);
    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        pscope2 ps(pgc,c->a,a0);
        momentum_rk *m = c->pmom;
        c->ppress->amr_finish(c->pp,c->a,pgc,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s),alpha);
    }

    fill_vel_after(s);

    // the faces at the patch edges (as pr_correct with alpha 1)
    for(const cfd_amr_cf &F : cf)
    {
        const double xf = pvec(F.pid,-1)[F.qf];
        const double xc = pr_xs(F,-1);
        const double rof = 1.0/(F.kf*F.dist*F.dxnf);
        field &ff = CP(F.pid)->pmom->amr_velout(CP(F.pid)->a,F.d,s);
        ff(F.ff[0],F.ff[1],F.ff[2]) = F.ustar - double(F.dir)*alpha*p->dt*(xc-xf)/(F.dist*rof);
    }

    mom0->amr_project_after(p,a0,pgc,u0,v0,w0);
    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        pscope2 ps(pgc,c->a,a0);
        momentum_rk *m = c->pmom;
        patch_relax_vel(id,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s));
        m->amr_project_after(c->pp,c->a,pgc,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s));
    }

    restrict_vel(s);

    for(int id=0; id<npat; ++id)
    CP(id)->pp->dt = dts[id];

    t_sync += MPI_Wtime()-t0;
    tm[2] += MPI_Wtime()-t0;
}
