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
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<set>
#include<iostream>

//  G 40: the placement of the CFD patches (reefamr3d: units, shadows).  A level-1 patch reads its
//  coarser level from its shadow, the level-0 cells under it and around it on its rank.  The
//  shadows are filled from level 0 when the patches read them (down: the sources of the fills of
//  level 1, the state of new patches, the coarse data of the faces at the patch edges, the
//  correction of the level-0 V-cycle), and what the patches write onto level 0 goes back to the
//  ranks that hold the cells (up):
//    - the restriction onto the cells and faces under the level-1 patches (set; a face shared by
//      two level-1 patches from the patch below it)
//    - the coarse faces at the patch edges, the mean of the fine faces (pr_predict, the
//      synchronisation of G 7 1; set)
//    - the coarse cells at the patch edges: the fluxes of the fine faces in the composite operator
//      (pr_apply) and the refluxing of heat and concentration (added)
//    - the coarse cell and side of every face at the patch edges: the level-0 row folds its coupling
//      to the covered cell into its diagonal (once per assembly), the level-0 multigrid takes the
//      couplings of the fine faces (pr_matrix)
//  A face is set on the rank that holds the cell below it and on the one that holds the cell above
//  (its copy in the halo).  Every exchange is collective; the shadows hold no state.

namespace
{
inline int fdiv(int a, int r) { return (a>=0) ? a/r : -((-a+r-1)/r); }

inline int cix(const lexer *q, int i, int j, int k)
{
    return (i-q->imin)*q->jmax*q->kmax + (j-q->jmin)*q->kmax + k-q->kmin;
}
}

void cfd_amr::sh_run(r3cplan &X, const vector<std::function<field&(int)>> &f, bool add)
{
    if(!place || f.empty())
    return;
    // (the arrays resolved once per grid, not per item)
    const int ng = (int)P.size()+1;
    vector<double*> tab(f.size()*(size_t)ng,nullptr);
    auto F = [&](int m, int g) -> double*
    {
        double *&t = tab[(size_t)m*ng+g+1];
        if(t==nullptr)
        t = f[m](g).V;
        return t;
    };
    cplan_run(X,(int)f.size(),F,F,add);
}

void cfd_amr::sh_fill_down(const vector<fspec> &fs)
{
    vector<std::function<field&(int)>> f;
    for(const fspec &F : fs)
    f.push_back(F.f);
    sh_run(sh_down,f,false);
}

void cfd_amr::sh_zero(r3cplan &X, const std::function<field&(int)> &f)
{
    for(size_t s=0; s<X.sp.size(); ++s)
    for(size_t q=0; q<X.sq[s].size(); ++q)
    f(X.sg[s][q]).V[X.sq[s][q]] = 0.0;
}

// the up plans of the shadows (after pr_setup; collective)
void cfd_amr::sh_plans()
{
    if(!place)
    return;

    auto sidx = [&](int g, const int *I) -> int
    {
        const r3patch *c = P[g];
        return cix(c->pp,I[0]-c->lo[0],I[1]-c->lo[1],I[2]-c->lo[2]);
    };
    auto push = [&](vector<vector<int>> &it, int r, int g, int q, const int *I, int tag)
    {
        const int e[6] = {g,q,I[0],I[1],I[2],tag};
        it[r].insert(it[r].end(),e,e+6);
    };
    // a face (the + face of level-0 cell I along d): on the rank of I and on the rank of the cell
    // above (its - face, in the halo there)
    auto push_face = [&](vector<vector<int>> &it, int g, int q, const int *I, int d)
    {
        int J[3] = {I[0],I[1],I[2]};
        J[d] += 1;
        int r1 = -1;
        if(in_domain(0,I))
        {
            r1 = owner0(I);
            push(it,r1,g,q,I,0);
        }
        if(in_domain(0,J) && owner0(J)!=r1)
        push(it,owner0(J),g,q,I,0);
    };

    // the restriction: the cells under the level-1 patches; the faces under them, the - face at the
    // low edge unless a level-1 patch covers the cell below it (that patch gives it, its + face)
    vector<vector<int>> ic(nranks), ifc[3];
    for(int d=0; d<3; ++d)
    ifc[d].assign(nranks,vector<int>());

    for(int id : lev[1])
    {
        const r3patch *c = P[id];
        const int g = c->sh;
        int a[3], b[3];
        for(int d=0; d<3; ++d)
        {
            a[d] = fdiv(c->lo[d],rr[d]);
            b[d] = fdiv(c->hi[d],rr[d]);
        }
        int I[3];
        for(I[0]=a[0]; I[0]<=b[0]; ++I[0])
        for(I[1]=a[1]; I[1]<=b[1]; ++I[1])
        for(I[2]=a[2]; I[2]<=b[2]; ++I[2])
        push(ic,owner0(I),g,sidx(g,I),I,0);

        for(int d=0; d<3; ++d)
        {
            if(rr[d]!=2)
            continue;
            int lo[3] = {a[0],a[1],a[2]};
            lo[d] -= 1;
            for(I[0]=lo[0]; I[0]<=b[0]; ++I[0])
            for(I[1]=lo[1]; I[1]<=b[1]; ++I[1])
            for(I[2]=lo[2]; I[2]<=b[2]; ++I[2])
            {
                if(I[d]==a[d]-1 && in_domain(0,I) && covered(0,I))
                continue;
                push_face(ifc[d],g,sidx(g,I),I,d);
            }
        }
    }
    cplan_build(pu_cell,false,ic);
    for(int d=0; d<3; ++d)
    cplan_build(pu_face[d],false,ifc[d]);

    // the faces at the edges of the level-1 patches: the coarse faces, the coarse cells, the keys
    vector<vector<int>> icf[3], icc(nranks), ik(nranks);
    for(int d=0; d<3; ++d)
    icf[d].assign(nranks,vector<int>());
    std::set<std::tuple<int,int,int,int,int>> sf, scc;
    kidx.clear();
    kval.assign(P.size()-npat,vector<double>());

    for(const cfd_amr_cf &F : cf)
    {
        if(P[F.pid]->lev!=1)
        continue;
        const int g = F.gc;
        const int *og = P[g]->lo;

        const int Cf[3] = {F.cf[0]+og[0],F.cf[1]+og[1],F.cf[2]+og[2]};
        if(sf.insert(std::make_tuple(g,F.d,Cf[0],Cf[1],Cf[2])).second)
        push_face(icf[F.d],g,cix(P[g]->pp,F.cf[0],F.cf[1],F.cf[2]),Cf,F.d);

        const int C[3] = {F.ci[0]+og[0],F.ci[1]+og[1],F.ci[2]+og[2]};
        if(scc.insert(std::make_tuple(g,C[0],C[1],C[2],0)).second)
        push(icc,owner0(C),g,F.qc,C,0);

        // the key: to the rank of the coarse cell (its row, its multigrid row), and to the rank of
        // the covered cell across the face when that is another one (its multigrid row); tag
        // 6 d + 3 (dir > 0) + h, h 0: both rows here, 1: the coarse cell's, 2: the covered cell's
        const auto key = std::make_tuple(g,C[0],C[1],C[2],F.d,F.dir);
        if(kidx.find(key)==kidx.end())
        {
            const int n = (int)kval[g-npat].size();
            kidx[key] = n;
            kval[g-npat].push_back(0.0);
            int Pc[3] = {C[0],C[1],C[2]};
            Pc[F.d] -= F.dir;
            const int rc = owner0(C), rp = owner0(Pc);
            const int t = 6*F.d+3*(F.dir>0 ? 1 : 0);
            push(ik,rc,g,n,C,t+(rp==rc ? 0 : 1));
            if(rp!=rc)
            push(ik,rp,g,n,Pc,t+2);
        }
    }
    for(int d=0; d<3; ++d)
    cplan_build(pu_cf[d],false,icf[d]);
    cplan_build(pu_cc,false,icc);
    cplan_build(pu_key,false,ik);

    // down for the solver vectors (sh_down1): only the cells of the shadows the operator and the
    // fills next to the level-1 patches read, the coarse cells at their edges with their transverse
    // neighbours (the stencils of pr_weights) and the parents of the cells next to the patches (the
    // fills up to depth 1, injection); sh_down sends every cell of the shadow arrays
    {
        vector<vector<char>> mark(P.size());
        auto need = [&](int g, int q)
        {
            if(g<npat || q<0)
            return;
            vector<char> &m = mark[g];
            if(m.empty())
            m.assign((size_t)P[g]->pp->imax*P[g]->pp->jmax*P[g]->pp->kmax,0);
            m[q] = 1;
        };
        for(int id : lev[1])
        {
            const r3patch *c = P[id];
            for(const r3fill &f : c->fill)
            if(f.kind==1 && f.depth<=1)
            need(f.g,cix(glex(f.g),f.s[0],f.s[1],f.s[2]));
        }
        for(const cfd_amr_cf &F : cf)
        {
            if(glevel(F.pid)!=1)
            continue;
            need(F.gc,F.qc);
            for(int t=0; t<3; ++t)
            {
                need(F.gc,F.tq[t][0]);
                need(F.gc,F.tq[t][1]);
            }
        }

        vector<vector<int>> it(nranks);
        for(int g=npat; g<(int)P.size(); ++g)
        {
            if(mark[g].empty())
            continue;
            const r3patch *c = P[g];
            const lexer *q = c->pp;
            const int sI = q->jmax*q->kmax;
            for(int n=0; n<(int)mark[g].size(); ++n)
            if(mark[g][n])
            {
                int I[3] = {c->lo[0] + n/sI + q->imin, c->lo[1] + (n/q->kmax)%q->jmax + q->jmin, c->lo[2] + n%q->kmax + q->kmin};
                wrap0(I);
                push(it,owner0(I),g,n,I,0);
            }
        }
        cplan_build(sh_down1,true,it);
    }

    // the porosity of the shadows (it does not change)
    vector<std::function<field&(int)>> por;
    por.push_back([this](int g) -> field& { return gfd(g)->porosity; });
    sh_down_f(por);
}
