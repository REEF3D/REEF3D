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

    for(int g=-1; g<(int)P.size(); ++g)
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
            if(covered(l,I))
            cq[g+1].push_back(qq);
            else
            {
                lq[g+1].push_back(qq);
                lr[g+1].push_back(row[qq]);
            }
        }
    }

    // the local grid of level lv whose interior holds the global cell I (-1: level 0)
    auto grid_of = [&](int lv, const int *I, int *loc) -> int
    {
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

    for(int id=0; id<(int)P.size(); ++id)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        const int l = c->lev;

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
                        ++bad;
                        continue;
                    }
                    lexer *qp = glex(F.gp);
                    lexer *qc = glex(F.gc);
                    F.qp = cix(qp,pl[0],pl[1],pl[2]);
                    F.qc = cix(qc,cl[0],cl[1],cl[2]);
                    const vector<int> &rowc = (F.gc<0) ? row0 : CP(F.gc)->row;
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

    bad = pgc0->globalimax(bad);
    if(bad>0)
    {
        if(myrank==0)
        cout<<"CFD AMR: a patch edge lies on a partition edge of level 0 (the coarse cell across it is on another rank); "
              "move the refinement box (G 15) or change the partition"<<endl;
        MPI_Abort(MPI_COMM_WORLD,-3120);
    }

    // vectors
    for(int k=0; k<NVEC+1; ++k)
    kv0.push_back(new field4(p0));
    for(int id=0; id<(int)P.size(); ++id)
    for(int k=0; k<NVEC+1; ++k)
    CP(id)->kv.push_back(new field4(CP(id)->pp));

    // multigrids
    mg0 = new reefmg_core;
    mg0->set_precision(p0->N14==32 ? 32 : 64);
    mg0->set_fine_operator(&pr_fineop_none);
    if(!mg0->setup(pgc0->cart(),p0->knox,p0->knoy,p0->knoz,p0->gknox,p0->gknoy,p0->N13,p0->DXN+marge-1,p0->DYN+marge-1))
    {
        if(p0->mpirank==0)
        cout<<"CFD AMR: level-0 multigrid "<<mg0->err()<<endl;
        MPI_Abort(MPI_COMM_WORLD,-3122);
    }

    for(int id=0; id<(int)P.size(); ++id)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        c->mg = new reefmg_core;
        c->mg->set_precision(pp->N14==32 ? 32 : 64);
        c->mg->set_fine_operator(&pr_fineop_none);
        if(!c->mg->setup(MPI_COMM_SELF,pp->knox,pp->knoy,pp->knoz,pp->knox,pp->knoy,0,pp->DXN+marge-1,pp->DYN+marge-1))
        cout<<"CFD AMR: patch multigrid "<<c->mg->err()<<endl;
    }

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

        auto leaf = [&](const int *cl) -> bool
        {
            const int G[3] = {cl[0]+og[0],cl[1]+og[1],cl[2]+og[2]};
            if(!in_domain(lc,G))
            return false;
            if(lc>0 && !refined(lc,G))
            return false;
            if(covered(lc,G))
            return false;
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
        F.rof = c->pp->W1*H + c->pp->W3*(1.0-H);
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
void cfd_amr::pr_hydro()
{
    const bool zgrav = (p0->W20==0.0 && p0->W21==0.0 && p0->W22!=0.0);
    const double g = fabs(p0->W22);


    // grid gg, column (i,j): p_h of cell kc, from the plane on top of cell ktop-1 (kc <= ktop)
    auto column = [&](int gg, int i, int j, int ktop, int kc, double &val) -> bool
    {
        lexer *q = glex(gg);
        fdm *a = gfd(gg);
        const int m = q->margin;
        const int k0 = MIN(kc,ktop-1);
        if(k0<-m || ktop>q->knoz+m-1 || i<-m || i>q->knox+m-1 || j<-m || j>q->knoy+m-1)
        return false;
        // densities as the grid has them: the local interface width (interface_width.h)
        auto rho = [&](double ph, double e) { const double H = heaviside(ph,e); return q->W1*H + q->W3*(1.0-H); };
        for(int k=k0; k<=MIN(kc+1,ktop); ++k)
        if(q->flag4[cix(q,i,j,k)]<=0 && k<ktop)
        return false;

        double h = g*rho(a->phi(i,j,ktop-1),interface_width(q,a->phi,i,j,ktop-1))*0.5*q->DZN[ktop-1+marge];   // cell ktop-1
        if(kc==ktop)
        {
            if(q->flag4[cix(q,i,j,ktop)]<=0)
            return false;
            val = h - g*rho(0.5*(a->phi(i,j,ktop-1)+a->phi(i,j,ktop)),interface_width_face(q,a->phi,i,j,ktop-1,0,0,1))*q->DZP[ktop-1+marge];
            return true;
        }
        for(int k=ktop-2; k>=kc; --k)
        {
            if(q->flag4[cix(q,i,j,k)]<=0)
            return false;
            h += g*rho(0.5*(a->phi(i,j,k)+a->phi(i,j,k+1)),interface_width_face(q,a->phi,i,j,k,0,0,1))*q->DZP[k+marge];
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
            const double roc = qc->W1*Hc + qc->W3*(1.0-Hc);
            const double rof = F.rof;

            const double *dp[3] = {qc->DXP,qc->DYP,qc->DZP};
            const int Lc = (F.dir>0) ? F.pi[F.d] : F.ci[F.d];
            const double Gc = double(F.dir)*(pC-pP)/(dp[F.d][Lc+marge]*roc);
            const double Gf = double(F.dir)*(pCs-pf)/(F.dist*rof);

            v += bs*p0->dt*(Gc-Gf);

            // a face normal to gravity: the hydrostatic part of the pressure gradient balances the
            // gravity of the coarse stage velocity, both out
            if(F.hok && F.d==2)
            v -= bs*p0->dt*p0->W22;
        }

        F.ustar = v;
        field &ff = CP(F.pid)->pmom->amr_velout(CP(F.pid)->a,F.d,s);
        ff(F.ff[0],F.ff[1],F.ff[2]) = v;
    }

    // the coarse faces: the mean of their fine faces (the right-hand sides of both sides see the
    // same flux)
    for(cfd_amr_cf &F : cf)
    gmom(F.gface)->amr_velout(gfd(F.gface),F.d,s)(F.cf[0],F.cf[1],F.cf[2]) = 0.0;
    for(cfd_amr_cf &F : cf)
    gmom(F.gface)->amr_velout(gfd(F.gface),F.d,s)(F.cf[0],F.cf[1],F.cf[2]) += F.area*F.ustar;
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
        mg0->coarsen();
    }

    // patches: their own matrix (before the couplings across the patch edges below), the couplings
    // to the cells around the patch dropped
    for(int id=0; id<(int)P.size(); ++id)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        sc_level &L = c->mg->fine();
        pr_clear(L);
        const matrix_diag &M = c->a->M;
        const int nx=pp->knox, ny=pp->knoy, nz=pp->knoz;

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
            pr_coef(c->mg,L,q,M.p[r],
                    (i+1<nx) ? M.n[r] : 0.0, (i>0) ? M.s[r] : 0.0,
                    (j+1<ny) ? M.w[r] : 0.0, (j>0) ? M.e[r] : 0.0,
                    (k+1<nz) ? M.t[r] : 0.0, (k>0) ? M.b[r] : 0.0);
            L.act[q] = 1;
        }
        c->mg->coarsen();
    }

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

        matrix_diag &Mc = gfd(F.gc)->M;
        double &ccn = mcoef(Mc,F.nc,F.d,-F.dir);
        Mc.p[F.nc] += ccn;
        ccn = 0.0;
    }

}

// the cells around the patches: the vector of the patch of the same level, or the value of the
// coarser cell (the patch rows couple to the coarse cell across a patch edge); the level-0 halo
void cfd_amr::pr_sync(int k)
{
    if(p0->mpi_size>1)
    pgc0->gcparax(p0,pfield(-1,k),4);

    vector<fspec> fs;
    fs.push_back({CELL,2,[this,k](int g) -> field& { return pfield(g,k); }});
    fill_all(fs);
}

// y = A x on the leaf cells
void cfd_amr::pr_apply(int kx, int ky)
{
    pr_sync(kx);

    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        const matrix_diag &M = gfd(g)->M;
        const double *x = pvec(g,kx);
        double *y = pvec(g,ky);
        const int sI = q->jmax*q->kmax;
        const int sJ = q->kmax;
        const vector<int> &LQ = lq[g+1], &LR = lr[g+1];

        const double *Mp=M.p.data(), *Mn=M.n.data(), *Ms=M.s.data(), *Mw=M.w.data(), *Me=M.e.data(), *Mt=M.t.data(), *Mb=M.b.data();

        for(size_t n=0; n<LQ.size(); ++n)
        {
            const int qq = LQ[n], r = LR[n];
            y[qq] = Mp[r]*x[qq] + Mn[r]*x[qq+sI] + Ms[r]*x[qq-sI] + Mw[r]*x[qq+sJ] + Me[r]*x[qq-sJ] + Mt[r]*x[qq+1] + Mb[r]*x[qq-1];
        }

        for(int qq : cq[g+1])
        y[qq] = 0.0;
    }

    for(const cfd_amr_cf &F : cf)
    {
        const double xs = pr_xs(F,kx);
        const double xf = pvec(F.pid,kx)[F.qf];
        pvec(F.pid,ky)[F.qf] -= F.kf*xs;
        pvec(F.gc,ky)[F.qc] += F.kc*(xs-xf);
    }
}

double cfd_amr::pr_dot(int ka, int kb)
{
    double s=0.0;
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *a = pvec(g,ka), *b = pvec(g,kb);
        for(int qq : lq[g+1])
        s += a[qq]*b[qq];
    }
    MPI_Allreduce(MPI_IN_PLACE,&s,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return s;
}

void cfd_amr::pr_start()
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vb = pvec(g,NS);
        double *vr = pvec(g,NR), *vrh = pvec(g,NRH), *vp = pvec(g,NPV), *vv = pvec(g,NVV);
        for(int qq : lq[g+1])
        {
            double r = vb[qq] - vv[qq];
            vr[qq] = r; vrh[qq] = r; vp[qq] = 0.0; vv[qq] = 0.0;
        }
    }
}

void cfd_amr::pr_p(double beta, double om)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vr = pvec(g,NR), *vv = pvec(g,NVV);
        double *vp = pvec(g,NPV);
        for(int qq : lq[g+1])
        vp[qq] = vr[qq] + beta*(vp[qq] - om*vv[qq]);
    }
}

void cfd_amr::pr_s(double alp)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vr = pvec(g,NR), *vv = pvec(g,NVV);
        double *vs = pvec(g,NS);
        for(int qq : lq[g+1])
        vs[qq] = vr[qq] - alp*vv[qq];
    }
}

void cfd_amr::pr_x(double alp, double om)
{
    for(int g=-1; g<(int)P.size(); ++g)
    {
        const double *vph = pvec(g,NPH), *vsh = pvec(g,NSH), *vs = pvec(g,NS), *vt = pvec(g,NT);
        double *x = pvec(g,-1), *vr = pvec(g,NR);
        for(int qq : lq[g+1])
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

// z = M^-1 r: one FAC V-cycle over the levels.  Down: on every level, finest first, a patch-local
// V-cycle on the composite residual (the cells under a finer patch take the residual of the finer
// level); level 0: a V-cycle on the residual left; up: the correction of the coarser level prolonged
// into the patches, a patch-local V-cycle on the residual left.
void cfd_amr::pr_prec(int kr, int kz)
{
    const int KT = REEFAMR_NTMP;        // residual
    const int KC = REEFAMR_NVEC;        // correction of a level, prolonged to the next finer one

    for(int g=-1; g<(int)P.size(); ++g)
    {
        lexer *q = glex(g);
        const size_t n = (size_t)q->imax*q->jmax*q->kmax;
        std::fill(pvec(g,kz),pvec(g,kz)+n,0.0);
        std::fill(pvec(g,KC),pvec(g,KC)+n,0.0);
    }

    vector<fspec> fc;
    fc.push_back({CELL,2,[this](int g) -> field& { return pfield(g,REEFAMR_NVEC); }});
    const fspec cl = {CELL,0,[this](int g) -> field& { return pfield(g,REEFAMR_NVEC); }};

    // down
    for(int l=maxlev; l>=1; --l)
    {
        if(l==maxlev)
        {
            // z = 0: the residual is r
            for(int id : lev[l])
            pr_patch_solve(id,kr,kz,-1);
        }
        else
        {
            pr_residual(kr,kz,KT,l);
            for(int id : lev[l])
            pr_patch_solve(id,KT,kz,-1);
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

        if(p0->mpi_size>1)
        pgc0->gcparax(p0,pfield(-1,KC),4);
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
                        const double v = prolong(cl,B.g,sp,o);
                        t(i,j,k) = v;
                        z(i,j,k) += v;
                    }
                }
            }
        }

        pr_residual(kr,kz,KT,l);
        for(int id : lev[l])
        pr_patch_solve(id,KT,kz,-1);

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

    std::fill(L.u.begin(),L.u.end(),0.0);
    std::fill(L.f.begin(),L.f.end(),0.0);
    for(int i=0; i<pp->knox; ++i)
    for(int j=0; j<pp->knoy; ++j)
    for(int k=0; k<pp->knoz; ++k)
    {
        const int qq = cix(pp,i,j,k);
        if(c->row[qq]>=0)
        L.f[L.idx(i,j,k)] = t[qq];
    }

    c->mg->vcycle(0,2,2);

    for(int i=0; i<pp->knox; ++i)
    for(int j=0; j<pp->knoy; ++j)
    for(int k=0; k<pp->knoz; ++k)
    {
        const int qq = cix(pp,i,j,k);
        if(c->row[qq]<0)
        continue;
        const double u = L.u[L.idx(i,j,k)];
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

    for(int id=0; id<(int)P.size(); ++id)
    {
        cfd_amr_patch *c = CP(id);
        pscope2 ps(pgc,c->a,a0);
        momentum_rk *m = c->pmom;
        c->pflow->pressure_io(c->pp,c->a,pgc);
        c->ppress->amr_prepare(c->pp,c->a,c->ppois,pgc,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s),alpha);
    }

    pr_matrix(s);

    // right-hand side into NS, A x0 into NVV (x0 = 0: pcorr cleared by pjm_corr::rhs)
    for(int g=-1; g<(int)P.size(); ++g)
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

    for(int id=0; id<(int)P.size(); ++id)
    {
        cfd_amr_patch *c = CP(id);
        pscope2 ps(pgc,c->a,a0);
        momentum_rk *m = c->pmom;
        c->ppress->amr_finish(c->pp,c->a,pgc,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s),alpha);
    }

    fill_vel_after(s);
    pr_correct(s);

    mom0->amr_project_after(p,a0,pgc,u0,v0,w0);
    for(int id=0; id<(int)P.size(); ++id)
    {
        cfd_amr_patch *c = CP(id);
        pscope2 ps(pgc,c->a,a0);
        momentum_rk *m = c->pmom;
        m->amr_project_after(c->pp,c->a,pgc,m->amr_velout(c->a,0,s),m->amr_velout(c->a,1,s),m->amr_velout(c->a,2,s));
    }

    restrict_vel(s);

    p->poissontime = MPI_Wtime()-t0;
}
