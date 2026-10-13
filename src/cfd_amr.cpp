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
#include"poisson_pcorr.h"
#include"convection.h"
#include"fou.h"
#include"weno_flux_nug.h"
#include"weno_hj_nug.h"
#include"weno_hj_df_nug.h"
#include"diff_void.h"
#include"ediff2.h"
#include"ediff2_2D.h"
#include"idiff2_FS.h"
#include"idiff2_FS_2D.h"
#include"bicgstab_ijk.h"
#include"bicgstab_ijk_2D.h"
#include"solver_void.h"
#include"kepsilon_void.h"
#include"kepsilon_IM1.h"
#include"komega_IM1.h"
#include"rans_io.h"
#include"6DOF_cfd.h"
#include"cfd_amr_6dof.h"
#include"convection_void.h"
#include"ifou.h"
#include"iweno_hj_df_nug.h"
#include"idiff2.h"
#include"heat_void.h"
#include"concentration_void.h"
#include"reini_RK3.h"
#include"patchBC_void.h"
#include"ioflow_void.h"
#include"vrans_v.h"
#include"initialize.h"
#include"reefmg_core.h"
#include<cmath>
#include<iostream>
#include<iomanip>
#include<algorithm>
#include<limits>

namespace
{
// floor(a/r) for negative a as well
inline int fdiv(int a, int r) { return (a>=0) ? a/r : -((-a+r-1)/r); }

// index of a cell in the arrays of lexer q (field layout)
inline int cix(const lexer *q, int i, int j, int k)
{
    return (i-q->imin)*q->jmax*q->kmax + (j-q->jmin)*q->kmax + k-q->kmin;
}

// The implicit diffusion of a patch: the cells around the patch are Dirichlet values (the stage
// velocity of the coarser level).  bicgstab_ijk keeps its own vector and exchanges the cells
// around a rank grid with the partition neighbours; a patch has no neighbours, so the couplings
// to those cells are taken into the right-hand side before the solve.  The faces at the high
// patch edges (the last face of a row of the patch, its own unknown) are Dirichlet values too:
// the coarse face velocity (edge, NaN elsewhere); their momentum predictor would carry the
// pressure jump across the edge into the patch.
class amr_dsolver : public solver
{
public:
    amr_dsolver(solver *s, const vector<signed char> *mk, const vector<double> *ed) : in(s), mkind(mk), edge(ed) {}

    void start(lexer *p, fdm *a, ghostcell *pgc, field &f, vec &rhsvec, int var) override
    {
        int *flag = p->flag4;
        int ul=0, vl=0, wl=0;
        if(var==1) { flag = p->flag1; ul = p->ulast; }
        if(var==2) { flag = p->flag2; vl = p->vlast; }
        if(var==3) { flag = p->flag3; wl = p->wlast; }

        const int sI = p->jmax*p->kmax, sJ = p->kmax;
        const vector<signed char> &mk = *mkind;
        matrix_diag &M = a->M;

        auto fold = [&](double &c, int q, int n)
        {
            if(c!=0.0 && mk[q]>=0 && mk[q]!=3)
            {
                rhsvec.V[n] -= c*f.V[q];
                c = 0.0;
            }
        };

        const vector<double> *ed = (var>=1 && var<=3 && !edge[var-1].empty()) ? &edge[var-1] : nullptr;

        int n=0;
        for(int i=0; i<p->knox-ul; ++i)
        for(int j=0; j<p->knoy-vl; ++j)
        for(int k=0; k<p->knoz-wl; ++k)
        {
            const int q = cix(p,i,j,k);
            if(flag[q]<=0)
            continue;
            if(ed!=nullptr && !std::isnan((*ed)[q]))
            {
                const double v = (*ed)[q];
                M.p[n] = 1.0;
                M.n[n] = M.s[n] = M.w[n] = M.e[n] = M.t[n] = M.b[n] = 0.0;
                rhsvec.V[n] = v;
                f.V[q] = v;
                ++n;
                continue;
            }
            fold(M.n[n],q+sI,n);
            fold(M.s[n],q-sI,n);
            if(p->j_dir==1)
            {
            fold(M.w[n],q+sJ,n);
            fold(M.e[n],q-sJ,n);
            }
            fold(M.t[n],q+1,n);
            fold(M.b[n],q-1,n);
            ++n;
        }

        in->start(p,a,pgc,f,rhsvec,var);
    }
    void startf(lexer *p, ghostcell *pgc, field &f, vec &r, matrix_diag &M, int var) override { in->startf(p,pgc,f,r,M,var); }
    void startF(lexer *p, ghostcell *pgc, double *f, vec &r, matrix_diag &M, int var) override { in->startF(p,pgc,f,r,M,var); }
    void startV(lexer *p, ghostcell *pgc, double *f, vec &r, matrix_diag &M, int var) override { in->startV(p,pgc,f,r,M,var); }
    void startM(lexer *p, ghostcell *pgc, double *x, double *r, double *M, int var) override { in->startM(p,pgc,x,r,M,var); }

private:
    solver *in;
    const vector<signed char> *mkind;
    const vector<double> *edge;
};

// patch work: no exchange, rank-local reductions, the fdm of the patch in the ghostcell object
struct pscope
{
    ghostcell *g;
    fdm *back;
    bool oldc, oldl;
    pscope(ghostcell *gg, fdm *a, fdm *a0) : g(gg), back(a0)
    {
        oldc = g->set_comms(false);
        oldl = g->set_local(true);
        g->fdm_update(a);
    }
    ~pscope()
    {
        g->fdm_update(back);
        g->set_comms(oldc);
        g->set_local(oldl);
    }
};
}

cfd_amr_patch::~cfd_amr_patch()
{
}

// --------------------------------------------------------------------- scope
bool cfd_amr::scope(lexer *p, ghostcell *pgc)
{
    if(p->G1<=0 || p->A10!=6)
    return false;

    const char *why = nullptr;

    if(!(p->N40==2 || p->N40==3))
    why = "N 40 2 or 3 (SSP Runge-Kutta, level set in the stages)";
    else if(p->F80!=0 || p->F300!=0)
    why = "the level set (F 30) or single phase flow, no VOF (F 80) or multiphase (F 300)";
    else if(!(p->T10==0 || p->T10==1 || p->T10==2 || p->T10==21 || p->T10==22))
    why = "laminar flow or a RANS model (T 10 0, 1, 2, 21, 22)";
    else if(p->H10!=0 || p->C10!=0 || p->W30!=0 || p->W90!=0)
    why = "no heat, concentration, compressibility or rheology (H 10, C 10, W 30, W 90)";
    else if(p->S10!=0 || p->Q10!=0)
    why = "no sediment (S 10, Q 10)";
    else if(p->X10>1 || p->Z10!=0 || p->Z20!=0 || p->Z30!=0)
    why = "floating bodies X 10 0 or 1, no structures (Z 10, Z 20, Z 30)";
    else if(p->B200!=0)
    why = "no porous media (B 200)";
    else if((p->solidread!=0 || p->toporead!=0) && p->G505!=0)
    why = "solids and topography with direct forcing (G 505 0)";
    else if(p->D30<1 || p->D30>3)
    why = "the pressure projection D 30 1-3";
    else if(p->B30!=0)
    why = "no pressure reference point (B 30 0)";
    else if(!(p->D10==1 || p->D10==4 || p->D10==5))
    why = "convection D 10 1, 4 or 5";
    else if(p->F30>0 && !(p->F35==1 || p->F35==4 || p->F35==5))
    why = "level-set convection F 35 1, 4 or 5";
    else if(!(p->D20==0 || p->D20==1 || p->D20==2))
    why = "diffusion D 20 0, 1 or 2";
    else if(p->F30>0 && p->F40!=3)
    why = "the reinitialisation F 40 3";
    else if(p->mz>1 || p->periodic1>0 || p->periodic2>0 || p->periodic3>0)
    why = "no partition in z (M 10 with the vertical not split) and no periodic boundaries";
    else if(p->B440>0 || p->B441>0 || p->B442>0)
    why = "no patch boundary conditions (B 440-442); inflow, outflow and waves on level 0 (the patches keep away from them)";

    int bad = (why!=nullptr) ? 1 : 0;
    bad = pgc->globalimax(bad);

    if(bad)
    {
        if(p->mpirank==0)
        cout<<"CFD AMR (G 1 "<<p->G1<<"): needs "<<(why ? why : "a supported set-up")<<"; running without mesh refinement"<<endl;
        return false;
    }
    return true;
}

// --------------------------------------------------------------------- set-up
cfd_amr::cfd_amr(lexer *p, fdm *a, ghostcell *pgc, momentum *pmom, pressure *ppress, poisson *ppois, ioflow *pflow, vrans *pvrans,
                 fsi *pfsi, initialize *pini, turbulence *pturb, sixdof *p6dof) : reefamr3d(p,pgc), p0(p), a0(a), pois0(ppois), pflow0(pflow), pvrans0(pvrans),
                 pfsi0(pfsi), pini0(pini)
{
    // RANS and solids on the patches
    pturb0 = pturb;
    turb = (p->T10==1 || p->T10==2 || p->T10==21 || p->T10==22) && dynamic_cast<rans_io*>(pturb)!=nullptr ? 1 : 0;
    nts = turb ? 8 : 6;
    solids = (p->solidread!=0 || p->toporead!=0) ? 1 : 0;

    // floating bodies: moved on level 0, the hull and its forcing on every grid
    b6 = (p->X10==1) ? dynamic_cast<sixdof_cfd*>(p6dof) : nullptr;
    bodies = (b6!=nullptr && b6->objects()>0) ? 1 : 0;

    mom0 = dynamic_cast<momentum_rk*>(pmom);
    press0 = dynamic_cast<pjm_corr*>(ppress);

    if(mom0==nullptr || press0==nullptr)
    {
        if(p->mpirank==0)
        cout<<"CFD AMR: needs momentum_rk and pjm_corr"<<endl;
        MPI_Abort(MPI_COMM_WORLD,-3110);
    }

    gcval_phi = 51;
    if(p->F50==2) gcval_phi = 52;
    if(p->F50==3) gcval_phi = 53;
    if(p->F50==4) gcval_phi = 54;

    // subcycling (G 7 1): the stage times
    sub = (p->G7==1) ? 1 : 0;
    th_end.assign(maxlev+1,0.0);
    if(mom0->amr_stages()==3)
    {
        tin[0]=0.0; tin[1]=1.0; tin[2]=0.5;
        tout[0]=1.0; tout[1]=0.5; tout[2]=1.0;
    }
    else
    {
        tin[0]=0.0; tin[1]=1.0; tin[2]=1.0;
        tout[0]=1.0; tout[1]=1.0; tout[2]=1.0;
    }

    // the hierarchy: static boxes G 15, adaptive with the flags G 24-26 of the initial level 0
    adapt = (p->G24>0.0 || p->G25>0 || p->G26>0.0) ? 1 : 0;
    if(adapt)
    {
        setup_base(p,pgc);
        vector<vector<unsigned char>> seed;
        flags(seed,true);
        setup(p,pgc,&seed);
    }
    else
    setup(p,pgc);

    hierarchy(true);
}

// after the setup or a regrid
void cfd_amr::hierarchy(bool first)
{
    lexer *p = p0;
    ghostcell *pgc = pgc0;

    // the patches have no inflow or outflow (the wave generation stays on level 0)
    int io = 0;
    for(auto q : P)
    for(int n=0; n<q->pp->gcb4_count; ++n)
    {
        const int g = q->pp->gcb4[n][4];
        if(g==1 || g==2 || g==6 || g==7 || g==8)
        ++io;
    }
    io = pgc->globalimax(io);
    if(io>0)
    {
        if(p->mpirank==0)
        cout<<"CFD AMR: a refined region (G 15, G 24-26) reaches an inflow or outflow boundary; keep the boxes away from them"<<endl;
        MPI_Abort(MPI_COMM_WORLD,-3111);
    }

    int np = (int)GP.size();
    if(p->mpirank==0)
    {
        if(first)
        cout<<"CFD AMR: "<<maxlev<<" refined levels, "<<np<<" patches";
        else
        cout<<"CFD AMR: regrid at step "<<p->count<<": "<<np<<" patches";
        for(int l=1; l<=maxlev; ++l)
        {
            long cells = 0;
            for(int g : glev[l])
            cells += (long)(GP[g].hi[0]-GP[g].lo[0]+1)*(GP[g].hi[1]-GP[g].lo[1]+1)*(GP[g].hi[2]-GP[g].lo[2]+1);
            cout<<"; level "<<l<<": "<<glev[l].size()<<" patches, "<<cells<<" cells";
        }
        cout<<endl;
    }

    // kinds of the cells around the patches
    for(int id=0; id<(int)P.size(); ++id)
    {
        cfd_amr_patch *c = CP(id);
        c->mkind.assign(c->pp->imax*c->pp->jmax*c->pp->kmax,-1);
        for(const r3fill &f : c->fill)
        c->mkind[cix(c->pp,f.d[0],f.d[1],f.d[2])] = (signed char)f.kind;
    }

    pr_setup();

    // subcycling: the state of the grids that are the coarser level of a patch
    if(sub)
    {
        if(maxlev>=1 && ts0==nullptr)
        ts0 = new cfd_amr_tlev(p,nts);
        for(int id=0; id<(int)P.size(); ++id)
        if(P[id]->lev<maxlev && CP(id)->ts==nullptr)
        CP(id)->ts = new cfd_amr_tlev(CP(id)->pp,nts);
    }
}

cfd_amr::~cfd_amr()
{
}

namespace
{
template<class T, class B> bool del_as(B *&b)
{
    T *t = dynamic_cast<T*>(b);
    if(t==nullptr)
    return false;
    delete t;
    b = nullptr;
    return true;
}

// the cell arrays of a patch lexer (build_lexer, flagini); the lexer has no destructor for them
void free_lexer_cells(lexer *q)
{
    int **ia[] = {&q->flag1,&q->flag2,&q->flag3,&q->flag4,&q->flag5,&q->flagsf1,&q->flagsf2,&q->flagsf3,&q->flagsf4,
                  &q->IO,&q->DF,&q->DF1,&q->DF2,&q->DF3,&q->DFF};
    for(int **a : ia)
    {
        delete [] *a;
        *a = nullptr;
    }
    delete [] q->flag_solid;
    delete [] q->flag_topo;
    q->flag_solid = q->flag_topo = nullptr;
    q->geometry_free();     // the horizontal geometry layer of build_lexer
}
}

// a patch of an old hierarchy (regrid): its objects, fields and the cell arrays of its lexer (the
// boundary lists and slices of the lexer are left)
void cfd_amr::free_patch(cfd_amr_patch *c)
{
    delete c->ts;
    for(field4 *f : c->kv)
    delete f;
    c->kv.clear();
    delete c->mg;

    delete c->pmom;
    delete c->p6dof;
    delete c->ppress;
    del_as<poisson_pcorr>(c->ppois);
    del_as<reini_RK3>(c->preini);
    del_as<amr_dsolver>(c->psolv);
    if(!del_as<bicgstab_ijk>(c->psolv_in))
    del_as<bicgstab_ijk_2D>(c->psolv_in);
    del_as<solver_void>(c->ppoissonsolv);
    del_as<heat_void>(c->pheat);
    del_as<concentration_void>(c->pconc);
    if(!del_as<kepsilon_void>(c->pturb) && !del_as<kepsilon_IM1>(c->pturb))
    del_as<komega_IM1>(c->pturb);
    if(!del_as<convection_void>(c->pturbdisc) && !del_as<ifou>(c->pturbdisc))
    del_as<iweno_hj_df_nug>(c->pturbdisc);
    if(!del_as<diff_void>(c->pturbdiff))
    del_as<idiff2>(c->pturbdiff);
    delete c->pBC;
    del_as<ioflow_v>(c->pflow);
    del_as<vrans_v>(c->pvrans);
    if(!del_as<diff_void>(c->pdiff) && !del_as<ediff2>(c->pdiff) && !del_as<ediff2_2D>(c->pdiff) && !del_as<idiff2_FS>(c->pdiff))
    del_as<idiff2_FS_2D>(c->pdiff);
    if(!del_as<fou>(c->pconvec) && !del_as<weno_flux_nug>(c->pconvec))
    del_as<weno_hj_nug>(c->pconvec);
    if(!del_as<fou>(c->pfsfdisc) && !del_as<weno_flux_nug>(c->pfsfdisc))
    del_as<weno_hj_df_nug>(c->pfsfdisc);

    delete c->a;
    free_lexer_cells(c->pp);
    delete c->pp;
    delete c;
}

r3patch* cfd_amr::patch_new()
{
    return new cfd_amr_patch;
}

// the objects of a CFD step on the patch, as driver::logic_cfd builds them for level 0
void cfd_amr::patch_objects(r3patch *rp)
{
    cfd_amr_patch *c = static_cast<cfd_amr_patch*>(rp);
    lexer *pp = c->pp;

    // the patch takes no volume correction of its own (level 0 does it over the whole domain)
    pp->F46 = 0;

    c->a = new fdm(pp);
    pscope ps(pgc0,c->a,a0);
    ghostcell *pgc = pgc0;
    fdm *a = c->a;

    pgc->sizeM_update(pp,a);

    // no solids: direct forcing flags
    for(int n=0; n<pp->imax*pp->jmax*pp->kmax; ++n)
    {
        pp->DF[n] = 1;
        pp->DF1[n] = pp->DF2[n] = pp->DF3[n] = 1;
    }
    pp->gcdf1_count = pp->gcdf2_count = pp->gcdf3_count = pp->gcdf4_count = 0;

    // cells around the patch: fill kind (-1 interior)
    c->mkind.assign(pp->imax*pp->jmax*pp->kmax,-1);

    // convection
    if(pp->D10==1)
    c->pconvec = new fou(pp);
    if(pp->D10==4)
    {
        weno_flux_nug *w = new weno_flux_nug(pp);
        w->set_weno_weights(pp->D12,pp->D13);
        c->pconvec = w;
    }
    if(pp->D10==5)
    {
        weno_hj_nug *w = new weno_hj_nug(pp);
        w->set_weno_weights(pp->D12,pp->D13);
        c->pconvec = w;
    }

    // level-set convection
    if(pp->F35==1)
    c->pfsfdisc = new fou(pp);
    if(pp->F35==4)
    {
        weno_flux_nug *w = new weno_flux_nug(pp);
        w->set_weno_weights(pp->F37,pp->F38);
        c->pfsfdisc = w;
    }
    if(pp->F35==5)
    {
        weno_hj_df_nug *w = new weno_hj_df_nug(pp);
        w->set_weno_weights(pp->F37,pp->F38);
        c->pfsfdisc = w;
    }
    if(c->pfsfdisc==nullptr)
    c->pfsfdisc = new fou(pp);

    // diffusion
    if(pp->D20==0)
    c->pdiff = new diff_void;
    if(pp->D20==1 && pp->j_dir==1)
    c->pdiff = new ediff2(pp);
    if(pp->D20==1 && pp->j_dir==0)
    c->pdiff = new ediff2_2D(pp);
    if(pp->D20==2 && pp->j_dir==1)
    c->pdiff = new idiff2_FS(pp);
    if(pp->D20==2 && pp->j_dir==0)
    c->pdiff = new idiff2_FS_2D(pp);

    if(pp->j_dir==0)
    c->psolv_in = new bicgstab_ijk_2D(pp,a,pgc);
    else
    c->psolv_in = new bicgstab_ijk(pp,a,pgc);
    c->psolv = new amr_dsolver(c->psolv_in,&c->mkind,c->edge);
    c->ppoissonsolv = new solver_void(pp,a,pgc);

    c->pheat = new heat_void(pp,a,pgc);
    c->pconc = new concentration_void(pp,a,pgc);

    // turbulence: the model of level 0, its convection (T 12) and diffusion
    if(turb && (pp->T10==1 || pp->T10==21))
    c->pturb = new kepsilon_IM1(pp,a,pgc);
    else if(turb && (pp->T10==2 || pp->T10==22))
    c->pturb = new komega_IM1(pp,a,pgc);
    else
    c->pturb = new kepsilon_void(pp,a,pgc);

    if(pp->T12==1)
    c->pturbdisc = new ifou(pp);
    else if(pp->T12==5 || pp->T12==55)
    c->pturbdisc = new iweno_hj_df_nug(pp);
    else
    c->pturbdisc = new convection_void(pp);

    if(turb && pp->D20>0)
    c->pturbdiff = new idiff2(pp);
    else
    c->pturbdiff = new diff_void;

    c->ppress = new pjm_corr(pp,a,pgc,c->pheat,c->pconc);
    c->ppois = new poisson_pcorr(pp,c->pheat,c->pconc);

    c->preini = new reini_RK3(pp,1);
    c->pBC = new patchBC_void(pp);
    c->pflow = new ioflow_v(pp,pgc,c->pBC);
    c->pvrans = new vrans_v(pp,pgc);

    c->pmom = new momentum_rk(pp,a,pgc,c->pconvec,c->pfsfdisc,c->pdiff,c->ppress,c->ppois,c->pturb,c->psolv,c->ppoissonsolv,
                              c->pflow,c->pheat,c->pconc,c->preini,pfsi0);

    if(bodies)
    c->p6dof = new cfd_amr_6dof(pp,b6);
}

// --------------------------------------------------------------------- fills
double cfd_amr::cell_slope(field &f, int g, const int *s, int d, int mode)
{
    if(mode==2 || rr[d]==1)
    return 0.0;

    int a[3] = {s[0],s[1],s[2]}, b[3] = {s[0],s[1],s[2]};
    a[d] -= 1;
    b[d] += 1;
    const double c0 = f(s[0],s[1],s[2]);
    const double dl = c0 - f(a[0],a[1],a[2]);
    const double dr = f(b[0],b[1],b[2]) - c0;

    if(mode==0)
    return 0.5*(dl+dr);

    // MC limiter
    if(dl*dr<=0.0)
    return 0.0;
    const double sg = dl>0.0 ? 1.0 : -1.0;
    return sg*MIN(MIN(2.0*fabs(dl),2.0*fabs(dr)),0.5*fabs(dl+dr));
}

// value of a cell (or the face of the cell in direction type) of level l from the parent cell s of
// grid g, the child at offset o
double cfd_amr::prolong(const fspec &F, int g, const int *s, const int *o)
{
    field &f = F.f(g);

    if(F.type==CELL)
    {
        double v = f(s[0],s[1],s[2]);
        for(int d=0; d<3; ++d)
        v += cell_slope(f,g,s,d,F.mode)*(o[d]==0 ? -0.25 : 0.25);
        return v;
    }

    const int d0 = F.type;

    // the coarse face value of cell c, linear in the transverse directions
    auto T = [&](const int *cc)
    {
        double v = f(cc[0],cc[1],cc[2]);
        for(int t=0; t<3; ++t)
        if(t!=d0)
        v += cell_slope(f,g,cc,t,F.mode)*(o[t]==0 ? -0.25 : 0.25);
        return v;
    };

    if(rr[d0]==2 && o[d0]==0)
    {
        int sm[3] = {s[0],s[1],s[2]};
        sm[d0] -= 1;
        return 0.5*(T(sm) + T(s));
    }
    return T(s);
}

void cfd_amr::fill(int l, vector<fspec> &fs, bool siblings_only)
{
    const int nv = (int)fs.size();

    fill_run(l,nv,
             [&](int kind, int g, const int *s, const int *o, double *v)
             {
                 for(int m=0; m<nv; ++m)
                 v[m] = (kind==0) ? fs[m].f(g)(s[0],s[1],s[2]) : (siblings_only ? 0.0 : prolong(fs[m],g,s,o));
             },
             [&](r3patch *c, int id, const r3fill &f, const double *v)
             {
                 if(siblings_only && f.kind==1)
                 return;
                 for(int m=0; m<nv; ++m)
                 (fs[m].fd ? fs[m].fd : fs[m].f)(id)(f.d[0],f.d[1],f.d[2]) = v[m];
             });
}

void cfd_amr::fill_all(vector<fspec> &fs)
{
    for(int l=1; l<=maxlev; ++l)
    fill(l,fs);
}

// density and viscosity of the cells around patch id from its level set (as fluid_update_fsf, with
// the local interface width)
void cfd_amr::margin_rovisc(int id)
{
    cfd_amr_patch *c = CP(id);
    lexer *pp = c->pp;
    fdm *a = c->a;

    for(const r3fill &f : c->fill)
    {
        if(f.kind==3)
        continue;
        const int i=f.d[0], j=f.d[1], k=f.d[2];
        // the width from the nearest cell whose neighbours are in the field (the outermost ghost
        // layer has no neighbour outside)
        const int m = pp->margin;
        const int iw = std::min(std::max(i,1-m),pp->knox+m-2);
        const int jw = pp->j_dir==1 ? std::min(std::max(j,1-m),pp->knoy+m-2) : j;
        const int kw = std::min(std::max(k,1-m),pp->knoz+m-2);
        const double H = heaviside(a->phi(i,j,k),interface_width(pp,a->phi,iw,jw,kw));
        a->ro(i,j,k) = pp->W1*H + pp->W3*(1.0-H);
        a->visc(i,j,k) = pp->W2*H + pp->W4*(1.0-H);
    }
}

void cfd_amr::bc_state(int id)
{
    cfd_amr_patch *c = CP(id);
    lexer *pp = c->pp;
    fdm *fa = c->a;
    pscope ps(pgc0,fa,a0);

    pgc0->start1(pp,fa->u,10);
    pgc0->start2(pp,fa->v,11);
    pgc0->start3(pp,fa->w,12);
    pgc0->start4(pp,fa->phi,gcval_phi);
    pgc0->start4(pp,fa->press,40);

    if(turb)
    {
        pgc0->start4(pp,fa->eddyv,24);
        c->pturb->gcupdate(pp,fa,pgc0);
    }
    if(solids)
    {
        pgc0->start4a(pp,fa->solid,150);
        pgc0->start4a(pp,fa->topo,150);
    }
}

// --------------------------------------------------------------------- restriction
// cells and faces of the grids of level l-1 under the patches of level l: the mean of the children
void cfd_amr::restrict_level(int l, const vector<std::function<field&(int)>> &cells, const vector<std::function<field&(int)>> &faces)
{
    const double wc = 1.0/double(rr[0]*rr[1]*rr[2]);

    for(int id : lev[l])
    {
        cfd_amr_patch *c = CP(id);

        for(const r3patch::pblock &B : c->par)
        {
            int og[3];
            goff(B.g,og);

            int P[3];
            for(P[0]=B.lo[0]; P[0]<=B.hi[0]; ++P[0])
            for(P[1]=B.lo[1]; P[1]<=B.hi[1]; ++P[1])
            for(P[2]=B.lo[2]; P[2]<=B.hi[2]; ++P[2])
            {
                const int pl[3] = {P[0]-og[0], P[1]-og[1], P[2]-og[2]};
                int f0[3];
                for(int d=0; d<3; ++d)
                f0[d] = P[d]*rr[d] - c->lo[d];

                for(auto &sel : cells)
                {
                    field &src = sel(id);
                    double s = 0.0;
                    for(int a=0; a<rr[0]; ++a)
                    for(int b=0; b<rr[1]; ++b)
                    for(int e=0; e<rr[2]; ++e)
                    s += src(f0[0]+a,f0[1]+b,f0[2]+e);
                    sel(B.g)(pl[0],pl[1],pl[2]) = wc*s;
                }

                for(int d0=0; d0<3; ++d0)
                {
                    if(rr[d0]!=2 || (faces.size()<3 && d0>=(int)faces.size()))
                    continue;
                    field &src = faces[d0](id);
                    field &dst = faces[d0](B.g);
                    const int nt = (rr[0]*rr[1]*rr[2])/2;

                    // the + face (index P) from the children at the high side, the - face (index
                    // P-1) from the children at the low side (their low face: index f0-1)
                    for(int side=0; side<2; ++side)
                    {
                        double s = 0.0;
                        int ff[3];
                        for(int a=0; a<rr[0]; ++a)
                        for(int b=0; b<rr[1]; ++b)
                        for(int e=0; e<rr[2]; ++e)
                        {
                            const int o[3] = {a,b,e};
                            if(o[d0]!=0)
                            continue;
                            for(int d=0; d<3; ++d)
                            ff[d] = f0[d]+o[d];
                            ff[d0] = (side==0) ? f0[d0]+1 : f0[d0]-1;
                            s += src(ff[0],ff[1],ff[2]);
                        }
                        int dl[3] = {pl[0],pl[1],pl[2]};
                        if(side==1)
                        dl[d0] -= 1;
                        dst(dl[0],dl[1],dl[2]) = s/double(nt);

                        // the - face at the low edge of the parent's grid belongs to the grid of
                        // the cell below it (its + face): the coarse leaf cell there takes its
                        // flux from that grid (the copies between the patches of a level go from
                        // the owner to the cells around the other patch)
                        if(side==1 && dl[d0]<0)
                        {
                            int G[3] = {P[0],P[1],P[2]};
                            G[d0] -= 1;
                            int go = -2, lo[3];
                            if(l-1==0)
                            {
                                go = -1;
                                for(int d=0; d<3; ++d)
                                {
                                    lo[d] = G[d]-org[d];
                                    if(lo[d]<0 || lo[d]>=kn0[d])
                                    go = -2;
                                }
                            }
                            else
                            for(int g : lev[l-1])
                            {
                                r3patch *q = this->P[g];
                                if(G[0]>=q->lo[0] && G[0]<=q->hi[0] && G[1]>=q->lo[1] && G[1]<=q->hi[1] && G[2]>=q->lo[2] && G[2]<=q->hi[2])
                                {
                                    go = g;
                                    for(int d=0; d<3; ++d)
                                    lo[d] = G[d]-q->lo[d];
                                    break;
                                }
                            }
                            if(go>=-1 && go!=B.g)
                            faces[d0](go)(lo[0],lo[1],lo[2]) = s/double(nt);
                        }
                    }
                }
            }
        }
    }
}

// --------------------------------------------------------------------- initialisation
// the patches from the initialised level 0: geometry fields of initialize, then velocities, level
// set, pressure prolonged into the patch and its surrounding cells, density and viscosity
// the turbulence model on the patches of level l (the cells around them filled)
void cfd_amr::turb_step(int l)
{
    for(int id : lev[l])
    {
        cfd_amr_patch *c = CP(id);
        pscope ps(pgc0,c->a,a0);
        c->pturb->start(c->a,c->pp,c->pturbdisc,c->pturbdiff,c->psolv,pgc0,c->pflow,c->pvrans);
    }
}

// --------------------------------------------------------------------- state fields
field& cfd_amr::kfield(int g, int m)
{
    rans_io *r = dynamic_cast<rans_io*>(g<0 ? pturb0 : CP(g)->pturb);
    return (m==0) ? r->kin : r->eps;
}

vector<int> cfd_amr::state_ids()
{
    vector<int> ids = {0,1,2,3,4,5};
    if(turb)
    {
        ids.push_back(6);
        ids.push_back(7);
        ids.push_back(10);
    }
    if(solids)
    {
        ids.push_back(8);
        ids.push_back(9);
    }
    return ids;
}

field& cfd_amr::sfield(fdm *a, turbulence *t, int m)
{
    if(m==0) return a->phi;
    if(m==1) return a->press;
    if(m==2) return a->eddyv;
    if(m==3) return a->u;
    if(m==4) return a->v;
    if(m==5) return a->w;
    if(m==6) return dynamic_cast<rans_io*>(t)->kin;
    if(m==7) return dynamic_cast<rans_io*>(t)->eps;
    if(m==8) return a->solid;
    if(m==9) return a->topo;
    return a->walld;
}

vector<cfd_amr::fspec> cfd_amr::state_specs()
{
    vector<fspec> fs;
    for(int m : state_ids())
    {
        const int type = (m>=3 && m<=5) ? m-3 : CELL;
        const int mode = (m==2 || (m>=3 && m<=7)) ? 1 : 0;
        fs.push_back({type,mode,[this,m](int g) -> field& { return sfield(gfd(g),g<0 ? pturb0 : CP(g)->pturb,m); }});
    }
    return fs;
}

// the fields of patch id before its data: initialize, ioflow, no porous media, solids or bodies
void cfd_amr::patch_fields_ini(int id)
{
    lexer *p = p0;
    ghostcell *pgc = pgc0;
    cfd_amr_patch *c = CP(id);
    lexer *pp = c->pp;
    pscope ps(pgc,c->a,a0);
    const double psi = pp->psi, psi0 = pp->psi0;
    pini0->start(c->a,pp,pgc);
    c->pflow->ini(pp,c->a,pgc);
    pp->psi = psi;      // (initialize sets it)
    pp->psi0 = psi0;
    if(turb)
    c->pturb->ini(pp,c->a,pgc);
    if(solids)
    {
        pp->solidread = p->solidread;
        pp->toporead = p->toporead;
        pp->topoforcing = p->topoforcing;
    }

    // no porous media, solids or bodies: also in the cells around the patch (initialize sets
    // the interior; PORVAL at the patch edge reads the cell outside)
    const int na = pp->imax*pp->jmax*pp->kmax;
    for(int n=0; n<na; ++n)
    {
        c->a->porosity.V[n] = 1.0;
        c->a->fb.V[n] = 1.0;
        c->a->topo.V[n] = 1.0;
        c->a->solid.V[n] = 1.0e8;
    }

    pp->maxlength = p->maxlength;
    pp->xcoormax = p->xcoormax; pp->xcoormin = p->xcoormin;
    pp->ycoormax = p->ycoormax; pp->ycoormin = p->ycoormin;
    pp->zcoormax = p->zcoormax; pp->zcoormin = p->zcoormin;
    pp->phimean = p->phimean;
    pp->count = p->count;
}

void cfd_amr::ini(lexer *p, fdm *a, ghostcell *pgc)
{
    for(int id=0; id<(int)P.size(); ++id)
    patch_fields_ini(id);

    // interior of every patch from its parent, level by level (with the geometry: solid, topo,
    // wall distance)
    vector<fspec> fs = state_specs();

    for(int l=1; l<=maxlev; ++l)
    {
        for(int id : lev[l])
        {
            cfd_amr_patch *c = CP(id);
            for(const r3patch::pblock &B : c->par)
            {
                int og[3];
                goff(B.g,og);
                int I[3];
                for(I[0]=c->lo[0]; I[0]<=c->hi[0]; ++I[0])
                for(I[1]=c->lo[1]; I[1]<=c->hi[1]; ++I[1])
                for(I[2]=c->lo[2]; I[2]<=c->hi[2]; ++I[2])
                {
                    int Pg[3], s[3], o[3];
                    bool in = true;
                    for(int d=0; d<3; ++d)
                    {
                        Pg[d] = fdiv(I[d],rr[d]);
                        o[d] = I[d]-Pg[d]*rr[d];
                        s[d] = Pg[d]-og[d];
                        if(Pg[d]<B.lo[d] || Pg[d]>B.hi[d])
                        in = false;
                    }
                    if(!in)
                    continue;
                    for(auto &F : fs)
                    F.f(id)(I[0]-c->lo[0],I[1]-c->lo[1],I[2]-c->lo[2]) = prolong(F,B.g,s,o);
                }
            }
        }
        fill(l,fs);

        // the boundary values: the cells outside the domain are the neighbours of the slopes of
        // the prolongation into the next level
        for(int id : lev[l])
        bc_state(id);
    }

    // density, viscosity; boundary values
    for(int id=0; id<(int)P.size(); ++id)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        fdm *fa = c->a;
        pscope ps(pgc,fa,a0);

        for(int i=0; i<pp->knox; ++i)
        for(int j=0; j<pp->knoy; ++j)
        for(int k=0; k<pp->knoz; ++k)
        {
            const double H = heaviside(fa->phi(i,j,k),interface_width(pp,fa->phi,i,j,k));
            fa->ro(i,j,k) = pp->W1*H + pp->W3*(1.0-H);
            fa->visc(i,j,k) = pp->W2*H + pp->W4*(1.0-H);
        }
        margin_rovisc(id);

        pgc->start1(pp,fa->u,10);
        pgc->start2(pp,fa->v,11);
        pgc->start3(pp,fa->w,12);
        pgc->start4(pp,fa->phi,gcval_phi);
        pgc->start4(pp,fa->press,40);
        pgc->start4(pp,fa->ro,1);
        pgc->start4(pp,fa->visc,1);

        c->ppress->ini(pp,fa,pgc);

        // solids and topography: the direct forcing fields and flags of the patch
        if(solids)
        pgc->solid_forcing_ini(pp,fa);
        if(bodies)
        c->p6dof->body(pp,fa,pgc);
        if(turb)
        {
            c->pturb->gcupdate(pp,fa,pgc);
            // the old time level of kin, eps (IM1: the time term) from the initial state, as
            // level 0 at the end of driver_ini_cfd
            c->pturb->ktimesave(pp,fa,pgc);
            c->pturb->etimesave(pp,fa,pgc);
        }
    }

    timestep(p,a,pgc);

    // the top and the bottom of a refined region should stay out of the interface band (the
    // patches across the free surface at their sides are fine): hint at the start
    {
        double zb = 1.0e20;
        for(const cfd_amr_cf &F : cf)
        if(F.d==2 && fabs(CP(F.pid)->a->phi.V[F.qf])<1.5*p->psi)
        zb = MIN(zb,CP(F.pid)->pp->ZP[F.fi[2]+marge]);
        zb = pgc->globalmin(zb);
        if(p->mpirank==0 && zb<1.0e19)
        cout<<"CFD AMR: the free surface is close to the top or the bottom of a refined box (z = "<<zb
            <<"); keep them more than the interface thickness ("<<p->psi<<") plus the wave amplitude away from it"<<endl;
    }

    // the output of the patches at the start (printer_CFD has written level 0)
    print(p,true);
}

// --------------------------------------------------------------------- time step
// the step size of the finest level (all levels step together): the criterion of ietimestep on
// every patch, the minimum with level 0
void cfd_amr::timestep(lexer *p, fdm *a, ghostcell *pgc)
{
    double cu = 1.0e10;
    double vpat[5] = {0.0,0.0,0.0,0.0,0.0};     // largest velocity on the patches: value, x, y, z, level

    for(int id=0; id<(int)P.size(); ++id)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        fdm *fa = c->a;
        double umax=0.0, vmax=0.0, wmax=0.0;
        double cp = 1.0e10;

        // the velocities of the patch and of the layer of cells around it (filled from the coarser
        // level or a sibling: what can enter the patch within a step)
        const int jl = (pp->j_dir==1) ? 1 : 0;
        for(int i=-1; i<=pp->knox; ++i)
        for(int j=-jl; j<pp->knoy+jl; ++j)
        for(int k=-1; k<=pp->knoz; ++k)
        {
            const int q = cix(pp,i,j,k);
            const bool in = (i>=0 && i<pp->knox && j>=0 && j<pp->knoy && k>=0 && k<pp->knoz);
            const double au = (pp->flag1[q]>0) ? fabs(fa->u(i,j,k)) : 0.0;
            const double av = (pp->flag2[q]>0) ? fabs(fa->v(i,j,k)) : 0.0;
            const double aw = (pp->flag3[q]>0) ? fabs(fa->w(i,j,k)) : 0.0;
            umax = MAX(umax,au);
            vmax = MAX(vmax,av);
            wmax = MAX(wmax,aw);
            if(in && MAX3(au,av,aw)>vpat[0])
            {
                vpat[0] = MAX3(au,av,aw);
                vpat[1] = pp->XP[i+marge]; vpat[2] = pp->YP[j+marge]; vpat[3] = pp->ZP[k+marge];
                vpat[4] = c->lev;
            }
        }


        const double mF = MAX3(fa->maxF,fa->maxG,fa->maxH);
        double vel = sqrt(umax*umax + vmax*vmax + wmax*wmax);

        // the first step (from rest): the velocity a fall over the height of the domain gives (the
        // flow starts from a pressure jump, which the finer grids feel within the first step)
        if(p->count==0)
        vel = MAX(vel,sqrt(2.0*sqrt(p->W20*p->W20+p->W21*p->W21+p->W22*p->W22)*(p->zcoormax-p->zcoormin)));

        for(int i=0; i<pp->knox; ++i)
        for(int j=0; j<pp->knoy; ++j)
        for(int k=0; k<pp->knoz; ++k)
        {
            if(pp->flag4[cix(pp,i,j,k)]<=0)
            continue;
            const double dx = MIN3(pp->DXN[i+marge],pp->DYN[j+marge],pp->DZN[k+marge]);
            if(p->N50==2)
            {
                cp = MIN(cp,2.0/(umax/pp->DXN[i+marge] + sqrt(4.0*fabs(fa->maxF)/pp->DXN[i+marge])));
                cp = MIN(cp,2.0/(wmax/pp->DZN[k+marge] + sqrt(4.0*fabs(fa->maxH)/pp->DZN[k+marge])));
                if(pp->j_dir==1)
                cp = MIN(cp,2.0/(vmax/pp->DYN[j+marge] + sqrt(4.0*fabs(fa->maxG)/pp->DYN[j+marge])));
            }
            else
            cp = MIN(cp,2.0/(vel/dx + sqrt(4.0*fabs(mF)/dx)));
        }

        // subcycling: the patch takes 2^l steps per step of level 0
        if(sub)
        cp *= double(1<<c->lev);
        cu = MIN(cu,cp);

        fa->maxF = fa->maxG = fa->maxH = 0.0;
    }

    double g = -cu;
    MPI_Allreduce(MPI_IN_PLACE,&g,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cu = -g;

    // the largest velocity of the level-0 leaf cells (diagnostics)
    double v0[5] = {0.0,0.0,0.0,0.0,0.0};
    for(int i=0; i<p->knox; ++i)
    for(int j=0; j<p->knoy; ++j)
    for(int k=0; k<p->knoz; ++k)
    {
        const int q = cix(p,i,j,k);
        const int I[3] = {i+org[0],j+org[1],k+org[2]};
        if(p->flag4[q]<=0 || covered(0,I))
        continue;
        const double vv = MAX3(fabs(a->u(i,j,k)),fabs(a->v(i,j,k)),fabs(a->w(i,j,k)));
        if(vv>v0[0])
        {
            v0[0] = vv;
            v0[1] = p->XP[i+marge]; v0[2] = p->YP[j+marge]; v0[3] = p->ZP[k+marge];
        }
    }
    {
        struct { double v; int r; } in = {v0[0],myrank}, out;
        MPI_Allreduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
        MPI_Bcast(v0,5,MPI_DOUBLE,out.r,MPI_COMM_WORLD);
        if(myrank==0 && (p->count%p->P12==0))
        cout<<"CFD AMR: largest velocity on level 0 (leaf cells) "<<setprecision(3)<<v0[0]<<" at ("<<v0[1]<<", "<<v0[2]<<", "<<v0[3]<<")"<<endl;
    }

    {
        struct { double v; int r; } in = {vpat[0],myrank}, out;
        MPI_Allreduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
        MPI_Bcast(vpat,5,MPI_DOUBLE,out.r,MPI_COMM_WORLD);
        if(myrank==0 && (p->count%p->P12==0))
        cout<<"CFD AMR: largest velocity on the patches "<<setprecision(3)<<vpat[0]<<" at ("<<vpat[1]<<", "<<vpat[2]<<", "<<vpat[3]<<"), level "<<int(vpat[4])<<endl;
    }

    p->dt = MIN(p->dt,p->N47*cu);

    for(auto q : P)
    {
        q->pp->dt = p->dt;
        q->pp->dt_old = p->dt_old;
    }
}

// --------------------------------------------------------------------- the step
// the stage input of the patches, level by level: velocities, level set, pressure, eddy viscosity
// in the cells around the patches; density and viscosity there from the level set
void cfd_amr::fill_stage_in(int l, int s)
{
    vector<fspec> fs;
    fs.push_back({0,1,[this,s](int g) -> field& { return gmom(g)->amr_vel(gfd(g),0,s); }});
    fs.push_back({1,1,[this,s](int g) -> field& { return gmom(g)->amr_vel(gfd(g),1,s); }});
    fs.push_back({2,1,[this,s](int g) -> field& { return gmom(g)->amr_vel(gfd(g),2,s); }});
    fs.push_back({CELL,0,[this](int g) -> field& { return gfd(g)->phi; }});
    fs.push_back({CELL,0,[this](int g) -> field& { return gfd(g)->press; }});
    fs.push_back({CELL,1,[this](int g) -> field& { return gfd(g)->eddyv; }});

    // the level-set field of the stage input: in the first stage momentum_rk copies phi into it
    // (the patch interior), the cells around the patch get the phi of the coarser level
    if(s==0)
    fs.push_back({CELL,0,[this](int g) -> field& { return gfd(g)->phi; },
                         [this](int g) -> field& { return gmom(g)->amr_phi_in(0); }});
    else
    fs.push_back({CELL,0,[this,s](int g) -> field& { return gmom(g)->amr_phi_in(s); }});

    fill(l,fs);

    for(int id : lev[l])
    margin_rovisc(id);
}

// the level set of the stage output around the patches of level l (the coarser level has finished
// its reinitialisation of the stage)
void cfd_amr::fill_phi_out(int l, int s)
{
    vector<fspec> fs;
    fs.push_back({CELL,0,[this,s](int g) -> field& { return gmom(g)->amr_phi_out(s); }});
    fill(l,fs);
}

// the Dirichlet values of the implicit diffusion around the patches of level l: the stage velocity
// of the coarser level (after its momentum, before the projection); next to a patch of the same
// level, whose momentum of the stage is not done yet, its stage input
void cfd_amr::fill_velout(int l, int s)
{
    vector<fspec> fs;
    for(int c=0; c<3; ++c)
    fs.push_back({c,1,[this,s,c,l](int g) -> field& { return (g>=0 && P[g]->lev==l) ? gmom(g)->amr_vel(gfd(g),c,s) : gmom(g)->amr_velout(gfd(g),c,s); },
                      [this,s,c](int g) -> field& { return gmom(g)->amr_velout(gfd(g),c,s); }});
    fill(l,fs);
}

// the projected stage velocity around the patches of all levels
void cfd_amr::fill_vel_after(int s)
{
    vector<fspec> fs;
    fs.push_back({0,1,[this,s](int g) -> field& { return gmom(g)->amr_velout(gfd(g),0,s); }});
    fs.push_back({1,1,[this,s](int g) -> field& { return gmom(g)->amr_velout(gfd(g),1,s); }});
    fs.push_back({2,1,[this,s](int g) -> field& { return gmom(g)->amr_velout(gfd(g),2,s); }});
    fill_all(fs);
}

// the Dirichlet values of the implicit diffusion at the high patch edges of level l: the coarse
// face velocity, linear in the transverse directions (G 7 0: the stage velocity of the coarser
// level after its momentum; G 7 1: the coarser level at the output time of the stage)
void cfd_amr::edge_faces(int l, int s)
{
    for(int id : lev[l])
    {
        cfd_amr_patch *c = CP(id);
        const size_t na = (size_t)c->pp->imax*c->pp->jmax*c->pp->kmax;
        for(int d=0; d<3; ++d)
        c->edge[d].assign(na,std::numeric_limits<double>::quiet_NaN());
    }

    for(cfd_amr_cf &F : cf)
    {
        cfd_amr_patch *c = CP(F.pid);
        if(c->lev!=l || F.dir<0)
        continue;
        field &fc = sub ? tsf(F.gface,F.d) : gmom(F.gface)->amr_velout(gfd(F.gface),F.d,s);
        double v = fc(F.cf[0],F.cf[1],F.cf[2]);
        for(int t=0; t<3; ++t)
        if(t!=F.d)
        v += cell_slope(fc,F.gface,F.cf,t,1)*(F.o[t]==0 ? -0.25 : 0.25);
        c->edge[F.d][cix(c->pp,F.ff[0],F.ff[1],F.ff[2])] = v;
    }
}

// restriction of the projected stage velocity (finest level first; the faces around the patches
// of a level are filled again before that level is restricted)
void cfd_amr::restrict_vel(int s)
{
    vector<std::function<field&(int)>> none;
    vector<std::function<field&(int)>> faces;
    for(int c=0; c<3; ++c)
    faces.push_back([this,s,c](int g) -> field& { return gmom(g)->amr_velout(gfd(g),c,s); });

    vector<fspec> fs;
    for(int c=0; c<3; ++c)
    fs.push_back({c,1,[this,s,c](int g) -> field& { return gmom(g)->amr_velout(gfd(g),c,s); }});

    // after the restriction onto a level of patches, their copies of each other again; the faces
    // at their edges to the coarser level (in the cells around them) keep the projected values
    for(int l=maxlev; l>=1; --l)
    {
        restrict_level(l,none,faces);
        if(l>=2)
        fill(l-1,fs,true);
    }

    pgc0->start1(p0,mom0->amr_velout(a0,0,s),10);
    pgc0->start2(p0,mom0->amr_velout(a0,1,s),11);
    pgc0->start3(p0,mom0->amr_velout(a0,2,s),12);
}

// restriction of the level set, pressure, density, viscosity at the end of the stage
void cfd_amr::restrict_scalars(int s)
{
    vector<std::function<field&(int)>> cells, none;
    cells.push_back([this](int g) -> field& { return gfd(g)->phi; });
    cells.push_back([this,s](int g) -> field& { return gmom(g)->amr_phi_out(s); });
    cells.push_back([this](int g) -> field& { return gfd(g)->press; });
    cells.push_back([this](int g) -> field& { return gfd(g)->ro; });
    cells.push_back([this](int g) -> field& { return gfd(g)->visc; });
    if(turb)
    {
        cells.push_back([this](int g) -> field& { return gfd(g)->eddyv; });
        cells.push_back([this](int g) -> field& { return kfield(g,0); });
        cells.push_back([this](int g) -> field& { return kfield(g,1); });
    }

    for(int l=maxlev; l>=1; --l)
    restrict_level(l,cells,none);

    if(turb)
    pturb0->gcupdate(p0,a0,pgc0);

    pgc0->start4(p0,a0->phi,gcval_phi);
    pgc0->start4(p0,mom0->amr_phi_out(s),gcval_phi);
    pgc0->start4(p0,a0->press,40);
    pgc0->start4(p0,a0->ro,1);
    pgc0->start4(p0,a0->visc,1);
}

void cfd_amr::step(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof)
{
    if(sub)
    {
        step_sub(p,a,pgc,pvrans,p6dof);
        return;
    }

    const double t0 = MPI_Wtime();

    // composite volume of the level set at the start of the step (F 46 2, 3)
    if(p->F46==2 || p->F46==3)
    {
        double V, dV;
        volume(V,dV);
        vol_step = V;
        if(!vol_ini)
        {
            vol_net = V;
            vol_ini = true;
        }
    }

    for(auto q : P)
    {
        lexer *pp = q->pp;
        pp->dt = p->dt;
        pp->dt_old = p->dt_old;
        pp->simtime = p->simtime;
        pp->count = p->count;
    }

    // the turbulence model on the patches (level 0: the driver), the cells around them from the
    // coarser level
    if(turb)
    for(int l=1; l<=maxlev; ++l)
    {
        vector<fspec> fs;
        fs.push_back({CELL,1,[this](int g) -> field& { return gfd(g)->eddyv; }});
        fs.push_back({CELL,1,[this](int g) -> field& { return kfield(g,0); }});
        fs.push_back({CELL,1,[this](int g) -> field& { return kfield(g,1); }});
        fill(l,fs);
        turb_step(l);
    }

    mom0->amr_step_begin(p,a,pgc);
    for(int id=0; id<(int)P.size(); ++id)
    {
        pscope ps(pgc,CP(id)->a,a0);
        CP(id)->pmom->amr_step_begin(CP(id)->pp,CP(id)->a,pgc);
    }

    const int stages = mom0->amr_stages();

    for(int s=0; s<stages; ++s)
    {
        double ta = MPI_Wtime();

        // level set
        mom0->amr_ls_transport(p,a,pgc,s);
        mom0->amr_ls_finish(p,a,pgc,s,false,-1);

        for(int l=1; l<=maxlev; ++l)
        {
            fill_stage_in(l,s);

            for(int id : lev[l])
            {
                pscope ps(pgc,CP(id)->a,a0);
                CP(id)->pmom->amr_ls_transport(CP(id)->pp,CP(id)->a,pgc,s);
            }

            fill_phi_out(l,s);

            // reinitialisation one iteration at a time, the cells around the patches filled in
            // between (patches of the same level next to each other)
            const int iters = mom0->amr_reini_iters(p,s);
            for(int it=0; it<iters; ++it)
            {
                if(it>0)
                fill_phi_out(l,s);

                for(int id : lev[l])
                {
                    pscope ps(pgc,CP(id)->a,a0);
                    if(it==0)
                    CP(id)->pmom->amr_ls_finish(CP(id)->pp,CP(id)->a,pgc,s,false,1);
                    else
                    CP(id)->pmom->amr_ls_reini(CP(id)->pp,CP(id)->a,pgc,s,1);
                }
            }
        }

        double tb = MPI_Wtime();
        tm[0] += tb-ta;

        // momentum
        mom0->amr_momentum(p,a,pgc,pvrans,p6dof,s,1);

        for(int l=1; l<=maxlev; ++l)
        {
            fill_velout(l,s);
            edge_faces(l,s);

            for(int id : lev[l])
            {
                pscope ps(pgc,CP(id)->a,a0);
                CP(id)->pmom->amr_momentum(CP(id)->pp,CP(id)->a,pgc,CP(id)->pvrans,(sixdof*)CP(id)->p6dof,s,patch_forcing());
            }
        }

        double tc = MPI_Wtime();
        tm[1] += tc-tb;

        // composite projection
        pr_project(p,pgc,s);

        double td = MPI_Wtime();
        tm[2] += td-tc;

        // end of the stage
        mom0->amr_stage_end(p,a,pgc,s);
        for(int id=0; id<(int)P.size(); ++id)
        {
            pscope ps(pgc,CP(id)->a,a0);
            CP(id)->pmom->amr_stage_end(CP(id)->pp,CP(id)->a,pgc,s);
        }

        restrict_scalars(s);

        tm[3] += MPI_Wtime()-td;
    }

    // the volume correction over the leaf cells of all grids (F 46 2, 3)
    volume_correct(p,pgc);

    if(turb)
    for(int id=0; id<(int)P.size(); ++id)
    {
        pscope ps(pgc,CP(id)->a,a0);
        CP(id)->pturb->ktimesave(CP(id)->pp,CP(id)->a,pgc);
        CP(id)->pturb->etimesave(CP(id)->pp,CP(id)->a,pgc);
    }

    // adaptive refinement: regrid every G 2 steps
    if(adapt && p->G2>0 && p->count%p->G2==0)
    regrid(p,pgc);

    tm[7] += MPI_Wtime()-t0;


    if(p->count%p->P12==0)
    {
        double V, dV;
        volume(V,dV);
        if(p->mpirank==0)
        {
            const std::streamsize pr = cout.precision();
            cout<<"CFD AMR: volume of the level set over the leaf cells "<<setprecision(8)<<V<<endl;
            cout.precision(pr);
        }
    }

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"CFD AMR: pressure iterations "<<pr_it_last<<"  res "<<setprecision(3)<<pr_res_last
        <<"  time: level set "<<tm[0]<<" momentum "<<tm[1]<<" projection "<<tm[2]<<" end "<<tm[3]<<" total "<<tm[7]<<" (solver: prec "<<tm[5]<<" apply "<<tm[6]<<" patch cycles "<<tp[0]<<" level-0 cycles "<<tp[1]<<")"<<endl;
}
