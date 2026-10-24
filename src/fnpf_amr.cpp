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

#include"reefckpt.h"
#include"fnpf_amr.h"
#include"lexer.h"
#include"bc_zone.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"slice4.h"
#include"fnpf_fsfbc.h"
#include"fnpf_fsfbc_wd.h"
#include"fnpf_sigma.h"
#include"fnpf_fsf_update.h"
#include"fnpf_bed_update.h"
#include"fnpf_laplace_cds2.h"
#include"reefmg_core.h"
#include"fnpf_amr_fill.h"
#include"fnpf_body.h"
#include"fnpf_6DOF.h"
#include"6DOF_obj.h"
#include<cmath>
#include<algorithm>
#include<mpi.h>
#include<iomanip>
#include<cstdio>
#include<sys/stat.h>
#include<sys/types.h>

namespace
{
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

// inside a patch kernel the global reductions return the local value (the patch lies on one rank)
struct fnpf_amr_local
{
    ghostcell *g;
    bool old;
    fnpf_amr_local(ghostcell *gg) : g(gg) { old = g->set_local(true); }
    ~fnpf_amr_local() { g->set_local(old); }
};

typedef reefamr_comms_off comms_off;

// the Laplace solver seen by fnpf_RK3 with mesh refinement: level 0 is assembled by its own
// solver (assemble_only), then fnpf_amr solves all grids together
class fnpf_laplace_amr : public fnpf_laplace
{
public:
    fnpf_laplace_amr(fnpf_amr *a, fnpf_laplace *outer, fnpf_laplace_cds2 *plap0) : a(a), outer(outer), plap0(plap0) {}

    void start(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf, double *f, slice &Fifsf) override
    {
        // without patches, or subcycled (G 7 1): level 0 alone with its own solver
        if(!a->active() || a->sub_on())
        {
            outer->start(p,c,pgc,psolv,pf,f,Fifsf);
            return;
        }

        a->lap_solve(p,c,pgc,psolv,pf,f,Fifsf);
    }

private:
    fnpf_amr *a;
    fnpf_laplace *outer;
    fnpf_laplace_cds2 *plap0;
};
}

fnpf_amr::fnpf_amr(lexer *p, fdm_fnpf *c, ghostcell *pgc) : reefamr(p,pgc)
{
    c0 = c;

    reefamr_param q;
    q.name = "FNPF AMR";
    q.maxlev = p->G1;
    q.regrid = 0;
    q.nbuf = MAX(p->G3,0);
    q.tile = MAX(p->G4,4);
    q.tile += q.tile%2;
    q.keep = MIN(MAX(p->G5,0),250);

    // cells computed beyond the patch box: the free-surface derivatives (WENO5, CDS4) reach
    // three cells; the patch arrays reach EXT + margin fine cells beyond the patch, which the
    // coarser level has to cover
    q.ext = 3;
    q.nest = MAX(3,(q.ext+p->margin+1)/2);

    // vertical refinement: every refined level doubles the sigma layers (G 6 1)
    vref = (p->G6==1) ? 2 : 1;
    q.vref.assign(q.maxlev+1,vref);

    for(int k=0; k<p->G10; ++k)
    {
        q.rbox.push_back(p->G10_xs[k]); q.rbox.push_back(p->G10_xe[k]);
        q.rbox.push_back(p->G10_ys[k]); q.rbox.push_back(p->G10_ye[k]);
    }
    for(int k=0; k<p->G11; ++k)
    {
        q.fbox.push_back(p->G11_xs[k]); q.fbox.push_back(p->G11_xe[k]);
        q.fbox.push_back(p->G11_ys[k]); q.fbox.push_back(p->G11_ye[k]);
    }

    // no refinement in the relaxation zones of the wave generation (B 98 2) and the numerical
    // beach (B 99 1, 2), measured from the ends of the domain in x, and in the B 520 zones (bc_zone)
    bc_zone_set::from_legacy(p,pgc).norefine_boxes(p,q.fbox);
    q.ioband = 4;

    // refinement around the resolved body (X 10 1): margin G 12 around the wetted hull,
    // rebuilt every G 2 steps so that the zone follows the body (G 2 0: static)
    q.zones = (p->X10==1 && p->G12>0);
    q.zr = p->G12_r;
    q.zL = p->G13_L;
    q.za = p->G13_a;

    // the FNPF bodies are moored or oscillate: the zone is aligned with x and y (the heading of
    // the SFLOW ship zone swings with every reversal of the motion), and the layout is kept as
    // long as it covers the flagged tiles with at most 50 % excess
    q.zalign = true;
    q.lazy = 1.5;

    // subcycling (G 7 1, fnpf_amr_sub.cpp): level l takes 2^l steps per level-0 step; bodies with
    // the zone around the hull (G 12: the finest level advances the body), its margin at least
    // 12 level-0 cells (the coarse levels step first with their own coarse picture of the body)
    // wetting and drying (A 343 1-3): the patches are fully wet. No patch within q.dryband
    // level-0 cells of a dry cell, a thin-film cell (WL < 10 A 344) or the coastline damping band
    // (3 times A 341 dx or A 342, as fnpf_fsfbc_wd), checked at every regrid (cell_unfit); a static
    // box (G 10) loses the tiles there.  The band covers the fills around a patch (EXT fine cells
    // and the bicubic stencil); with run-up and run-down (A 343 2, 3) it is doubled
    // G 30 1: the patches may cover the coastline (wd, dryband off): the static coastline (A 343 1)
    // or the dynamic wetting and drying (A 343 2, 3: dyn)
    wd = (p->A343>=1 && p->G30==1);
    dyn = (wd && p->A343>=2);
    wfill = wd;

    // dyn: the cells around a patch and the covered cells of a coarser grid keep the flag of their
    // source cells through the wetting and drying of their grid (lexer::wetfix)
    if(dyn)
    {
        wfix0.assign(p->imax*p->jmax,-1);
        p->wetfix = wfix0.data();
    }

    if(p->A343>=1 && !wd)
    {
        q.dryband = (p->A343>=2) ? 8 : 4;

        double d3 = 0.0;
        if(p->A341>0.0)
        d3 = p->A341*0.5*(p->DXD+p->DYD);
        if(p->A342>0.0)
        d3 = p->A342;
        dist_coast = 3.0*d3;
    }

    // wave breaking (A 350 1 Baquet, 2 Kennedy, 3 Romero): on every grid (detection and the implicit
    // damping with the breaking viscosity), the cells around a patch take the viscosity of their source
    // grid.  The detection runs on every grid with its own state (Kennedy: onset time and flag,
    // Romero: local wavenumber and onset time), a new patch starts without breaking
    brk = (p->A350>=1 && p->A350<=3);

    // solution-adaptive refinement (as nhflow_amr, tag): a wet cell of level l flags level l+1
    // where the surface jumps by more than G 20 to a neighbour or its second difference along x or
    // y exceeds G 21 (m); where it breaks (G 23 1 with A 350 1: the breaking flag of the grid, not
    // the coastline viscosity); within G 22 cells (1-3) of a dry cell (with G 30 1: the patches
    // cover the static shoreline).  Without G 30 1 the dryband keeps them off the coast as before
    tol_eta = MAX(p->G20,0.0);
    tol_curv = MAX(p->G21,0.0);
    adaptive = (tol_eta>0.0 || tol_curv>0.0);

    if(wd && p->G22>0)
    {
    nshore = MIN(p->G22,3);
    adaptive = true;
    }

    if(brk && p->G23==1)
    {
    flagbreak = true;
    adaptive = true;
    }

    // the zone follows the body, the patches follow the flags: new patches every G 2 steps
    // (G 2 0: static)
    if(q.zones || adaptive)
    q.regrid = MAX(p->G2,0);

    sub = (p->G7==1 && (p->X10==0 || q.zones)) ? 1 : 0;
    zr_user = q.zr;
    if(sub==1 && q.zones)
    q.zr = MAX(q.zr,12.0*p->DXM);

    // several ranks (G 40 1): the patches are placed for the load of the ranks; parents and old
    // patches on other ranks are reached through the block plans and old_run of the core, the
    // hull triangles are taken by the rank of the finest grid at their centroid (owns_point)
    q.place = (p->G40>=1 && p->mpi_size>1) ? MIN(p->G40,2) : 0;
    q.rebalance = 0.1;
    q.place_whole_zones = true;

    configure(q);

    gcval_eta = 55;
    gcval_fifsf = 60;

    Se0 = &c0->eta;
    Sf0 = &c0->Fifsf;
    stg = 2;

    printtime_amr = 0.0;
    printcount_amr = 0;
    lap_it_total = lap_solves = 0;
    lap_it_last = 0;
    lap_res_last = 0.0;
    for(int k=0; k<9; ++k)
    tm[k] = 0.0;
}

fnpf_amr::~fnpf_amr()
{
    free_patches();

    for(auto v : kv0)
    delete [] v;
    delete mg0;
}

fnpf_laplace* fnpf_amr::laplace(fnpf_laplace *outer, fnpf_laplace_cds2 *plap0)
{
    if(maxlev<1)
    return outer;

    fnpf_amr::plap0 = plap0;

    return new fnpf_laplace_amr(this,outer,plap0);
}

// the run settings with which ini keeps the refinement on
bool fnpf_amr::keys_ok(lexer *p)
{
    int ok=1;
    if(p->A310!=3) ok=0;
    if(p->A350<0 || p->A350>3) ok=0;
    if(p->X10!=0 && p->X10!=1) ok=0;
    if(p->A380!=0) ok=0;
    if(p->A311==3) ok=0;
    if(p->j_dir!=1 && (p->geo_mb>0 || p->geo_curv==1)) ok=0;     // 2D grid (one cell in y): Cartesian

    return ok==1;
}

void fnpf_amr::attach_body(fnpf_body *b)
{
    body = (b!=nullptr && b->present()) ? b : nullptr;

    // the hull triangles (X 185) are sized for the finest level (sixdof_obj::amr_hfac, as
    // nhflow_amr and the CFD mesh refinement): every triangle is integrated on the finest grid
    // that holds its centroid, so on a patch the loads are summed over patch-sized triangles.
    // Set here, before fnpf_6DOF::initialize builds them (inidisc_step2, before ini); only when
    // ini keeps the refinement on.  With G 6 1 the finest level also has 2^G1 times the sigma
    // layers: the vertical spacing is scaled as well (amr_vfac), as on a uniform grid with those layers.
    if(body!=nullptr && maxlev>0 && keys_ok(p0))
    {
        vector<sixdof_obj*> obj;
        body->amr_bodies(obj);

        for(auto o : obj)
        {
        o->amr_hfac = 1.0/double(1<<maxlev);
        o->amr_vfac = (vref==2) ? 1.0/double(1<<maxlev) : 1.0;
        }
    }
}

void fnpf_amr::zone_bodies(vector<sixdof_obj*> &obj)
{
    if(body!=nullptr)
    body->amr_bodies(obj);
}

fnpf_fsf* fnpf_amr::patch_fsf(int n)
{
    return FP(n)->pf;
}

slice& fnpf_amr::patch_tendency(int n, int m)
{
    return (m==0) ? static_cast<slice&>(*FP(n)->ek) : static_cast<slice&>(*FP(n)->fk);
}

void fnpf_amr::prolong_cols(double **f, const vector<char> &need)
{
    auto sel = [&](int g) -> double* { return f[g+1]; };
    vector<reefamr_patch*> nd;
    for(int n=0; n<(int)P.size(); ++n)
    if(need[n])
    nd.push_back(P[n]);
    auto need_fn = [&](reefamr_patch *c) { return std::find(nd.begin(),nd.end(),c)!=nd.end(); };
    for(int l=1; l<=maxlev; ++l)
    prolong_interior_col(l,need_fn,sel,sel);
}

int fnpf_amr::klev(int l) const
{
    return p0->knoz*((vref==2) ? (1<<l) : 1);
}

bool fnpf_amr::owns_point(double x, double y, int id)
{
    // placed patches: the rank of the finest grid at (x,y), then its local grid
    if(par.place>0)
    return finest_rank(x,y)==p0->mpirank && finest_at(x,y)==id;

    if(!(x >= p0->originx && x < p0->endx))
    return false;
    if(p0->j_dir==1 && !(y >= p0->originy && y < p0->endy))
    return false;
    return finest_at(x,y)==id;
}

// the finest local grid whose interior holds (x,y)
int fnpf_amr::finest_at(double x, double y)
{
    int g=-1, l=0;
    for(int n=0; n<(int)P.size(); ++n)
    {
        reefamr_patch *c = P[n];
        lexer *pp = c->pp;
        if(c->lev<=l)
        continue;
        if(x>=pp->XN[EXT+marge] && x<pp->XN[EXT+c->nx+marge] && (sy==0 || (y>=pp->YN[EY+marge] && y<pp->YN[EY+c->ny+marge])))
        {
            g=n;
            l=c->lev;
        }
    }
    return g;
}

// --------------------------------------------------------------------- grids
int fnpf_amr::gknoz(int g)
{
    return glex(g)->knoz;
}

slice& fnpf_amr::sval(int g, int s, int m)
{
    if(g<0)
    return (m==0) ? *Se0 : *Sf0;

    fnpf_amr_patch *c = FP(g);
    if(s==0) return (m==0) ? *c->erk1 : *c->frk1;
    if(s==1) return (m==0) ? *c->erk2 : *c->frk2;
    return (m==0) ? c->c->eta : c->c->Fifsf;
}

// --------------------------------------------------------------------- setup
void fnpf_amr::ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    if(!keys_ok(p))
    {
        if(p->mpirank==0)
        cout<<"FNPF AMR (G 1): only for A 310 3, A 311 other than 3, A 350 0-3, X 10 0/1, A 380 0 and on a 2D grid Cartesian -- refinement switched off"<<endl;
        maxlev=0;
        return;
    }

    setup(p,pgc);

    lap_mix = (p->A328==1 && p->N10==1);

    mkdir("./REEF3D_FNPF_AMR",0777);

    if(p->mpirank==0)
    {
        if(reefckpt_open(logout,"./REEF3D_FNPF_AMR/REEF3D_FNPF_AMR_log.dat",p->I45))   // I 45: continued
        logout<<"# count \t simtime \t dt \t patches \t cells \t Laplace iterations (last phi) \t residual \t max phi iterations \t max psi iterations \t solves \t iterations \t layout \t regrids skipped"<<endl;
    }

    Se0 = &c0->eta;
    Sf0 = &c0->Fifsf;
    stg = 2;

    for(int it=0; it<maxlev; ++it)
    regrid(p,pgc,true);

    if(body!=nullptr)
    body->amr_grids(p,pgc);

    // load: the patches stay with the rank of the level-0 cells below them
    const double cmax = pgc->globalmax(double(cells_local));
    const double c0n = pgc->globalmax(double(p->knox*p->knoy));

    if(p->mpirank==0)
    {
        cout<<"FNPF AMR: "<<maxlev<<" level(s), "<<patches_total<<" patch(es), "<<cells_total<<" refined columns";
        if(vref==2)
        cout<<", sigma layers doubled on every level (G 6)";
        const char *follows = (par.zones && adaptive) ? ", the zone follows the body, the patches the flags"
                            : (adaptive ? ", the patches follow the flags" : ", the zone follows the body");
        if(regrid_int>0 && sub==0)
        cout<<", regrid every "<<regrid_int<<" steps (G 2)"<<follows;
        if(regrid_int>0 && sub==1)
        cout<<", regrid every "<<MAX(regrid_int>>maxlev,1)<<" level-0 steps (G 2)"<<follows;
        if(sub==1)
        cout<<", subcycled (G 7 1: level l takes 2^l steps per level-0 step)";
        cout<<endl;
        if(p->G7==1 && sub==0)
        cout<<"FNPF AMR: G 7 1 (subcycling) with a body needs the zone around the hull (G 12) -- one time step for all levels"<<endl;
        if(wd)
        cout<<"FNPF AMR: the patches may cover the "<<(dyn ? "moving" : "static")<<" coastline (A 343 "<<p->A343<<", G 30 1)"<<endl;
        else if(p->A343>=1)
        cout<<"FNPF AMR: the patches stay "<<par.dryband<<" level-0 cells off the coastline (A 343 "<<p->A343<<")"<<endl;
        if(brk)
        cout<<"FNPF AMR: wave breaking on every grid (A 350 "<<p->A350<<")"<<endl;
        if(adaptive)
        {
            cout<<"FNPF AMR: flags:";
            if(tol_eta>0.0)
            cout<<" surface jump > "<<tol_eta<<" m (G 20)";
            if(tol_curv>0.0)
            cout<<" second difference > "<<tol_curv<<" m (G 21)";
            if(nshore>0)
            cout<<" shoreline, "<<nshore<<" cell(s) (G 22)";
            if(flagbreak)
            cout<<" breaking (G 23)";
            cout<<endl;
            if(regrid_int==0)
            cout<<"FNPF AMR: G 2 0 -- the flags are taken once, at the start"<<endl;
        }
        if(p->G22>0 && !wd)
        cout<<"FNPF AMR: G 22 needs G 30 1 (patches over the coastline) -- ignored"<<endl;
        if(p->G23==1 && !brk)
        cout<<"FNPF AMR: G 23 needs A 350 1-3 -- ignored"<<endl;
        if(sub==1 && par.zones && par.zr>zr_user)
        cout<<"FNPF AMR: G 7 1 with a body: zone margin "<<par.zr<<" m (12 cells of level 0) instead of G 12 "<<zr_user<<endl;
        cout<<"FNPF AMR: refined columns per rank max "<<(long)cmax<<", mean "<<setprecision(4)<<double(cells_total)/p->mpi_size
            <<" (level-0 columns per rank "<<(long)c0n<<")"<<endl;
    }
}

// --------------------------------------------------------------------- patches
reefamr_patch* fnpf_amr::patch_new()
{
    fnpf_amr_patch *c = new fnpf_amr_patch;
    c->serial = serial_next++;
    return c;
}

// sigma grid and 3D flags of the patch, as driver::makegrid_sigma: a fluid column is fluid
// from the bed node k=0 to the free-surface node k=knoz
void fnpf_amr::build_lexer3D(fnpf_amr_patch &c)
{
    lexer *pp = c.pp;
    const int n3 = pp->imax*pp->jmax*pp->kmax;
    const int n7 = pp->imax*pp->jmax*(pp->kmax+2);

    pp->Iarray(pp->flag4,n3);
    pp->Iarray(pp->flag7,n7);
    for(int q=0; q<n7; ++q)
    pp->flag7[q] = -10;

    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        const bool fl = pp->flagslice4[lij(pp,ii,jj)]>0;

        for(int kk=pp->kmin; kk<pp->kmin+pp->kmax; ++kk)
        pp->flag4[(ii-pp->imin)*pp->jmax*pp->kmax + (jj-pp->jmin)*pp->kmax + kk-pp->kmin] = (fl && kk>=0 && kk<pp->knoz) ? 1 : -10;

        if(fl)
        for(int kk=0; kk<=pp->knoz; ++kk)
        pp->flag7[fidx(pp,ii,jj,kk)] = 1;
    }

    pp->Darray(pp->ZSN,pp->imax*pp->jmax*(pp->kmax+1));
    pp->Darray(pp->ZSP,n3);
    pp->Darray(pp->bed,pp->imax*pp->jmax);

    // the free surface spans the whole patch: no inflow clamp of velcalc_sig, no domain ends
    pp->origin_i = pp->origin_j = 100000000;
}

void fnpf_amr::free_lexer3D(fnpf_amr_patch &c)
{
    lexer *pp = c.pp;
    const int n3 = pp->imax*pp->jmax*pp->kmax;
    const int n7 = pp->imax*pp->jmax*(pp->kmax+2);

    pp->del_Iarray(pp->flag4,n3);
    pp->del_Iarray(pp->flag7,n7);
    pp->del_Darray(pp->ZSN,pp->imax*pp->jmax*(pp->kmax+1));
    pp->del_Darray(pp->ZSP,n3);
    pp->del_Darray(pp->bed,pp->imax*pp->jmax);
    pp->del_Darray(pp->sig,n7);
    pp->del_Darray(pp->sigx,n7);
    pp->del_Darray(pp->sigy,n7);
    pp->del_Darray(pp->sigz,pp->imax*pp->jmax);
    pp->del_Darray(pp->sigxx,n7);
}

// columns of the patch with a solid neighbour in x or y (wall ghost nodes of Fi, as gc_fivec)
void fnpf_amr::build_walls(fnpf_amr_patch &c)
{
    lexer *pp = c.pp;
    c.wall.clear();

    for(int ii=0; ii<pp->knox; ++ii)
    for(int jj=0; jj<pp->knoy; ++jj)
    {
        if(pp->flagslice4[lij(pp,ii,jj)]<0)
        continue;

        if(pp->flagslice4[lij(pp,ii-1,jj)]<0 || pp->flagslice4[lij(pp,ii+1,jj)]<0
        || pp->flagslice4[lij(pp,ii,jj-1)]<0 || pp->flagslice4[lij(pp,ii,jj+1)]<0)
        {
            c.wall.push_back(ii);
            c.wall.push_back(jj);
        }
    }
}

// wall ghost nodes of a Fi-layout array: mirrored, as ghostcell::fivec (no in/outflow on a patch)
void fnpf_amr::walls_fi(fnpf_amr_patch &c, double *f)
{
    lexer *pp = c.pp;
    const int sI = pp->jmax*pp->kmaxF;
    const int sJ = pp->kmaxF;
    const int *fl = pp->flag7;

    for(size_t n=0; n<c.wall.size(); n+=2)
    for(int kk=0; kk<=pp->knoz; ++kk)
    {
        const int q = fidx(pp,c.wall[n],c.wall[n+1],kk);

        if(fl[q-sI]<0) { f[q-sI] = f[q]; f[q-2*sI] = f[q]; f[q-3*sI] = f[q]; }
        if(fl[q+sI]<0) { f[q+sI] = f[q]; f[q+2*sI] = f[q]; f[q+3*sI] = f[q]; }
        if(fl[q-sJ]<0) { f[q-sJ] = f[q]; f[q-2*sJ] = f[q]; f[q-3*sJ] = f[q]; }
        if(fl[q+sJ]<0) { f[q+sJ] = f[q]; f[q+2*sJ] = f[q]; f[q+3*sJ] = f[q]; }
    }
}

void fnpf_amr::walls_sl(fnpf_amr_patch &c, slice &f, int gcv)
{
    comms_off guard(pgc0);
    pgc0->gcsl_start4(c.pp,f,gcv);
}

// FNPF objects of a new patch (comms off); bed and state are set by regrid_state
void fnpf_amr::patch_objects(reefamr_patch *q, ghostcell *pgc)
{
    fnpf_amr_patch *c = FP(q);

    build_bc2D(pgc,*c);
    build_lexer3D(*c);

    lexer *pp = c->pp;

    c->c = new fdm_fnpf(pp);
    if(wd)
    c->pf = new fnpf_fsfbc_wd(pp,c->c,pgc);
    else
    c->pf = new fnpf_fsfbc(pp,c->c,pgc);
    c->psig = new fnpf_sigma(pp,c->c,pgc);
    c->pfu = new fnpf_fsf_update(pp,c->c,pgc);
    c->pbu = new fnpf_bed_update(pp);
    c->plap = new fnpf_laplace_cds2(pp);
    c->plap->assemble_only = true;

    c->erk1 = new slice4(pp); c->erk2 = new slice4(pp);
    c->frk1 = new slice4(pp); c->frk2 = new slice4(pp);
    c->ek = new slice4(pp); c->fk = new slice4(pp);
    c->en = new slice4(pp);

    // the patches are fully wet and away from the coastline (A 343: cell_unfit): no coastline
    // band (the coastline viscosity of A 337 in fnpf_breaking reads c->coastline)
    fdm_fnpf *cc = c->c;
    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        cc->bc(ii,jj) = 0;
        cc->etaloc(ii,jj) = pp->knoz;
        cc->coastline(ii,jj) = 1.0e20;
        pp->wet[lij(pp,ii,jj)] = 1;
    }

    if(dyn)
    {
        c->wfix.assign(pp->imax*pp->jmax,-1);
        pp->wetfix = c->wfix.data();
    }

    build_walls(*c);
}

void fnpf_amr::patch_delete(reefamr_patch *q)
{
    fnpf_amr_patch *c = FP(q);

    for(auto v : c->kvw)
    delete [] v;
    for(auto v : c->kv)
    delete [] v;
    delete c->mg;

    delete c->erk1; delete c->erk2; delete c->frk1; delete c->frk2; delete c->ek; delete c->fk; delete c->en;
    delete c->plap;
    delete c->pbu;
    delete c->pfu;
    delete c->psig;
    delete c->pf;

    // fdm_fnpf has no destructor (on level 0 it lives for the whole run): its 3D arrays are
    // freed here, otherwise every regrid that replaces a patch leaks them
    {
        fdm_fnpf *cc = c->c;
        lexer *pp = c->pp;
        const int n7 = pp->imax*pp->jmax*(pp->kmax+2);

        pp->del_Darray(cc->U,n7);
        pp->del_Darray(cc->V,n7);
        pp->del_Darray(cc->W,n7);
        pp->del_Darray(cc->Fi,n7);
        pp->del_Darray(cc->Uin,n7);

        if(pp->X10>0)
        {
        pp->del_Darray(cc->FBF,n7);
        pp->del_Darray(cc->FBu,n7);
        pp->del_Darray(cc->FBv,n7);
        pp->del_Darray(cc->FBw,n7);
        }
    }
    delete c->c;

    free_lexer3D(*c);
}

void fnpf_amr::tag(int l, vector<unsigned char> &M)
{
    // the boxes (G 10) and the body zone (G 12) are marked by the core; the solution-adaptive
    // flags on the end-of-step state of the level-l grids (level 0 under the patches: the
    // restricted surface and its own breaking detection)
    if(!adaptive)
    return;

    const int lf = l+1;

    auto test = [&](lexer *q, fdm_fnpf *d, int ii, int jj)
    {
        // a wet cell above the thin film (as the dryband, cell_unfit)
        auto ok = [&](int a, int b)
        {
            const int n = lij(q,a,b);
            return q->flagslice4[n]>0 && q->wet[n]==1 && d->WL(a,b)>=10.0*q->A344;
        };

        // breaking (G 23)
        if(flagbreak && q->flagslice4[lij(q,ii,jj)]>0 && d->breaking(ii,jj)>0)
        return true;

        // shoreline (G 22): a fluid cell with a cell of the other wet state within nshore cells
        if(nshore>0 && q->flagslice4[lij(q,ii,jj)]>0)
        {
            const int w0 = q->wet[lij(q,ii,jj)];
            for(int a=-nshore; a<=nshore; ++a)
            for(int b=-nshore*q->j_dir; b<=nshore*q->j_dir; ++b)
            {
                const int n = lij(q,ii+a,jj+b);
                if(q->flagslice4[n]>0 && q->wet[n]!=w0)
                return true;
            }
        }

        if(!ok(ii,jj))
        return false;

        const double ec = d->eta(ii,jj);

        if(tol_eta>0.0)
        {
            const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};
            for(int k=0; k<4; ++k)
            if((dj[k]==0 || q->j_dir==1) && ok(ii+di[k],jj+dj[k]) && fabs(d->eta(ii+di[k],jj+dj[k])-ec)>tol_eta)
            return true;
        }

        if(tol_curv>0.0)
        {
            if(ok(ii+1,jj) && ok(ii-1,jj) && fabs(d->eta(ii+1,jj)-2.0*ec+d->eta(ii-1,jj))>tol_curv)
            return true;
            if(q->j_dir==1 && ok(ii,jj+1) && ok(ii,jj-1) && fabs(d->eta(ii,jj+1)-2.0*ec+d->eta(ii,jj-1))>tol_curv)
            return true;
        }
        return false;
    };

    if(l==0)
    {
        for(int ii=0; ii<NX0; ++ii)
        for(int jj=0; jj<NY0; ++jj)
        if(test(p0,c0,ii,jj))
        tag_cell(lf,ii+O0i,jj+O0j,nbuf,M);
    }

    for(int n : lev[l])
    {
        fnpf_amr_patch *pc = FP(n);
        for(int ii=EXT; ii<EXT+pc->nx; ++ii)
        for(int jj=EY; jj<EY+pc->ny; ++jj)
        if(test(pc->pp,pc->c,ii,jj))
        tag_cell(lf,ii-EXT+pc->I0,jj-EY+pc->J0,nbuf,M);
    }
}

// --------------------------------------------------------------------- interpolation
// interpolation from coarse cell (ic,jc) to the centre of its child (ox,oy), a quarter coarse
// cell away in x and y: weights w[(di+2)*5+(dj+2)] of the coarse cells (ic+di,jc+dj).
// pord 4: bicubic (tensor product, the 4 cells around the child centre in each direction),
// where all 16 cells are fluid; otherwise biquadratic on the 3x3 cells, slopes and
// curvatures switched off next to solids
void fnpf_amr::pweights(lexer *q, int ic, int jc, int ox, int oy, double *w)
{
    // wfill (coastline in the patches): dry cells of the coarse grid as solids
    auto fl = [&](int a, int b) { const int n = lij(q,a,b); return q->flagslice4[n]>0 && (!wfill || q->wet[n]==1); };

    for(int m=0; m<25; ++m)
    w[m] = 0.0;

    // a dry coarse cell: the mean of its wet neighbours (a wet child at the coastline takes the
    // water of the coarse grid), its own value where it has none
    if(wfill && q->flagslice4[lij(q,ic,jc)]>0 && q->wet[lij(q,ic,jc)]==0)
    {
        int nw = 0;
        for(int a=-1; a<=1; ++a)
        for(int b=-sy; b<=sy; ++b)
        if((a!=0 || b!=0) && fl(ic+a,jc+b))
        ++nw;

        if(nw==0)
        {
        w[2*5+2] = 1.0;
        return;
        }

        for(int a=-1; a<=1; ++a)
        for(int b=-sy; b<=sy; ++b)
        if((a!=0 || b!=0) && fl(ic+a,jc+b))
        w[(a+2)*5+(b+2)] = 1.0/double(nw);
        return;
    }

    // 2D grid (oy 0): cubic in x on the 4 cells around the child centre, else quadratic in x
    if(oy==0)
    {
        const int a0 = (ox>0) ? -1 : -2;
        bool all = (pord>=4);
        for(int a=a0; a<a0+4 && all; ++a)
        if(!fl(ic+a,jc))
        all=false;
        if(all)
        {
            const double cp[4] = {-0.0546875, 0.8203125, 0.2734375, -0.0390625};
            for(int a=0; a<4; ++a)
            w[(a0+a+2)*5+2] = (ox>0) ? cp[a] : cp[3-a];
            return;
        }
        const double x = 0.25*ox;
        w[2*5+2] = 1.0;
        if(fl(ic+1,jc) && fl(ic-1,jc))
        {
            w[3*5+2] += 0.5*x + 0.5*x*x;
            w[1*5+2] += -0.5*x + 0.5*x*x;
            w[2*5+2] -= x*x;
        }
        return;
    }

    if(pord>=4)
    {
        bool all=true;
        const int a0 = (ox>0) ? -1 : -2, b0 = (oy>0) ? -1 : -2;
        for(int a=a0; a<a0+4 && all; ++a)
        for(int b=b0; b<b0+4 && all; ++b)
        if(!fl(ic+a,jc+b))
        all=false;

        if(all)
        {
            // Lagrange weights at +1/4 on the nodes -1,0,1,2 (mirrored for -1/4)
            const double cp[4] = {-0.0546875, 0.8203125, 0.2734375, -0.0390625};
            double wx[4], wy[4];
            for(int m=0; m<4; ++m)
            {
                wx[m] = (ox>0) ? cp[m] : cp[3-m];
                wy[m] = (oy>0) ? cp[m] : cp[3-m];
            }
            for(int a=0; a<4; ++a)
            for(int b=0; b<4; ++b)
            w[(a0+a+2)*5+(b0+b+2)] = wx[a]*wy[b];
            return;
        }
    }

    const bool ex = fl(ic+1,jc) && fl(ic-1,jc);
    const bool ey = fl(ic,jc+1) && fl(ic,jc-1);
    const bool exy = ex && ey && fl(ic+1,jc+1) && fl(ic-1,jc+1) && fl(ic+1,jc-1) && fl(ic-1,jc-1);

    const double x = 0.25*ox, y = 0.25*oy;
    auto W = [&](int di, int dj) -> double& { return w[(di+2)*5+(dj+2)]; };

    W(0,0) = 1.0;

    if(ex)
    {
        W(1,0) += 0.5*x + 0.5*x*x;
        W(-1,0) += -0.5*x + 0.5*x*x;
        W(0,0) -= x*x;
    }
    if(ey)
    {
        W(0,1) += 0.5*y + 0.5*y*y;
        W(0,-1) += -0.5*y + 0.5*y*y;
        W(0,0) -= y*y;
    }
    if(exy)
    {
        const double c = 0.25*x*y;
        W(1,1) += c; W(-1,1) -= c; W(1,-1) -= c; W(-1,-1) += c;
    }
}

// slice value at the child centre
double fnpf_amr::pq(slice &f, lexer *q, int ic, int jc, int ox, int oy)
{
    double w[25];
    pweights(q,ic,jc,ox,oy,w);

    double r = 0.0;
    for(int di=-2; di<=2; ++di)
    for(int dj=-2; dj<=2; ++dj)
    {
        const double c = w[(di+2)*5+(dj+2)];
        if(c!=0.0)
        r += c*f(ic+di,jc+dj);
    }
    return r;
}

// the cached weights of pcol (see pwcache): their number, or -1 where not cached
int fnpf_amr::pw_get(int g, lexer *q, int ic, int jc, int ox, int oy, const int *&offs, const double *&ws)
{
    if(wfill || body!=nullptr)
    return -1;

    if((int)pwc.size()<(int)P.size()+1)
    pwc.resize(P.size()+1);
    pwcache &C = pwc[g+1];
    const int stamp = regrids;
    if(C.stamp!=stamp || C.q!=q)
    {
        C.stamp = stamp;
        C.q = q;
        C.slot.assign((size_t)q->imax*q->jmax*8,-1);
        C.num.clear(); C.off.clear(); C.wt.clear();
    }

    // (pord: the preconditioner prolongs biquadratically, the fills bicubically)
    const int qd = ((ox>0) ? 1 : 0) + ((oy>0) ? 2 : 0) + ((pord>=4) ? 4 : 0);
    const long key = (long)lij(q,ic,jc)*8 + qd;
    if(key<0 || key>=(long)C.slot.size())
    return -1;

    if(C.slot[key]==-1)
    {
        double w[25];
        pweights(q,ic,jc,ox,oy,w);
        const int sI = q->jmax*q->kmaxF;
        const int sJ = q->kmaxF;
        const int at = (int)C.off.size();
        int nw = 0;
        for(int di=-2; di<=2; ++di)
        for(int dj=-2; dj<=2; ++dj)
        {
            const double c = w[(di+2)*5+(dj+2)];
            if(c==0.0)
            continue;
            C.off.push_back(di*sI + dj*sJ);
            C.wt.push_back(c);
            ++nw;
        }
        C.slot[key] = (int)C.num.size();
        C.num.push_back(at);
        C.num.push_back(nw);
    }

    const int e = C.slot[key];
    const int at = C.num[e], nw = C.num[e+1];
    offs = C.off.data() + at;
    ws = C.wt.data() + at;
    return nw;
}

// column of a child of coarse cell (ic,jc) of grid g: nodes 0..knf of the fine column;
// horizontally biquadratic on the coarse nodes, vertically linear in sigma (the fine nodes
// subdivide the coarse ones evenly)
void fnpf_amr::pcol(int g, int ic, int jc, int ox, int oy, const double *src, int knf, double *v)
{
    lexer *q = glex(g);
    const int knc = q->knoz;
    const int fz = knf/knc;
    const int sI = q->jmax*q->kmaxF;
    const int sJ = q->kmaxF;

    const double *s[25];
    int nw=0;
    double ww[25];
    const int n0 = fidx(q,ic,jc,0);

    const int *co;
    const double *cw;
    const int nc = pw_get(g,q,ic,jc,ox,oy,co,cw);
    if(nc>=0)
    {
        // (the same order of the weights as below)
        for(int m=0; m<nc; ++m)
        {
            s[m] = src + n0 + co[m];
            ww[m] = cw[m];
        }
        nw = nc;
    }
    else
    {
    double w[25];
    pweights(q,ic,jc,ox,oy,w);
    for(int di=-2; di<=2; ++di)
    for(int dj=-2; dj<=2; ++dj)
    {
        const double c = w[(di+2)*5+(dj+2)];
        if(c==0.0)
        continue;
        s[nw] = src + n0 + di*sI + dj*sJ;
        ww[nw] = c;
        ++nw;
    }
    }

    auto node = [&](int K)
    {
        double r = 0.0;
        for(int m=0; m<nw; ++m)
        r += ww[m]*s[m][K];
        return r;
    };

    if(fz==1)
    {
        for(int kk=0; kk<=knf; ++kk)
        v[kk] = node(kk);
        return;
    }

    double lo = node(0);
    for(int K=0; K<knc; ++K)
    {
        const double hi = node(K+1);
        for(int r=0; r<fz; ++r)
        v[K*fz+r] = lo + (hi-lo)*(double(r)/double(fz));
        lo = hi;
    }
    v[knf] = lo;
}

// --------------------------------------------------------------------- regridding
void fnpf_amr::regrid_prepare(ghostcell *pgc)
{
    Se0 = &c0->eta;
    Sf0 = &c0->Fifsf;

    auto sel = [&](int g, int m) -> slice& { return sval(g,2,m); };
    for(int l=1; l<=maxlev; ++l)
    fill_sl(l,2,7500+l,sel);
}

void fnpf_amr::regrid_static(ghostcell *pgc)
{
}

// state of the new patches, coarse to fine: bed, depth and its derivatives, eta and Fifsf,
// sigma grid, Fi and Fz.  Where an old patch of the same level was (the zone moved with the
// body), its values are taken over; elsewhere they are interpolated from the parent level.
// The bed is always prolonged from the parent (static, the same values as before).
void fnpf_amr::regrid_state(ghostcell *pgc, vector<reefamr_patch*> &oldP)
{
    Se0 = &c0->eta;
    Sf0 = &c0->Fifsf;

    // the patch layout changed: rows, multigrids and the body grids are rebuilt
    {
        bool changed = false;
        for(auto q : P)
        if(q->fresh)
        changed = true;
        for(auto q : oldP)
        if(std::find(P.begin(),P.end(),q)==P.end())
        changed = true;
        if(changed)
        ++layout_id;
    }

    auto sbed = [&](int g, int m) -> slice& { return gfd(g)->bed; };
    auto sfz  = [&](int g, int m) -> slice& { return gfd(g)->Fz; };
    auto sst  = [&](int g, int m) -> slice& { return sval(g,2,m); };
    auto sfi  = [&](int g) -> double* { return gfd(g)->Fi; };
    auto scst = [&](int g, int m) -> slice& { return gfd(g)->coastline; };
    auto svb  = [&](int g, int m) -> slice& { return gfd(g)->vb; };
    auto fresh = [&](reefamr_patch *c) { return c->fresh; };

    for(int l=1; l<=maxlev; ++l)
    {
        // bed (geometry: all cells)
        wfill = false;
        prolong_interior_sl(l,[&](reefamr_patch *c) { return c->fresh; },1,sbed);
        fill_sl(l,1,7510+l,sbed);
        wfill = wd;

        for(int id : lev[l])
        if(P[id]->fresh)
        {
            fnpf_amr_patch *c = FP(id);
            lexer *p = c->pp;
            fdm_fnpf *cc = c->c;

            walls_sl(*c,cc->bed,50);

            comms_off guard(pgc);

            SLICEBASELOOP
            {
            p->bed[IJ] = cc->bed(i,j);
            cc->depth(i,j) = p->wd - cc->bed(i,j);
            }
            pgc->gcsl_start4(p,cc->depth,50);

            for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
            for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
            p->bed[lij(p,ii,jj)] = cc->bed(ii,jj);
        }

        // G 30 1: the coastline distance of the parent grid (geometry: all cells), for the cells of
        // a new patch whose coastline may lie outside the patch (patch_coast_dist)
        if(wd)
        {
            wfill = false;
            prolong_interior_sl(l,fresh,1,scst);
            wfill = wd;
        }

        // surface: prolonged, then the cells of the old patches (on any rank)
        prolong_interior_sl(l,[&](reefamr_patch *c) { return c->fresh; },2,sst);
        prolong_interior_sl(l,[&](reefamr_patch *c) { return c->fresh; },1,sfz);
        old_run(l,3,7590+l,oldP,
                [&](reefamr_patch *o, int io, int jo, double *v)
                {
                    fdm_fnpf *oc = FP(o)->c;
                    v[0] = oc->eta(io,jo);
                    v[1] = oc->Fifsf(io,jo);
                    v[2] = oc->Fz(io,jo);
                },
                [&](reefamr_patch *q, int id, int ii, int jj, const double *v)
                {
                    fdm_fnpf *cc = FP(q)->c;
                    cc->eta(ii,jj) = v[0];
                    cc->Fifsf(ii,jj) = v[1];
                    cc->Fz(ii,jj) = v[2];
                });
        fill_sl(l,2,7520+l,sst);
        fill_sl(l,1,7530+l,sfz);

        // dyn: the wet flags of the parent grids prolonged (all cells), a new cell is wet where this
        // exceeds 1/2 and its water level the criterion plus the rewetting margin
        vector<slice4*> wsl(P.size()+1,nullptr);
        if(dyn)
        {
            auto make = [&](int g)
            {
                lexer *q = glex(g);
                wsl[g+1] = new slice4(q);
                for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
                for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
                (*wsl[g+1])(ii,jj) = double(q->wet[lij(q,ii,jj)]);
            };
            if(l==1)
            make(-1);
            for(int g=0; g<(int)P.size(); ++g)
            if(P[g]->lev==l-1 || (P[g]->lev==l && P[g]->fresh))
            make(g);

            wfill = false;
            prolong_interior_sl(l,fresh,1,[&](int g, int m) -> slice& { return *wsl[g+1]; });
            wfill = wd;
        }

        // breaking (A 350 1): the viscosity of the last stage, which damps the first stage of the
        // next step (the detection comes after the damping): prolonged, then that of the old patches
        if(brk)
        {
            prolong_interior_sl(l,fresh,1,svb);
            old_run(l,1,7650+l,oldP,
                    [&](reefamr_patch *o, int io, int jo, double *v) { v[0] = FP(o)->c->vb(io,jo); },
                    [&](reefamr_patch *q, int id, int ii, int jj, const double *v) { FP(q)->c->vb(ii,jj) = v[0]; });
        }

        for(int id : lev[l])
        if(P[id]->fresh)
        {
            fnpf_amr_patch *c = FP(id);
            lexer *p = c->pp;
            fdm_fnpf *cc = c->c;

            walls_sl(*c,cc->eta,gcval_eta);
            walls_sl(*c,cc->Fifsf,gcval_fifsf);

            comms_off guard(pgc);

            // as fnpf_RK3::inidisc_step1 (it sets every cell wet), before the wet mask below
            c->pf->fsfdisc_ini(p,cc,pgc,cc->eta,cc->Fifsf);

            // coastline in the patch: wet mask and coastline distance from the fine bed, the
            // distance of the parent where the coastline may lie outside the patch
            if(wd)
            {
                vector<double> cpar, esave;
                cpar.reserve((size_t)c->nx*c->ny);
                for(int ii=EXT; ii<EXT+c->nx; ++ii)
                for(int jj=EY; jj<EY+c->ny; ++jj)
                {
                cpar.push_back(cc->coastline(ii,jj));
                esave.push_back(cc->eta(ii,jj));
                }

                patch_wd_ini(*c,pgc);
                patch_coast_dist(*c,cpar);

                // dyn: the wet mask of the water now (the flags of the old patches below), not of
                // the still water level of patch_wd_ini; wet above the criterion plus the rewetting
                // margin A 331 (a dry parent gives the dry level, its interpolation sits on the
                // criterion)
                if(dyn)
                {
                    slice4 &wi = *wsl[id+1];
                    size_t q = 0;
                    for(int ii=EXT; ii<EXT+c->nx; ++ii)
                    for(int jj=EY; jj<EY+c->ny; ++jj, ++q)
                    {
                        cc->eta(ii,jj) = esave[q];
                        const int n = lij(p,ii,jj);
                        if(p->flagslice4[n]>0)
                        p->wet[n] = (wi(ii,jj)>0.5 && esave[q] + cc->depth(ii,jj) > (1.0+p->A331)*cc->wd_criterion) ? 1 : 0;
                    }
                }
            }

            if(brk)
            for(int ii=EXT; ii<EXT+c->nx; ++ii)
            for(int jj=EY; jj<EY+c->ny; ++jj)
            cc->vb(ii,jj) = MAX(cc->vb(ii,jj),0.0);
        }

        for(auto w : wsl)
        delete w;

        // dyn: the flags of the old patches (on any rank), those of the cells around the patches
        // from their source grids, then the state of the wetting and drying from these flags
        if(dyn)
        {
            old_run(l,1,7740+l,oldP,
                    [&](reefamr_patch *o, int io, int jo, double *v) { v[0] = o->pp->wet[lij(o->pp,io,jo)]; },
                    [&](reefamr_patch *q, int id, int ii, int jj, const double *v) { q->pp->wet[lij(q->pp,ii,jj)] = (int)v[0]; });
            fill_wet(l);
        }

        for(int id : lev[l])
        if(P[id]->fresh)
        {
            fnpf_amr_patch *c = FP(id);
            lexer *p = c->pp;
            fdm_fnpf *cc = c->c;

            comms_off guard(pgc);

            if(dyn)
            {
                fnpf_amr_local loc(pgc);
                static_cast<fnpf_fsfbc_wd*>(c->pf)->wd_reset(p,cc,pgc,cc->eta,cc->Fifsf);
            }

            // sigma grid as fnpf_sigma::sigma_ini (the bed is prolonged, not smoothed again)
            FLOOP
            p->sig[FIJK] = p->ZN[KP];

            SLICELOOP4
            {
                k=0;
                p->sig[FIJKm1] = p->ZN[KM1];
                p->sig[FIJKm2] = p->ZN[KM2];
                p->sig[FIJKm3] = p->ZN[KM3];
                k=p->knoz;
                p->sig[FIJKp1] = p->ZN[KP1];
                p->sig[FIJKp2] = p->ZN[KP2];
                p->sig[FIJKp3] = p->ZN[KP3];
            }

            c->pf->fsfdisc(p,cc,pgc,cc->eta,cc->Fifsf);
            c->psig->sigma_update(p,cc,pgc,c->pf,cc->eta);
        }

        // Fi: prolonged, then the columns of the old patches
        prolong_interior_col(l,[&](reefamr_patch *c) { return c->fresh; },sfi,sfi);
        {
        const int knf = klev(l);
        old_run(l,knf+1,7595+l,oldP,
                [&](reefamr_patch *o, int io, int jo, double *v)
                {
                    for(int kk=0; kk<=knf; ++kk)
                    v[kk] = FP(o)->c->Fi[fidx(o->pp,io,jo,kk)];
                },
                [&](reefamr_patch *q, int id, int ii, int jj, const double *v)
                {
                    for(int kk=0; kk<=knf; ++kk)
                    FP(q)->c->Fi[fidx(q->pp,ii,jj,kk)] = v[kk];
                });
        }
        fill_col(l,7540+l,sfi);

        for(int id : lev[l])
        if(P[id]->fresh)
        {
            fnpf_amr_patch *c = FP(id);
            comms_off guard(pgc);
            walls_fi(*c,c->c->Fi);
            c->pfu->fsfbc_sig(c->pp,c->c,pgc,c->c->Fifsf,c->c->Fi);
            c->pbu->bedbc_sig(c->pp,c->c,pgc,c->c->Fi,c->pf);
            c->pfu->velcalc_sig(c->pp,c->c,pgc,c->c->Fi);
        }
    }
}

// level 0 consistent with the patches
void fnpf_amr::regrid_finish(ghostcell *pgc, int old_total)
{
    Se0 = &c0->eta;
    Sf0 = &c0->Fifsf;

    restrict_sl(2,[&](int g, int m) -> slice& { return sval(g,2,m); });
    restrict_col([&](int g) -> double* { return gfd(g)->Fi; });
    restrict_flags(pgc);

    if(patches_total>0 || old_total>0)
    {
        pgc->gcsl_start4(p0,c0->eta,gcval_eta);
        pgc->gcsl_start4(p0,c0->Fifsf,gcval_fifsf);
        pgc->start7V(p0,c0->Fi,c0->bc,250);
    }
}

// --------------------------------------------------------------------- time stepping
// the tendencies of the patches (kinematic and dynamic FSBC), as fnpf_RK3: before the body
// loads, which need them for the psi_0 data
void fnpf_amr::stage_tendency(lexer *p, fdm_fnpf *c, ghostcell *pgc, int s)
{
    if(!active())
    return;

    // G 7 1: level 0 steps alone; at its first stage (all grids at the start of the step) the
    // parent columns of level 1, the body saved before its predicted stages, and the patch
    // tendencies for the loads of all grids
    if(sub_on())
    {
        if(s>0)
        return;
        sub_snapshot(-1);
        fnpf_6DOF *b = dynamic_cast<fnpf_6DOF*>(body);
        if(b!=nullptr)
        b->amr_save();
        if(b==nullptr)
        return;
    }

    double t0 = MPI_Wtime();

    comms_off guard(pgc);

    for(auto q : P)
    {
        fnpf_amr_patch *pc = FP(q);
        lexer *pp = pc->pp;
        fdm_fnpf *cc = pc->c;

        pp->dt = p->dt;
        pp->dt_old = p->dt_old;
        pp->simtime = p->simtime;
        pp->count = p->count;
        pp->amr_follow(*p);

        slice4 &ek = *pc->ek, &fk = *pc->fk;
        slice4 &erk1 = *pc->erk1, &erk2 = *pc->erk2;
        slice &Ein = (s==0) ? static_cast<slice&>(cc->eta) : (s==1) ? static_cast<slice&>(erk1) : static_cast<slice&>(erk2);

        lexer *p = pp;
        fdm_fnpf *c = cc;

        pc->pf->kfsfbc(p,c,pgc);
        SLICELOOP4
        ek(i,j) = c->K(i,j);

        pc->pf->dfsfbc(p,c,pgc,Ein);
        SLICELOOP4
        fk(i,j) = c->K(i,j);
    }

    tm[0] += MPI_Wtime()-t0;
}

// the patches form their stage values (as fnpf_RK3) and the footprint of the body, the fine
// values go to the coarser levels, then the cells around the patches are filled
void fnpf_amr::stage_surface(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &Se, slice &Sf, int s)
{
    Se0 = &Se;
    Sf0 = &Sf;
    stg = s;

    if(!active() || sub_on())
    return;

    double t0 = MPI_Wtime();

    for(int id=0; id<(int)P.size(); ++id)
    {
        fnpf_amr_patch *pc = FP(id);

        {
        comms_off guard(pgc);

        fdm_fnpf *cc = pc->c;
        slice4 &ek = *pc->ek, &fk = *pc->fk;
        slice4 &erk1 = *pc->erk1, &erk2 = *pc->erk2, &frk1 = *pc->frk1, &frk2 = *pc->frk2;

        lexer *p = pc->pp;
        fdm_fnpf *c = cc;

        if(s==0 && brk)
        {
            slice4 &en = *pc->en;
            SLICELOOP4
            en(i,j) = c->eta(i,j);
        }

        if(s==0)
        {
            SLICELOOP4
            erk1(i,j) = c->eta(i,j) + p->dt*ek(i,j);
            SLICELOOP4
            frk1(i,j) = c->Fifsf(i,j) + p->dt*fk(i,j);
        }
        if(s==1)
        {
            SLICELOOP4
            erk2(i,j) = 0.75*c->eta(i,j) + 0.25*erk1(i,j) + 0.25*p->dt*ek(i,j);
            SLICELOOP4
            frk2(i,j) = 0.75*c->Fifsf(i,j) + 0.25*frk1(i,j) + 0.25*p->dt*fk(i,j);
        }
        if(s==2)
        {
            SLICELOOP4
            c->eta(i,j) = (1.0/3.0)*c->eta(i,j) + (2.0/3.0)*erk2(i,j) + (2.0/3.0)*p->dt*ek(i,j);
            SLICELOOP4
            c->Fifsf(i,j) = (1.0/3.0)*c->Fifsf(i,j) + (2.0/3.0)*frk2(i,j) + (2.0/3.0)*p->dt*fk(i,j);
        }

        }

        if(body!=nullptr && !stage_ring())
        body->amr_surface(pc->pp,pgc,id,sval(id,s,0),sval(id,s,1));
    }

    auto sel = [&](int g, int m) -> slice& { return sval(g,s,m); };

    // breaking (implicit damping) and the coastline in the patches read the cells around a patch:
    // these first take the stage values of the coarser grid (damped already, level by level), then
    // the damping and the coastline of the level, as fnpf_RK3
    if(stage_ring())
    for(int l=1; l<=maxlev; ++l)
    {
        fill_sl(l,2,7640+l,sel);

        if(dyn)
        fill_wet(l);

        for(int id : lev[l])
        {
            walls_sl(*FP(id),sel(id,0),gcval_eta);
            walls_sl(*FP(id),sel(id,1),gcval_fifsf);
        }

        stage_damp_coast(l,pgc,s);
    }

    restrict_sl(2,sel);
    restrict_flags(pgc);

    pgc->gcsl_start4(p,Se,gcval_eta);
    pgc->gcsl_start4(p,Sf,gcval_fifsf);

    for(int l=1; l<=maxlev; ++l)
    {
        fill_sl(l,2,7500+l,sel);

        if(dyn)
        for(int id : lev[l])
        ring_dry(*FP(id),sel(id,0));

        for(int id : lev[l])
        {
            walls_sl(*FP(id),sel(id,0),gcval_eta);
            walls_sl(*FP(id),sel(id,1),gcval_fifsf);
        }
    }

    tm[0] += MPI_Wtime()-t0;
}

void fnpf_amr::step_end(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    // G 7 1: the finer levels after level 0
    if(sub_on())
    sub_step(p,c,pgc);

    {
    comms_off guard(pgc);
    for(auto q : P)
    {
        fnpf_amr_patch *pc = FP(q);
        pc->pbu->bedbc_sig(pc->pp,pc->c,pgc,pc->c->Fi,pc->pf);
        pc->pfu->velcalc_sig(pc->pp,pc->c,pgc,pc->c->Fi);

        if(wd)
        {
        fnpf_amr_local loc(pgc);
        pc->pf->coastline_vel(pc->pp,pc->c,pgc,pc->c->U);
        pc->pf->coastline_vel(pc->pp,pc->c,pgc,pc->c->V);
        pc->pf->coastline_vel(pc->pp,pc->c,pgc,pc->c->W);
        }
    }
    }

    // A 343: the patches are fully wet; a patch that the run-up or the drawdown reaches (a static
    // box, or between two regrids) has no wetting and drying of its own
    if(p->A343>=1 && !wd && !thin_warned)
    {
        int thin = 0;
        for(auto q : P)
        {
            lexer *pp = q->pp;
            fdm_fnpf *cc = FP(q)->c;
            for(int ii=EXT; ii<EXT+q->nx; ++ii)
            for(int jj=EY; jj<EY+q->ny; ++jj)
            if(pp->flagslice4[lij(pp,ii,jj)]>0 && cc->eta(ii,jj)+pp->wd-cc->bed(ii,jj) < 10.0*p->A344)
            thin = 1;
        }
        thin = pgc->globalimax(thin);

        if(thin>0)
        {
            thin_warned = true;
            if(p->mpirank==0)
            cout<<"FNPF AMR WARNING - a patch has a dry or thin-film cell at t = "<<p->simtime
                <<" (the patches have no wetting and drying: keep the refinement away from the coast)"<<endl;
        }
    }

    // the zone around the body moves: new patches every G 2 steps (G 7 1: G 2 counts steps of the
    // finest level)
    const int rint = (sub==1) ? MAX(regrid_int>>maxlev,1) : regrid_int;
    // (replace_due: G 41 0 with a static refinement, the placement once more with the measured model)
    if((regrid_int>0 && p->count%rint==0) || replace_due(p))
    {
        double t0 = MPI_Wtime();
        Se0 = &c0->eta;
        Sf0 = &c0->Fifsf;
        stg = 2;
        regrid(p,pgc,false);
        tm[4] += MPI_Wtime()-t0;
    }

    // load statistics of the ranks (G 43)
    load_step(p);
}

// time step of the finest grid, same CFL as fnpf_timestep
void fnpf_amr::timestep(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    if(p->N48==0)
    return;

    // G 7 1: a patch cell of level l limits the level-0 step as a cell 2^l times larger
    auto sc = [&](reefamr_patch *q) { return (sub==1) ? double(1<<q->lev) : 1.0; };

    double r=1.0;
    for(auto q : P)
    r = MIN(r, sc(q)*q->pp->DXM/p->DXM);
    r = pgc->globalmin(r);

    if(p->count==0)
    {
        p->dt *= r;
        p->dt_old = p->dt;
        return;
    }

    double depthmax=0.0;
    SLICELOOP4
    depthmax = MAX(depthmax,c->depth(i,j));
    depthmax = pgc->globalmax(depthmax);

    double cu = 1.0e10;
    for(auto q : P)
    {
        lexer *pp = q->pp;
        for(int ii=EXT; ii<EXT+q->nx; ++ii)
        for(int jj=EY; jj<EY+q->ny; ++jj)
        if(pp->flagslice4[lij(pp,ii,jj)]>0)
        {
            double dx = sc(q)*((pp->j_dir==1) ? MIN(pp->DXN[ii+marge],pp->DYN[jj+marge]) : pp->DXN[ii+marge]);
            cu = MIN(cu, 1.0/((fabs(MAX(p->umax, sqrt(9.81*depthmax)))/dx)));
            cu = MIN(cu, 1.0/((fabs(MAX(p->vmax, sqrt(9.81*depthmax)))/dx)));
        }
    }

    double dtp = pgc->globalmin(p->N47*cu);
    p->dt = MIN(p->dt,dtp);
}

// --------------------------------------------------------------------- wetting and drying, breaking
// A 343 1-3: a level-0 cell that no patch may cover (the band of q.dryband cells around it is kept
// free as well): dry, thin film (WL < 10 A 344, the cells the coastline viscosity of A 337 and the
// wet-dry clamp act on) or inside the coastline damping band
bool fnpf_amr::cell_unfit(int ii, int jj)
{
    if(p0->A343==0 || wd)
    return false;

    const int n = lij(p0,ii,jj);

    if(p0->flagslice4[n]<=0)
    return false;

    if(p0->wet[n]==0 || c0->WL(ii,jj)<10.0*p0->A344)
    return true;

    if(dist_coast>0.0 && c0->coastline(ii,jj)>=0.0 && c0->coastline(ii,jj)<dist_coast)
    return true;

    return false;
}

// A 350 1 and A 343 1 with G 30 1: the implicit damping of the stage values (alpha 1, 1/4, 2/3), then
// the clamp and the coastline relaxation, in the order of fnpf_RK3; the cells around the patch hold
// the stage values of the coarser grid (stage_ring)
void fnpf_amr::stage_damp_coast(int l, ghostcell *pgc, int s)
{
    if(brk)
    {
        const double al = (s==0) ? 1.0 : (s==1) ? 0.25 : 2.0/3.0;
        level_damping(l,pgc,s,0,gcval_eta,al);
        level_damping(l,pgc,s,1,gcval_fifsf,al);
    }

    if(wd)
    for(int id : lev[l])
    {
        comms_off guard(pgc);
        patch_coast(*FP(id),pgc,sval(id,s,0),sval(id,s,1));
    }

    if(body!=nullptr)
    for(int id : lev[l])
    body->amr_surface(FP(id)->pp,pgc,id,sval(id,s,0),sval(id,s,1));
}

// A 350 1: the implicit damping of the stage value m (0 eta, 1 Fifsf) of stage s with the breaking
// viscosity on the patches of level l, as fnpf_fsfbc_wd::damping on level 0: one system over all
// patches of the level (on any rank), the cells around them from a patch of the same level are
// unknowns of the system (their values through the fills of the level), those from the coarser
// grid hold its damped stage value (stage_ring); solids, the edge of the domain and dry neighbours
// (wd) take the value of the cell at the start of the solve, dry cells stay unchanged (as
// fnpf_fsfbc_wd::damping).  Jacobi-preconditioned
// BiCGStab with the reductions over all ranks; a patch cut into pieces (ranks, G 40 1, several
// patches) is solved as the whole patch.
void fnpf_amr::level_damping(int l, ghostcell *pgc, int s, int m, int gcval, double alpha)
{
    bool on = false;
    double alpha_eta = 1.0;

    if((gcval==60 || gcval==160) && (p0->A357==1 || p0->A357==2))
    on = true;

    if((gcval==55 || gcval==155) && (p0->A357==1 || p0->A357==3))
    {
    on = true;
    alpha_eta = p0->A366;
    }

    if(!on)
    return;

    double vbmax = 0.0;
    for(int id : lev[l])
    {
        fnpf_amr_patch *pc = FP(id);
        for(int ii=EXT; ii<EXT+pc->nx; ++ii)
        for(int jj=EY; jj<EY+pc->ny; ++jj)
        vbmax = MAX(vbmax,fabs(pc->c->vb(ii,jj)));
    }
    vbmax = pgc->globalmax(vbmax);

    if(vbmax<=0.0)
    return;

    // the system of every patch: interior cells q = (i-EXT)*ny + (j-EXT); a neighbour in the box
    // (nb), or a cell around the patch that a patch of the same level fills (ri, rj: its value in
    // the work slices Y, Z)
    struct sysp
    {
        int nx, ny;
        vector<double> ap, b, x, r, r0, pv, v, sv, tv;
        vector<double> ac[4];
        vector<int> nb[4], ri[4], rj[4];
        slice4 *Y = nullptr, *Z = nullptr;
    };
    vector<sysp> S(P.size());
    vector<slice4*> W(P.size(),nullptr);

    for(int id : lev[l])
    {
        fnpf_amr_patch *pc = FP(id);
        S[id].Y = new slice4(pc->pp);
        S[id].Z = new slice4(pc->pp);
        W[id] = S[id].Y;
    }

    // which cells around a patch come from a patch of the same level (1) or the coarser grid (0)
    for(int id : lev[l])
    {
        lexer *pp = FP(id)->pp;
        slice4 &Y = *S[id].Y;
        for(int ii=0; ii<pp->knox; ++ii)
        for(int jj=0; jj<pp->knoy; ++jj)
        Y(ii,jj) = 0.0;
    }
    fill_run(l,1,7660+l,
             [&](const reefamr_fill &f, double *v) { v[0] = (f.kind==0 && f.g>=0 && P[f.g]->lev==l) ? 1.0 : 0.0; },
             [&](reefamr_patch *c, int id, const reefamr_fill &f, const double *w) { (*S[id].Y)(f.di,f.dj) = w[0]; });

    for(int id : lev[l])
    {
        fnpf_amr_patch *pc = FP(id);
        lexer *p = pc->pp;
        fdm_fnpf *c = pc->c;
        slice &f = sval(id,s,m);
        slice4 &sib = *S[id].Y;
        sysp &Q = S[id];
        comms_off guard(pgc);

        const int i0 = EXT, i1 = EXT+pc->nx, j0 = EY, j1 = EY+pc->ny;
        const int nj = pc->ny;
        const int N = pc->nx*pc->ny;
        Q.nx = pc->nx;
        Q.ny = pc->ny;
        Q.ap.assign(N,1.0);
        Q.b.assign(N,0.0);
        Q.x.assign(N,0.0);
        for(int d=0; d<4; ++d)
        {
            Q.ac[d].assign(N,0.0);
            Q.nb[d].assign(N,-1);
            Q.ri[d].assign(N,-1);
            Q.rj[d].assign(N,-1);
        }

        const double idt = 1.0/(alpha*p->dt);

        for(i=i0; i<i1; ++i)
        for(j=j0; j<j1; ++j)
        {
            const int q = (i-i0)*nj + (j-j0);
            Q.x[q] = f(i,j);

            if(p->flagslice4[lij(p,i,j)]<0 || (wd && p->wet[lij(p,i,j)]==0))
            {
            Q.ap[q] = 1.0;
            Q.b[q] = f(i,j);
            continue;
            }

            const double visc = alpha_eta*c->vb(i,j);
            const double cf[4] = {visc/(p->DXP[IM1]*p->DXN[IP])*p->x_dir,
                                  visc/(p->DXP[IP]*p->DXN[IP])*p->x_dir,
                                  visc/(p->DYP[JM1]*p->DYN[JP])*p->y_dir,
                                  visc/(p->DYP[JP]*p->DYN[JP])*p->y_dir};
            const int di[4] = {-1,1,0,0}, dj[4] = {0,0,-1,1};

            Q.ap[q] = cf[0]+cf[1]+cf[2]+cf[3] + idt;
            Q.b[q] = f(i,j)*idt;

            for(int d=0; d<4; ++d)
            {
                if(cf[d]==0.0)
                continue;

                const int a = i+di[d], e = j+dj[d];
                const bool inbox = (a>=i0 && a<i1 && e>=j0 && e<j1);

                // solid, the edge of the domain and a dry neighbour (wd): the value of the cell at the
                // start of the solve, as fnpf_fsfbc_wd::damping on level 0
                if(p->flagslice4[lij(p,a,e)]<0 || (wd && p->wet[lij(p,a,e)]==0))
                Q.b[q] += cf[d]*f(i,j);

                else if(inbox)
                {
                Q.ac[d][q] = -cf[d];
                Q.nb[d][q] = (a-i0)*nj + (e-j0);
                }

                // around the patch, from a patch of the same level: an unknown of the system
                else if(sib(a,e)>0.5)
                {
                Q.ac[d][q] = -cf[d];
                Q.ri[d][q] = a;
                Q.rj[d][q] = e;
                }

                // around the patch, from the coarser grid: its damped stage value
                else
                Q.b[q] += cf[d]*f(a,e);
            }
        }
    }

    // y = A u: u in the interior of the work slices U (W -> U), the cells around the patches from
    // the patches of the same level (0 from the coarser grid: its part is in b)
    auto apply = [&](bool z, const vector<vector<double>*> &u, vector<vector<double>*> &y)
    {
        for(int id : lev[l])
        {
            slice4 &U = z ? *S[id].Z : *S[id].Y;
            W[id] = &U;
            const vector<double> &uu = *u[id];
            const int nj = S[id].ny;
            for(int ii=0; ii<S[id].nx; ++ii)
            for(int jj=0; jj<nj; ++jj)
            U(ii+EXT,jj+EY) = uu[ii*nj+jj];
        }

        fill_run(l,1,7670+l,
                 [&](const reefamr_fill &f, double *v) { v[0] = (f.kind==0 && f.g>=0 && P[f.g]->lev==l) ? (*W[f.g])(f.si,f.sj) : 0.0; },
                 [&](reefamr_patch *c, int id, const reefamr_fill &f, const double *w) { (*W[id])(f.di,f.dj) = w[0]; });

        for(int id : lev[l])
        {
            sysp &Q = S[id];
            slice4 &U = *W[id];
            const vector<double> &uu = *u[id];
            vector<double> &yy = *y[id];
            const int N = Q.nx*Q.ny;
            for(int q=0; q<N; ++q)
            {
                double t = Q.ap[q]*uu[q];
                for(int d=0; d<4; ++d)
                {
                    if(Q.nb[d][q]>=0)
                    t += Q.ac[d][q]*uu[Q.nb[d][q]];
                    else if(Q.ri[d][q]>=0)
                    t += Q.ac[d][q]*U(Q.ri[d][q],Q.rj[d][q]);
                }
                yy[q] = t;
            }
        }
    };

    auto gdot = [&](vector<double> sysp::*a, vector<double> sysp::*b)
    {
        double t = 0.0;
        for(int id : lev[l])
        {
            const vector<double> &u = S[id].*a, &w = S[id].*b;
            for(size_t q=0; q<u.size(); ++q)
            t += u[q]*w[q];
        }
        return pgc->globalsum(t);
    };

    vector<vector<double>> yv(P.size()), zv(P.size());
    vector<vector<double>*> pu(P.size(),nullptr), py(P.size(),nullptr);

    for(int id : lev[l])
    {
        sysp &Q = S[id];
        const size_t N = Q.ap.size();
        Q.r.assign(N,0.0); Q.pv.assign(N,0.0); Q.v.assign(N,0.0); Q.sv.assign(N,0.0); Q.tv.assign(N,0.0);
        yv[id].assign(N,0.0);
        zv[id].assign(N,0.0);
    }

    // r = b - A x
    for(int id : lev[l])
    {
        pu[id] = &S[id].x;
        py[id] = &S[id].r;
    }
    apply(false,pu,py);
    for(int id : lev[l])
    {
        sysp &Q = S[id];
        for(size_t q=0; q<Q.r.size(); ++q)
        Q.r[q] = Q.b[q]-Q.r[q];
        Q.r0 = Q.r;
    }

    const double bn = sqrt(gdot(&sysp::b,&sysp::b));
    const double tol = 1.0e-13*bn + 1.0e-300;
    double rho=1.0, al=1.0, om=1.0;

    for(int it=0; it<500 && sqrt(gdot(&sysp::r,&sysp::r))>tol; ++it)
    {
        const double rho1 = gdot(&sysp::r0,&sysp::r);
        if(fabs(rho1)<1.0e-300)
        break;

        const double beta = (rho1/rho)*(al/om);
        for(int id : lev[l])
        {
            sysp &Q = S[id];
            for(size_t q=0; q<Q.r.size(); ++q)
            {
            Q.pv[q] = Q.r[q] + beta*(Q.pv[q] - om*Q.v[q]);
            yv[id][q] = Q.pv[q]/Q.ap[q];
            }
            pu[id] = &yv[id];
            py[id] = &Q.v;
        }
        apply(false,pu,py);

        al = rho1/gdot(&sysp::r0,&sysp::v);
        for(int id : lev[l])
        {
            sysp &Q = S[id];
            for(size_t q=0; q<Q.r.size(); ++q)
            {
            Q.sv[q] = Q.r[q] - al*Q.v[q];
            zv[id][q] = Q.sv[q]/Q.ap[q];
            }
            pu[id] = &zv[id];
            py[id] = &Q.tv;
        }
        apply(true,pu,py);

        const double tt = gdot(&sysp::tv,&sysp::tv);
        om = (tt>0.0) ? gdot(&sysp::tv,&sysp::sv)/tt : 0.0;

        for(int id : lev[l])
        {
            sysp &Q = S[id];
            for(size_t q=0; q<Q.r.size(); ++q)
            {
            Q.x[q] += al*yv[id][q] + om*zv[id][q];
            Q.r[q] = Q.sv[q] - om*Q.tv[q];
            }
        }

        rho = rho1;

        if(om==0.0)
        break;
    }

    for(int id : lev[l])
    {
        fnpf_amr_patch *pc = FP(id);
        slice &f = sval(id,s,m);
        const int nj = pc->ny;
        for(int ii=EXT; ii<EXT+pc->nx; ++ii)
        for(int jj=EY; jj<EY+pc->ny; ++jj)
        f(ii,jj) = S[id].x[(ii-EXT)*nj + (jj-EY)];

        delete S[id].Y;
        delete S[id].Z;
    }
}

// A 350 1: breaking detection on the stage values of the patches (all levels: l -1, or level l),
// as fnpf_RK3 does on level 0; the cells around the patches then take the breaking viscosity of
// their source grid: their surface is the coarse one (interpolated), their own detection would
// mark breaking there and spread it into the patch (NHFLOW AMR step 7)
void fnpf_amr::patch_breaking(int l, ghostcell *pgc, int s, double alpha)
{
    const int la = (l<0) ? 1 : l, lb = (l<0) ? maxlev : l;

    // the surface slopes of the cells around the patches from their source grid: the detection
    // spreads a breaking flag by two cells, so it reads the cells around a patch, whose own slopes
    // there reach beyond the cells the patch computes
    for(int L=la; L<=lb; ++L)
    fill_sl(L,2,7680+L,[&](int g, int m) -> slice& { return (m==0) ? static_cast<slice&>(gfd(g)->Ex) : static_cast<slice&>(gfd(g)->Ey); });

    {
    comms_off guard(pgc);
    fnpf_amr_local loc(pgc);

    for(int id=0; id<(int)P.size(); ++id)
    {
        if(l>=0 && P[id]->lev!=l)
        continue;

        fnpf_amr_patch *pc = FP(id);
        lexer *pp = pc->pp;

        // the breaking count is printed by level 0
        const int p12 = pp->P12;
        pp->P12 = 1<<30;
        pc->pf->breaking(pp,pc->c,pgc,sval(id,s,0),*pc->en,sval(id,s,1),alpha);
        pp->P12 = p12;
    }
    }

    for(int L=la; L<=lb; ++L)
    fill_run(L,1,7570+L,
             [&](const reefamr_fill &f, double *v)
             {
                 v[0] = (f.kind==0 || f.kind==1) ? gfd(f.g)->vb(f.si,f.sj) : 0.0;
             },
             [&](reefamr_patch *c, int id, const reefamr_fill &f, const double *w)
             {
                 FP(id)->c->vb(f.di,f.dj) = w[0];
             });
}

// A 343 1 with G 30 1: the wet mask, the clamp state and the coastline distance of a new patch from
// its fine bed (fnpf_fsfbc_wd::wetdry with count 0, as at the start of the run on level 0)
void fnpf_amr::patch_wd_ini(fnpf_amr_patch &pc, ghostcell *pgc)
{
    lexer *pp = pc.pp;
    fdm_fnpf *cc = pc.c;
    fnpf_amr_local loc(pgc);

    // isolated wet pockets (A 338) are found on level 0 only: a patch holds a part of the water
    // body, whose connection may run outside it
    const int cnt = pp->count, a338 = pp->A338;
    pp->count = 0;
    pp->A338 = 0;
    pc.pf->wetdry(pp,cc,pgc,cc->eta,cc->Fifsf);
    pp->count = cnt;
    pp->A338 = a338;
}

// A 343 2, 3 with G 30 1: the wet flags of the cells around the level-l patches from their source
// cells (a patch of the same level, or the parent cell), kept through the wetting and drying of the
// patch (lexer::wetfix)
void fnpf_amr::fill_wet(int l)
{
    fill_run(l,1,7690+l,
             [&](const reefamr_fill &f, double *v)
             {
                 lexer *q = glex(f.g);
                 v[0] = (f.kind==0 || f.kind==1) ? double(q->wet[lij(q,f.si,f.sj)]) : 1.0;
             },
             [&](reefamr_patch *c, int id, const reefamr_fill &f, const double *w)
             {
                 lexer *pp = c->pp;
                 const int n = lij(pp,f.di,f.dj);
                 pp->wet[n] = (int)w[0];
                 FP(id)->wfix[n] = (int)w[0];
             });
}

// A 343 2, 3 with G 30 1: a dry cell around a patch holds the water level of a dry cell (the
// criterion A 344 on the bed of the patch), not the interpolated surface of the wet coarse cells
void fnpf_amr::ring_dry(fnpf_amr_patch &pc, slice &Se)
{
    lexer *pp = pc.pp;
    fdm_fnpf *cc = pc.c;
    for(int ii=0; ii<pp->knox; ++ii)
    for(int jj=0; jj<pp->knoy; ++jj)
    {
        const int n = lij(pp,ii,jj);
        if(pc.wfix[n]==0 && pp->flagslice4[n]>0)
        Se(ii,jj) = cc->wd_criterion - cc->depth(ii,jj);
    }
}

// A 343 2, 3 with G 30 1: the wet flags of the covered cells of the coarser grids after a
// restriction: wet only if all four children are wet (as nhflow_amr::restrict_flags: a coarse cell
// over the fine shoreline would take the mean of its wet children and rewet its coarse neighbours
// with that false slope); they keep these flags through the wetting and drying of their grid
// (lexer::wetfix), the other cells of the coarser grids follow their own wetting and drying
void fnpf_amr::restrict_flags(ghostcell *pgc)
{
    if(!dyn)
    return;

    // the covered cells of the window's coarser grids are set below; the others are free again
    // (the cells around a patch keep their source flag)
    if(wlo==0)
    for(int &v : wfix0)
    v = -1;

    for(auto q : P)
    if(q->lev>=MAX(wlo,1) && q->lev<wtop())
    {
        fnpf_amr_patch *c = FP(q);
        lexer *pp = c->pp;
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EY; jj<EY+c->ny; ++jj)
        c->wfix[lij(pp,ii,jj)] = -1;
    }

    for(int l=wtop(); l>=wlo+1; --l)
    block_up(l,1,7730+l,
             [&](reefamr_patch *c, int id, int k, double *v)
             {
                 lexer *pp = c->pp;
                 const int nby = c->ny>>sy;
                 const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                 int wt = 1;
                 for(int a=0; a<2; ++a)
                 for(int b=0; b<=sy; ++b)
                 if(pp->wet[lij(pp,i0+a,j0+b)]!=1)
                 wt = 0;
                 v[0] = wt;
             },
             [&](const reefamr_block &B, int key, const double *v)
             {
                 lexer *q = glex(B.g);
                 const int n = lij(q,B.ic,B.jc);
                 q->wet[n] = (int)v[0];
                 if(B.g<0)
                 wfix0[n] = (int)v[0];
                 else
                 FP(B.g)->wfix[n] = (int)v[0];
             });

    if(wlo==0)
    {
        pgc->gcsl_start4Vint(p0,p0->wet,50);
        pgc->gcsl_start4Vint(p0,p0->wetfix,50);
    }
}

// the coastline distance of a new patch: the patch sees only the coastline inside it and the EXT
// cells around it (wet mask from the fine bed), so where a cell is farther from the coastline than
// from the edge of these cells inside the domain, the coastline may lie outside: there the smaller
// of both distances, the parent's interpolated (cpar, the interior cells in i-major order); the
// sign (wet or dry) is that of the fine bed.  The coastline bands (a few cells) are exact.
void fnpf_amr::patch_coast_dist(fnpf_amr_patch &pc, const vector<double> &cpar)
{
    lexer *pp = pc.pp;
    fdm_fnpf *cc = pc.c;
    const int l = pc.lev;
    const int nxl = GNX<<l, nyl = GNY<<ysh(l);
    const double dx = pp->DXM;
    const double ext = EXT*dx;

    size_t q = 0;
    for(int ii=EXT; ii<EXT+pc.nx; ++ii)
    for(int jj=EY; jj<EY+pc.ny; ++jj, ++q)
    {
        const int I = ii-EXT+pc.I0, J = jj-EY+pc.J0;
        double e = 1.0e20;
        if(pc.I0>0)
        e = MIN(e,(I-pc.I0+0.5)*dx+ext);
        if(pc.I1<nxl-1)
        e = MIN(e,(pc.I1-I+0.5)*dx+ext);
        if(pp->j_dir==1 && pc.J0>0)
        e = MIN(e,(J-pc.J0+0.5)*dx+ext);
        if(pp->j_dir==1 && pc.J1<nyl-1)
        e = MIN(e,(pc.J1-J+0.5)*dx+ext);

        const double d = cc->coastline(ii,jj);
        if(fabs(d)>e)
        cc->coastline(ii,jj) = (d>=0.0 ? 1.0 : -1.0)*MIN(fabs(d),fabs(cpar[q]));
    }
}

// A 343 1 with G 30 1: the clamp of the thin cells and the coastline relaxation of the stage values
// on a patch, in the order of fnpf_RK3 (after the damping)
void fnpf_amr::patch_coast(fnpf_amr_patch &pc, ghostcell *pgc, slice &Se, slice &Sf)
{
    lexer *pp = pc.pp;
    fdm_fnpf *cc = pc.c;
    fnpf_amr_local loc(pgc);

    pc.pf->wetdry(pp,cc,pgc,Se,Sf);
    pc.pf->coastline_eta(pp,cc,pgc,Se);
    pc.pf->coastline_fi(pp,cc,pgc,Sf);
}
