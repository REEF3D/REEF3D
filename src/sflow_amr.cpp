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

#include"sflow_amr.h"
#include"sflow_boussinesq.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"lagoon_output.h"
#include"slice.h"
#include"sflow_HLL.h"
#include"sflow_signal_speed.h"
#include"sflow_reconstruct_hires.h"
#include"sflow_reconstruct_weno.h"
#include"sflow_diffusion_void.h"
#include"sflow_hydrostatic.h"
#include"sflow_eta.h"
#include"sflow_forcing.h"
#include"sflow_momentum_RK3.h"
#include"sflow_pjm_lin.h"
#include"sflow_pjm_quad.h"
#include"sflow_ediff.h"
#include"sflow_amr_ship.h"
#include"sflow_amr_linesolve.h"
#include"6DOF_sflow.h"
#include"reefmg_core.h"
#include"reefmg2D.h"
#include"vec2D.h"
#include"ioflow_void.h"
#include<cmath>
#include<mpi.h>
#include<iomanip>
#include<cstdio>
#include<sys/stat.h>
#include<sys/types.h>

namespace
{
// floor(a/2^l) for negative a as well
inline int fsh(int a, int l) { return (a>=0) ? (a>>l) : -(((-a)-1)>>l)-1; }

inline double mmod(double a, double b)
{
    if(a*b<=0.0)
    return 0.0;
    return fabs(a)<fabs(b)?a:b;
}

// index of cell (i,j) in the lexer 2D arrays (wet, flagslice4)
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

// switches the MPI exchange of the ghostcell class off while patch kernels run
typedef reefamr_comms_off comms_off;
}

sflow_amr::sflow_amr(lexer *p, fdm2D *b, ghostcell *pgc, patchBC_interface *ppBC, sixdof *pp6dof,
                     sflow_HLL *pphll, sflow_momentum_RK3 *ppmom) : reefamr(p,pgc), eps(1.0e-6)
{
    b0 = b;
    pBC = ppBC;
    p6dof = pp6dof;
    phll0 = pphll;
    pmom0 = ppmom;

    reefamr_param q;
    q.name = "SFLOW AMR";
    q.maxlev = p->G1;
    q.regrid = p->G2;
    q.nbuf = MAX(p->G3,0);
    tol_eta = p->G20;
    shore = p->G22;
    q.tile = MAX(p->G4,4);
    q.tile += q.tile%2;
    q.nest = 3;
    q.keep = MIN(MAX(p->G5,0),250);

    // cells computed beyond the patch box: the wet-dry step of a stage needs the new state
    // two cells further out, the discharge limiter (B 60) three
    // (non-hydrostatic without the shallow-water switch A 221: the deep flag needs the
    // wet state three cells further out; the dispersion correction of A 220 3 two more:
    // without them a patch cut at a partition edge sees a 4e-4 relative difference to the serial run)
    // (Boussinesq A 220 4: the dispersive operators reach three cells, the flux M five through the
    // reconstruction; u_a is solved on the leaf cells of all levels together, bous_solve)
    q.ext = (p->B60>=1) ? 3 : 2;
    if(p->A220>=1 && p->A220<=3 && p->A221==0)
    q.ext = 4;
    if(p->A220==3 && p->A224>1.0)
    q.ext += 2;
    if(p->A220==4)
    q.ext = 6;

    nh = (p->A220>=1 && p->A220<=3) ? 1 : 0;
    bous = (p->A220==4) ? 1 : 0;
    if(bous==1)
    pmom0->bous_amr = this;
    NV = bous==1 ? 14 : 10;

    // the patch arrays reach EXT + margin fine cells beyond the patch: the coarser level has to
    // cover them
    if(q.ext>4)
    q.nest = MAX(q.nest,(q.ext+p->margin+1)/2);
    nh_it_total = nh_solves = 0;
    nh_it_last = 0;
    sub = 0;
    hstage = 0;
    tint = -1.0;
    nh_rebuild0 = true;

    // moving body (X 10 2: direct forcing, X 10 3: pressure of the ship)
    shipmode = 0;
    ship6 = nullptr;
    if(p->X10==2 || p->X10==3)
    {
        ship6 = dynamic_cast<sixdof_sflow*>(pp6dof);
        if(ship6!=nullptr)
        shipmode = p->X10;
    }

    // boxes only: static patches
    if(tol_eta<=0.0 && shore==0 && (p->G12==0 || shipmode==0))
    q.regrid = 0;

    // refinement boxes (G 10), no refinement in the G 11 boxes and next to in- and outflow
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
    q.ioband = 4;

    // refinement around the moving body (G 12 margin, G 13 wake wedge)
    q.zones = (shipmode>0 && p->G12>0);
    q.zr = p->G12_r;
    q.zL = p->G13_L;
    q.za = p->G13_a;

    // several ranks (G 40 1): the patches are placed for the load of the ranks; parents and old
    // patches on other ranks are reached through the block plans and old_run of the core.  Not
    // with the moving body (X 10 2/3): its fields on the patches come from the level-0 cells
    // below the patch (fs0_at)
    q.place = (p->G40>=1 && p->mpi_size>1 && shipmode==0) ? MIN(p->G40,2) : 0;
    q.rebalance = 0.1;
    if(p->G40>=1 && shipmode>0 && p->mpirank==0)
    cout<<"SFLOW AMR: G 40 is not available with the moving body (X 10 2/3), the patches stay on the rank below"<<endl;

    configure(q);

    printcount_amr = 0;
    printtime_amr = 0.0;
    m0 = 0.0;
    for(int k=0; k<10; ++k)
    tm[k]=0.0;

    pflow_void = new ioflow_v(p,pgc,pBC);

    phll0->amr = this;
    phll0->amr_id = 0;
}

sflow_amr::~sflow_amr()
{
    free_patches();

    for(auto v : nh0_v)
    delete v;
    told_free(told0);
    delete nhmg0;
    delete nhr0;
}

// --------------------------------------------------------------------- grids
sflow_amr::gh sflow_amr::grid(int g)
{
    gh r;
    if(g<0)
    {
        r.q = p0; r.b = b0; r.m = pmom0;
        r.oi = O0i; r.oj = O0j;
        return r;
    }
    sflow_amr_patch *c = SP(g);
    r.q = c->pp; r.b = c->b; r.m = c->pmom;
    r.oi = c->I0-EXT; r.oj = c->J0-EXT;
    return r;
}

// --------------------------------------------------------------------- setup
void sflow_amr::ini(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    // scope
    int ok=1;
    if(p->A220<0 || p->A220>4) ok=0;
    if(p->A210==2) ok=0;
    if(p->A260!=0) ok=0;
    if(p->S10!=0) ok=0;
    if(p->X10>3) ok=0;
    if(p->W90!=0) ok=0;
    if(p->j_dir!=1) ok=0;

    if(ok==0)
    {
        if(p->mpirank==0)
        cout<<"SFLOW AMR (G 1): only for A 220 0-4, A 210 3, A 260 0, S 10 0, X 10 0-3, W 90 0 and 2D grids -- refinement switched off"<<endl;
        maxlev=0;
        return;
    }

    // subcycling (G 7 1, sflow_amr_sub.cpp): hydrostatic, non-hydrostatic (A 220 1-3: level solves
    // and synchronisation projections) and Boussinesq (A 220 4: u_a per level), sflow_amr_subnh.cpp;
    // not with the moving body
    sub = (p->G7==1) ? 1 : 0;
    if(sub==1 && shipmode>0)
    {
        if(p->mpirank==0)
        cout<<"SFLOW AMR: G 7 1 (subcycling) not with a moving body -- one time step for all levels"<<endl;
        sub = 0;
    }

    // hierarchy: rank boxes, tiles, level-0 flags, no-refinement cells
    setup(p,pgc);

    // bed and still water depth in the level-0 halo (the fine bed is prolonged from it)
    pgc->gcslparax(p,b->bed,4);
    pgc->gcslparax(p,b->depth,4);

    mkdir("./REEF3D_SFLOW_AMR",0777);

    if(p->mpirank==0)
    {
        logout.open("./REEF3D_SFLOW_AMR/REEF3D_SFLOW_AMR_log.dat");
        logout<<"# count \t simtime \t dt \t patches \t cells \t mass \t rel. mass change"<<endl;
    }

    // initial hierarchy: every pass can add one level
    for(int it=0; it<maxlev; ++it)
    regrid(p,pgc,true);

    m0 = mass(p,b,pgc);

    if(p->mpirank==0)
    {
        cout<<"SFLOW AMR: "<<maxlev<<" level(s), "<<patches_total<<" patch(es), "<<cells_total<<" refined cells";
        if(regrid_int>0 && sub==0)
        cout<<", regrid every "<<regrid_int<<" steps";
        if(regrid_int>0 && sub==1)
        cout<<", regrid every "<<MAX(regrid_int>>maxlev,1)<<" level-0 steps";
        if(sub==1)
        cout<<", subcycled (G 7 1: level l takes 2^l steps per level-0 step)";
        cout<<endl;
        if(shipmode>0)
        cout<<"SFLOW AMR: moving body X 10 "<<shipmode<<" on the patches"<<(p->G12>0?", refinement around the body (G 12)":"")<<endl;
        if(p->G12>0 && shipmode==0)
        cout<<"SFLOW AMR: G 12 needs a moving body (X 10 2/3) -- ignored"<<endl;
    }
}

// --------------------------------------------------------------------- REEFAMR hooks
reefamr_patch* sflow_amr::patch_new()
{
    return new sflow_amr_patch;
}

// wall lists and SFLOW objects of a new patch (comms off); the bed and the state are set later
void sflow_amr::patch_objects(reefamr_patch *q, ghostcell *pgc)
{
    sflow_amr_patch *c = SP(q);
    lexer *p = p0;

    build_bc2D(pgc,*c);

    lexer *pp = c->pp;

    c->b = new fdm2D(pp);
    c->phll = new sflow_HLL(pp,pgc,pBC);
    c->phll->amr = this;
    c->pss = new sflow_signal_speed(pp);

    if(p->A211<=3)
    c->precon = new sflow_reconstruct_hires(pp,pBC);
    if(p->A211>=4)
    c->precon = new sflow_reconstruct_weno(pp,pBC,1);

    // diffusion: explicit on the patches (A 212 2 is implicit on level 0 only)
    if(p->A212>=1)
    c->pdiff = new sflow_ediff(pp);
    else
    c->pdiff = new sflow_diffusion_void(pp);
    if(nh==1)
    {
    if(p->A220==1)
    c->pnh = new sflow_pjm_lin(pp,c->b,pBC);
    else
    c->pnh = new sflow_pjm_quad(pp,c->b,pgc,pBC);
    c->ppress = c->pnh;
    }
    else
    c->ppress = new sflow_hydrostatic(pp,c->b,pBC);
    c->pfsf = new sflow_eta(pp,c->b,pgc,pBC);
    c->psfdf = new sflow_forcing(pp);
    sixdof *p6 = p6dof;
    if(shipmode>0)
    {
    c->pship = new sflow_amr_ship(pp);
    p6 = c->pship;
    }
    if(bous==1)
    c->psolv = new sflow_amr_linesolve;
    c->pmom = new sflow_momentum_RK3(pp,c->b,pgc,c->phll,c->pss,c->precon,c->pdiff,c->ppress,c->psolv,nullptr,pflow_void,c->pfsf,c->psfdf,p6);
    c->pmom->nh_defer = (nh==1);
    c->pmom->bous_defer = (bous==1);

    for(int ip=0; ip<5; ++ip)
    {
        c->rec[ip][0].assign(c->ny,0.0);
        c->rec[ip][1].assign(c->ny,0.0);
        c->rec[ip][2].assign(c->nx,0.0);
        c->rec[ip][3].assign(c->nx,0.0);
    }

    if(sub==1)
    for(int ip=0; ip<5; ++ip)
    {
        c->freg[ip][0].assign(c->ny,0.0);
        c->freg[ip][1].assign(c->ny,0.0);
        c->freg[ip][2].assign(c->nx,0.0);
        c->freg[ip][3].assign(c->nx,0.0);
    }
}

void sflow_amr::patch_delete(reefamr_patch *q)
{
    sflow_amr_patch *c = SP(q);

    for(auto v : c->nv)
    delete v;
    delete c->mg;
    told_free(c->told);

    delete c->pmom;
    delete c->pship;
    delete c->psolv;
    delete c->psfdf;
    delete c->pfsf;
    delete c->ppress;
    delete c->pdiff;
    delete c->precon;
    delete c->pss;
    delete c->phll;
    delete c->b;
}

void sflow_amr::zone_bodies(vector<sixdof_obj*> &obj)
{
    if(ship6!=nullptr)
    for(int nb=0; nb<ship6->objects(); ++nb)
    obj.push_back(ship6->object(nb));
}

double sflow_amr::fill_aux(reefamr_patch *c, int ii, int jj)
{
    return SP(c)->b->depth(ii,jj);
}

double sflow_amr::serve_aux(int l, int I, int J)
{
    return p0->wd - bed_at(l,I,J);
}

void sflow_amr::regrid_ids()
{
    for(int n=0; n<(int)P.size(); ++n)
    SP(n)->phll->amr_id = n+1;
}

// current hierarchy: final state with filled ghost cells
void sflow_amr::regrid_prepare(ghostcell *pgc)
{
    cache_stage(0);
    for(int l=1; l<=maxlev; ++l)
    fill_level(pgc,l,0);
}

// bed of the new patches: prolonged from level 0, cells of other ranks from their owner
void sflow_amr::regrid_static(ghostcell *pgc)
{
    lexer *p = p0;
    const int me = p->mpirank;

    {
    vector<vector<int>> req(p->mpi_size), srv;
    vector<vector<double*>> dst(p->mpi_size);

    for(auto q : P)
    if(q->fresh)
    {
        sflow_amr_patch *c = SP(q);
        lexer *pp = c->pp;
        const int l = c->lev;
        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            int I = ii-EXT+c->I0, J = jj-EXT+c->J0;
            int o = owner(fsh(I,l),fsh(J,l));

            if(o<0 || o==me)
            c->b->bed(ii,jj) = bed_at(l,I,J);

            if(o>=0 && o!=me)
            {
                req[o].push_back(l); req[o].push_back(I); req[o].push_back(J);
                dst[o].push_back(&c->b->bed(ii,jj));
            }
        }
    }

    xsetup(req,srv,3);

    reefamr_xplan X;
    for(int r=0; r<p->mpi_size; ++r)
    {
        if(!srv[r].empty())
        {
            X.speer.push_back(r);
            vector<double> v;
            for(size_t k=0; k<srv[r].size(); k+=3)
            v.push_back(bed_at(srv[r][k],srv[r][k+1],srv[r][k+2]));
            X.sbuf.push_back(v);
        }
        if(!req[r].empty())
        {
            X.rpeer.push_back(r);
            X.rcount.push_back((int)req[r].size()/3);
        }
    }
    xrun(X,1,7001);

    for(size_t k=0; k<X.rpeer.size(); ++k)
    {
        int r = X.rpeer[k];
        for(size_t m=0; m<dst[r].size(); ++m)
        *dst[r][m] = X.rbuf[k][m];
    }

    for(auto q : P)
    if(q->fresh)
    {
        sflow_amr_patch *c = SP(q);
        lexer *pp = c->pp;
        fdm2D *pb = c->b;
        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            pb->bed0(ii,jj) = pb->bed(ii,jj);
            pb->topobed(ii,jj) = pb->bed(ii,jj);
            pb->solidbed(ii,jj) = pb->bed(ii,jj);
            pb->depth(ii,jj) = p->wd - pb->bed(ii,jj);
            pb->ks(ii,jj) = p->B50;
        }

        // wall ghost cells as sflow_eta::depth_update sets them at the end of every step (mirrored);
        // bed_at gives them the bed of the nearest level-0 cell instead.  A patch created at a regrid
        // would otherwise run its first step with other ghost depths than a kept patch covering the
        // same cells, and the result would depend on how the level is split into patches and ranks
        {
        comms_off guard(pgc);
        pgc->gcsl_start4(pp,pb->depth,50);
        }

        // face depth as sflow_reconstruct::reconstruct_WL
        for(int ii=pp->imin; ii<pp->imin+pp->imax-1; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        pb->dfx(ii,jj) = 0.5*(pb->depth(ii+1,jj)+pb->depth(ii,jj));

        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax-1; ++jj)
        pb->dfy(ii,jj) = 0.5*(pb->depth(ii,jj+1)+pb->depth(ii,jj));
    }
    }
}

// state of the new patches, coarse to fine: the cells an old patch of the same level held (on
// any rank), the other blocks prolonged from their parent (on its rank)
void sflow_amr::regrid_state(ghostcell *pgc, vector<reefamr_patch*> &oldP)
{
    cache_stage(0);

    for(int l=1; l<=maxlev; ++l)
    {
        {
        comms_off guard(pgc);
        for(int n : lev[l])
        if(P[n]->fresh)
        ini_patch_defaults(*SP(n));
        }

        ini_patch_old(l,oldP);
        ini_patch_prolong(l);

        {
        comms_off guard(pgc);
        for(int n : lev[l])
        if(P[n]->fresh)
        ini_patch_finish(p0,*SP(n));
        }

        fill_level(pgc,l,0);
    }
}

// level 0 consistent with the patches
void sflow_amr::regrid_finish(ghostcell *pgc, int old_total)
{
    restrict_levels(p0,2);

    if(patches_total>0 || old_total>0)
    exchange_level0(p0,b0,pgc,2);

    // non-hydrostatic pressure outside the patch interiors, as after a pressure solve: the rows
    // of the next solve take it in the computed cells beyond the patch (Uest, bed acceleration),
    // where a new patch otherwise has 0 and a kept one the values of the last solve
    if(nh==1)
    for(int l=1; l<=maxlev; ++l)
    nh_qfill(l,-1);
}

// fine bed: limited-linear prolongation, level by level from level 0 (children average = parent);
// solid cells and cells next to solids take the parent value without slopes
double sflow_amr::bed_at(int l, int I, int J)
{
    if(l==0)
    {
        // cells outside the domain take the bed of the nearest domain cell
        I = MAX(MIN(I,GNX-1),0);
        J = MAX(MIN(J,GNY-1),0);
        int ii = MAX(MIN(I-O0i, p0->imin+p0->imax-1), p0->imin);
        int jj = MAX(MIN(J-O0j, p0->jmin+p0->jmax-1), p0->jmin);
        return b0->bed(ii,jj);
    }

    const uint64_t key = (uint64_t(l)<<58) ^ (uint64_t(uint32_t(I+(1<<28)))<<29) ^ uint64_t(uint32_t(J+(1<<28)));
    auto it = bedmemo.find(key);
    if(it!=bedmemo.end())
    return it->second;

    int Ic = fsh(I,1), Jc = fsh(J,1);
    double bc = bed_at(l-1,Ic,Jc);

    if(flag0(fsh(I,l),fsh(J,l))<0)
    {
        bedmemo[key]=bc;
        return bc;
    }

    const int lc = l-1;
    auto fl = [&](int a, int d) { return flag0(fsh(a,lc),fsh(d,lc)); };

    double sx=0.0, sy=0.0;
    if(fl(Ic+1,Jc)>0 && fl(Ic-1,Jc)>0)
    sx = mmod(bed_at(lc,Ic+1,Jc)-bc, bc-bed_at(lc,Ic-1,Jc));
    if(fl(Ic,Jc+1)>0 && fl(Ic,Jc-1)>0)
    sy = mmod(bed_at(lc,Ic,Jc+1)-bc, bc-bed_at(lc,Ic,Jc-1));

    double ox = (I-2*Ic==0) ? -0.25 : 0.25;
    double oy = (J-2*Jc==0) ? -0.25 : 0.25;
    double r = bc + sx*ox + sy*oy;
    bedmemo[key]=r;
    return r;
}

// --------------------------------------------------------------------- refinement flags
// cells of level l on this rank that need level l+1 (evaluated on the level-l grids)
void sflow_amr::tag(int l, vector<unsigned char> &M)
{
    const int lf = l+1;
    const double hfilm = 3.0*p0->A244;

    auto test = [&](lexer *q, fdm2D *bb, slice &WL, int ii, int jj)
    {
        if(q->flagslice4[lij(q,ii,jj)]<0)
        return false;

        int wc = q->wet[lij(q,ii,jj)];
        double ec = WL(ii,jj)-bb->depth(ii,jj);
        const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};

        for(int k=0; k<4; ++k)
        {
            int a=ii+di[k], d=jj+dj[k];
            if(q->flagslice4[lij(q,a,d)]<0)
            continue;

            int wn = q->wet[lij(q,a,d)];

            if(tol_eta>0.0 && wc==1 && wn==1)
            if(fabs(ec-(WL(a,d)-bb->depth(a,d)))>tol_eta)
            return true;

            if(shore>0 && wc==1 && wn==0 && WL(ii,jj)>hfilm)
            return true;
        }
        return false;
    };

    if(l==0)
    {
        for(int ii=0; ii<NX0; ++ii)
        for(int jj=0; jj<NY0; ++jj)
        if(test(p0,b0,b0->WL,ii,jj))
        tag_cell(lf,ii+O0i,jj+O0j,nbuf,M);
        return;
    }

    for(int n : lev[l])
    {
        sflow_amr_patch *c = SP(n);
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
        if(test(c->pp,c->b,c->b->WL,ii,jj))
        tag_cell(lf,ii-EXT+c->I0,jj-EXT+c->J0,nbuf,M);
    }
}

// interior of a new patch: old patches of the same level where they existed, otherwise
// prolonged from the new coarser level (exact mass and momentum of fully wet parents)
void sflow_amr::ini_patch_defaults(sflow_amr_patch &c)
{
    lexer *pp = c.pp;
    fdm2D *b = c.b;
    const double wd = p0->A244;

    // defaults for all cells (cells outside the interior are filled afterwards)
    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        b->WL(ii,jj) = wd;
        b->UH(ii,jj) = b->VH(ii,jj) = b->WH(ii,jj) = 0.0;
        b->eta(ii,jj) = wd - b->depth(ii,jj);
        b->U(ii,jj) = b->V(ii,jj) = b->W(ii,jj) = 0.0;
        b->press(ii,jj) = 0.0;
        b->UA(ii,jj) = b->VA(ii,jj) = b->MX(ii,jj) = b->MY(ii,jj) = 0.0;
        pp->wet[lij(pp,ii,jj)] = 0;
        pp->deep[lij(pp,ii,jj)] = 0;
    }

    c.blkold.assign((size_t)(c.nx/2)*(c.ny/2),0);
}

// a solid block (its first cell) keeps the defaults
bool sflow_amr::block_solid(sflow_amr_patch &c, int bi, int bj)
{
    return c.pp->flagslice4[lij(c.pp,EXT+2*bi,EXT+2*bj)]<0;
}

// the cells of the fresh level-l patches that an old patch held (all four cells of a block lie in
// the same old patch: patch boxes are tile aligned)
void sflow_amr::ini_patch_old(int l, vector<reefamr_patch*> &oldP)
{
    old_run(l,15,7140+l,oldP,
            [&](reefamr_patch *q, int io, int jo, double *v)
            {
                sflow_amr_patch *o = SP(q);
                lexer *qp = o->pp;
                fdm2D *ob = o->b;
                v[0] = ob->WL(io,jo);  v[1] = ob->UH(io,jo);  v[2] = ob->VH(io,jo);
                v[3] = ob->eta(io,jo); v[4] = ob->U(io,jo);   v[5] = ob->V(io,jo);
                v[6] = ob->WH(io,jo);  v[7] = ob->W(io,jo);   v[8] = ob->press(io,jo);
                v[9] = ob->UA(io,jo);  v[10] = ob->VA(io,jo); v[11] = ob->MX(io,jo); v[12] = ob->MY(io,jo);
                v[13] = qp->wet[lij(qp,io,jo)];
                v[14] = qp->deep[lij(qp,io,jo)];
            },
            [&](reefamr_patch *q, int id, int ii, int jj, const double *v)
            {
                sflow_amr_patch &c = *SP(q);
                const int bi = (ii-EXT)/2, bj = (jj-EXT)/2;
                if(block_solid(c,bi,bj))
                return;
                lexer *pp = c.pp;
                fdm2D *b = c.b;
                b->WL(ii,jj) = v[0];  b->UH(ii,jj) = v[1];  b->VH(ii,jj) = v[2];
                b->eta(ii,jj) = v[3]; b->U(ii,jj) = v[4];   b->V(ii,jj) = v[5];
                b->WH(ii,jj) = v[6];  b->W(ii,jj) = v[7];   b->press(ii,jj) = v[8];
                b->UA(ii,jj) = v[9];  b->VA(ii,jj) = v[10]; b->MX(ii,jj) = v[11]; b->MY(ii,jj) = v[12];
                pp->wet[lij(pp,ii,jj)] = (int)v[13];
                pp->deep[lij(pp,ii,jj)] = (int)v[14];
                c.blkold[(size_t)bi*(c.ny/2)+bj] = 1;
            });
}

// the other blocks of the fresh level-l patches: prolonged from the parent (prolong_parts of the
// four children and the parent's values on its rank), finished with the child depths, then the
// mass and momentum of a fully wet parent restored exactly
void sflow_amr::ini_patch_prolong(int l)
{
    const int nv = 4*NPR+6;

    block_down_if(l,nv,7160+l,
                  [&](reefamr_patch *q) { return q->fresh; },
                  [&](const reefamr_block &B, int key, double *v)
                  {
                      gh G = grid(B.g);
                      const int ic = B.ic, jc = B.jc;
                      for(int a=0;a<2;++a)
                      for(int d=0;d<2;++d)
                      prolong_parts(B.g,ic,jc,a==0?-1:1,d==0?-1:1,&v[(2*a+d)*NPR]);
                      double *w = &v[4*NPR];
                      w[0] = G.q->wet[lij(G.q,ic,jc)];
                      w[1] = (*sWL[B.g+1])(ic,jc);
                      w[2] = (*sUH[B.g+1])(ic,jc);
                      w[3] = (*sVH[B.g+1])(ic,jc);
                      w[4] = (*sWH[B.g+1])(ic,jc);
                      w[5] = G.b->press(ic,jc);
                  },
                  [&](reefamr_patch *q, int id, int kb, const double *v)
                  {
                      sflow_amr_patch &c = *SP(q);
                      const int nby = c.ny/2;
                      const int bi = kb/nby, bj = kb%nby;
                      if(block_solid(c,bi,bj) || c.blkold[kb])
                      return;
                      ini_block(c,bi,bj,v);
                  });
}

void sflow_amr::ini_block(sflow_amr_patch &c, int bi, int bj, const double *pv)
{
    lexer *pp = c.pp;
    fdm2D *b = c.b;
    const double wd = p0->A244;
    const int i0 = EXT+2*bi, j0 = EXT+2*bj;
    const double *w = &pv[4*NPR];
    double v[NVMAX];

    for(int a=0;a<2;++a)
    for(int d=0;d<2;++d)
    {
        int ii = i0+a, jj = j0+d;
        prolong_finish(&pv[(2*a+d)*NPR],b->depth(ii,jj),v);
        b->WL(ii,jj)=v[0]; b->UH(ii,jj)=v[1]; b->VH(ii,jj)=v[2]; b->WH(ii,jj)=v[8];
        b->press(ii,jj) = w[5];
        if(bous==1)
        {
        b->UA(ii,jj)=v[10]; b->VA(ii,jj)=v[11]; b->MX(ii,jj)=v[12]; b->MY(ii,jj)=v[13];
        }
        pp->wet[lij(pp,ii,jj)]=(int)v[6];
    }

    int nw = pp->wet[lij(pp,i0,j0)] + pp->wet[lij(pp,i0+1,j0)] + pp->wet[lij(pp,i0,j0+1)] + pp->wet[lij(pp,i0+1,j0+1)];

    if(nw==4 && (int)w[0]==1)
    {
        double sw = b->WL(i0,j0)+b->WL(i0+1,j0)+b->WL(i0,j0+1)+b->WL(i0+1,j0+1);
        double su = b->UH(i0,j0)+b->UH(i0+1,j0)+b->UH(i0,j0+1)+b->UH(i0+1,j0+1);
        double sv = b->VH(i0,j0)+b->VH(i0+1,j0)+b->VH(i0,j0+1)+b->VH(i0+1,j0+1);
        double sh = b->WH(i0,j0)+b->WH(i0+1,j0)+b->WH(i0,j0+1)+b->WH(i0+1,j0+1);
        double du = 4.0*w[2] - su;
        double dv = 4.0*w[3] - sv;
        double dw = 4.0*w[4] - sh;
        double fw = 4.0*w[1]/sw;

        for(int a=0;a<2;++a)
        for(int d=0;d<2;++d)
        {
            b->UH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*du;
            b->VH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*dv;
            b->WH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*dw;
            b->WL(i0+a,j0+d) *= fw;
        }
    }

    for(int a=0;a<2;++a)
    for(int d=0;d<2;++d)
    {
        int ii = i0+a, jj = j0+d;
        int wt = pp->wet[lij(pp,ii,jj)];
        double wlvl = fabs(b->WL(ii,jj))>wd ? b->WL(ii,jj) : 1.0e20;
        b->eta(ii,jj) = b->WL(ii,jj) - b->depth(ii,jj);
        b->U(ii,jj) = wt==1 ? b->UH(ii,jj)/wlvl : 0.0;
        b->V(ii,jj) = wt==1 ? b->VH(ii,jj)/wlvl*p0->y_dir : 0.0;
        b->W(ii,jj) = wt==1 ? b->WH(ii,jj)/wlvl : 0.0;
        if(bous==1)
        {
        b->U(ii,jj) = wt==1 ? b->MX(ii,jj)/wlvl : 0.0;
        b->V(ii,jj) = wt==1 ? b->MY(ii,jj)/wlvl*p0->y_dir : 0.0;
        }
        pp->deep[lij(pp,ii,jj)] = wt;
    }
}

void sflow_amr::ini_patch_finish(lexer *p, sflow_amr_patch &c)
{
    lexer *pp = c.pp;
    fdm2D *b = c.b;

    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        b->eta_n(ii,jj) = b->eta(ii,jj);
        b->hp(ii,jj) = b->WL(ii,jj);
        b->breaking(ii,jj) = 0;
        pp->wet_n[lij(pp,ii,jj)] = pp->wet[lij(pp,ii,jj)];
    }

    pp->dt = p->dt;
    pp->dt_old = p->dt_old;
    pp->simtime = p->simtime;
    pp->count = p->count;
}

// --------------------------------------------------------------------- stage data
void sflow_amr::cache_stage(int s)
{
    sWL.assign(P.size()+1,nullptr);
    sUH.assign(P.size()+1,nullptr);
    sVH.assign(P.size()+1,nullptr);
    sWH.assign(P.size()+1,nullptr);

    slice *WLo,*UHo,*VHo,*WHo;
    pmom0->stage_io(s,b0,sWL[0],sUH[0],sVH[0],WLo,UHo,VHo);
    pmom0->stage_io_w(s,b0,sWH[0],WHo);

    for(size_t n=0; n<P.size(); ++n)
    {
    SP(n)->pmom->stage_io(s,SP(n)->b,sWL[n+1],sUH[n+1],sVH[n+1],WLo,UHo,VHo);
    SP(n)->pmom->stage_io_w(s,SP(n)->b,sWH[n+1],WHo);
    }
}

// well balanced prolongation from grid g (coarse cell ic,jc): the surface eta is interpolated
// with minmod slopes (switched off next to dry or solid cells), the fine depth follows from the
// fine bed; momentum through u.  Thin films (h_c < 3 A 244) copy the parent depth.
// v: WL, UH, VH, eta, U, V, wet, deep, WH, W
// prolongation of the coarse cell (ic,jc) of grid g into its child in quadrant (ox,oy), in two
// parts: prolong_parts on the rank of the parent (interpolated surface and velocities, without
// the child), prolong_finish with the child depth (on the rank of the child: placed patches)
void sflow_amr::prolong(int g, int ic, int jc, int ox, int oy, double depthf, double *v)
{
    double r[NPR];
    prolong_parts(g,ic,jc,ox,oy,r);
    prolong_finish(r,depthf,v);
}

// r[0] 0 dry parent, 1 shallow parent (its values), 2 interpolated; r[1] surface (2) or water
// level (1), r[2..4] velocities (2) or UH, VH, WH (1), r[5..8] Boussinesq u_a, v_a, M/H
void sflow_amr::prolong_parts(int g, int ic, int jc, int ox, int oy, double *r)
{
    gh G = grid(g);
    fdm2D *pb = G.b;
    lexer *q = G.q;
    slice &WLp = *sWL[g+1], &UHp = *sUH[g+1], &VHp = *sVH[g+1], &WHp = *sWH[g+1];
    const double wd = p0->A244;

    for(int m=0; m<NPR; ++m)
    r[m] = 0.0;

    // parent state: the stage input, or with G 7 1 (tint >= 0) linear in time between the start
    // of the parent step (told) and its end (the current state); wet where either was wet and the
    // interpolated water level is above the wet-dry depth
    const bool ti = (tint>=0.0);
    sflow_amr_told *O = ti ? &told(g) : nullptr;
    auto tv = [&](int k, slice &N, int a, int bb) { return (1.0-tint)*(*O->f[k])(a,bb) + tint*N(a,bb); };
    auto wlat = [&](int a, int bb) { return ti ? tv(0,pb->WL,a,bb) : WLp(a,bb); };
    auto uhat = [&](int a, int bb) { return ti ? tv(1,pb->UH,a,bb) : UHp(a,bb); };
    auto vhat = [&](int a, int bb) { return ti ? tv(2,pb->VH,a,bb) : VHp(a,bb); };
    auto what = [&](int a, int bb) { return ti ? tv(3,pb->WH,a,bb) : WHp(a,bb); };
    auto wetat = [&](int a, int bb, double wlc)
    {
        const int k = lij(q,a,bb);
        if(!ti || tint==1.0)
        return q->wet[k];
        if(tint==0.0)
        return O->wet[k];
        return ((O->wet[k]==1 || q->wet[k]==1) && wlc>wd+eps) ? 1 : 0;
    };

    auto get = [&](int a, int bb, double &e, double &u, double &vv, double &ww, int &wt)
    {
        double wlc = wlat(a,bb);
        e = wlc - pb->depth(a,bb);
        wt = wetat(a,bb,wlc);
        if(q->flagslice4[lij(q,a,bb)]<0)
        wt = -1;
        double wlvl = fabs(wlc)>wd ? wlc : 1.0e20;
        u = wt==1 ? uhat(a,bb)/wlvl : 0.0;
        vv = wt==1 ? vhat(a,bb)/wlvl : 0.0;
        ww = wt==1 ? what(a,bb)/wlvl : 0.0;
    };

    // Boussinesq: u_a and M/H of the coarse cell (the parent state at the start of the stage; G 7 1:
    // linear in time as above)
    auto getb = [&](int a, int bb, double *w)
    {
        if(ti)
        {
            double wlc = wlat(a,bb);
            double wlvl = fabs(wlc)>wd ? wlc : 1.0e20;
            int wt = wetat(a,bb,wlc)==1 && q->flagslice4[lij(q,a,bb)]>0;
            w[0] = wt ? tv(5,pb->UA,a,bb) : 0.0;
            w[1] = wt ? tv(6,pb->VA,a,bb) : 0.0;
            w[2] = wt ? tv(7,pb->MX,a,bb)/wlvl : 0.0;
            w[3] = wt ? tv(8,pb->MY,a,bb)/wlvl : 0.0;
            return;
        }
        double wlc = WLp(a,bb);
        double wlvl = fabs(wlc)>wd ? wlc : 1.0e20;
        int wt = q->wet[lij(q,a,bb)]==1 && q->flagslice4[lij(q,a,bb)]>0;
        w[0] = wt ? pb->UA(a,bb) : 0.0;
        w[1] = wt ? pb->VA(a,bb) : 0.0;
        w[2] = wt ? pb->MX(a,bb)/wlvl : 0.0;
        w[3] = wt ? pb->MY(a,bb)/wlvl : 0.0;
    };

    // Boussinesq: limited slopes of u_a and M/H
    auto partb = [&](double sx, double sy)
    {
        double c0[4],cE[4],cW[4],cN[4],cS[4];
        getb(ic,jc,c0);
        getb(ic+1,jc,cE); getb(ic-1,jc,cW); getb(ic,jc+1,cN); getb(ic,jc-1,cS);
        for(int k=0; k<4; ++k)
        r[5+k] = c0[k] + sx*mmod(cE[k]-c0[k],c0[k]-cW[k]) + sy*mmod(cN[k]-c0[k],c0[k]-cS[k]);
    };

    double e0,u0,v0,w0v; int w0;
    get(ic,jc,e0,u0,v0,w0v,w0);

    if(w0!=1)
    return;

    double hc = wlat(ic,jc);
    if(hc<3.0*wd)
    {
        r[0] = 1.0;
        r[1] = hc; r[2] = uhat(ic,jc); r[3] = vhat(ic,jc); r[4] = what(ic,jc);
        if(bous==1)
        partb(0.0,0.0);
        return;
    }

    double eE,uE,vE,wEv,eW,uW,vW,wWv,eN,uN,vN,wNv,eS,uS,vS,wSv;
    int wE,wW,wN,wS;
    get(ic+1,jc,eE,uE,vE,wEv,wE);
    get(ic-1,jc,eW,uW,vW,wWv,wW);
    get(ic,jc+1,eN,uN,vN,wNv,wN);
    get(ic,jc-1,eS,uS,vS,wSv,wS);

    double sxe=0.0,sye=0.0,sxu=0.0,syu=0.0,sxv=0.0,syv=0.0,sxw=0.0,syw=0.0;
    if(wE==1 && wW==1 && wN==1 && wS==1)
    {
        sxe = mmod(eE-e0,e0-eW); sye = mmod(eN-e0,e0-eS);
        sxu = mmod(uE-u0,u0-uW); syu = mmod(uN-u0,u0-uS);
        sxv = mmod(vE-v0,v0-vW); syv = mmod(vN-v0,v0-vS);
        sxw = mmod(wEv-w0v,w0v-wWv); syw = mmod(wNv-w0v,w0v-wSv);
    }

    const double fx = 0.25*ox, fy = 0.25*oy;
    r[0] = 2.0;
    r[1] = e0 + sxe*fx + sye*fy;
    r[2] = u0 + sxu*fx + syu*fy;
    r[3] = v0 + sxv*fx + syv*fy;
    r[4] = w0v + sxw*fx + syw*fy;

    if(bous==1)
    {
        bool full = wE==1 && wW==1 && wN==1 && wS==1;
        partb(full ? fx : 0.0, full ? fy : 0.0);
    }
}

void sflow_amr::prolong_finish(const double *r, double depthf, double *v)
{
    const double wd = p0->A244;

    auto dry = [&]()
    {
        v[0]=wd; v[1]=v[2]=0.0; v[3]=wd-depthf; v[4]=v[5]=0.0; v[6]=0.0; v[7]=0.0; v[8]=v[9]=0.0;
        if(bous==1)
        v[10]=v[11]=v[12]=v[13]=0.0;
    };

    auto fin = [&](double wl, double uh, double vh, double wh)
    {
        double wlvl = fabs(wl)>wd ? wl : 1.0e20;
        v[0]=wl; v[1]=uh; v[2]=vh;
        v[3]=wl-depthf;
        v[4]=uh/wlvl;
        v[5]=vh/wlvl*p0->y_dir;
        v[6]=1.0; v[7]=1.0;
        v[8]=wh; v[9]=wh/wlvl;
    };

    // Boussinesq: u_a and M, U = M/H
    auto finb = [&](double hf)
    {
        for(int k=0; k<4; ++k)
        {
            if(k<2) v[10+k] = r[5+k];
            else    v[10+k] = hf*r[5+k];
        }
        v[4] = v[12]/(fabs(hf)>wd ? hf : 1.0e20);
        v[5] = v[13]/(fabs(hf)>wd ? hf : 1.0e20)*p0->y_dir;
    };

    if(r[0]==0.0)
    {
        dry();
        return;
    }

    if(r[0]==1.0)
    {
        fin(r[1],r[2],r[3],r[4]);
        if(bous==1)
        finb(r[1]);
        return;
    }

    double hf = r[1] + depthf;

    if(hf<=wd+eps)
    {
        dry();
        return;
    }

    fin(hf, hf*r[2], hf*r[3], hf*r[4]);

    if(bous==1)
    finb(hf);
}

void sflow_amr::eval_fill(const reefamr_fill &f, double *v)
{
    if(f.kind==1)
    {
        prolong(f.g,f.si,f.sj,f.ox,f.oy,f.aux,v);
        return;
    }

    if(f.kind==0)
    {
        gh G = grid(f.g);
        v[0] = (*sWL[f.g+1])(f.si,f.sj);
        v[1] = (*sUH[f.g+1])(f.si,f.sj);
        v[2] = (*sVH[f.g+1])(f.si,f.sj);
        v[3] = G.b->eta(f.si,f.sj);
        v[4] = G.b->U(f.si,f.sj);
        v[5] = G.b->V(f.si,f.sj);
        v[6] = G.q->wet[lij(G.q,f.si,f.sj)];
        v[7] = G.q->deep[lij(G.q,f.si,f.sj)];
        v[8] = (*sWH[f.g+1])(f.si,f.sj);
        v[9] = G.b->W(f.si,f.sj);
        if(bous==1)
        {
        v[10] = G.b->UA(f.si,f.sj);
        v[11] = G.b->VA(f.si,f.sj);
        v[12] = G.b->MX(f.si,f.sj);
        v[13] = G.b->MY(f.si,f.sj);
        }
        return;
    }

    const double wd = p0->A244;
    v[0]=wd; v[1]=v[2]=0.0; v[3]=wd-f.aux; v[4]=v[5]=0.0; v[6]=v[7]=0.0; v[8]=v[9]=0.0;
    if(bous==1)
    v[10]=v[11]=v[12]=v[13]=0.0;
}

// values of a filled cell into patch c (grid id n)
void sflow_amr::store_fill(sflow_amr_patch *c, int n, int ii, int jj, const double *w)
{
    lexer *pp = c->pp;
    fdm2D *b = c->b;
    (*sWL[n+1])(ii,jj) = w[0];
    (*sUH[n+1])(ii,jj) = w[1];
    (*sVH[n+1])(ii,jj) = w[2];
    b->eta(ii,jj) = w[3];
    b->U(ii,jj) = w[4];
    b->V(ii,jj) = w[5];
    b->hp(ii,jj) = w[0];
    pp->wet[lij(pp,ii,jj)] = (int)w[6];
    pp->deep[lij(pp,ii,jj)] = (int)w[7];
    (*sWH[n+1])(ii,jj) = w[8];
    b->W(ii,jj) = w[9];
    if(bous==1)
    {
    b->UA(ii,jj) = w[10];
    b->VA(ii,jj) = w[11];
    b->MX(ii,jj) = w[12];
    b->MY(ii,jj) = w[13];
    }
}

// cells of the level-l patches outside their interior, then the wall conditions
void sflow_amr::fill_level(ghostcell *pgc, int l, int s)
{
    fill_run(l,NV,7100+l,
             [&](const reefamr_fill &f, double *v) { eval_fill(f,v); },
             [&](reefamr_patch *c, int n, const reefamr_fill &f, const double *w) { store_fill(SP(c),n,f.di,f.dj,w); });

    comms_off guard(pgc);

    for(int n : lev[l])
    apply_bc(pgc,*SP(n),s);
}

// wall ghost cells of a patch from its (filled) fluid cells, as at the end of a level-0 stage
void sflow_amr::apply_bc(ghostcell *pgc, sflow_amr_patch &c, int s)
{
    lexer *pp = c.pp;
    fdm2D *b = c.b;
    slice *WLi,*UHi,*VHi,*WLo,*UHo,*VHo;
    c.pmom->stage_io(s,b,WLi,UHi,VHi,WLo,UHo,VHo);

    slice *WHi,*WHo;
    c.pmom->stage_io_w(s,b,WHi,WHo);

    int gcval_eta = 50+p0->F50;

    pgc->gcsl_start4(pp,b->eta,gcval_eta);
    pgc->gcsl_start4Vint(pp,pp->wet,50);
    pgc->gcsl_start4Vint(pp,pp->deep,50);
    c.pmom->ghostcells(pp,b,pgc,*UHi,*VHi,*WHi,*WLi);
}

// every patch into the covered cells of its parents (conservative averages of the 2x2 children),
// finest level first; the parent cells may lie on another rank (block plans)
void sflow_amr::restrict_levels(lexer *p, int s)
{
    for(int l=maxlev; l>=1; --l)
    restrict_level(p,l,s);
}

void sflow_amr::restrict_level(lexer *p, int l, int s)
{
    const int nv = (bous==1) ? 8 : 4;

    block_up(l,nv,7020+l,
             [&](reefamr_patch *q, int id, int k, double *v)
             {
                 sflow_amr_patch &c = *SP(q);
                 fdm2D *b = c.b;
                 slice *WLi,*UHi,*VHi,*WLo,*UHo,*VHo;
                 c.pmom->stage_io(s,b,WLi,UHi,VHi,WLo,UHo,VHo);
                 slice *WHi,*WHo;
                 c.pmom->stage_io_w(s,b,WHi,WHo);

                 const int nby = c.ny/2;
                 const int i0 = EXT+2*(k/nby), j0 = EXT+2*(k%nby);

                 v[0] = 0.25*((*WLo)(i0,j0)+(*WLo)(i0+1,j0)+(*WLo)(i0,j0+1)+(*WLo)(i0+1,j0+1));
                 v[1] = 0.25*((*UHo)(i0,j0)+(*UHo)(i0+1,j0)+(*UHo)(i0,j0+1)+(*UHo)(i0+1,j0+1));
                 v[2] = 0.25*((*VHo)(i0,j0)+(*VHo)(i0+1,j0)+(*VHo)(i0,j0+1)+(*VHo)(i0+1,j0+1));
                 v[3] = 0.25*((*WHo)(i0,j0)+(*WHo)(i0+1,j0)+(*WHo)(i0,j0+1)+(*WHo)(i0+1,j0+1));

                 // Boussinesq: u_a and M
                 if(bous==1)
                 {
                 auto avg = [&](slice &f) { return 0.25*(f(i0,j0)+f(i0+1,j0)+f(i0,j0+1)+f(i0+1,j0+1)); };
                 v[4] = avg(b->UA);
                 v[5] = avg(b->VA);
                 v[6] = avg(b->MX);
                 v[7] = avg(b->MY);
                 }
             },
             [&](const reefamr_block &B, int key, const double *v)
             {
                 gh G = grid(B.g);
                 slice *pWLi,*pUHi,*pVHi,*pWLo,*pUHo,*pVHo;
                 G.m->stage_io(s,G.b,pWLi,pUHi,pVHi,pWLo,pUHo,pVHo);
                 slice *pWHi,*pWHo;
                 G.m->stage_io_w(s,G.b,pWHi,pWHo);

                 const int ic = B.ic, jc = B.jc;
                 const double wl = v[0], uh = v[1], vh = v[2], wh = v[3];

                 (*pWLo)(ic,jc)=wl; (*pUHo)(ic,jc)=uh; (*pVHo)(ic,jc)=vh; (*pWHo)(ic,jc)=wh;

                 int w = wl>p->A244+eps ? 1 : 0;
                 G.q->wet[lij(G.q,ic,jc)] = w;

                 double wlvl = fabs(wl)>p->A244 ? wl : 1.0e20;
                 G.b->eta(ic,jc) = wl - G.b->depth(ic,jc);
                 G.b->U(ic,jc) = w==1 ? uh/wlvl : 0.0;
                 G.b->V(ic,jc) = w==1 ? vh/wlvl : 0.0;
                 G.b->W(ic,jc) = w==1 ? wh/wlvl : 0.0;
                 G.b->hp(ic,jc) = wl;

                 // Boussinesq: u_a and M with V, U = M/H
                 if(bous==1)
                 {
                 G.b->UA(ic,jc) = v[4];
                 G.b->VA(ic,jc) = v[5];
                 G.b->MX(ic,jc) = v[6];
                 G.b->MY(ic,jc) = v[7];
                 G.b->U(ic,jc) = w==1 ? G.b->MX(ic,jc)/wlvl : 0.0;
                 G.b->V(ic,jc) = w==1 ? G.b->MY(ic,jc)/wlvl : 0.0;
                 }
             });
}

// level-0 halo after the restriction: the neighbour ranks see the restricted values
void sflow_amr::exchange_level0(lexer *p, fdm2D *b, ghostcell *pgc, int s)
{
    if(p->mpi_size<=1)
    return;

    slice *WLi,*UHi,*VHi,*WLo,*UHo,*VHo;
    pmom0->stage_io(s,b,WLi,UHi,VHi,WLo,UHo,VHo);

    slice *WHi,*WHo;
    pmom0->stage_io_w(s,b,WHi,WHo);

    slice *f[13] = {WLo,UHo,VHo,WHo,&b->eta,&b->U,&b->V,&b->W,&b->hp,&b->UA,&b->VA,&b->MX,&b->MY};
    for(int k=0; k<(bous==1?13:9); ++k)
    {
        pgc->gcslparax(p,*f[k],4);
        pgc->gcslparacox(p,*f[k],10);
    }

    pgc->gcslparaxV_int(p,p->wet,4);
    pgc->gcslparacoxV_int(p,p->wet,50);
}

// fine face values of the level-l patches on the partition edges, to the rank of the coarse cell
void sflow_amr::exchange_fluxes(int l)
{
    face_run(l,5,7200+l,
             [&](reefamr_patch *q, int side, int r, double *sb)
             {
                 sflow_amr_patch *c = SP(q);
                 for(int ip=0; ip<5; ++ip)
                 sb[ip] = 0.5*(c->rec[ip][side][r]+c->rec[ip][side][r+1]);
             },
             [&](reefamr_match &m, const double *v)
             {
                 for(int ip=0; ip<5; ++ip)
                 m.val[ip] = v[ip];
             });
}

// --------------------------------------------------------------------- stages
void sflow_amr::step_begin(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // G 7 1: level 0 steps alone, the patches follow in step_end
    if(sub==1)
    {
        if(maxlev>=1 && patches_total>0)
        sub_begin(p,b,pgc);
        return;
    }

    comms_off guard(pgc);

    for(auto q : P)
    {
        sflow_amr_patch *c = SP(q);
        c->pp->dt = p->dt;
        c->pp->dt_old = p->dt_old;
        c->pp->simtime = p->simtime;
        c->pp->count = p->count;
        c->pmom->inflow(c->pp,c->b,pgc,pflow_void);
    }

    // body on the patches: X 10 3 updates the ship pressure once per step (as level 0),
    // X 10 2 moves the body in every stage
    ship_patches(shipmode==3);
}

void sflow_amr::stage_begin(lexer *p, fdm2D *b, ghostcell *pgc, int s)
{
    if(maxlev<1 || patches_total==0)
    return;

    // G 7 1: level 0 alone, its coarse face fluxes go into the flux registers (hll_hook)
    if(sub==1)
    {
        hstage = s;
        return;
    }

    double t0 = MPI_Wtime();
    cache_stage(s);

    for(int l=1; l<=maxlev; ++l)
    fill_level(pgc,l,s);

    if(shipmode==2 && s>0)
    ship_patches(false);
    tm[0] += MPI_Wtime()-t0;

    // finest first: the coarser grids take the recorded fine fluxes
    for(int l=maxlev; l>=1; --l)
    {
        t0 = MPI_Wtime();
        {
        comms_off guard(pgc);
        for(int n : lev[l])
        SP(n)->pmom->rk_stage(SP(n)->pp,SP(n)->b,pgc,s);
        }
        tm[1] += MPI_Wtime()-t0;

        t0 = MPI_Wtime();
        exchange_fluxes(l);
        tm[2] += MPI_Wtime()-t0;
    }
}

void sflow_amr::stage_end(lexer *p, fdm2D *b, ghostcell *pgc, int s)
{
    if(maxlev<1 || patches_total==0 || sub==1)
    return;

    double t0 = MPI_Wtime();
    restrict_levels(p,s);

    exchange_level0(p,b,pgc,s);
    tm[3] += MPI_Wtime()-t0;
}

void sflow_amr::step_end(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    // G 7 1: the patches take their steps now (two per step of the parent), refluxing, restriction
    if(sub==1)
    {
        if(patches_total>0)
        sub_end(p,b,pgc);
    }
    else
    // end of the step as for level 0 in sflow_f::mainloop: breaking flags cleared (A 248 0),
    // water depth and wet-dry state updated
    {
    comms_off guard(pgc);
    for(auto q : P)
    {
        sflow_amr_patch *c = SP(q);
        lexer *pp = c->pp;
        c->pmom->rk_finish(pp,c->b,pgc);

        if(p->A248==0)
        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        c->b->breaking(ii,jj)=0;

        c->pfsf->depth_update(pp,c->b,pgc,c->b->WL);
    }
    }

    // G 7 1: G 2 counts steps of the finest level, as without subcycling (a level-0 step is 2^G1 of them)
    const int rint = (sub==1) ? MAX(regrid_int>>maxlev,1) : regrid_int;

    if(regrid_int>0 && p->count%rint==0)
    {
    double t0 = MPI_Wtime();
    regrid(p,pgc,false);
    tm[4] += MPI_Wtime()-t0;
    }
}

// --------------------------------------------------------------------- flux matching
void sflow_amr::hll_hook(lexer *p, fdm2D *b, int ipol, int id)
{
    slice &Fx = (ipol==4) ? b->FEx : b->Fx;
    slice &Fy = (ipol==4) ? b->FEy : b->Fy;

    // a patch records its boundary faces
    if(id>0)
    {
        sflow_amr_patch &c = *SP(id-1);
        const int il = EXT-1, ih = EXT+c.nx-1;
        const int jl = EXT-1, jh = EXT+c.ny-1;

        for(int r=0; r<c.ny; ++r)
        {
            c.rec[ipol][0][r] = Fx(il,EXT+r);
            c.rec[ipol][1][r] = Fx(ih,EXT+r);
        }
        for(int r=0; r<c.nx; ++r)
        {
            c.rec[ipol][2][r] = Fy(EXT+r,jl);
            c.rec[ipol][3][r] = Fy(EXT+r,jh);
        }

        if(ipol==4)
        {
            for(int r=0; r<c.ny; ++r)
            {
                c.rec[0][0][r] = b->dfx(il,EXT+r);
                c.rec[0][1][r] = b->dfx(ih,EXT+r);
            }
            for(int r=0; r<c.nx; ++r)
            {
                c.rec[0][2][r] = b->dfy(EXT+r,jl);
                c.rec[0][3][r] = b->dfy(EXT+r,jh);
            }
        }

        // G 7 1: fine register, the face fluxes of this stage with its weight in the step
        if(sub==1)
        {
            const double wdt = (hstage==2 ? 2.0/3.0 : 1.0/6.0)*p->dt;
            for(int side=0; side<4; ++side)
            for(size_t r=0; r<c.rec[ipol][side].size(); ++r)
            c.freg[ipol][side][r] += wdt*c.rec[ipol][side][r];
        }
    }

    if(id>=(int)match.size())
    return;

    // G 7 1: the coarse faces keep their own flux (the fine one is not there yet, refluxing at the
    // end of the fine steps) and add it to the coarse register; the face depth is the fine one
    // as in the synchronous coupling (well balanced: the same face depth in the flux and in the
    // bed slope term of the coarse cell, the face depths of the bed are static)
    if(sub==1)
    {
        const double wdt = (hstage==2 ? 2.0/3.0 : 1.0/6.0)*p->dt;

        for(size_t k=0; k<match[id].size(); ++k)
        {
            reefamr_match &m = match[id][k];
            if(ipol==4)
            {
                sflow_amr_patch &c = *SP(m.child);
                double df = 0.5*(c.rec[0][m.side][m.r]+c.rec[0][m.side][m.r+1]);
                if(m.dir==0)
                b->dfx(m.fi,m.fj) = df;
                else
                b->dfy(m.fi,m.fj) = df;
            }
            cregL[id][5*k+ipol] += wdt*(m.dir==0 ? Fx(m.fi,m.fj) : Fy(m.fi,m.fj));
        }

        for(size_t k=0; k<rmatch[id].size(); ++k)
        {
            reefamr_match &m = rmatch[id][k];
            if(ipol==4)
            {
                if(m.dir==0)
                b->dfx(m.fi,m.fj) = m.val[0];
                else
                b->dfy(m.fi,m.fj) = m.val[0];
            }
            cregR[id][5*k+ipol] += wdt*(m.dir==0 ? Fx(m.fi,m.fj) : Fy(m.fi,m.fj));
        }
        return;
    }

    // coarse faces next to the patches of this rank: mean of the two fine faces
    for(auto &m : match[id])
    {
        sflow_amr_patch &c = *SP(m.child);
        double val = 0.5*(c.rec[ipol][m.side][m.r]+c.rec[ipol][m.side][m.r+1]);

        if(m.dir==0)
        Fx(m.fi,m.fj) = val;
        else
        Fy(m.fi,m.fj) = val;

        if(ipol==4)
        {
            double df = 0.5*(c.rec[0][m.side][m.r]+c.rec[0][m.side][m.r+1]);
            if(m.dir==0)
            b->dfx(m.fi,m.fj) = df;
            else
            b->dfy(m.fi,m.fj) = df;
        }
    }

    // coarse faces next to patches of other ranks
    for(auto &m : rmatch[id])
    {
        if(m.dir==0)
        Fx(m.fi,m.fj) = m.val[ipol];
        else
        Fy(m.fi,m.fj) = m.val[ipol];

        if(ipol==4)
        {
            if(m.dir==0)
            b->dfx(m.fi,m.fj) = m.val[0];
            else
            b->dfy(m.fi,m.fj) = m.val[0];
        }
    }
}

// --------------------------------------------------------------------- time step
void sflow_amr::timestep(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    // initial time step (sflow_etimestep::ini): linear in DXM, so scale with the finest patch
    if(p->count==0)
    {
        // G 7 1: every level at its own CFL number, the step of level 0 as without patches
        if(sub==1)
        return;

        double r=1.0;
        for(auto c : P)
        r = MIN(r, c->pp->DXM/p->DXM);

        r = pgc->globalmin(r);
        p->dt *= r;
        p->dt_old = p->dt;
        return;
    }

    // same CFL as sflow_etimestep (unsplit 2D: x and y Courant numbers add up), over the wet real patch cells
    const double g = fabs(p->W22);
    double cmin = 1.0e20;
    vector<double> cl(maxlev+1,1.0e20);     // G 7 1: per level

    for(auto q : P)
    {
        sflow_amr_patch *c = SP(q);
        lexer *pp = c->pp;
        fdm2D *pb = c->b;
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
        if(pp->wet[lij(pp,ii,jj)]==1 && pp->flagslice4[lij(pp,ii,jj)]>0)
        {
            double cc = sqrt(g*MAX(pb->WL(ii,jj),p->A244));
            double sigma = (fabs(pb->U(ii,jj))+cc)/pp->DXN[ii+marge];

            if(p->j_dir==1)
            sigma += (fabs(pb->V(ii,jj))+cc)/pp->DYN[jj+marge];

            double cv = 1.0/sigma;

            if(p->A219==2)
            cv = MIN(cv, pp->DXN[ii+marge]/(fabs(pb->U(ii,jj))>1.0e-20?fabs(pb->U(ii,jj)):1.0e-20));

            cmin = MIN(cmin, cv);
            cl[c->lev] = MIN(cl[c->lev], cv);
        }
    }

    // explicit dispersion correction of A 220 3, as sflow_etimestep (G 7 1: per level)
    double dtd = 1.0e20;
    vector<double> dl(maxlev+1,1.0e20);
    if(p->A220==3 && p->A224>1.0)
    {
        const double B = (p->A224-1.0)/3.0;
        for(auto q : P)
        {
            sflow_amr_patch *c = SP(q);
            lexer *pp = c->pp;
            fdm2D *pb = c->b;
            for(int ii=EXT; ii<EXT+c->nx; ++ii)
            for(int jj=EXT; jj<EXT+c->ny; ++jj)
            if(pp->wet[lij(pp,ii,jj)]==1 && pp->flagslice4[lij(pp,ii,jj)]>0)
            {
                double hh = MAX(pb->WL(ii,jj),p->A244);
                double dx = pp->DXN[ii+marge], dy = pp->DYN[jj+marge];
                const double dv = 1.2/((1.0/(dx*dx) + p->y_dir/(dy*dy))*sqrt(g*B*hh*hh*hh));
                dtd = MIN(dtd, dv);
                dl[c->lev] = MIN(dl[c->lev], dv);
            }
        }
    }

    // Boussinesq (A 220 4): Courant number 0.5 at most, as sflow_etimestep
    const double cfl = (p->A220==4) ? MIN(p->N47,0.25) : p->N47;

    double dtp = MIN(cfl*2.0*cmin, dtd);

    // G 7 1: level l takes 2^l steps per step of level 0
    if(sub==1)
    {
        dtp = 1.0e20;
        for(int l=1; l<=maxlev; ++l)
        dtp = MIN(dtp, MIN(cfl*2.0*cl[l], dl[l])*double(1<<l));
    }

    dtp = pgc->globalmin(dtp);

    if(p->N48==1)
    p->dt = MIN(p->dt,dtp);
}

// --------------------------------------------------------------------- output
double sflow_amr::mass(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double m=0.0;
    SLICELOOP4
    m += b->WL(i,j)*p->DXN[IP]*p->DYN[JP];

    return pgc->globalsum(m);
}

void sflow_amr::print(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    bool doprint=false;

    if((p->count%p->P181==0 && p->P182<0.0 && p->P10==1 && p->P181>0) || (p->count==0 && p->P182<0.0))
    doprint=true;

    if((p->simtime>printtime_amr && p->P182>0.0 && p->P10==1) || (p->count==0 && p->P182>0.0))
    {
        doprint=true;
        printtime_amr += p->P182;
    }

    double m = mass(p,b,pgc);

    if(p->P51>0)
    gauges(p,b,pgc);

    if(p->mpirank==0 && (p->count%p->P12==0 || doprint))
    logout<<p->count<<" \t "<<setprecision(10)<<p->simtime<<" \t "<<p->dt<<" \t "<<patches_total<<" \t "<<cells_total<<" \t "
          <<setprecision(15)<<m<<" \t "<<setprecision(6)<<(m-m0)/(fabs(m0)>1.0e-20?m0:1.0)<<endl;

    if(p->mpirank==0 && doprint)
    cout<<"SFLOW AMR: "<<patches_total<<" patches, "<<cells_total<<" cells; time in patch stages "<<setprecision(4)<<tm[1]<<" s, fill "<<tm[0]<<" s, restriction "<<tm[3]<<" s, regrid "<<tm[4]<<" s";
    if(p->mpirank==0 && doprint && shipmode>0)
    cout<<", body "<<tm[8]<<" s";
    if(p->mpirank==0 && doprint && bous==1)
    cout<<"; u_a "<<tm[9]<<" s, mean iterations "<<(bq_solves>0 ? double(bq_it_total)/bq_solves : 0.0);
    if(p->mpirank==0 && doprint && nh==1)
    cout<<"; pressure "<<tm[5]<<" s (preconditioner "<<tm[6]<<" s, operator "<<tm[7]<<" s), mean iterations "<<(nh_solves>0 ? double(nh_it_total)/nh_solves : 0.0);
    if(p->mpirank==0 && doprint && nh==1 && sub==1)
    cout<<" (level solves "<<(sub_lv_n>0 ? double(sub_lv_it)/sub_lv_n : 0.0)<<", synchronisation "<<(sub_sy_n>0 ? double(sub_sy_it)/sub_sy_n : 0.0)<<")";
    if(p->mpirank==0 && doprint)
    cout<<endl;

    if(!doprint)
    return;

    // P 18: the grids in the LAGOON store; P 18 1: instead of the .vtr and .vtm files
    bool stored = false;
    if(p->P18>0)
    stored = print_lagoon(p,b,pgc);

    if(!lagoon_amr_output::files_needed(p,stored))
    {
        ++printcount_amr;
        return;
    }

    write_vtr0(p,b);

    for(int n=0; n<(int)P.size(); ++n)
    write_vtr(p,*SP(n),n);

    // multiblock index: level 0 of every rank + all patches
    int np = (int)P.size();
    vector<int> all(p->mpi_size,0);
    MPI_Allgather(&np,1,MPI_INT,&all[0],1,MPI_INT,MPI_COMM_WORLD);

    if(p->mpirank==0)
    {
        char name[256];
        snprintf(name,sizeof(name),"./REEF3D_SFLOW_AMR/REEF3D-SFLOW-AMR-%08i.vtm",printcount_amr);
        ofstream out(name);
        out<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">\n<vtkMultiBlockDataSet>\n";
        out<<"<Block index=\"0\" name=\"level 0\">\n";
        for(int r=0; r<p->mpi_size; ++r)
        out<<"<DataSet index=\""<<r<<"\" file=\"REEF3D-SFLOW-AMR-L0-"<<setw(8)<<setfill('0')<<printcount_amr<<"-"<<setw(4)<<r+1<<".vtr\"/>\n";
        out<<setfill(' ')<<"</Block>\n<Block index=\"1\" name=\"patches\">\n";
        int idx=0;
        for(int r=0; r<p->mpi_size; ++r)
        for(int q=0; q<all[r]; ++q)
        {
        out<<"<DataSet index=\""<<idx<<"\" file=\"REEF3D-SFLOW-AMR-"<<setw(8)<<setfill('0')<<printcount_amr<<"-"<<setw(4)<<r+1<<"-"<<setw(4)<<q+1<<".vtr\"/>\n";
        out<<setfill(' ');
        ++idx;
        }
        out<<"</Block>\n</vtkMultiBlockDataSet>\n</VTKFile>\n";
        out.close();
    }

    ++printcount_amr;
}

// surface elevation at the P 51 gauges from the finest grid that holds them
void sflow_amr::gauges(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const int ng = p->P51;
    vector<double> lv(ng,-1.0), val(ng,0.0);

    auto find = [&](lexer *q, int i0, int i1, int j0, int j1, double x, double y, int &ii, int &jj)
    {
        ii=-1; jj=-1;
        for(int a=i0; a<=i1; ++a)
        if(x>=q->XN[a+marge] && x<q->XN[a+1+marge]) { ii=a; break; }
        for(int a=j0; a<=j1; ++a)
        if(y>=q->YN[a+marge] && y<q->YN[a+1+marge]) { jj=a; break; }
        return ii>=0 && jj>=0;
    };

    for(int k=0; k<ng; ++k)
    {
        int ii,jj;
        if(find(p,0,NX0-1,0,NY0-1,p->P51_x[k],p->P51_y[k],ii,jj))
        {
            lv[k]=0.0;
            val[k]=b->eta(ii,jj);
        }
        for(auto c : P)
        if(c->lev>lv[k] && find(c->pp,EXT,EXT+c->nx-1,EXT,EXT+c->ny-1,p->P51_x[k],p->P51_y[k],ii,jj))
        {
            lv[k]=c->lev;
            val[k]=SP(c)->b->eta(ii,jj);
        }
    }

    vector<double> lmax(ng);
    MPI_Allreduce(&lv[0],&lmax[0],ng,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    for(int k=0; k<ng; ++k)
    if(lv[k]!=lmax[k])
    val[k]=0.0;
    vector<double> vs(ng);
    MPI_Allreduce(&val[0],&vs[0],ng,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    if(p->mpirank==0)
    {
        if(!gaugeout.is_open())
        {
            gaugeout.open("./REEF3D_SFLOW_AMR/REEF3D_SFLOW_AMR_gauges.dat");
            gaugeout<<"# simtime";
            for(int k=0; k<ng; ++k)
            gaugeout<<" \t eta("<<p->P51_x[k]<<","<<p->P51_y[k]<<")";
            gaugeout<<" \t levels"<<endl;
        }
        gaugeout<<setprecision(10)<<p->simtime;
        for(int k=0; k<ng; ++k)
        gaugeout<<" \t "<<setprecision(10)<<vs[k];
        for(int k=0; k<ng; ++k)
        gaugeout<<" \t "<<int(lmax[k]);
        gaugeout<<endl;
    }
}

// level 0 of this rank, cell centred, same layout as the patches
void sflow_amr::write_vtr0(lexer *p, fdm2D *b)
{
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_SFLOW_AMR/REEF3D-SFLOW-AMR-L0-%08i-%04i.vtr",printcount_amr,p->mpirank+1);

    const int nx=p->knox, ny=p->knoy, m=marge;
    ofstream out(name);
    out<<"<?xml version=\"1.0\"?>\n";
    out<<"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out<<"<RectilinearGrid WholeExtent=\"0 "<<nx<<" 0 "<<ny<<" 0 0\">\n";
    out<<"<FieldData><DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\">"<<p->simtime<<"</DataArray>";
    out<<"<DataArray type=\"Int32\" Name=\"level\" NumberOfTuples=\"1\">0</DataArray></FieldData>\n";
    out<<"<Piece Extent=\"0 "<<nx<<" 0 "<<ny<<" 0 0\">\n<CellData Scalars=\"eta\">\n";

    auto field = [&](const char *nm, slice &f, double shift)
    {
        out<<"<DataArray type=\"Float64\" Name=\""<<nm<<"\" format=\"ascii\">\n";
        for(int jj=0; jj<ny; ++jj)
        {
            for(int ii=0; ii<nx; ++ii)
            out<<setprecision(17)<<f(ii,jj)+shift<<" ";
            out<<"\n";
        }
        out<<"</DataArray>\n";
    };
    field("eta",b->eta,0.0);
    field("elevation",b->eta,p->wd);
    field("WL",b->WL,0.0);
    field("u",b->U,0.0);
    field("v",b->V,0.0);
    field("bed",b->bed,0.0);
    field("press",b->press,0.0);
    field("w",b->W,0.0);

    out<<"<DataArray type=\"Int32\" Name=\"wetdry\" format=\"ascii\">\n";
    for(int jj=0; jj<ny; ++jj)
    {
        for(int ii=0; ii<nx; ++ii)
        out<<p->wet[lij(p,ii,jj)]<<" ";
        out<<"\n";
    }
    out<<"</DataArray>\n</CellData>\n<Coordinates>\n";
    out<<"<DataArray type=\"Float64\" Name=\"x\" format=\"ascii\">";
    for(int ii=0; ii<=nx; ++ii) out<<setprecision(12)<<p->XN[ii+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"y\" format=\"ascii\">";
    for(int jj=0; jj<=ny; ++jj) out<<setprecision(12)<<p->YN[jj+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"z\" format=\"ascii\">0</DataArray>\n";
    out<<"</Coordinates>\n</Piece>\n</RectilinearGrid>\n</VTKFile>\n";
    out.close();
}

void sflow_amr::write_vtr(lexer *p, sflow_amr_patch &c, int id)
{
    lexer *pp = c.pp;
    fdm2D *pb = c.b;
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_SFLOW_AMR/REEF3D-SFLOW-AMR-%08i-%04i-%04i.vtr",printcount_amr,p->mpirank+1,id+1);

    ofstream out(name);
    out<<"<?xml version=\"1.0\"?>\n";
    out<<"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out<<"<RectilinearGrid WholeExtent=\"0 "<<c.nx<<" 0 "<<c.ny<<" 0 0\">\n";
    out<<"<FieldData><DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\">"<<p->simtime<<"</DataArray>";
    out<<"<DataArray type=\"Int32\" Name=\"level\" NumberOfTuples=\"1\">"<<c.lev<<"</DataArray></FieldData>\n";
    out<<"<Piece Extent=\"0 "<<c.nx<<" 0 "<<c.ny<<" 0 0\">\n";
    out<<"<CellData Scalars=\"eta\">\n";

    auto field = [&](const char *nm, slice &f, double shift)
    {
        out<<"<DataArray type=\"Float64\" Name=\""<<nm<<"\" format=\"ascii\">\n";
        for(int jj=EXT; jj<EXT+c.ny; ++jj)
        {
            for(int ii=EXT; ii<EXT+c.nx; ++ii)
            out<<setprecision(17)<<f(ii,jj)+shift<<" ";
            out<<"\n";
        }
        out<<"</DataArray>\n";
    };
    field("eta",pb->eta,0.0);
    field("elevation",pb->eta,p->wd);
    field("WL",pb->WL,0.0);
    field("u",pb->U,0.0);
    field("v",pb->V,0.0);
    field("bed",pb->bed,0.0);
    field("press",pb->press,0.0);
    field("w",pb->W,0.0);

    out<<"<DataArray type=\"Int32\" Name=\"wetdry\" format=\"ascii\">\n";
    for(int jj=EXT; jj<EXT+c.ny; ++jj)
    {
        for(int ii=EXT; ii<EXT+c.nx; ++ii)
        out<<pp->wet[lij(pp,ii,jj)]<<" ";
        out<<"\n";
    }
    out<<"</DataArray>\n";
    out<<"</CellData>\n<Coordinates>\n";

    const int m = marge;
    out<<"<DataArray type=\"Float64\" Name=\"x\" format=\"ascii\">";
    for(int ii=EXT; ii<=EXT+c.nx; ++ii) out<<setprecision(12)<<pp->XN[ii+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"y\" format=\"ascii\">";
    for(int jj=EXT; jj<=EXT+c.ny; ++jj) out<<setprecision(12)<<pp->YN[jj+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"z\" format=\"ascii\">0</DataArray>\n";
    out<<"</Coordinates>\n</Piece>\n</RectilinearGrid>\n</VTKFile>\n";
    out.close();
}

// P 18: this rank's grids as its .vtr files have them (the fields of write_vtr0 and
// write_vtr, wetdry as int32), gathered into the LAGOON store; true when the output is there
bool sflow_amr::print_lagoon(lexer *p, fdm2D *b, ghostcell *pgc)
{
    static lagoon_amr_output *writer = nullptr;
    if(writer==nullptr)
    writer = new lagoon_amr_output("SFLOW", {{"eta",false}, {"elevation",false}, {"WL",false}, {"u",false}, {"v",false},
                                             {"bed",false}, {"press",false}, {"w",false}, {"wetdry",true}});

    const int m = marge;
    vector<lagoon_amr::grid> grids;
    auto take = [&](lexer *q, fdm2D *f, int i0, int j0, int nx, int ny, int level)
    {
        lagoon_amr::grid g;
        g.level = level;
        g.nx = nx;
        g.ny = ny;
        for(int ii=i0; ii<=i0+nx; ++ii) g.x.push_back(q->XN[ii+m]);
        for(int jj=j0; jj<=j0+ny; ++jj) g.y.push_back(q->YN[jj+m]);
        auto field = [&](slice &s, double shift)
        {
            for(int jj=j0; jj<j0+ny; ++jj)
            for(int ii=i0; ii<i0+nx; ++ii)
            g.values.push_back(s(ii,jj)+shift);
        };
        field(f->eta,0.0);
        field(f->eta,p->wd);
        field(f->WL,0.0);
        field(f->U,0.0);
        field(f->V,0.0);
        field(f->bed,0.0);
        field(f->press,0.0);
        field(f->W,0.0);
        for(int jj=j0; jj<j0+ny; ++jj)
        for(int ii=i0; ii<i0+nx; ++ii)
        g.values.push_back(double(q->wet[lij(q,ii,jj)]));
        grids.push_back(std::move(g));
    };
    take(p,b,0,0,p->knox,p->knoy,0);
    for(int n=0; n<(int)P.size(); ++n)
    take(SP(n)->pp,SP(n)->b,EXT,EXT,SP(n)->nx,SP(n)->ny,SP(n)->lev);

    return writer->write(p,pgc,grids,printcount_amr);
}
