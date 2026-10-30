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

#include"sflow_curv.h"
#include"reef3d_mpmd.h"
#include"sflow_amr.h"
#include"gcsl_yskip.h"
#include"reefckpt.h"
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
#include"sflow_momentum_SI.h"
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
#include"ioflow_amr_block.h"
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

sflow_amr::sflow_amr(lexer *p, fdm2D *b, ghostcell *pgc, patchBC_interface *ppBC, sixdof *pp6dof, ioflow *ppflow,
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
    // second difference of the surface (G 21), breaking (G 23 1 with A 246 1/2): as fnpf_amr
    tol_curv = MAX(p->G21,0.0);
    flagbreak = (p->G23==1 && p->A246>=1);
    q.tile = MAX(p->G4,4);
    q.tile += q.tile%2;
    q.nest = 3;
    q.keep = MIN(MAX(p->G5,0),250);
    // the coarse faces next to a finer patch also on the grids that compute them around their box: the
    // wet-dry update of a stage reads the new state of the cells around the box, which otherwise depended
    // on which grid holds the face (another partition, another layout of the patches)
    q.ring = true;

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
    NV = bous==1 ? 16 : 12;
    ISF = NV-2;
    ISC = NV-1;

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
    if(tol_eta<=0.0 && tol_curv<=0.0 && !flagbreak && shore==0 && (p->G12==0 || shipmode==0))
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

    // no refinement in the relaxation zones of the wave generation and the numerical beach (B 98 2,
    // B 99 1/2, B 107, B 108, B 520): the patches run without the in- and outflow (pflow_of), so their
    // cells are not relaxed, and the restriction would replace the relaxed level-0 cells by them
    if(ppflow!=nullptr)
    {
        const size_t n0 = q.fbox.size();
        ppflow->relax_boxes(p,q.fbox);
        nrelax = int((q.fbox.size()-n0)/4);
    }

    // refinement around the moving body (G 12 margin, G 13 wake wedge)
    q.zones = (shipmode>0 && p->G12>0);
    q.zr = p->G12_r;
    q.zL = p->G13_L;
    q.za = p->G13_a;

    // several ranks (G 40 1): the patches are placed for the load of the ranks; parents and old
    // patches on other ranks are reached through the block plans and old_run of the core.  With the
    // moving body (X 10 2/3) the level-0 level set below the patches comes from the ranks that hold it
    // (ship_fs_exchange); the draft and the pressure are evaluated on the patch itself
    q.place = (p->G40>=1 && p->mpi_size>1) ? MIN(p->G40,2) : 0;
    q.rebalance = 0.1;

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
    delete sfs0;
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
    r.oi = c->I0-EXT; r.oj = c->J0-EY;
    return r;
}

// --------------------------------------------------------------------- setup
void sflow_amr::ini(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(maxlev<1 && p->geo_mb==0)
    return;

    // scope
    int ok=1;
    if(p->A220<0 || p->A220>4) ok=0;
    if(p->A210==2 && si==0) ok=0;
    // semi-implicit free surface (A 226 1): hydrostatic, one global time step, no moving bodies (sflow_amr_si.cpp)
    if(si==1 && (p->A220!=0 || p->G7==1 || p->X10!=0)) ok=0;
    if(p->A260!=0 && p->A260!=3 && p->A260!=4) ok=0;   // the algebraic eddy viscosity of every grid (sflow_ediff)
    // sediment transport (S 10 1): bed load, Exner (S 31 1/2, S 32 1/2), sand slide, no fractions (sflow_amr_sed.cpp)
    if(p->S10!=0 && (p->S10!=1 || p->S51>0 || p->S31==3 || p->S32>=3)) ok=0;
    if(p->X10>3) ok=0;
    if(p->W90!=0) ok=0;
    if(p->j_dir!=1 && (p->geo_mb>0 || p->geo_curv==1)) ok=0;     // 2D grid (one cell in y): Cartesian

    // multi-block grid (G 50 1, DIVEMesh R 20 1): hydrostatic, without subcycling; moving bodies X 10 2/3 (P5c: the
    // body on the blocks 1.. as on the patches, sflow_amr_ship.cpp)
    if(p->geo_mb>0)
    {
        if(p->A220!=0) ok=0;
        if(p->G7==1) ok=0;
        if(p->X10!=0 && p->X10!=2 && p->X10!=3) ok=0;
        if(ok==0)
        {
            if(p->mpirank==0)
            cout<<"!!! SFLOW multi-block grid: only for A 220 0, A 210 3 (A 226 1: A 210 2/3), A 260 0/3/4, S 10 0/1 (S 51 0, S 31 1/2, S 32 1/2), X 10 0/2/3 (A 226 1: X 10 0), W 90 0, without G 7 1 !!!"<<endl;
            exit(1);
        }
    }

    if(ok==0)
    {
        if(p->mpirank==0)
        cout<<"SFLOW AMR (G 1): only for A 220 0-4, A 210 3 (A 226 1: A 210 2/3, A 220 0, X 10 0, without G 7 1), A 260 0/3/4, S 10 0/1 (S 51 0, S 31 1/2, S 32 1/2), X 10 0-3, W 90 0, and with one cell in y Cartesian -- refinement switched off"<<endl;
        maxlev=0;
        return;
    }

    // subcycling (G 7 1, sflow_amr_sub.cpp): hydrostatic, non-hydrostatic (A 220 1-3: level solves
    // and synchronisation projections) and Boussinesq (A 220 4: u_a per level), sflow_amr_subnh.cpp;
    // the moving body (X 10 2/3) on the patches linear in time within the level-0 step
    // (sflow_amr_ship.cpp)
    sub = (p->G7==1) ? 1 : 0;

    // hierarchy: rank boxes, tiles, level-0 flags, no-refinement cells
    setup(p,pgc);

    // bed and still water depth in the level-0 halo (the fine bed is prolonged from it)
    pgc->gcslparax(p,b->bed,4);
    pgc->gcslparax(p,b->depth,4);

    // shoreline: restricted water surface of the covered level-0 cells
    sfs0 = new slice4(p);
    std::fill(sfs0->V,sfs0->V+(size_t)p->imax*p->jmax,NAN);

    // multi-block grid: blocks 1.. as level-0 patches, their initial state, the fill plans
    if(p->geo_mb>0)
    {
        blocks_setup(p,pgc);
        blocks_ini(p,pgc);
        blocks_plans(pgc);
        cache_stage(0);
        fill_level0(pgc,0);

        // S 10: the sediment of the blocks starts from their bed
        if(p->S10==1)
        for(int id : lev[0])
        sed_patch_ini(pgc,*SP(id));
    }

    mkdir("./REEF3D_SFLOW_AMR",0777);

    if(p->mpirank==0)
    {
        if(reefckpt_open(logout,"./REEF3D_SFLOW_AMR/REEF3D_SFLOW_AMR_log.dat",p->I45))   // I 45: continued
        logout<<"# count \t simtime \t dt \t patches \t cells \t mass \t rel. mass change"<<endl;
    }

    // initial hierarchy: every pass can add one level.  The initial state of the fine cells: as the level-0
    // initial state is a water surface (F 60, relaxation zones), a fine cell below the surface of the wet
    // neighbours of its dry parent is wet (dryfill, regrid_state), so a lake at rest starts at rest at a
    // shoreline the coarser grid does not resolve
    ini_pass = true;
    for(int it=0; it<maxlev; ++it)
    regrid(p,pgc,true);
    ini_pass = false;

    if(p->geo_mb>0 && maxlev>=1)
    {
        cache_stage(0);
        fill_level0(pgc,0);
    }

    m0 = mass(p,b,pgc);

    if(p->mpirank==0 && maxlev>=1)
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
        if(nrelax>0)
        cout<<"SFLOW AMR: no refinement in "<<nrelax<<" relaxation zone(s) of the wave generation and the numerical beach"<<endl;
        if(tol_eta>0.0 || tol_curv>0.0 || shore>0 || flagbreak)
        {
            cout<<"SFLOW AMR: flags:";
            if(tol_eta>0.0)
            cout<<" surface jump > "<<tol_eta<<" m (G 20)";
            if(tol_curv>0.0)
            cout<<" second difference > "<<tol_curv<<" m (G 21)";
            if(shore>0)
            cout<<" shoreline (G 22)";
            if(flagbreak)
            cout<<" breaking (G 23)";
            cout<<endl;
        }
        if(p->G23==1 && !flagbreak)
        cout<<"SFLOW AMR: G 23 needs the breaking detection A 246 1/2 -- ignored"<<endl;
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

    // S 10: the sediment of the patch (sflow_amr_sed.cpp), its state set by regrid_state
    if(p->S10==1)
    {
    pp->gcin_count = pp->gcout_count = 0;   // no 3D in- and outflow lists on a patch (sediment_f::ini_guard)
    c->sed = new sediment_f(pp,pgc,nullptr,pBC,1);
    }
    c->phll = new sflow_HLL(pp,pgc,pBC);
    c->phll->amr = this;
    c->pss = new sflow_signal_speed(pp);

    if(p->A211<=3)
    c->precon = new sflow_reconstruct_hires(pp,pBC);
    if(p->A211>=4)
    c->precon = new sflow_reconstruct_weno(pp,pBC,1);

    // diffusion: explicit on the patches (A 212 2 is implicit on level 0 only)
    if(p->A212>=1)
    {
    sflow_ediff *ed = new sflow_ediff(pp);
    ed->own_eddyv = 1;
    c->pdiff = ed;
    }
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

    // multi-block: the open faces of a block 1.. (G 54) and its edge zones (B 522; the lists from build_bc2D); the
    // iowave edges of level 0 (open_xm .. open_yp) are not the block's
    const int zin = (c->blk>0 && c->lev==0) ? bzone[0][c->blk] : -1;
    const int zout = (c->blk>0 && c->lev==0) ? bzone[1][c->blk] : -1;
    // relaxation and beach zones (B 522 methods 1, 2)
    const int rin = (c->blk>0 && c->lev==0) ? bzr[0][c->blk] : -1;
    const int rout = (c->blk>0 && c->lev==0) ? bzr[1][c->blk] : -1;
    if(c->blk>0)
    {
        pp->open_xm = (zin>=0) ? 1 : 0;
        pp->open_xp = (zout>=0) ? 1 : 0;
        pp->open_ym = pp->open_yp = 0;
    }
    if(c->blk>0 && c->lev==0 && (bopen[0][c->blk]>=0 || bopen[1][c->blk]>=0 || zin>=0 || zout>=0 || rin>=0 || rout>=0))
    {
        c->pio = new ioflow_amr_block(pp,pgc,pBC);
        pp->amr_own_io = 1;
        c->pio->open_faces(pp,bopen[0][c->blk],bopen[1][c->blk]);
        c->pio->open_edges(pp,zin>=0 ? &bzones[zin] : nullptr,zout>=0 ? &bzones[zout] : nullptr,bzbg);
        // the zones (the faces closed: walls, as the native grid with B 98 2); the waves of the relaxation zones with
        // the B 92 waves (level 0 has no wave generation on a curvilinear grid, sflow_curv::check)
        c->pio->open_zones(pp,rin>=0 ? &bzones[rin] : nullptr,rout>=0 ? &bzones[rout] : nullptr,bzbg,p0->B92>0,
                           EXT,blk[c->blk].ni,EXT,blk[c->blk].nj);
        if(bopen[0][c->blk]>=0 || bopen[1][c->blk]>=0 || zin>=0 || zout>=0)
        c->phll->open_io();
    }

    c->pmom = new sflow_momentum_RK3(pp,c->b,pgc,c->phll,c->pss,c->precon,c->pdiff,c->ppress,c->psolv,nullptr,pflow_of(c),c->pfsf,c->psfdf,p6);
    c->pmom->nh_defer = (nh==1);
    c->pmom->bous_defer = (bous==1);

    // semi-implicit free surface (A 226 1): the step of the grid (sflow_amr_si.cpp); the RK3 object of the patch serves
    // the fills and the restriction (its stage arrays of stage 0 in and stage 2 out are the fields of the grid). A
    // refined patch has no open boundaries (no refinement next to them)
    if(si==1)
    {
        if(c->lev>=1)
        pp->open_xm = pp->open_xp = pp->open_ym = pp->open_yp = 0;
        c->psim = new sflow_momentum_SI(pp,c->b,pgc,c->phll,c->pss,c->precon,c->pdiff,c->ppress,c->psolv,nullptr,pflow_of(c),c->pfsf,c->psfdf,p6,true);
        c->psim->amrgrid = true;
    }

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

    c->sfc.assign((size_t)pp->imax*pp->jmax,NAN);
    c->gsf.assign((size_t)pp->imax*pp->jmax,NAN);
}

ioflow* sflow_amr::pflow_of(sflow_amr_patch *c)
{
    return (c->pio!=nullptr) ? static_cast<ioflow*>(c->pio) : pflow_void;
}

void sflow_amr::patch_delete(reefamr_patch *q)
{
    sflow_amr_patch *c = SP(q);

    for(auto v : c->nv)
    delete v;
    delete c->mg;
    told_free(c->told);

    if(c->sed!=nullptr)
    {
    delete c->sed;
    if(c->pp->DFBED!=nullptr)
    c->pp->del_Iarray(c->pp->DFBED,c->pp->imax*c->pp->jmax);
    }
    delete c->pmom;
    delete c->psim;
    delete c->pio;
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

    // multi-block: the cells around the blocks and the ring of level 0 hold the restricted values of
    // their blocks (tags, slopes of the prolongation next to the interfaces)
    if(nblk>0)
    fill_level0(pgc,0);

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
            int I = ii-EXT+c->I0, J = jj-EY+c->J0;

            // multi-block: the cell an index around a block stands for
            if(nblk>0)
            {
                int I2,J2;
                if(vcanon(l,I,J,I2,J2))
                {
                    I = I2;
                    J = J2;
                }
            }
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
        gcsl_yrows(pp,pb->depth);       // (one cell in y: the y ghost rows, read below)
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
        // S 10: the bed of the fresh patches from the current bed (it evolved)
        if(p0->S10==1)
        sed_bed_regrid(pgc,l,oldP);

        {
        comms_off guard(pgc);
        for(int n : lev[l])
        if(P[n]->fresh)
        ini_patch_defaults(*SP(n));
        }

        ini_patch_old(l,oldP);

        // shoreline: in the initial hierarchy the children of a dry parent whose bed lies (about) at or above
        // the surface of its wet neighbours are wet below that surface (the initial state is a surface);
        // at a later regrid the children of a dry parent stay dry (its mass is what it is)
        dryfill = ini_pass;
        ini_patch_prolong(l);
        dryfill = false;

        {
        comms_off guard(pgc);
        for(int n : lev[l])
        if(P[n]->fresh)
        ini_patch_finish(p0,*SP(n));
        }

        fill_level(pgc,l,0);

        // S 10: the sediment of the fresh patches starts from their bed
        for(int n : lev[l])
        if(P[n]->fresh)
        sed_patch_ini(pgc,*SP(n));
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

// fine bed: limited-linear prolongation, level by level from level 0 (children average = parent,
// area weighted on a curvilinear grid);
// solid cells and cells next to solids take the parent value without slopes
double sflow_amr::bed_at(int l, int I, int J)
{
    // multi-block: the cell an index around a block stands for
    if(nblk>0)
    {
        int I2,J2;
        if(vcanon(l,I,J,I2,J2))
        {
            I = I2;
            J = J2;
        }
    }

    if(l==0 && nblk>0)
    {
        // blocks 1..: from their nodes (as blocks_ini), beyond a wall the nearest cell of the block
        const int bn = block_near(I,J);
        if(bn>=1)
        {
            int i,j;
            block_own(0,bn,I,J,i,j);
            return block_bed(bn,MAX(MIN(i,blk[bn].ni-1),0),MAX(MIN(j,blk[bn].nj-1),0));
        }

        // native cells beyond the rank arrays (a block's patch reaches them from another rank)
        const int ii = I-O0i, jj = J-O0j;
        if(I>=0 && I<GNX && J>=0 && J<GNY && (ii<p0->imin || ii>=p0->imin+p0->imax || jj<p0->jmin || jj>=p0->jmin+p0->jmax))
        return native_bed(I,J);
    }

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

    int Ic = fsh(I,1), Jc = fsh(J,sy);
    double bc = bed_at(l-1,Ic,Jc);

    if(flag0(fsh(I,l),fsh(J,ysh(l)))<0)
    {
        bedmemo[key]=bc;
        return bc;
    }

    const int lc = l-1;
    auto fl = [&](int a, int d) { return flag0(fsh(a,lc),fsh(d,ysh(lc))); };

    double sx=0.0, sy=0.0;
    if(fl(Ic+1,Jc)>0 && fl(Ic-1,Jc)>0)
    sx = mmod(bed_at(lc,Ic+1,Jc)-bc, bc-bed_at(lc,Ic-1,Jc));
    if(fl(Ic,Jc+1)>0 && fl(Ic,Jc-1)>0)
    sy = mmod(bed_at(lc,Ic,Jc+1)-bc, bc-bed_at(lc,Ic,Jc-1));

    double ox = (I-2*Ic==0) ? -0.25 : 0.25;
    double oy = (this->sy==0) ? 0.0 : ((J-2*Jc==0) ? -0.25 : 0.25);     // (2D grid: no offset in y; sy here is the slope)
    double r = bc + sx*ox + sy*oy;

    // curvilinear grid: the area-weighted mean of the children is the parent bed (the restriction
    // of WL is area weighted; with the plain mean a lake at rest got a surface step at the patches)
    if(p0->geo_curv==1 && (sx!=0.0 || sy!=0.0))
    {
        double sa=0.0, so=0.0;
        for(int a=0; a<2; ++a)
        for(int d=0; d<2; ++d)
        {
            const double A = area_at(l,2*Ic+a,2*Jc+d);
            sa += A;
            so += A*(sx*(a==0 ? -0.25 : 0.25) + sy*(d==0 ? -0.25 : 0.25));
        }
        if(sa>0.0)
        r -= so/sa;
    }

    bedmemo[key]=r;
    return r;
}

// curvilinear grid: area of cell (I,J) of level l, its nodes subdividing the level-0 nodes
// bilinearly (as the patch lexers, reefamr::build_lexer)
double sflow_amr::area_at(int l, int I, int J)
{
    const double rl = double(1<<l);
    double x[2][2], y[2][2];

    // multi-block: the nodes of the block that holds the cell
    const int bk = (nblk>0) ? block_at(l,I,J) : 0;
    auto node = [&](int K, int L, double &xx, double &yy)
    {
        double z;
        if(bk>=1)
        block_node_any(bk,K-blk[bk].vi,L-blk[bk].vj,xx,yy,z);
        else
        p0->curv_node(K,L,xx,yy,z);
    };

    for(int a=0; a<2; ++a)
    for(int d=0; d<2; ++d)
    {
        const int NI = I+a, NJ = J+d;
        const int K = fsh(NI,l), L = fsh(NJ,l);
        const double fr = double(NI-(K<<l))/rl, fs = double(NJ-(L<<l))/rl;
        double x00,y00,x10,y10,x01,y01,x11,y11;
        node(K,L,x00,y00); node(K+1,L,x10,y10);
        node(K,L+1,x01,y01); node(K+1,L+1,x11,y11);
        x[a][d] = x00 + (x10-x00)*fr + (x01-x00)*fs + (x11-x10-x01+x00)*fr*fs;
        y[a][d] = y00 + (y10-y00)*fr + (y01-y00)*fs + (y11-y10-y01+y00)*fr*fs;
    }

    return geo2d::area(x[0][0],y[0][0],x[1][0],y[1][0],x[0][1],y[0][1],x[1][1],y[1][1]);
}

// multi-block: bed of cell (i,j) of block b >= 1, the mean of its node levels (as blocks_ini)
double sflow_amr::block_bed(int bk, int i, int j)
{
    double x,y,z00,z10,z01,z11;
    block_node_any(bk,i,j,x,y,z00);
    block_node_any(bk,i+1,j,x,y,z10);
    block_node_any(bk,i,j+1,x,y,z01);
    block_node_any(bk,i+1,j+1,x,y,z11);
    return MAX(p0->global_zmin,0.25*(z00 + z10 + z01 + z11));
}

// multi-block: bed of a native cell beyond the rank arrays, from the curvilinear nodes
double sflow_amr::native_bed(int I, int J)
{
    double x,y,z00,z10,z01,z11;
    p0->curv_node(I,J,x,y,z00);
    p0->curv_node(I+1,J,x,y,z10);
    p0->curv_node(I,J+1,x,y,z01);
    p0->curv_node(I+1,J+1,x,y,z11);
    return MAX(p0->global_zmin,0.25*(z00 + z10 + z01 + z11));
}

// --------------------------------------------------------------------- refinement flags
// cells of level l on this rank that need level l+1 (evaluated on the level-l grids)
void sflow_amr::tag(int l, vector<unsigned char> &M)
{
    const int lf = l+1;

    // shoreline (G 22): the cells on both sides of a change between wet and dry (the wet flag of the
    // solver, A 244), however thin the film at the change, if water deeper than 3 A 244 lies within
    // two cells (films left on a dry slope are not flagged).  Their buffer covers the distance the
    // shoreline can travel until the next regrid, in level-l cells: per step of the finest level at
    // most 2 N 47 of its cells (sflow_amr::timestep; 2: one cell in y), over the steps between two
    // regrids (G 2, under G 7 1 in whole level-0 steps).  At least G 3.
    const double hfilm = 3.0*p0->A244;
    int nshore = nbuf;
    if(shore>0)
    {
        const double cfl = (p0->A220==4) ? MIN(p0->N47,0.25) : p0->N47;
        const int nsteps = (sub==1) ? (MAX(regrid_int>>maxlev,1)<<maxlev) : MAX(regrid_int,1);
        const double cells = 2.0*cfl*double(nsteps)/double(1<<(maxlev-l));
        nshore = MAX(nbuf, int(ceil(cells-1.0e-9)));
    }

    // water deeper than hfilm within two cells of (ii,jj)
    auto deep_near = [&](lexer *q, slice &WL, int ii, int jj)
    {
        for(int a=ii-2; a<=ii+2; ++a)
        for(int d=jj-2*sy; d<=jj+2*sy; ++d)
        if(q->flagslice4[lij(q,a,d)]>=0 && q->wet[lij(q,a,d)]==1 && WL(a,d)>hfilm)
        return true;
        return false;
    };

    // a wet cell (G 21: the three cells of the second difference)
    auto wetc = [&](lexer *q, int a, int d) { return q->flagslice4[lij(q,a,d)]>=0 && q->wet[lij(q,a,d)]==1; };

    // 0: no flag, 1: surface jump (G 20), second difference (G 21), breaking (G 23), 2: shoreline (G 22)
    auto test = [&](lexer *q, fdm2D *bb, slice &WL, int ii, int jj)
    {
        if(q->flagslice4[lij(q,ii,jj)]<0)
        return 0;

        int wc = q->wet[lij(q,ii,jj)];
        double ec = WL(ii,jj)-bb->depth(ii,jj);
        const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};
        int r = 0;

        // breaking (G 23): the flag of the grid in the step (breaking_print: A 248 0 clears breaking at
        // the end of the step, before the regrid of the patches)
        if(flagbreak && wc==1 && (bb->breaking(ii,jj)>0 || bb->breaking_print(ii,jj)>0.0))
        r = 1;

        // second difference of the surface along x or y (G 21), between wet cells
        if(tol_curv>0.0 && wc==1)
        {
            auto es = [&](int a, int d) { return WL(a,d)-bb->depth(a,d); };
            if(wetc(q,ii+1,jj) && wetc(q,ii-1,jj) && fabs(es(ii+1,jj)-2.0*ec+es(ii-1,jj))>tol_curv)
            r = 1;
            if(sy==1 && wetc(q,ii,jj+1) && wetc(q,ii,jj-1) && fabs(es(ii,jj+1)-2.0*ec+es(ii,jj-1))>tol_curv)
            r = 1;
        }

        for(int k=0; k<2+2*sy; ++k)     // (2D grid: the x neighbours)
        {
            int a=ii+di[k], d=jj+dj[k];
            if(q->flagslice4[lij(q,a,d)]<0)
            continue;

            int wn = q->wet[lij(q,a,d)];

            if(shore>0 && wc!=wn && deep_near(q,WL,ii,jj))
            return 2;

            if(tol_eta>0.0 && wc==1 && wn==1)
            if(fabs(ec-(WL(a,d)-bb->depth(a,d)))>tol_eta)
            r = 1;
        }
        return r;
    };

    auto flag = [&](int r, int I, int J)
    {
        if(r>0)
        tag_cell(lf,I,J,(r==2)?nshore:nbuf,M);
    };

    if(l==0)
    {
        for(int ii=0; ii<NX0; ++ii)
        for(int jj=0; jj<NY0; ++jj)
        if(nblk==0 || nring[lij(p0,ii,jj)]==0)
        flag(test(p0,b0,b0->WL,ii,jj),ii+O0i,jj+O0j);

        // multi-block: the blocks 1.. (level-0 patches)
        for(int n : lev[0])
        {
            sflow_amr_patch *c = SP(n);
            for(int ii=EXT; ii<EXT+c->nx; ++ii)
            for(int jj=EY; jj<EY+c->ny; ++jj)
            flag(test(c->pp,c->b,c->b->WL,ii,jj),ii-EXT+c->I0,jj-EY+c->J0);
        }
        return;
    }

    for(int n : lev[l])
    {
        sflow_amr_patch *c = SP(n);
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EY; jj<EY+c->ny; ++jj)
        flag(test(c->pp,c->b,c->b->WL,ii,jj),ii-EXT+c->I0,jj-EY+c->J0);
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
        b->breaking(ii,jj) = 0;
    }

    c.blkold.assign((size_t)(c.nx/2)*(c.ny>>sy),0);
}

// a solid block (its first cell) keeps the defaults
bool sflow_amr::block_solid(sflow_amr_patch &c, int bi, int bj)
{
    return c.pp->flagslice4[lij(c.pp,EXT+2*bi,EY+ychild(bj))]<0;
}

// the cells of the fresh level-l patches that an old patch held (all four cells of a block lie in
// the same old patch: patch boxes are tile aligned)
void sflow_amr::ini_patch_old(int l, vector<reefamr_patch*> &oldP)
{
    old_run(l,17,7140+l,oldP,
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
                // shoreline offset of a covered cell, as a patch that is kept holds it (prolongation of
                // the fresh children)
                v[15] = o->sfc[lij(qp,io,jo)];
                // breaking with persistence (A 248 1): the flag carries over to the next step
                v[16] = ob->breaking(io,jo);
            },
            [&](reefamr_patch *q, int id, int ii, int jj, const double *v)
            {
                sflow_amr_patch &c = *SP(q);
                const int bi = (ii-EXT)/2, bj = (jj-EY)>>sy;
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
                c.sfc[lij(pp,ii,jj)] = v[15];
                if(p0->A248==1)
                b->breaking(ii,jj) = (int)v[16];
                c.blkold[(size_t)bi*(c.ny>>sy)+bj] = 1;
            });
}

// the other blocks of the fresh level-l patches: prolonged from the parent (prolong_parts of the
// four children and the parent's values on its rank), finished with the child depths, then the
// mass and momentum of a fully wet parent restored exactly
void sflow_amr::ini_patch_prolong(int l)
{
    const int nv = 4*NPR+7;

    block_down_if(l,nv,7160+l,
                  [&](reefamr_patch *q) { return q->fresh; },
                  [&](const reefamr_block &B, int key, double *v)
                  {
                      gh G = grid(B.g);
                      const int ic = B.ic, jc = B.jc;
                      for(int a=0;a<2;++a)
                      for(int d=0;d<=sy;++d)
                      prolong_parts(B.g,ic,jc,a==0?-1:1,(sy==0) ? 0 : (d==0?-1:1),&v[(2*a+d)*NPR]);
                      double *w = &v[4*NPR];
                      w[0] = G.q->wet[lij(G.q,ic,jc)];
                      w[1] = (*sWL[B.g+1])(ic,jc);
                      w[2] = (*sUH[B.g+1])(ic,jc);
                      w[3] = (*sVH[B.g+1])(ic,jc);
                      w[4] = (*sWH[B.g+1])(ic,jc);
                      w[5] = G.b->press(ic,jc);
                      w[6] = G.b->breaking(ic,jc);
                  },
                  [&](reefamr_patch *q, int id, int kb, const double *v)
                  {
                      sflow_amr_patch &c = *SP(q);
                      const int nby = c.ny>>sy;
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
    const int i0 = EXT+2*bi, j0 = EY+ychild(bj);
    const double *w = &pv[4*NPR];
    double v[NVMAX];

    for(int a=0;a<2;++a)
    for(int d=0;d<=sy;++d)
    {
        int ii = i0+a, jj = j0+d;
        prolong_finish(&pv[(2*a+d)*NPR],b->depth(ii,jj),v);
        b->WL(ii,jj)=v[0]; b->UH(ii,jj)=v[1]; b->VH(ii,jj)=v[2]; b->WH(ii,jj)=v[8];
        b->press(ii,jj) = w[5];
        // breaking with persistence (A 248 1): the children take the flag of the parent
        if(p0->A248==1)
        b->breaking(ii,jj) = (int)w[6];
        if(bous==1)
        {
        b->UA(ii,jj)=v[10]; b->VA(ii,jj)=v[11]; b->MX(ii,jj)=v[12]; b->MY(ii,jj)=v[13];
        }
        pp->wet[lij(pp,ii,jj)]=(int)v[6];
    }

    // (2D grid: the 2x1 children)
    const int nch = 2<<sy;
    int nw = (sy==0) ? pp->wet[lij(pp,i0,j0)] + pp->wet[lij(pp,i0+1,j0)]
                     : pp->wet[lij(pp,i0,j0)] + pp->wet[lij(pp,i0+1,j0)] + pp->wet[lij(pp,i0,j0+1)] + pp->wet[lij(pp,i0+1,j0+1)];

    // curvilinear grid: the volumes of the children (area weighted) sum to the parent's
    if(nw==nch && (int)w[0]==1 && pp->geo_curv==1)
    {
        double ar[2][2], at=0.0;
        for(int a=0;a<2;++a)
        for(int d=0;d<=sy;++d)
        {
        ar[a][d] = pp->AREA2D[lij(pp,i0+a,j0+d)];
        at += ar[a][d];
        }
        double sw=0.0, su=0.0, sv=0.0, sh=0.0;
        for(int a=0;a<2;++a)
        for(int d=0;d<=sy;++d)
        {
        sw += ar[a][d]*b->WL(i0+a,j0+d);
        su += ar[a][d]*b->UH(i0+a,j0+d);
        sv += ar[a][d]*b->VH(i0+a,j0+d);
        sh += ar[a][d]*b->WH(i0+a,j0+d);
        }
        const double du = at*w[2] - su, dv = at*w[3] - sv, dw = at*w[4] - sh;
        const double fw = at*w[1]/sw;

        for(int a=0;a<2;++a)
        for(int d=0;d<=sy;++d)
        {
            b->UH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*du;
            b->VH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*dv;
            b->WH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*dw;
            b->WL(i0+a,j0+d) *= fw;
        }
    }
    else
    if(nw==nch && (int)w[0]==1 && sy==0)
    {
        double sw = b->WL(i0,j0)+b->WL(i0+1,j0);
        double su = b->UH(i0,j0)+b->UH(i0+1,j0);
        double sv = b->VH(i0,j0)+b->VH(i0+1,j0);
        double sh = b->WH(i0,j0)+b->WH(i0+1,j0);
        double du = 2.0*w[2] - su;
        double dv = 2.0*w[3] - sv;
        double dw = 2.0*w[4] - sh;
        double fw = 2.0*w[1]/sw;

        for(int a=0;a<2;++a)
        {
            b->UH(i0+a,j0) += b->WL(i0+a,j0)/sw*du;
            b->VH(i0+a,j0) += b->WL(i0+a,j0)/sw*dv;
            b->WH(i0+a,j0) += b->WL(i0+a,j0)/sw*dw;
            b->WL(i0+a,j0) *= fw;
        }
    }
    else
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
        for(int d=0;d<=sy;++d)
        {
            b->UH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*du;
            b->VH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*dv;
            b->WH(i0+a,j0+d) += b->WL(i0+a,j0+d)/sw*dw;
            b->WL(i0+a,j0+d) *= fw;
        }
    }

    for(int a=0;a<2;++a)
    for(int d=0;d<=sy;++d)
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
        pp->wet_n[lij(pp,ii,jj)] = pp->wet[lij(pp,ii,jj)];
    }

    pp->dt = p->dt;
    pp->dt_old = p->dt_old;
    pp->simtime = p->simtime;
    pp->count = p->count;
}

// --------------------------------------------------------------------- multi-block level 0
// blocks 1..: bed from the block nodes (mean of the four node bed levels, as lexer::geometry_bed),
// state as sflow_f::ini_fsf (F 60, F 72 boxes by the cell centroid), wet flags as level 0
void sflow_amr::blocks_ini(lexer *p, ghostcell *pgc)
{
    // the fine bed of the patches takes level-0 beds beyond the rank arrays from the nodes
    // (native_bed): check that the native bed of block 0 is the mean of the node levels
    if(maxlev>=1)
    {
        int bad=0;
        for(int ii=0; ii<p->knox; ++ii)
        for(int jj=0; jj<p->knoy; ++jj)
        if(nring[lij(p,ii,jj)]==0 && p->flagslice4[lij(p,ii,jj)]>0 && b0->bed(ii,jj)!=native_bed(ii+O0i,jj+O0j))
        ++bad;
        bad = pgc->globalisum(bad);
        if(bad>0 && p->mpirank==0)
        cout<<"SFLOW AMR: multi-block grid, "<<bad<<" level-0 cells with a bed other than the mean of their nodes (bit for bit);"
            <<" the fine bed next to the rank boundaries may differ"<<endl;
    }

    comms_off guard(pgc);

    for(int id : lev[0])
    {
        sflow_amr_patch *c = SP(id);
        lexer *pp = c->pp;
        fdm2D *pb = c->b;
        c->phll->amr_id = id+1;

        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            const int i = ii-EXT, j = jj-EXT;
            // as lexer::geometry_bed (the same order of the sum)
            double x,y,z00,z10,z01,z11;
            block_node_any(c->blk,i,j,x,y,z00);
            block_node_any(c->blk,i+1,j,x,y,z10);
            block_node_any(c->blk,i,j+1,x,y,z01);
            block_node_any(c->blk,i+1,j+1,x,y,z11);
            const double zb = 0.25*(z00 + z10 + z01 + z11);
            pb->bed(ii,jj) = pb->bed0(ii,jj) = pb->topobed(ii,jj) = pb->solidbed(ii,jj) = MAX(zb,p->global_zmin);
            pb->depth(ii,jj) = p->wd - pb->bed(ii,jj);
            pb->ks(ii,jj) = p->B50;
        }

        for(int ii=pp->imin; ii<pp->imin+pp->imax-1; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        pb->dfx(ii,jj) = 0.5*(pb->depth(ii+1,jj)+pb->depth(ii,jj));

        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax-1; ++jj)
        pb->dfy(ii,jj) = 0.5*(pb->depth(ii,jj+1)+pb->depth(ii,jj));

        ini_patch_defaults(*c);

        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            const int q = lij(pp,ii,jj);
            double eta = 0.0;
            for(int qn=0; qn<p->F72; ++qn)
            if(pp->XC2D[q]>=p->F72_xs[qn] && pp->XC2D[q]<p->F72_xe[qn] && pp->YC2D[q]>=p->F72_ys[qn] && pp->YC2D[q]<p->F72_ye[qn])
            eta = p->F72_h[qn]-p->wd;
            pb->eta(ii,jj) = eta;
            pb->WL(ii,jj) = MAX(eta + pb->depth(ii,jj),0.0);
        }

        // X 10 3 (X 207 0): the draft of the ship depresses the water level at t = 0, as on level 0
        // (sixdof_obj_2D::ini_shipwave)
        if(shipmode==3 && p->X207==0)
        {
            zone_setup(p0);
            ship_fields(*c,false);
            for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
            for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
            {
                i=ii; j=jj;
                const double H = c->pship->Hsolid(pp,pb->fs);
                pb->eta(ii,jj) = -H*c->pship->draft(ii,jj);
                pb->WL(ii,jj) = MAX(pb->eta(ii,jj) + pb->depth(ii,jj),0.0);
            }
        }

        for(int ii=0; ii<pp->knox; ++ii)
        for(int jj=0; jj<pp->knoy; ++jj)
        pp->wet[lij(pp,ii,jj)]=1;

        // B 522 relaxation and beach zones: the state at the start as sflow_f::ini
        if(c->pio!=nullptr)
        c->pio->zones_ini2D(pp,pb,pgc);

        c->pfsf->depth_update(pp,pb,pgc,pb->WL);

        ini_patch_finish(p,*c);
    }
}

// the cells around the blocks and the ring of level 0, from the stage input of their holders
void sflow_amr::fill_level0(ghostcell *pgc, int s)
{
    fill_run0(NV,7099,
              [&](const reefamr_fill &f, double *v) { eval_fill(f,v); },
              [&](reefamr_patch *c, int n, const reefamr_fill &f, const double *w) { store_fill(SP(c),n,f.di,f.dj,w); },
              [&](const reefamr_fill &f, const double *w) { store_native(f.di,f.dj,w); });

    {
    comms_off guard(pgc);
    for(int n : lev[0])
    apply_bc(pgc,*SP(n),s);
    }

    // the ghost copies of the ring cells on the other ranks
    exchange_level0(p0,b0,pgc,s,true);
}

// values of a ring cell into level 0 (stage input, as store_fill)
void sflow_amr::store_native(int ii, int jj, const double *w)
{
    lexer *p = p0;
    fdm2D *b = b0;
    (*sWL[0])(ii,jj) = w[0];
    (*sUH[0])(ii,jj) = w[1];
    (*sVH[0])(ii,jj) = w[2];
    b->eta(ii,jj) = w[3];
    b->U(ii,jj) = w[4];
    b->V(ii,jj) = w[5];
    b->hp(ii,jj) = w[0];
    p->wet[lij(p,ii,jj)] = (int)w[6];
    p->deep[lij(p,ii,jj)] = (int)w[7];
    (*sWH[0])(ii,jj) = w[8];
    b->W(ii,jj) = w[9];
    if(bous==1)
    {
    b->UA(ii,jj) = w[10];
    b->VA(ii,jj) = w[11];
    b->MX(ii,jj) = w[12];
    b->MY(ii,jj) = w[13];
    }
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

// shoreline: the restricted water surface of grid g (index id+1 as sWL), NaN where none
double *sflow_amr::sfc_of(int g)
{
    if(g<0)
    return sfs0->V;
    return SP(g)->sfc.data();
}

// G 7 1: the covered cells of the grids of level l at a shoreline after a stage of their own (the restriction
// follows after the steps of the finer level): their surface is the restricted surface of the wet children,
// moved with the water level.  The halo of level 0 holds the offsets and wet flags of the neighbour ranks
// (exchange_level0), so it is shifted alike without an exchange; cells without an offset keep their bits
void sflow_amr::sfc_stage(int l)
{
    auto shift = [&](lexer *q, fdm2D *b, const double *D)
    {
        for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
        for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
        {
            const int k = lij(q,ii,jj);
            if(!std::isnan(D[k]) && q->wet[k]==1)
            b->eta(ii,jj) += D[k];
        }
    };

    if(l==0)
    shift(p0,b0,sfs0->V);

    for(int n : lev[l])
    shift(SP(n)->pp,SP(n)->b,SP(n)->sfc.data());
}

// the grids of level l hold no restricted surface (before their cells are restricted again)
void sflow_amr::sfc_reset(int l)
{
    if(l==0)
    std::fill(sfs0->V,sfs0->V+(size_t)p0->imax*p0->jmax,NAN);

    for(int n : lev[l])
    std::fill(SP(n)->sfc.begin(),SP(n)->sfc.end(),NAN);
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

// r[0] 0 dry parent, 1 shallow parent (its values), 2 interpolated, 3 dry parent whose bed lies above
// the surface of its wet neighbours (dryfill), 4 shallow parent with a water surface (dryfill); r[1] surface (2, 3, 4)
// or water level (1), r[2..4] velocities (2) or UH, VH, WH (1, 4), r[5..8] Boussinesq u_a, v_a, M/H,
// r[9] the parent's water level (4)
//
// Shoreline: the surface of a covered cell is the restricted surface of its wet children (sfc: its offset), not its
// mean water level plus mean bed, which lies above the water where children are dry.  With dryfill (the
// fresh patches of the initial hierarchy: the initial state is a water surface) the children of a dry or
// shallow parent take the surface of the parent or of its wet neighbours: a lake at rest starts at rest
// where the coarser grid does not resolve the shoreline.  Not for the ghost cells (a wet ghost cell would
// exchange mass with a coarse cell that cannot take it) and not at a later regrid (the children of a dry
// parent stay dry: its mass is what it is)
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

    const double *SF = sfc_of(g);

    auto get = [&](int a, int bb, double &e, double &u, double &vv, double &ww, int &wt)
    {
        double wlc = wlat(a,bb);
        const double sf = SF[lij(q,a,bb)];
        e = wlc - pb->depth(a,bb) + (std::isnan(sf) ? 0.0 : sf);
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

    // shoreline: the water surface at the parent, from its own (restricted) surface or, dry, from its
    // wet neighbours (the lowest), NaN if there is none
    auto surface = [&]() -> double
    {
        if(w0==-1)
        return NAN;
        if(w0==1 || !std::isnan(SF[lij(q,ic,jc)]))
        return e0;
        double sm = NAN;
        for(int a=-1; a<=1; ++a)
        for(int d=-sy; d<=sy; ++d)
        if(a!=0 || d!=0)
        {
            double en,un,vn,wn; int wtn;
            get(ic+a,jc+d,en,un,vn,wn,wtn);
            if(wtn==1 || (wtn==0 && !std::isnan(SF[lij(q,ic+a,jc+d)])))
            sm = std::isnan(sm) ? en : MIN(sm,en);
        }
        return sm;
    };

    if(w0!=1)
    {
        if(!dryfill)
        return;
        // only where the surface lies below the parent's (mean) bed plus the wet-dry depth: the coarser grid
        // holds it dry because its bed reaches (about) to the water, not because a wet front has yet to reach it
        const double sm = surface();
        if(!std::isnan(sm) && sm + pb->depth(ic,jc) <= wd+eps)
        {
            r[0] = 3.0;
            r[1] = sm;
        }
        return;
    }

    double hc = wlat(ic,jc);
    if(hc<3.0*wd)
    {
        r[0] = 1.0;
        r[1] = hc; r[2] = uhat(ic,jc); r[3] = vhat(ic,jc); r[4] = what(ic,jc);
        if(dryfill)
        {
            r[0] = 4.0;
            r[1] = surface();
            r[9] = hc;
        }
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
    // (2D grid: the slopes in x only, the cells beyond the strip are no neighbours)
    const bool nbs = (sy==0) ? (wE==1 && wW==1) : (wE==1 && wW==1 && wN==1 && wS==1);
    if(sy==0 && nbs)
    {
        sxe = mmod(eE-e0,e0-eW);
        sxu = mmod(uE-u0,u0-uW);
        sxv = mmod(vE-v0,v0-vW);
        sxw = mmod(wEv-w0v,w0v-wWv);
    }
    else
    if(nbs)
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
        const bool full = nbs;
        partb(full ? fx : 0.0, full ? fy : 0.0);
    }
}

void sflow_amr::prolong_finish(const double *r, double depthf, double *v)
{
    const double wd = p0->A244;

    // shoreline: the water surface of the parent (ghost cells, sflow_amr::hll_hook), NaN for a dry one
    v[ISF] = (r[0]==2.0 || r[0]==4.0) ? r[1] : NAN;

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

    // shoreline (dryfill): wet below the surface, at rest (dry parent) or with the parent's flux per
    // width, at most its velocity (shallow parent)
    if(r[0]==3.0 || r[0]==4.0)
    {
        const double hs = r[1] + depthf;
        if(hs<=wd+eps)
        {
            dry();
            return;
        }
        const double f = (r[0]==4.0) ? MIN(1.0,hs/r[9]) : 0.0;
        fin(hs, f*r[2], f*r[3], f*r[4]);
        if(bous==1)
        {
            if(r[0]==3.0)
            v[10]=v[11]=v[12]=v[13]=0.0;
            else
            finb(hs);
        }
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
        v[ISC] = NAN;
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
        v[ISF] = NAN;
        // shoreline: the offset of a covered cell (restriction), so that the prolongation from this grid
        // sees the same surface around its box as the grid that holds the cell
        v[ISC] = sfc_of(f.g)[lij(G.q,f.si,f.sj)];
        return;
    }

    const double wd = p0->A244;
    v[0]=wd; v[1]=v[2]=0.0; v[3]=wd-f.aux; v[4]=v[5]=0.0; v[6]=v[7]=0.0; v[8]=v[9]=0.0;
    if(bous==1)
    v[10]=v[11]=v[12]=v[13]=0.0;
    v[ISF] = NAN;
    v[ISC] = NAN;
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
    c->gsf[lij(pp,ii,jj)] = w[ISF];
    c->sfc[lij(pp,ii,jj)] = w[ISC];
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

    // shoreline: the surface of the wet children (v[nv]); the parents of other cells hold none
    sfc_reset(l-1);

    block_up(l,nv+1,7020+l,
             [&](reefamr_patch *q, int id, int k, double *v)
             {
                 sflow_amr_patch &c = *SP(q);
                 fdm2D *b = c.b;
                 slice *WLi,*UHi,*VHi,*WLo,*UHo,*VHo;
                 c.pmom->stage_io(s,b,WLi,UHi,VHi,WLo,UHo,VHo);
                 slice *WHi,*WHo;
                 c.pmom->stage_io_w(s,b,WHi,WHo);

                 const int nby = c.ny>>sy;
                 const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);

                 // shoreline: the water surface of the wet children (area weighted on a curvilinear grid)
                 // where some children are dry (or have dry children themselves), NaN otherwise.  The
                 // parent's own surface (mean water level + mean bed) lies above it: the bed of the dry
                 // children reaches above the water
                 {
                 lexer *pp = c.pp;
                 double se=0.0, sa=0.0;
                 int nf=0, nw=0, np=0;
                 for(int a=0;a<2;++a)
                 for(int d=0;d<=sy;++d)
                 if(pp->flagslice4[lij(pp,i0+a,j0+d)]>0)
                 {
                 ++nf;
                 // a child covered by a finer patch at a shoreline: the surface of its own wet children
                 const double dc = c.sfc[lij(pp,i0+a,j0+d)];
                 if(!std::isnan(dc))
                 ++np;
                 if((*WLo)(i0+a,j0+d)>p->A244+eps)
                 {
                 ++nw;
                 const double A = (pp->geo_curv==1) ? pp->AREA2D[lij(pp,i0+a,j0+d)] : 1.0;
                 se += A*((*WLo)(i0+a,j0+d) - b->depth(i0+a,j0+d) + (std::isnan(dc) ? 0.0 : dc));
                 sa += A;
                 }
                 }
                 v[nv] = (sa>0.0 && (nw<nf || np>0)) ? se/sa : NAN;
                 }

                 // curvilinear grid: area weighted (the children of a cell differ in area)
                 if(c.pp->geo_curv==1)
                 {
                 lexer *pp = c.pp;
                 double w[4], ws=0.0;
                 for(int a=0;a<2;++a)
                 for(int d=0;d<=sy;++d)
                 {
                 w[2*a+d] = pp->AREA2D[lij(pp,i0+a,j0+d)];
                 ws += w[2*a+d];
                 }
                 auto avg = [&](slice &f) { return (w[0]*f(i0,j0)+w[1]*f(i0,j0+1)+w[2]*f(i0+1,j0)+w[3]*f(i0+1,j0+1))/ws; };
                 v[0] = avg(*WLo);
                 v[1] = avg(*UHo);
                 v[2] = avg(*VHo);
                 v[3] = avg(*WHo);
                 if(bous==1)
                 {
                 v[4] = avg(b->UA);
                 v[5] = avg(b->VA);
                 v[6] = avg(b->MX);
                 v[7] = avg(b->MY);
                 }
                 return;
                 }

                 v[0] = (sy==0) ? 0.5*((*WLo)(i0,j0)+(*WLo)(i0+1,j0)) : 0.25*((*WLo)(i0,j0)+(*WLo)(i0+1,j0)+(*WLo)(i0,j0+1)+(*WLo)(i0+1,j0+1));
                 v[1] = (sy==0) ? 0.5*((*UHo)(i0,j0)+(*UHo)(i0+1,j0)) : 0.25*((*UHo)(i0,j0)+(*UHo)(i0+1,j0)+(*UHo)(i0,j0+1)+(*UHo)(i0+1,j0+1));
                 v[2] = (sy==0) ? 0.5*((*VHo)(i0,j0)+(*VHo)(i0+1,j0)) : 0.25*((*VHo)(i0,j0)+(*VHo)(i0+1,j0)+(*VHo)(i0,j0+1)+(*VHo)(i0+1,j0+1));
                 v[3] = (sy==0) ? 0.5*((*WHo)(i0,j0)+(*WHo)(i0+1,j0)) : 0.25*((*WHo)(i0,j0)+(*WHo)(i0+1,j0)+(*WHo)(i0,j0+1)+(*WHo)(i0+1,j0+1));

                 // Boussinesq: u_a and M
                 if(bous==1)
                 {
                 auto avg = [&](slice &f) { return (sy==0) ? 0.5*(f(i0,j0)+f(i0+1,j0)) : 0.25*(f(i0,j0)+f(i0+1,j0)+f(i0,j0+1)+f(i0+1,j0+1)); };
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

                 // shoreline: the surface the uncovered neighbours see (reconstruction, pre-balanced source) is
                 // the surface of the wet children, not the mean water level plus the mean bed; kept as the
                 // offset from the latter (sfc), which the water level carries along (G 7 1: the coarse steps)
                 const double ds = std::isnan(v[nv]) ? NAN : v[nv] - G.b->eta(ic,jc);
                 if(w==1 && !std::isnan(ds))
                 G.b->eta(ic,jc) = v[nv];
                 G.b->U(ic,jc) = w==1 ? uh/wlvl : 0.0;
                 G.b->V(ic,jc) = w==1 ? vh/wlvl : 0.0;
                 G.b->W(ic,jc) = w==1 ? wh/wlvl : 0.0;
                 G.b->hp(ic,jc) = wl;

                 sfc_of(B.g)[lij(G.q,ic,jc)] = ds;

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
void sflow_amr::exchange_level0(lexer *p, fdm2D *b, ghostcell *pgc, int s, bool inputs)
{
    if(p->mpi_size<=1)
    return;

    slice *WLi,*UHi,*VHi,*WLo,*UHo,*VHo;
    pmom0->stage_io(s,b,WLi,UHi,VHi,WLo,UHo,VHo);

    slice *WHi,*WHo;
    pmom0->stage_io_w(s,b,WHi,WHo);

    // (multi-block: the stage input after the ring of level 0 was filled)
    if(inputs)
    {
        WLo=WLi; UHo=UHi; VHo=VHi; WHo=WHi;
    }

    slice *f[13] = {WLo,UHo,VHo,WHo,&b->eta,&b->U,&b->V,&b->W,&b->hp,&b->UA,&b->VA,&b->MX,&b->MY};
    for(int k=0; k<(bous==1?13:9); ++k)
    {
        pgc->gcslparax(p,*f[k],4);
        pgc->gcslparacox(p,*f[k],10);
    }

    pgc->gcslparaxV_int(p,p->wet,4);
    pgc->gcslparacoxV_int(p,p->wet,50);

    // shoreline: the restricted surface of the covered cells in the halo (prolongation next to the partition edge)
    pgc->gcslparax(p,*sfs0,4);
}

// fine face values of the level-l patches on the partition edges, to the rank of the coarse cell
void sflow_amr::exchange_fluxes(int l)
{
    face_run(l,5,7200+l,
             [&](reefamr_patch *q, int side, int r, double *sb)
             {
                 sflow_amr_patch *c = SP(q);
                 for(int ip=0; ip<5; ++ip)
                 sb[ip] = 0.5*(c->rec[ip][side][r]+c->rec[ip][side][r+fr1(side)]);
             },
             [&](reefamr_match &m, const double *v)
             {
                 for(int ip=0; ip<5; ++ip)
                 m.val[ip] = v[ip];
             });

    // the same faces for the grids of other ranks that compute them around their box (q.ring)
    face_run_ring(l,5,7330+l,
             [&](reefamr_patch *q, int side, int r, double *sb)
             {
                 sflow_amr_patch *c = SP(q);
                 for(int ip=0; ip<5; ++ip)
                 sb[ip] = 0.5*(c->rec[ip][side][r]+c->rec[ip][side][r+fr1(side)]);
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
        {
            // shoreline: the covered cells of level 0 start the step with the restricted surface (depth_update
            // at the end of the last step set it from the mean water level and the mean bed)
            sfc_stage(0);
            sub_begin(p,b,pgc);
        }
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
        c->pp->amr_follow(*p);
        c->pmom->inflow(c->pp,c->b,pgc,pflow_of(c));
    }

    // body on the patches: X 10 3 updates the ship pressure once per step (as level 0),
    // X 10 2 moves the body in every stage
    ship_patches(shipmode==3);
}

void sflow_amr::stage_begin(lexer *p, fdm2D *b, ghostcell *pgc, int s)
{
    // multi-block: the blocks and the ring of level 0 filled from the stage input, the patches
    // (finest first) and the blocks run their stages, level 0 follows
    if(nblk>0)
    {
        double t0 = MPI_Wtime();
        cache_stage(s);
        fill_level0(pgc,s);
        if(patches_total>0)
        for(int l=1; l<=maxlev; ++l)
        fill_level(pgc,l,s);

        // X 10 2: the body moves in every stage (on the blocks and the patches)
        if(shipmode==2 && s>0)
        ship_patches(false);
        tm[0] += MPI_Wtime()-t0;

        if(patches_total>0)
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

        t0 = MPI_Wtime();
        {
        comms_off guard(pgc);
        for(int n : lev[0])
        SP(n)->pmom->rk_stage(SP(n)->pp,SP(n)->b,pgc,s);
        }
        tm[1] += MPI_Wtime()-t0;
        return;
    }

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
    // G 7 1: level 0 alone, its covered cells at a shoreline keep the restricted surface
    if(maxlev>=1 && patches_total>0 && sub==1)
    {
        sfc_stage(0);
        return;
    }

    if(maxlev<1 || patches_total==0 || sub==1)
    return;

    double t0 = MPI_Wtime();
    restrict_levels(p,s);

    exchange_level0(p,b,pgc,s);
    tm[3] += MPI_Wtime()-t0;
}

void sflow_amr::step_end(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(!active())
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
    // the cells around the patches from their sources first: the wet-dry update of a cell next to
    // the box reads its neighbours there, which otherwise held the patch's own last stage, so that
    // the result depended on where the boxes lie (another partition, another layout of the patches)
    if(patches_total>0)
    {
        double t0 = MPI_Wtime();
        cache_stage(0);
        for(int l=1; l<=maxlev; ++l)
        fill_level(pgc,l,0);
        tm[0] += MPI_Wtime()-t0;
    }

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

    // (replace_due: G 41 0 with a static refinement, the placement once more with the measured model)
    if(maxlev>=1 && ((regrid_int>0 && p->count%rint==0) || replace_due(p)))
    {
    double t0 = MPI_Wtime();
    regrid(p,pgc,false);
    tm[4] += MPI_Wtime()-t0;
    }

    // multi-block: the ring of level 0 again from the blocks, for the time step, the output and the
    // non-staged fields of the next step (the outer ring cells drift from the wall of the native box
    // during the step; the CFL number of level 0 runs over the ring; after a regrid: the restricted
    // values under the new patches)
    if(nblk>0)
    {
    double t0 = MPI_Wtime();
    cache_stage(0);
    fill_level0(pgc,0);
    tm[0] += MPI_Wtime()-t0;
    }

    // load statistics of the ranks (G 43)
    load_step(p);
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
        const int jl = EY-1, jh = EY+c.ny-1;

        // shoreline: the boundary faces the fine fluxes of which replace the coarse ones (refluxing).
        //  - no mass flux between a dry cell of the patch and a wet cell beyond its boundary, as between two
        //    cells of a grid (sflow_HLL::aij_E zeroes the faces of the wet cells it visits, the ghost cells it
        //    does not)
        //  - a face between a dry cell of the patch and a dry ghost cell of a wet parent (the bed of both
        //    reaches above the surface of the parent) is a wall of the parent, as the face between a wet and
        //    a dry cell on the coarse grid: the pressure of the parent's surface with the face depth, no mass
        //    or momentum flux.  The flux between the two dry cells (their surfaces are bed + film) would
        //    move the water of the parent, which is at rest
        {
            lexer *pp = c.pp;
            const double g = fabs(p0->W22);
            auto drycell = [&](int a, int d) { return pp->wet[lij(pp,a,d)]==0 && pp->flagslice4[lij(pp,a,d)]>0; };
            auto wetcell = [&](int a, int d) { return pp->wet[lij(pp,a,d)]==1 && pp->flagslice4[lij(pp,a,d)]>0; };

            // the face (dir 1: i-face, 2: j-face) at (fi,fj) between the patch cell (a,d) and the ghost cell (a2,d2)
            auto face = [&](int dir, int fi, int fj, int a, int d, int a2, int d2)
            {
                slice &F = (dir==1) ? Fx : Fy;
                const double es = c.gsf[lij(pp,a2,d2)];
                if(ipol==4 && drycell(a,d) && wetcell(a2,d2))
                F(fi,fj) = 0.0;

                if(!(drycell(a,d) && drycell(a2,d2) && !std::isnan(es)))
                return;

                if(ipol==3 || ipol==4)
                {
                    F(fi,fj) = 0.0;
                    return;
                }

                const double df = (dir==1) ? c.b->dfx(fi,fj) : c.b->dfy(fi,fj);
                const double P = 0.5*g*es*es + g*es*df;
                double nx = (dir==1) ? 1.0 : 0.0, ny = (dir==1) ? 0.0 : 1.0, len;
                if(pp->geo_curv==1)
                sflow_curv::normal(pp,dir,fi,fj,nx,ny,len);
                F(fi,fj) = P*((ipol==1) ? nx : ny);
            };

            for(int r=0; r<c.ny; ++r)
            {
                face(1,il,EY+r,il+1,EY+r,il,EY+r);
                face(1,ih,EY+r,ih,EY+r,ih+1,EY+r);
            }
            if(sy==1)
            for(int r=0; r<c.nx; ++r)
            {
                face(2,EXT+r,jl,EXT+r,jl+1,EXT+r,jl);
                face(2,EXT+r,jh,EXT+r,jh,EXT+r,jh+1);
            }
        }

        for(int r=0; r<c.ny; ++r)
        {
            c.rec[ipol][0][r] = Fx(il,EY+r);
            c.rec[ipol][1][r] = Fx(ih,EY+r);
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
                c.rec[0][0][r] = b->dfx(il,EY+r);
                c.rec[0][1][r] = b->dfx(ih,EY+r);
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
                double df = 0.5*(c.rec[0][m.side][m.r]+c.rec[0][m.side][m.r+fr1(m.side)]);
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

        // the same faces around the box of this grid (q.ring): the fine face depth, no register
        if(ipol==4 && id<(int)xmatch.size())
        {
            for(auto &m : xmatch[id])
            {
                sflow_amr_patch &c = *SP(m.child);
                const double df = 0.5*(c.rec[0][m.side][m.r]+c.rec[0][m.side][m.r+fr1(m.side)]);
                if(m.dir==0)
                b->dfx(m.fi,m.fj) = df;
                else
                b->dfy(m.fi,m.fj) = df;
            }
            for(auto &m : xrmatch[id])
            {
                if(m.dir==0)
                b->dfx(m.fi,m.fj) = m.val[0];
                else
                b->dfy(m.fi,m.fj) = m.val[0];
            }
        }
        return;
    }

    // coarse faces next to the patches of this rank: mean of the two fine faces (multi-block: with
    // the sign of the coarse face, the face depth without)
    for(auto &m : match[id])
    {
        sflow_amr_patch &c = *SP(m.child);
        double val = 0.5*(c.rec[ipol][m.side][m.r]+c.rec[ipol][m.side][m.r+fr1(m.side)]);
        if(m.sgn<0)
        val = -val;

        if(m.dir==0)
        Fx(m.fi,m.fj) = val;
        else
        Fy(m.fi,m.fj) = val;

        if(ipol==4)
        {
            double df = 0.5*(c.rec[0][m.side][m.r]+c.rec[0][m.side][m.r+fr1(m.side)]);
            if(m.dir==0)
            b->dfx(m.fi,m.fj) = df;
            else
            b->dfy(m.fi,m.fj) = df;
        }
    }

    // coarse faces next to patches of other ranks
    for(auto &m : rmatch[id])
    {
        const double val = (m.sgn<0) ? -m.val[ipol] : m.val[ipol];
        if(m.dir==0)
        Fx(m.fi,m.fj) = val;
        else
        Fy(m.fi,m.fj) = val;

        if(ipol==4)
        {
            if(m.dir==0)
            b->dfx(m.fi,m.fj) = m.val[0];
            else
            b->dfy(m.fi,m.fj) = m.val[0];
        }
    }

    // the same coarse faces where this grid computes them around its box (q.ring)
    if(id<(int)xmatch.size())
    {
        for(auto &m : xmatch[id])
        {
            sflow_amr_patch &c = *SP(m.child);
            const double val = 0.5*(c.rec[ipol][m.side][m.r]+c.rec[ipol][m.side][m.r+fr1(m.side)]);
            if(m.dir==0)
            Fx(m.fi,m.fj) = val;
            else
            Fy(m.fi,m.fj) = val;

            if(ipol==4)
            {
                const double df = 0.5*(c.rec[0][m.side][m.r]+c.rec[0][m.side][m.r+fr1(m.side)]);
                if(m.dir==0)
                b->dfx(m.fi,m.fj) = df;
                else
                b->dfy(m.fi,m.fj) = df;
            }
        }

        for(auto &m : xrmatch[id])
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
}

// --------------------------------------------------------------------- time step
void sflow_amr::timestep(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(!active())
    return;

    // semi-implicit free surface: the advective limit over the cells of the blocks and the patches (sflow_amr_si.cpp)
    if(si==1)
    {
        if(p->count>0)
        si_timestep(p,pgc);
        return;
    }

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
        for(int jj=EY; jj<EY+c->ny; ++jj)
        if(pp->wet[lij(pp,ii,jj)]==1 && pp->flagslice4[lij(pp,ii,jj)]>0)
        {
            double cc = sqrt(g*MAX(pb->WL(ii,jj),p->A244));
            double sigma;

            if(pp->geo_curv==1)
            sigma = sflow_curv::sigma(pp,pb,ii,jj,cc);

            else
            {
            sigma = (fabs(pb->U(ii,jj))+cc)/pp->DXN[ii+marge];

            if(p->j_dir==1)
            sigma += (fabs(pb->V(ii,jj))+cc)/pp->DYN[jj+marge];
            }

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
            for(int jj=EY; jj<EY+c->ny; ++jj)
            if(pp->wet[lij(pp,ii,jj)]==1 && pp->flagslice4[lij(pp,ii,jj)]>0)
            {
                double hh = MAX(pb->WL(ii,jj),p->A244);
                double dx = pp->DXN[ii+marge], dy = pp->DYN[jj+marge];
                const double i2 = (pp->geo_curv==1) ? sflow_curv::inv2(pp,ii,jj) : 1.0/(dx*dx) + p->y_dir/(dy*dy);
                const double dv = 1.2/(i2*sqrt(g*B*hh*hh*hh));
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
    if(nblk==0 || nring[IJ]==0)
    m += b->WL(i,j)*((p->geo_curv==1) ? p->AREA2D[IJ] : p->DXN[IP]*p->DYN[JP]);

    // multi-block: the own cells of the blocks 1.. (level 0)
    for(int n : lev[0])
    {
        sflow_amr_patch *c = SP(n);
        lexer *pp = c->pp;
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EY; jj<EY+c->ny; ++jj)
        m += c->b->WL(ii,jj)*pp->AREA2D[lij(pp,ii,jj)];
    }

    return pgc->globalsum(m);
}

void sflow_amr::print(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(!active())
    return;

    dump_grids(p);

    bool doprint=false;

    if((p->count%p->P181==0 && p->P182<0.0 && p->P10>0 && p->P181>0) || (p->count==0 && p->P182<0.0))
    doprint=true;

    if((p->simtime>printtime_amr && p->P182>0.0 && p->P10>0) || (p->count==0 && p->P182>0.0))
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
    MPI_Allgather(&np,1,MPI_INT,&all[0],1,MPI_INT,reef3d_comm());

    if(p->mpirank==0)
    {
        char name[256];
        snprintf(name,sizeof(name),"./REEF3D_SFLOW_AMR/REEF3D-SFLOW-AMR-%08i.vtm",printcount_amr);
        ofstream out(name);
        out<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">\n<vtkMultiBlockDataSet>\n";
        out<<"<Block index=\"0\" name=\"level 0\">\n";
        for(int r=0; r<p->mpi_size; ++r)
        out<<"<DataSet index=\""<<r<<"\" file=\"REEF3D-SFLOW-AMR-L0-"<<setw(8)<<setfill('0')<<printcount_amr<<"-"<<setw(4)<<r+1<<(p->geo_curv==1?".vts":".vtr")<<"\"/>\n";
        out<<setfill(' ')<<"</Block>\n<Block index=\"1\" name=\"patches\">\n";
        int idx=0;
        for(int r=0; r<p->mpi_size; ++r)
        for(int q=0; q<all[r]; ++q)
        {
        out<<"<DataSet index=\""<<idx<<"\" file=\"REEF3D-SFLOW-AMR-"<<setw(8)<<setfill('0')<<printcount_amr<<"-"<<setw(4)<<r+1<<"-"<<setw(4)<<q+1<<(p->geo_curv==1?".vts":".vtr")<<"\"/>\n";
        out<<setfill(' ');
        ++idx;
        }
        out<<"</Block>\n</vtkMultiBlockDataSet>\n</VTKFile>\n";
        out.close();
    }

    ++printcount_amr;
}

// surface elevation at the P 51 gauges from the finest grid that holds them (interpolated between
// the cell centres, as the WSF file)
void sflow_amr::gauges(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const int ng = p->P51;
    vector<double> lv(ng,-1.0), val(ng,0.0);

    // curvilinear grid: the cell whose quadrilateral holds the point (the ring of a multi-block
    // native grid excluded)
    auto inquad = [&](lexer *q, int a, int d, double x, double y)
    {
        const int o[5][2] = {{0,0},{1,0},{1,1},{0,1},{0,0}};
        for(int k=0; k<4; ++k)
        {
            const int n0 = q->nij(a+o[k][0],d+o[k][1]), n1 = q->nij(a+o[k+1][0],d+o[k+1][1]);
            const double ex = q->XN2D[n1]-q->XN2D[n0], ey = q->YN2D[n1]-q->YN2D[n0];
            if(ex*(y-q->YN2D[n0]) - ey*(x-q->XN2D[n0]) < 0.0)
            return false;
        }
        return true;
    };

    auto find = [&](lexer *q, int i0, int i1, int j0, int j1, double x, double y, int &ii, int &jj)
    {
        ii=-1; jj=-1;
        if(q->geo_curv==1 && q->XN2D!=nullptr)
        {
            for(int a=i0; a<=i1; ++a)
            for(int d=j0; d<=j1; ++d)
            {
                if(q==p0 && nblk>0 && nring[lij(q,a,d)]==1)
                continue;
                if(inquad(q,a,d,x,y))
                {
                    ii=a; jj=d;
                    return true;
                }
            }
            return false;
        }
        for(int a=i0; a<=i1; ++a)
        if(x>=q->XN[a+marge] && x<q->XN[a+1+marge]) { ii=a; break; }
        if(sy==0)
        jj=j0;      // 2D grid: one cell in y
        else
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
            val[k]=p->ccslipol4(b->eta,p->P51_x[k],p->P51_y[k]);
        }
        for(auto c : P)
        if(c->lev>lv[k] && find(c->pp,EXT,EXT+c->nx-1,EY,EY+c->ny-1,p->P51_x[k],p->P51_y[k],ii,jj))
        {
            lv[k]=c->lev;
            val[k]=c->pp->ccslipol4(SP(c)->b->eta,p->P51_x[k],p->P51_y[k]);
        }
    }

    vector<double> lmax(ng);
    MPI_Allreduce(&lv[0],&lmax[0],ng,MPI_DOUBLE,MPI_MAX,reef3d_comm());
    for(int k=0; k<ng; ++k)
    if(lv[k]!=lmax[k])
    val[k]=0.0;
    vector<double> vs(ng);
    MPI_Allreduce(&val[0],&vs[0],ng,MPI_DOUBLE,MPI_SUM,reef3d_comm());

    if(p->mpirank==0)
    {
        if(!gaugeout.is_open())
        {
            if(reefckpt_open(gaugeout,"./REEF3D_SFLOW_AMR/REEF3D_SFLOW_AMR_gauges.dat",p->I45))   // I 45: continued
            {
            gaugeout<<"# simtime";
            for(int k=0; k<ng; ++k)
            gaugeout<<" \t eta("<<p->P51_x[k]<<","<<p->P51_y[k]<<")";
            gaugeout<<" \t levels"<<endl;
            }
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
    // curvilinear grid (G 50): structured grids (.vts) on the physical nodes
    const bool cv = (p->geo_curv==1);
    snprintf(name,sizeof(name),"./REEF3D_SFLOW_AMR/REEF3D-SFLOW-AMR-L0-%08i-%04i.%s",printcount_amr,p->mpirank+1,cv?"vts":"vtr");

    const int nx=p->knox, ny=p->knoy, m=marge;
    ofstream out(name);
    out<<"<?xml version=\"1.0\"?>\n";
    out<<"<VTKFile type=\""<<(cv?"StructuredGrid":"RectilinearGrid")<<"\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out<<"<"<<(cv?"StructuredGrid":"RectilinearGrid")<<" WholeExtent=\"0 "<<nx<<" 0 "<<ny<<" 0 0\">\n";
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
    if(cv)
    {
    out<<"</DataArray>\n</CellData>\n<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int jj=0; jj<=ny; ++jj)
    for(int ii=0; ii<=nx; ++ii)
    out<<setprecision(12)<<p->XN2D[p->nij(ii,jj)]<<" "<<p->YN2D[p->nij(ii,jj)]<<" 0\n";
    out<<"</DataArray>\n</Points>\n</Piece>\n</StructuredGrid>\n</VTKFile>\n";
    out.close();
    return;
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
    const bool cv = (p->geo_curv==1);
    snprintf(name,sizeof(name),"./REEF3D_SFLOW_AMR/REEF3D-SFLOW-AMR-%08i-%04i-%04i.%s",printcount_amr,p->mpirank+1,id+1,cv?"vts":"vtr");

    ofstream out(name);
    out<<"<?xml version=\"1.0\"?>\n";
    out<<"<VTKFile type=\""<<(cv?"StructuredGrid":"RectilinearGrid")<<"\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out<<"<"<<(cv?"StructuredGrid":"RectilinearGrid")<<" WholeExtent=\"0 "<<c.nx<<" 0 "<<c.ny<<" 0 0\">\n";
    out<<"<FieldData><DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\">"<<p->simtime<<"</DataArray>";
    out<<"<DataArray type=\"Int32\" Name=\"level\" NumberOfTuples=\"1\">"<<c.lev<<"</DataArray></FieldData>\n";
    out<<"<Piece Extent=\"0 "<<c.nx<<" 0 "<<c.ny<<" 0 0\">\n";
    out<<"<CellData Scalars=\"eta\">\n";

    auto field = [&](const char *nm, slice &f, double shift)
    {
        out<<"<DataArray type=\"Float64\" Name=\""<<nm<<"\" format=\"ascii\">\n";
        for(int jj=EY; jj<EY+c.ny; ++jj)
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
    for(int jj=EY; jj<EY+c.ny; ++jj)
    {
        for(int ii=EXT; ii<EXT+c.nx; ++ii)
        out<<pp->wet[lij(pp,ii,jj)]<<" ";
        out<<"\n";
    }
    out<<"</DataArray>\n";

    if(cv)
    {
    out<<"</CellData>\n<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int jj=EY; jj<=EY+c.ny; ++jj)
    for(int ii=EXT; ii<=EXT+c.nx; ++ii)
    out<<setprecision(12)<<pp->XN2D[pp->nij(ii,jj)]<<" "<<pp->YN2D[pp->nij(ii,jj)]<<" 0\n";
    out<<"</DataArray>\n</Points>\n</Piece>\n</StructuredGrid>\n</VTKFile>\n";
    out.close();
    return;
    }

    out<<"</CellData>\n<Coordinates>\n";

    const int m = marge;
    out<<"<DataArray type=\"Float64\" Name=\"x\" format=\"ascii\">";
    for(int ii=EXT; ii<=EXT+c.nx; ++ii) out<<setprecision(12)<<pp->XN[ii+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"y\" format=\"ascii\">";
    for(int jj=EY; jj<=EY+c.ny; ++jj) out<<setprecision(12)<<pp->YN[jj+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"z\" format=\"ascii\">0</DataArray>\n";
    out<<"</Coordinates>\n</Piece>\n</RectilinearGrid>\n</VTKFile>\n";
    out.close();
}

// P 18: this rank's grids as its .vtr files have them (the fields of write_vtr0 and
// write_vtr, wetdry as int32), gathered into the LAGOON store; true when the output is there
bool sflow_amr::print_lagoon(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // curvilinear grid: the LAGOON grids are rectilinear (x and y lines); the .vts files instead
    if(p->geo_curv==1)
    {
        static bool told_curv = false;
        if(p->mpirank==0 && !told_curv)
        cout<<"LAGOON: P 18 not yet for curvilinear grids (G 50), the .vts files are written"<<endl;
        told_curv = true;
        return false;
    }

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
    take(SP(n)->pp,SP(n)->b,EXT,EY,SP(n)->nx,SP(n)->ny,SP(n)->lev);

    return writer->write(p,pgc,grids,printcount_amr);
}
