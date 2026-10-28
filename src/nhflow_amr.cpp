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

#include"reefckpt.h"
#include"nhflow_amr.h"
#include"nhflow_momentum_SI.h"
#include"nhflow_amr_fill.h"
#include"lexer.h"
#include"bc_zone.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice4.h"
#include"nhflow_momentum_RK2.h"
#include"nhflow_momentum_RK3.h"
#include"nhflow_fsf_body.h"
#include"nhflow_signal_speed.h"
#include"nhflow_reconstruct_hires.h"
#include"nhflow_reconstruct_weno.h"
#include"nhflow_HLL.h"
#include"nhflow_diff_void.h"
#include"nhflow_ediff.h"
#include"nhflow_idiff.h"
#include"bicgstab_ijk.h"
#include"nhflow_pjm.h"
#include"nhflow_pjm_corr.h"
#include"nhflow_pjm_hs.h"
#include"nhflow_komega_void.h"
#include"vrans_nhflow_v.h"
#include"nhflow_fsf_f.h"
#include"nhflow_forcing.h"
#include"nhflow_force.h"
#include"nhflow_amr_6dof.h"
#include"nhflow_hullbody.h"
#include"nhflow_thinbody.h"
#include"6DOF_nhflow.h"
#include"6DOF_obj.h"
#include"sediment_void.h"
#include"ioflow_void.h"
#include"geo_mesh.h"
#include"ioflow_amr_block.h"
#include"patchBC_void.h"
#include"reefmg_core.h"
#include"nhflow_curv.h"
#include"definitions.h"
#include<cmath>
#include<mpi.h>
#include<iomanip>
#include<algorithm>
#include<sys/stat.h>

using namespace nhflow_amr_detail;

namespace
{
typedef reefamr_comms_off comms_off;

inline double mmod(double a, double b)
{
    if(a*b<=0.0)
    return 0.0;
    return fabs(a)<fabs(b) ? a : b;
}
}

nhflow_amr::nhflow_amr(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_momentum *pmom, nhflow_convection *pconv, nhflow_timestep *pstep,
                       sixdof *p6dof)
                      : reefamr(p,pgc)
{
    d0 = d;
    mom0 = dynamic_cast<nhflow_momentum_func*>(pmom);
    pconv0 = pconv;
    pstep0 = pstep;

    // floating bodies on the hierarchy (X 10 1 two-way, X 10 2 prescribed motion)
    if(p->X10==1 || p->X10==2)
    b6 = dynamic_cast<sixdof_nhflow*>(p6dof);

    // fixed solids (direct forcing, A 581-590) on every grid
    solids = (p->A581>0 || p->A583>0 || p->A584>0 || p->A585>0 || p->A586>0 || p->A587>0 || p->A588>0 || p->A589>0 || p->A590>0);

    reefamr_param q;
    q.name = "NHFLOW AMR";
    q.maxlev = p->G1;
    q.regrid = 0;
    q.nbuf = MAX(p->G3,0);
    q.tile = MAX(p->G4,4);
    q.tile += q.tile%2;
    q.keep = MIN(MAX(p->G5,0),250);

    // cells computed beyond the patch box: the WENO5 reconstruction reaches three cells; the patch
    // arrays reach EXT + margin fine cells beyond the patch, which the coarser level has to cover
    q.ext = 3;
    q.nest = MAX(3,(q.ext+p->margin+1)/2);

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

    // refinement around the floating body: margin G 12 around the wetted hull, a rectangle
    // aligned with x and y (moored or oscillating bodies), built at t = 0
    q.zones = (b6!=nullptr && p->G12>0);
    q.zr = p->G12_r;
    q.zalign = true;

    // optional wake wedge (G 13 L a, as SFLOW): the zone is oriented along the direction of
    // motion (the yaw angle at rest) and a wedge of half angle a reaches back from the bow by L
    // at most, where the bow has been; without G 13 the rectangle aligned with x and y
    if(p->G13_L>0.0)
    {
    q.zalign = false;
    q.zL = p->G13_L;
    q.za = p->G13_a;
    }

    // solution-adaptive flags (as sflow_amr): a cell of level l gets level l+1 where the surface
    // jumps by more than G 20 to a neighbour cell or its second difference along x or y exceeds
    // G 21 (both in m, so a finer level flags only steeper or shorter waves)
    tol_eta = MAX(p->G20,0.0);
    tol_curv = MAX(p->G21,0.0);
    adaptive = (tol_eta>0.0 || tol_curv>0.0);

    // wetting and drying: with G 30 1 the patches may cover dry and shallow cells (the NHFLOW
    // wetting and drying runs on every grid, the coupling interpolates from wet cells only);
    // otherwise they are fully wet: no patch within 4 level-0 cells of a dry or shallow cell,
    // checked at every regrid (cell_unfit; a static box over dry cells loses those tiles as well)
    shore = (p->G30==1);
    q.dryband = shore ? 0 : 4;

    // shoreline flag (G 22 n, with G 30 1): a cell flags the next level where a cell of the other
    // wet state lies within n cells (n <= 3: the halo of level 0 and the EXT cells of a patch),
    // so the patches follow run-up and run-down
    if(shore && p->G22>0)
    {
    nshore = MIN(p->G22,3);
    adaptive = true;
    }

    // breaking flag (G 23 1, with A 550 1): a cell with breaking viscosity flags the next level,
    // so the breaking zone stays on the finest grid as it moves
    if(p->G23==1 && p->A550==1)
    {
    flagbreak = true;
    adaptive = true;
    }

    // the zone follows the body, the patches follow the flags: new patches every G 2 steps
    // (G 2 0: static), the layout is kept as long as it covers the flagged tiles with at most
    // 50 % excess (as fnpf_amr)
    if(q.zones || adaptive)
    {
    q.regrid = MAX(p->G2,0);
    q.lazy = 1.5;
    }

    // vertical refinement (G 6 1): every level doubles the sigma layers of its parent, the
    // coarse nodes are nodes of the fine grid
    vr = (p->G6==1) ? 2 : 1;
    q.vref.assign(q.maxlev+1,vr);

    // several ranks (G 40 1): the patches are not cut at the rank boxes but placed for the load of
    // the ranks; parents and old patches on other ranks are reached through the block plans and
    // old_run of the core
    q.place = (p->G40>=1 && p->mpi_size>1) ? MIN(p->G40,2) : 0;
    q.rebalance = 0.1;

    // implicit diffusion (A 512 2): one composite solve over the leaf cells of all grids (G 31 1,
    // nhflow_amr_diff.cpp), or every grid its own (G 31 0)
    // (a turbulence model with A 512 0 diffuses implicitly as A 512 2, driver_logic_nhflow)
    cdiff = ((p->A512==2 || (p->A512==0 && p->A560>0)) && (p->G31==1 || p->geo_mb>0));
    brk_split = (cdiff && p->A550==1);
    brk0 = (p->geo_mb>0 && p->A550==1);

    // subcycling (G 7 1, nhflow_amr_sub.cpp): floating bodies with the zone around the hull (G 12:
    // the finest level advances the body); every grid solves its own implicit diffusion (G 31 1
    // needs all grids at the same time)
    sub = (p->G7==1 && (b6==nullptr || p->G12>0)) ? 1 : 0;

    // the hull as a sharp thin body (X 15 8, nhflow_hullbody): every grid holds the hull (blocked links, wall
    // velocity, the thin-body service of nhflow_thinbody), the loads are summed over the leaf grids.  The hull
    // with its band (the local fluxes within two cells, the footprint of the lower cells, the lid correction)
    // lies inside the finest level: the zone around the hull (G 12) with a margin of at least 4 level-0 cells,
    // one time step for all levels (G 7 0), on several ranks the zone placed whole (G 40 1)
    hb = (b6!=nullptr && p->X15==8);
    if(hb)
    sub = 0;
    if(hb)
    q.place_whole_zones = true;

    if(sub==1)
    cdiff = brk_split = false;

    // multi-block breaking: block 0 runs the parts of the breaking with the blocks (breaking_level0)
    if(brk0)
    if(nhflow_breaking *b = dynamic_cast<nhflow_breaking*>(mom0))
    b->brk_defer = true;

    // G 7 1 with bodies: the coarse levels step first with their own coarse picture of the body,
    // and its near field reaches the finest level through the fills and the edge pressure of the
    // level solves.  The zone margin is at least 12 cells of level 0: the towed box of 0020 with
    // two and three levels then within the error of one step for all levels (0.2 m = 4 cells:
    // 4.5x that error).  Following the finer levels in the covered cells during a coarse step was
    // tried and brought little (0.39 -> 0.35 N, the heave box worse)
    zr_user = q.zr;
    if(sub==1 && q.zones)
    q.zr = MAX(q.zr, 12.0*p->DXM);
    if(hb && q.zones)
    q.zr = MAX(q.zr, 4.0*p->DXM);

    // multi-block grid: passes of the block-local V-cycles of level 0 in the preconditioner of the
    // composite pressure (G 51), the corrections of the neighbour blocks filled in between
    npass0 = MAX(p->G51,1);
    csw = MAX(p->G53,0);

    configure(q);

    // the hull triangles (X 185) are sized for the finest level (before the 6DOF initialisation builds
    // them), as on the uniform fine grid: horizontally (sixdof_obj::amr_hfac) whenever the refinement stays
    // on, since every triangle is sampled on the finest grid at its centroid, whether the zone (G 12), a box
    // (G 10) or the adaptive flags put a patch over the hull; vertically (amr_vfac) with the sigma layers
    // doubled on every level (G 6 1), as fnpf_amr
    if(b6!=nullptr && p->G1>0 && keys_ok(p)==1)
    for(int nb=0; nb<b6->objects(); ++nb)
    {
    b6->object(nb)->amr_hfac = 1.0/double(1<<p->G1);
    b6->object(nb)->amr_vfac = (vr==2) ? 1.0/double(1<<p->G1) : 1.0;
    }

    if(p->F50==1) gcval_eta = 51;
    if(p->F50==2) gcval_eta = 52;
    if(p->F50==3) gcval_eta = 53;
    if(p->F50==4) gcval_eta = 54;

    bc5 = 21;
    bc6 = 3;

    printtime_amr = 0.0;
    printcount_amr = 0;
    for(int k=0; k<8; ++k)
    tm[k] = 0.0;

    pflowv = nullptr;
    pBCv = nullptr;
    psedv = nullptr;
}

nhflow_amr::~nhflow_amr()
{
    free_patches();
    df_free();

    for(auto v : kv0)
    delete [] v;
    delete mg0;
    cs_free();
}

nhflow_amr::pscope::pscope(ghostcell *gg, fdm_nhf *dp, fdm_nhf *d00) : g(gg), dl0(d00)
{
    old = g->set_comms(false);
    oldlocal = g->set_local(true);
    g->fdm_nhf_update(dp);
}

nhflow_amr::pscope::~pscope()
{
    g->fdm_nhf_update(dl0);
    g->set_local(oldlocal);
    g->set_comms(old);
}

// --------------------------------------------------------------------- grids
// stage input of grid g at stage s (s<=0: the state at the start of the step)
nhflow_amr::stg nhflow_amr::stage_in(int g, int s)
{
    fdm_nhf *d = gfd(g);
    if(s<=0)
    return {&d->WL,d->UH,d->VH,d->WH};

    nhflow_momentum_func *m = gmom(g);
    return {&m->stage_WL(d,s-1),m->stage_UH(d,s-1,0),m->stage_UH(d,s-1,1),m->stage_UH(d,s-1,2)};
}

// stage output of grid g at stage s
nhflow_amr::stg nhflow_amr::stage_out(int g, int s)
{
    fdm_nhf *d = gfd(g);
    nhflow_momentum_func *m = gmom(g);
    return {&m->stage_WL(d,s),m->stage_UH(d,s,0),m->stage_UH(d,s,1),m->stage_UH(d,s,2)};
}

// the keys the refinement runs with: 1 on, 0 off, -1 off for the sharp hull (X 15 8). The constructor sizes the
// hull triangles with it, before ini switches the refinement off with a message
int nhflow_amr::keys_ok(lexer *p)
{
    int ok=1;
    if(mom0==nullptr) ok=0;
    if(p->A510!=2 && p->A510!=3) ok=0;
    if(p->A511!=1 && p->A511!=2) ok=0;
    if(p->A520<0 || p->A520>2) ok=0;
    if(p->A512<0 || p->A512>2) ok=0;
    // turbulence (nhflow_amr_turb.cpp): k-epsilon and k-omega on every grid, one time step for all levels, without
    // the HLES split and QCR
    if(p->A560!=0 && p->A560!=1 && p->A560!=21 && p->A560!=2 && p->A560!=22) ok=0;
    if(p->A560!=0 && (p->G7==1 || p->A561!=0 || p->A563!=0.0)) ok=0;
    // the side-wall friction and wall functions (A 519 2) need the side walls in the boundary lists of the grids; the
    // patches and the blocks have the bed and the surface only
    if(p->A519==2) ok=0;
    // breaking (A 550 1, RK2 and RK3): it acts through the implicit diffusion (A 512 2), composite (G 31 1)
    // or per grid (G 31 0); on a multi-block grid the blocks break as one grid (breaking_level0). The eddy viscosity
    // (A 550 3) is not on the patches
    if(p->A550!=0 && (p->A550!=1 || (p->A510!=2 && p->A510!=3))) ok=0;
    if(p->B200!=0 || p->X330!=0 || p->A599!=0) ok=0;
    // sediment transport (S 10 1): bed load, Exner (S 31 1/2, S 32 1/2), sand slide; no suspended load, no
    // sediment fractions (nhflow_amr_sed.cpp)
    if(p->S10!=0 && (p->S10!=1 || p->S12!=0 || p->S51>0 || p->S31==3 || p->S32>=3)) ok=0;
    if(p->X10!=0 && (b6==nullptr || p->X60!=1 || p->X16!=0 || p->X320!=0 || p->A516==2 || p->A516==4)) ok=0;
    // fixed solids (A 581-590): every grid forces its own (solid_patch); not with subcycling, not with the
    // forcing of the water level at the solids (A 516 2/4: the slice lists of level 0)
    if(solids && (p->G7==1 || p->A516==2 || p->A516==4)) ok=0;
    if(p->A580==1) ok=0;
    if(p->E10>0 || p->L10>0 || p->Z20>0 || p->Z30>0) ok=0;
    if(p->j_dir!=1 && (p->geo_mb>0 || p->geo_curv==1)) ok=0;     // 2D grid (one cell in y): Cartesian

    // X 15 8 (hb): one body with the zone around it (G 12), the full projection (A 520 1), no multi-block grid, on
    // several ranks the placed patches (G 40 1, the zone whole on one rank)
    if(ok==1 && hb && (p->G12<=0 || b6->objects()!=1 || p->A520!=1 || p->geo_mb>0 || (p->mpi_size>1 && p->G40<1)))
    ok=-1;

    return ok;
}

// --------------------------------------------------------------------- setup
void nhflow_amr::ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(maxlev<1 && p->geo_mb==0)
    return;

    // multi-block grid (G 50 1, DIVEMesh R 20 1): the blocks run on the hierarchy, with or without
    // refinement
    if(p->geo_mb>0)
    blocks_check(p);

    // the keys (keys_ok): 0 refinement off, -1 off for the sharp hull (X 15 8)
    int ok = keys_ok(p);

    if(ok==0 && p->geo_mb>0)
    {
        if(p->mpirank==0)
        cout<<"!!! NHFLOW multi-block grid: only for A 510 2/3, A 511 1/2, A 520 0/1/2, A 512 0/1/2, A 560 0/1/2/21/22 (A 561 0, A 563 0), A 519 0/1, A 550 0/1, B 200 0, S 10 0/1 (S 12 0, S 51 0, S 31 1/2, S 32 1/2), no solids, membranes, nets, DEM, particles, rods or FEM structures !!!"<<endl;
        exit(1);
    }

    if(ok==-1)
    {
        if(p->mpirank==0)
        cout<<"NHFLOW AMR (G 1) with the sharp hull (X 15 8): only for one body with the zone around it (G 12), A 520 1, "
            <<"Cartesian grids, and on several ranks with placed patches (G 40 1) -- refinement switched off"<<endl;
        maxlev=0;
        return;
    }

    if(ok==0)
    {
        if(p->mpirank==0)
        cout<<"NHFLOW AMR (G 1): only for A 510 2/3, A 511 1/2, A 520 0/1/2, A 512 0/1/2, A 560 0/1/2/21/22 (not with G 7 1, A 561, A 563), A 519 0/1, A 550 0/1, B 200 0, S 10 0/1 (S 12 0, S 51 0, S 31 1/2, S 32 1/2), "
            <<"X 10 0/1/2 (X 60 1, X 16 0, A 516 0/1/3), fixed solids A 581-590 (not with G 7 1, A 516 2/4), no solids of the grid file (A 580 1), membranes, nets, DEM, particles, rods or FEM structures, and on a 2D grid Cartesian "
            <<"-- refinement switched off"<<endl;
        maxlev=0;
        return;
    }

    turb = (p->A560!=0) ? 1 : 0;
    if(turb && tb0==nullptr)
    {
        if(p->mpirank==0)
        cout<<"!!! NHFLOW AMR: the turbulence of level 0 is not attached (nhflow_amr::turb_attach) !!!"<<endl;
        exit(1);
    }

    // boundary types of the bed and the free surface, as level 0 has them
    for(int n=0; n<p->gcb4_count; ++n)
    {
        if(p->gcb4[n][3]==5) bc5 = p->gcb4[n][4];
        if(p->gcb4[n][3]==6) bc6 = p->gcb4[n][4];
    }
    bc5 = (int)pgc->globalmax(double(bc5));
    bc6 = (int)pgc->globalmax(double(bc6));

    // objects shared by all patches (their NHFLOW calls do nothing)
    pBCv = new patchBC_void(p);
    pflowv = new ioflow_v(p,pgc,pBCv);
    psedv = new sediment_void();

    // face entries of the flux exchange: 4 variables per layer of the coarse side, and dfx
    NF = 4*klev(MAX(maxlev-1,0)) + 1;

    // fixed solids: the refinement boxes around them (G 12)
    if(solids)
    solid_zone(p,pgc);

    setup(p,pgc);

    // multi-block grid: blocks 1.. as level-0 patches, their initial state, the fill plans, the bed
    // and the state of the cells around the blocks and of the ring
    if(p->geo_mb>0)
    {
        blocks_setup(p,pgc);
        blocks_ini(p,pgc);
        blocks_plans(pgc);
        fill_bed0(pgc);
        fill_level0(pgc,-1);
        blocks_sigma(pgc);
        fill_level0(pgc,-1);

        // A 560: the turbulence of the blocks
        turb_blocks_ini(pgc);

        // S 10: the sediment of the blocks starts from their bed
        for(int id : lev[0])
        sed_patch_ini(pgc,*NP(id));

        // I 11: the potential flow over the blocks (nhflow_potential_f::start leaves it to this)
        potential(p,pgc);
    }

    mkdir("./REEF3D_NHFLOW_AMR",0777);

    if(p->mpirank==0)
    {
        if(reefckpt_open(logout,"./REEF3D_NHFLOW_AMR/REEF3D_NHFLOW_AMR_log.dat",p->I45))   // I 45: continued
        logout<<"# count \t simtime \t dt \t patches \t cells \t pressure iterations \t residual \t water volume \t relative change \t layout"<<endl;
    }

    for(int it=0; it<maxlev; ++it)
    regrid(p,pgc,true);

    if(nblk>0 && maxlev>=1)
    fill_level0(pgc,-1);

    // the level-0 objects hand their time step, fluxes and time step limits to the hierarchy
    mom0->attach_runner(this);
    pconv0->set_hook(this,0);
    pstep0->phook = this;

    // initial time step with the patch cells
    pstep0->ini(p,d,pgc);

    // floating bodies: the loads of every hull triangle from the finest grid at its centroid,
    // the initial loads with the patches
    if(b6!=nullptr)
    {
        for(int nb=0; nb<b6->objects(); ++nb)
        b6->object(nb)->amr_grid_nhflow = [this](double x, double y)
        {
            const int g = finest_at(x,y);
            if(g<0)
            return sixdof_obj_nhflow::nhflow_grid{glex(-1),d0,(cur_stage<0) ? &d0->WL : stage_out(-1,cur_stage).WL};
            sixdof_obj_nhflow::nhflow_grid r{glex(g),NP(g)->d,(cur_stage<0) ? &NP(g)->d->WL : stage_out(g,cur_stage).WL};
            // multi-block grid: the blocks have their own logical spacing, the band by the level
            if(nblk>0)
            r.hr = 1.0/double(1<<P[g]->lev);
            return r;
        };

        // placed patches: the rank of the finest grid at the centroid takes the triangle
        if(par.place>0 && nblk==0)
        for(int nb=0; nb<b6->objects(); ++nb)
        b6->object(nb)->amr_owner_nhflow = [this](double x, double y)
        {
            return finest_rank(x,y)==p0->mpirank;
        };

        // multi-block grid (P5c): the rank of the native cell that holds the centroid (not the ring, it stands for
        // the cells of the neighbour blocks) or of the piece of a block 1.. that holds it; the patches of the
        // blocks are on the rank of the block piece below (G 40 0)
        if(nblk>0)
        for(int nb=0; nb<b6->objects(); ++nb)
        b6->object(nb)->amr_point_grid = true;

        if(nblk>0)
        for(int nb=0; nb<b6->objects(); ++nb)
        b6->object(nb)->amr_owner_nhflow = [this](double x, double y)
        {
            int a,b;
            if(p0->geo_cell(x,y,a,b) && !native_ring(a+O0i,b+O0j))
            return true;
            for(reefamr_patch *c : P)
            if(c->lev==0 && c->blk>0)
            if(c->pp->geo_cell(x,y,a,b) && a>=EXT && a<EXT+c->nx && b>=EXT && b<EXT+c->ny)
            return true;
            return false;
        };

        // X 15 8: the hull loads summed over the leaf part of every grid of this rank
        if(hb)
        b6->object(0)->amr_grids_nhflow = [this](const std::function<void(const sixdof_obj_nhflow::nhflow_grid&)> &f)
        {
            f(sixdof_obj_nhflow::nhflow_grid{glex(-1),d0,(cur_stage<0) ? &d0->WL : stage_out(-1,cur_stage).WL});
            for(int g=0; g<(int)P.size(); ++g)
            f(sixdof_obj_nhflow::nhflow_grid{glex(g),NP(g)->d,(cur_stage<0) ? &NP(g)->d->WL : stage_out(g,cur_stage).WL});
        };

        cur_stage = -1;
        body_loads(p,pgc);
    }

    m0 = mass(p,d,pgc);

    if(p->mpirank==0)
    {
    if(nblk>0)
    cout<<"NHFLOW multi-block grid: "<<nblk<<" blocks"<<(maxlev>=1 ? "; " : "\n");
    if(maxlev>=1)
    {
    cout<<"NHFLOW AMR: "<<maxlev<<" level(s), "<<patches_total<<" patch(es), "<<cells_total<<" refined columns, dt "<<p->dt;
    if(vr==2)
    cout<<", sigma layers doubled on every level (G 6)";
    if(shore)
    cout<<", wetting and drying in the patches (G 30)";
    if(adaptive)
    {
    cout<<", flags:";
    if(tol_eta>0.0)
    cout<<" surface jump "<<tol_eta<<" m (G 20)";
    if(tol_curv>0.0)
    cout<<" second difference "<<tol_curv<<" m (G 21)";
    if(nshore>0)
    cout<<" shoreline within "<<nshore<<" cells (G 22)";
    if(flagbreak)
    cout<<" breaking (G 23)";
    }
    if(regrid_int>0 && sub==0)
    cout<<", regrid every "<<regrid_int<<" steps (G 2)"<<(par.zones ? ", the zone follows the body" : "");
    if(regrid_int>0 && sub==1)
    cout<<", regrid every "<<MAX(regrid_int>>maxlev,1)<<" level-0 steps (G 2)";
    if(sub==1)
    cout<<", subcycled (G 7 1: level l takes 2^l steps per level-0 step)";
    cout<<endl;
    if(p->G7==1 && sub==0 && !hb)
    cout<<"NHFLOW AMR: G 7 1 (subcycling) with floating bodies needs the zone around the hull (G 12) -- one time step for all levels"<<endl;
    if(hb)
    {
    cout<<"NHFLOW AMR: the sharp hull (X 15 8) on every grid, the loads from the finest grid at every blocked face";
    if(p->G7==1)
    cout<<"; G 7 1 not with X 15 8 -- one time step for all levels";
    if(par.zr>zr_user)
    cout<<"; zone margin "<<par.zr<<" m (4 cells of level 0) instead of G 12 "<<zr_user;
    cout<<endl;
    }
    if(sub==1 && par.zones && par.zr>zr_user)
    cout<<"NHFLOW AMR: G 7 1 with bodies: zone margin "<<par.zr<<" m (12 cells of level 0) instead of G 12 "<<zr_user<<endl;
    if(p->G7==1 && sub==1 && p->A512==2 && p->G31==1)
    cout<<"NHFLOW AMR: G 7 1: every grid solves its own implicit diffusion (G 31 1 needs all grids at the same time)"<<endl;
    }
    }
}

// --------------------------------------------------------------------- patches
reefamr_patch* nhflow_amr::patch_new()
{
    return new nhflow_amr_patch;
}

// 3D flags, boundary list and the sigma-grid arrays of the patch, as driver::makegrid_sigma and
// the DIVEMesh boundary list for level 0: a fluid column is fluid from k=0 to k=knoz-1 (cells)
// and from the bed node to the free-surface node (F layout); bed and free surface are the only
// boundaries in the vertical
void nhflow_amr::build_lexer3D(nhflow_amr_patch &c)
{
    lexer *pp = c.pp;
    const int n3 = pp->imax*pp->jmax*pp->kmax;
    const int n7 = pp->imax*pp->jmax*(pp->kmax+2);
    const int nsl = pp->imax*pp->jmax;

    pp->Iarray(pp->flag1,n3);
    pp->Iarray(pp->flag2,n3);
    pp->Iarray(pp->flag3,n3);
    pp->Iarray(pp->flag4,n3);
    pp->Iarray(pp->flag5,n3);
    pp->Iarray(pp->flag7,n7);
    pp->Iarray(pp->IO,n3);
    pp->Iarray(pp->DF,n3);
    pp->Iarray(pp->DFF,pp->imax*pp->jmax*(pp->kmax+1));
    pp->Darray(pp->ZSN,pp->imax*pp->jmax*(pp->kmax+1));
    pp->Darray(pp->ZSP,n3);
    pp->Darray(pp->bed,nsl);
    pp->Darray(pp->WL,nsl);
    // sediment (S 10): the erodible-bed flag of the columns (sediment_f::active_nhflow sets it)
    pp->Iarray(pp->DFBED,nsl);
    for(int q=0; q<nsl; ++q)
    pp->DFBED[q] = 1;
    // no in- and outflow cells of the 3D lists on a patch (sediment_f::ini_guard reads them)
    pp->gcin_count = pp->gcout_count = 0;

    for(int q=0; q<n7; ++q)
    pp->flag7[q] = -10;

    auto fl = [&](int ii, int jj)
    {
        if(ii<pp->imin || ii>=pp->imin+pp->imax || jj<pp->jmin || jj>=pp->jmin+pp->jmax)
        return false;
        return pp->flagslice4[lij(pp,ii,jj)]>0;
    };

    int nfluid=0;
    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        const bool f = fl(ii,jj);

        for(int kk=pp->kmin; kk<pp->kmin+pp->kmax; ++kk)
        {
            const int q = cidx(pp,ii,jj,kk);
            const bool in = f && kk>=0 && kk<pp->knoz;
            pp->flag4[q] = in ? WATER_FLAG : -10;
            pp->flag1[q] = pp->flag4[q];
            pp->flag2[q] = pp->flag4[q];
            pp->flag3[q] = pp->flag4[q];
            pp->flag5[q] = pp->flag4[q];
            pp->DF[q] = 1;
            pp->IO[q] = 0;

            // faces on a wall (as ghostcell::flagfield with the boundary list)
            if(in && !fl(ii+1,jj))
            pp->flag1[q] = OBJ_FLAG;
            if(in && !fl(ii,jj+1))
            pp->flag2[q] = OBJ_FLAG;
            if(in && kk==pp->knoz-1)
            pp->flag3[q] = OBJ_FLAG;
        }

        if(f)
        for(int kk=0; kk<=pp->knoz; ++kk)
        pp->flag7[fidx(pp,ii,jj,kk)] = WATER_FLAG;

        if(f && ii>=0 && ii<pp->knox && jj>=0 && jj<pp->knoy)
        ++nfluid;
    }

    for(int q=0; q<pp->imax*pp->jmax*(pp->kmax+1); ++q)
    pp->DFF[q] = 1;

    // boundary list: bed and free surface of every fluid column of the computed range
    pp->gcb4_count = 2*nfluid;

    // cells of the computed range (the residual of bicgstab_ijk for the implicit diffusion)
    pp->cellnumtot = nfluid*pp->knoz;
    pp->Iarray(pp->gcb4,MAX(pp->gcb4_count,1),6);
    int n=0;
    for(int ii=0; ii<pp->knox; ++ii)
    for(int jj=0; jj<pp->knoy; ++jj)
    if(fl(ii,jj))
    {
        pp->gcb4[n][0]=ii; pp->gcb4[n][1]=jj; pp->gcb4[n][2]=0; pp->gcb4[n][3]=5; pp->gcb4[n][4]=bc5; pp->gcb4[n][5]=0;
        ++n;
        pp->gcb4[n][0]=ii; pp->gcb4[n][1]=jj; pp->gcb4[n][2]=pp->knoz-1; pp->gcb4[n][3]=6; pp->gcb4[n][4]=bc6; pp->gcb4[n][5]=0;
        ++n;
    }

    // no partition neighbours in 3D (the ghostcell exchange of the V-type arrays runs over these lists)
    pp->gcpara1_count = pp->gcpara2_count = pp->gcpara3_count = pp->gcpara4_count = pp->gcpara5_count = pp->gcpara6_count = 0;
    pp->gcparaco1_count = pp->gcparaco2_count = pp->gcparaco3_count = pp->gcparaco4_count = pp->gcparaco5_count = pp->gcparaco6_count = 0;

    // wet and deep everywhere (patches are fully wet), the domain ends never lie on a patch
    for(int q=0; q<nsl; ++q)
    {
        pp->wet[q] = pp->wet_n[q] = pp->deep[q] = (pp->flagslice4[q]>0) ? 1 : 0;
    }
}

// the 3D in- and outflow cells of the open faces of a block (G 54): the fluid cells of its first (-i, gcin)
// and last (+i, gcout) column of cells, the ghost cells beyond them marked IO 1 / 2 (as ioflow_gcio_lists and
// ioflow_gcio_marks with the DIVEMesh boundary groups 1 and 2 on the native grid)
void nhflow_amr::open_lists(nhflow_amr_patch &c, bool in, bool out)
{
    lexer *pp = c.pp;
    const int i0 = EXT, i1 = EXT+blk[c.blk].ni-1;
    auto fluid = [&](int ii, int jj) { return pp->flagslice4[lij(pp,ii,jj)]>0; };

    int nin=0, nout=0;
    for(int jj=0; jj<pp->knoy; ++jj)
    {
        if(in && fluid(i0,jj) && !fluid(i0-1,jj))
        nin += pp->knoz;
        if(out && fluid(i1,jj) && !fluid(i1+1,jj))
        nout += pp->knoz;
    }

    pp->gcin_count = nin;
    pp->gcout_count = nout;
    pp->Iarray(pp->gcin,nin,6);
    pp->Iarray(pp->gcout,nout,6);

    int a=0, b=0;
    for(int jj=0; jj<pp->knoy; ++jj)
    for(int kk=0; kk<pp->knoz; ++kk)
    {
        if(in && fluid(i0,jj) && !fluid(i0-1,jj))
        {
            pp->gcin[a][0]=i0; pp->gcin[a][1]=jj; pp->gcin[a][2]=kk; pp->gcin[a][3]=1; pp->gcin[a][4]=1; pp->gcin[a][5]=0;
            pp->IO[cidx(pp,i0-1,jj,kk)] = 1;
            ++a;
        }
        if(out && fluid(i1,jj) && !fluid(i1+1,jj))
        {
            pp->gcout[b][0]=i1; pp->gcout[b][1]=jj; pp->gcout[b][2]=kk; pp->gcout[b][3]=4; pp->gcout[b][4]=2; pp->gcout[b][5]=0;
            pp->IO[cidx(pp,i1+1,jj,kk)] = 2;
            ++b;
        }
    }
}

void nhflow_amr::free_lexer3D(nhflow_amr_patch &c)
{
    lexer *pp = c.pp;
    const int n3 = pp->imax*pp->jmax*pp->kmax;
    const int n7 = pp->imax*pp->jmax*(pp->kmax+2);
    const int nsl = pp->imax*pp->jmax;

    pp->del_Iarray(pp->flag1,n3);
    pp->del_Iarray(pp->flag2,n3);
    pp->del_Iarray(pp->flag3,n3);
    pp->del_Iarray(pp->flag4,n3);
    pp->del_Iarray(pp->flag5,n3);
    pp->del_Iarray(pp->flag7,n7);
    pp->del_Iarray(pp->IO,n3);
    pp->del_Iarray(pp->DF,n3);
    pp->del_Iarray(pp->DFF,pp->imax*pp->jmax*(pp->kmax+1));
    pp->del_Darray(pp->ZSN,pp->imax*pp->jmax*(pp->kmax+1));
    pp->del_Darray(pp->ZSP,n3);
    pp->del_Darray(pp->bed,nsl);
    pp->del_Darray(pp->WL,nsl);
    pp->del_Iarray(pp->DFBED,nsl);
    pp->del_Iarray(pp->gcb4,MAX(pp->gcb4_count,1),6);
    pp->del_Darray(pp->sig,n7);
    pp->del_Darray(pp->sigx,n7);
    pp->del_Darray(pp->sigy,n7);
    pp->del_Darray(pp->sigz,nsl);
    pp->del_Darray(pp->sigt,n7);
    pp->del_Darray(pp->sigxx,n7);
}

void nhflow_amr::regrid_ids()
{
    for(int n=0; n<(int)P.size(); ++n)
    NP(n)->pconv->set_hook(this,n+1);
}

// NHFLOW objects of a new patch (comms off); bed and state are set by regrid_state
void nhflow_amr::patch_objects(reefamr_patch *q, ghostcell *pgc)
{
    nhflow_amr_patch *c = NP(q);
    lexer *pp = c->pp;

    // multi-block: the faces of the blocks 1.. without an interface are walls; the open faces of a block (G 54)
    // are its in- and outflow as B 60 1 on the native grid (only these i-faces: river_xm, river_xp; a block with
    // open j faces is turned, geo_mesh::turn_blocks), its edge zones (B 522) its open edges as the iowave edges of
    // the native grid (open_xm, open_xp), its relaxation and beach zones (B 522) those of ioflow_amr_block; the
    // zones of iowave (B 98, B 99) belong to block 0
    const int oin = (c->blk>0 && c->lev==0) ? bopen[0][c->blk] : -1;
    const int oout = (c->blk>0 && c->lev==0) ? bopen[1][c->blk] : -1;
    const int zin = (c->blk>0 && c->lev==0) ? bzone[0][c->blk] : -1;
    const int zout = (c->blk>0 && c->lev==0) ? bzone[1][c->blk] : -1;
    const int rin = (c->blk>0 && c->lev==0) ? bzr[0][c->blk] : -1;
    const int rout = (c->blk>0 && c->lev==0) ? bzr[1][c->blk] : -1;
    const bool open = oin>=0 || oout>=0 || zin>=0 || zout>=0;
    const bool zoned = open || rin>=0 || rout>=0;
    if(c->blk>0)
    {
        pp->B60 = (oin>=0 || oout>=0) ? 1 : 0;
        pp->river_xm = (oin>=0) ? 1 : 0;
        pp->river_xp = (oout>=0) ? 1 : 0;
        pp->B98 = 0;
        pp->B99 = 0;
        pp->open_xm = (zin>=0) ? 1 : 0;
        pp->open_xp = (zout>=0) ? 1 : 0;
        pp->open_ym = pp->open_yp = 0;
        // an edge zone: the non-hydrostatic pressure at the -i edge 0 as at the iowave edges (nhflow_poisson, B 90);
        // a relaxation or beach zone: the relaxation of the semi-implicit step (nhflow_momentum_SI::relax)
        if(zin>=0 || zout>=0 || rin>=0 || rout>=0)
        pp->B90 = 1;
        // a relaxation zone on the -i face: that face is the generation side of the block as x- of the native grid
        // with B 98 2 (the velocities zero gradient in the ghost cells, start4V; not a G 52 wall, nhflow_curv::iwall)
        if(rin>=0 && bzones[rin].method==bc_method::relax)
        pp->B98 = 2;
        if(c->lev==0)
        pp->amr_block = p0->gridgeo->blk[c->blk].orig;   // the block of DIVEMesh (G 55: of the piece)
        // the blocks of DIVEMesh (the turns of the blocks: the faces of the B 521 edges, nhflow_momentum_SI)
        pp->gridgeo = p0->gridgeo;
    }

    build_bc2D(pgc,*c);
    build_lexer3D(*c);

    if(open)
    open_lists(*c,oin>=0 || zin>=0,oout>=0 || zout>=0);

    // vectors of the Poisson rows (nhflow_poisson fills them over FBASELOOP and LOOP)
    pp->veclength = pp->imax*pp->jmax*pp->kmaxF;

    c->d = new fdm_nhf(pp);
    fdm_nhf *d = c->d;
    pscope ps(pgc,d,d0);

    // G 30: the cells around the patch keep the flags of their source cell through the
    // wetting and drying of the patch (set by every fill), so that a fine face on the patch box
    // is open only where the coarse face is
    if(shore)
    {
        c->wfix.assign(pp->imax*pp->jmax,-1);
        pp->wetfix = c->wfix.data();
    }

    // A 550: the filled cells take their breaking from the source grid (nhflow_breaking does not
    // detect in them)
    if(pp->A550==1)
    {
        c->vbfill.assign(pp->imax*pp->jmax,-1.0);
        pp->amrvb = c->vbfill.data();
    }

    c->pBC = pBCv;
    c->pflow = pflowv;
    if(zoned)
    {
        c->pio = new ioflow_amr_block(pp,pgc,pBCv);
        pp->amr_own_io = 1;
        c->pio->open_faces(pp,oin,oout);
        c->pio->open_edges(pp,zin>=0 ? &bzones[zin] : nullptr,zout>=0 ? &bzones[zout] : nullptr,bzbg);
        // relaxation and beach zones (B 522 methods 1, 2), the waves with B 90 1 and B 98 2 on level 0
        c->pio->open_zones(pp,rin>=0 ? &bzones[rin] : nullptr,rout>=0 ? &bzones[rout] : nullptr,bzbg,p0->B90>0 && p0->B98==2,
                           EXT,blk[c->blk].ni,EXT,blk[c->blk].nj);
        c->pflow = c->pio;
    }
    // the floating bodies on the patch grid (all calls do nothing without bodies)
    // (multi-block grid: the blocks have their own logical spacing, the band by the level)
    c->p6dof = new nhflow_amr_6dof(pp,b6,(b6==nullptr) ? 0.0 : (nblk>0) ? b6->object(0)->nhflow_dsm()/double(1<<c->lev)
                                                          : b6->object(0)->nhflow_dsm()*pp->DXM/glex(-1)->DXM);
    c->psed = psedv;
    // S 10: the sediment of the patch (nhflow_amr_sed.cpp), its state set by regrid_state
    if(pp->S10==1)
    c->sed = new sediment_f(pp,pgc,nullptr,c->pBC,1);

    c->pss = new nhflow_signal_speed(pp);

    if(pp->A511==1)
    c->pconv = new nhflow_HLL(pp,pgc,c->pBC);
    if(pp->A511==2)
    c->pconv = new nhflow_HLL(pp,pgc,c->pBC);   // HLLC mode inside nhflow_HLL (A 511 2)

    // diffusion as on level 0 (A 512 1 explicit, 2 implicit with the breaking viscosity vb), the
    // implicit one with a solver of its own: bicgstab_ijk on the patch, its global sums stay on
    // the patch (pscope)
    // (a turbulence model with A 512 0: implicit, as driver_logic_nhflow on level 0)
    const int a512 = (pp->A512==0 && pp->A560>0) ? 2 : pp->A512;
    if(a512==1)
    c->pdiff = new nhflow_ediff(pp);
    else if(a512==2)
    c->pdiff = new nhflow_idiff(pp);
    else
    c->pdiff = new nhflow_diff_void(pp);
    if(a512==2)
    c->psolv = new bicgstab_ijk(pp,nullptr,pgc);

    // the breaking count of a patch is not printed (nhflow_breaking, P 12)
    pp->P12 = 1<<30;

    if(pp->A514<=3)
    c->precon = new nhflow_reconstruct_hires(pp,c->pBC);
    if(pp->A514==4 || pp->A514==5)
    c->precon = new nhflow_reconstruct_weno(pp,c->pBC);

    if(pp->A520==0)
    c->ppress = new nhflow_pjm_hs(pp,d,c->pBC);
    if(pp->A520==1)
    c->ppress = new nhflow_pjm(pp,d,pgc,c->pBC);
    if(pp->A520==2)
    c->ppress = new nhflow_pjm_corr(pp,d,pgc,c->pBC);

    if(turb)
    turb_objects(*c,pgc);
    else
    c->pturb = new nhflow_komega_func_void(pp,d,pgc);
    c->pvrans = new vrans_nhflow_v(pp,d,pgc);
    c->pfsf = new nhflow_fsf_f(pp,d,pgc,c->pflow,c->pBC);
    c->pdf = new nhflow_forcing(pp,d,pgc);
    // fixed solids: cast within the patch arrays (the patch lexer has the extent of the level-0 rank),
    // the force objects of the boxes the patch reaches into
    if(solids)
    {
    c->pdf->amr_ray_box(pp->XN[marge],pp->XN[pp->knox+marge],pp->YN[marge],pp->YN[pp->knoy+marge]);
    force_objects(c);
    }

    if(pp->A510==2 && pp->A534!=1)
    c->pmom = new nhflow_momentum_RK2(pp,d,pgc,c->p6dof,c->pvrans,c->pdf,c->psed);
    if(pp->A510==3 && pp->A534!=1)
    c->pmom = new nhflow_momentum_RK3(pp,d,pgc,c->p6dof,c->pvrans,c->pdf);
    if(pp->A534==1)
    c->pmom = new nhflow_momentum_SI(pp,d,pgc,c->p6dof,c->pvrans,c->pdf,true);
    
    // pierced-column water level (X 13): on the patch too, local to it
    if(c->pmom!=nullptr && c->pmom->pfsfbody!=nullptr)
    c->pmom->pfsfbody->set_patch();

    // the hull as a sharp thin body (X 15 8): the patch holds it too, with its thin-body service (stage() from the
    // forcing of the patch, the flux hook through this class, the projection terms through the patch's rows)
    if(hb)
    {
    d->hullbody = new nhflow_hullbody(pp,d,pgc,b6->object(0));
    d->hullbody->amr_patch = true;
    d->thinbody->amr_patch = true;
    }

    c->S = {c->pflow,c->pss,c->precon,c->pconv,c->pdiff,c->ppress,nullptr,c->psolv,nullptr,c->pfsf,c->pturb,c->pvrans};

    // breaking with G 31 1: the stage runner runs the parts of the breaking (breaking_patches); on a
    // multi-block grid the blocks 1.. with block 0 (breaking_level0)
    if(brk_split || (brk0 && c->blk>0 && c->lev==0))
    if(nhflow_breaking *b = dynamic_cast<nhflow_breaking*>(c->pmom))
    b->brk_defer = true;

    // material and porosity as driver_ini_nhflow and nhflow_f::ini on level 0
    const int n7 = pp->imax*pp->jmax*(pp->kmax+2);
    for(int n=0; n<n7; ++n)
    {
        d->POR[n] = 1.0;
        d->RO[n] = pp->W1;
        d->VISC[n] = pp->W2;
    }
    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        d->ks(ii,jj) = pp->B50;
        d->breaking(ii,jj) = 0;
    }

    // multi-block: porosity, density and viscosity of a block as level 0 sets them (nhflow_f::ini,
    // driver_ini_nhflow): in the fluid cells, the ghost cells by start5Vfull (porosity; the corner below
    // the bed and beyond a wall stays 0) and start4V (density, viscosity; below the bed and above the
    // surface they stay 0); the Poisson and diffusion rows next to the walls, the bed and the surface take them
    if(c->blk>0)
    {
        for(int n=0; n<n7; ++n)
        d->POR[n] = d->RO[n] = d->VISC[n] = 0.0;
        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        for(int kk=pp->kmin; kk<pp->kmin+pp->kmax; ++kk)
        if(pp->flag4[cidx(pp,ii,jj,kk)]>0)
        {
            d->POR[cidx(pp,ii,jj,kk)] = 1.0;
            d->RO[cidx(pp,ii,jj,kk)] = pp->W1;
            d->VISC[cidx(pp,ii,jj,kk)] = pp->W2;
        }
        pgc->start5Vfull(pp,d->POR,1);
        pgc->start4V(pp,d->RO,1);
        pgc->start4V(pp,d->VISC,1);
    }
}

void nhflow_amr::patch_delete(reefamr_patch *q)
{
    nhflow_amr_patch *c = NP(q);

    for(auto v : c->kv)
    delete [] v;
    delete c->mg;

    delete c->sed;
    delete c->pmom;
    delete c->pio;
    if(c->pio!=nullptr)
    {
        c->pp->del_Iarray(c->pp->gcin,c->pp->gcin_count,6);
        c->pp->del_Iarray(c->pp->gcout,c->pp->gcout_count,6);
    }
    delete c->p6dof;
    for(auto f : c->forces)
    delete f;
    c->forces.clear();
    delete c->pdf;
    delete c->pfsf;
    delete c->pvrans;
    turb_free(*c);
    delete c->pturb;
    delete c->ppress;
    delete c->precon;
    delete c->pdiff;
    // the solver interface has no virtual destructor: delete as the class it is
    delete static_cast<bicgstab_ijk*>(c->psolv);
    delete c->pconv;
    delete c->pss;

    // fdm_nhf has no destructor (on level 0 it lives for the whole run): its 3D arrays are
    // freed here, otherwise every regrid that replaces a patch leaks them
    {
        fdm_nhf *d = c->d;
        lexer *pp = c->pp;
        const int n7 = pp->imax*pp->jmax*(pp->kmax+2);

        // X 15 8: the hull of the patch, its thin-body service and the mobility fields of the projection
        delete d->hullbody;
        d->hullbody = nullptr;
        delete d->thinbody;
        d->thinbody = nullptr;
        double **mb[] = {&d->MBETA,&d->MCHI,&d->MRCX,&d->MRCY,&d->MBX,&d->MBY,&d->MBZ};
        for(auto f : mb)
        if(*f!=nullptr)
        pp->del_Darray(*f,n7);
        double **a[] = {&d->U,&d->V,&d->W,&d->omegaF,&d->UH,&d->VH,&d->WH,&d->P,&d->RO,&d->VISC,&d->EV,&d->EV0,
                        &d->F,&d->G,&d->H,&d->L,&d->Fext,&d->Gext,&d->Hext,&d->POR,&d->PORPART,&d->PORDEM,&d->test,
                        &d->KIN,&d->CONC,&d->SOLID,&d->FB,&d->FHB,&d->PORSTRUC,&d->Fx,&d->Fy,&d->Fz,&d->FEx,&d->FEy,
                        &d->FSW,&d->DWDT,&d->Fs,&d->Fn,&d->Fe,&d->Fw,&d->Ss,&d->Sn,&d->Se,&d->Sw,&d->SSx,&d->SSy,
                        &d->Us,&d->Un,&d->Ue,&d->Uw,&d->Ub,&d->Ut,&d->Vs,&d->Vn,&d->Ve,&d->Vw,&d->Vb,&d->Vt,
                        &d->Ws,&d->Wn,&d->We,&d->Ww,&d->Wb,&d->Wt,&d->UHs,&d->UHn,&d->UHe,&d->UHw,&d->UHb,&d->UHt,
                        &d->VHs,&d->VHn,&d->VHe,&d->VHw,&d->VHb,&d->VHt,&d->WHs,&d->WHn,&d->WHe,&d->WHw,&d->WHb,&d->WHt};
        for(auto f : a)
        pp->del_Darray(*f,n7);
        pp->del_Iarray(d->NODEVAL,pp->imax*pp->jmax*(pp->kmax+3));
    }
    delete c->d;

    free_lexer3D(*c);
    nhflow_curv::pair_forget(c->pp);
}

// cells of level l that need level l+1 (the boxes G 10 and the body zone are marked by the
// core): the end-of-step surface of the level-l grids, wet cells deeper than A 545 A 544 only
void nhflow_amr::tag(int l, vector<unsigned char> &M)
{
    if(!adaptive)
    return;

    const int lf = l+1;

    auto test = [&](lexer *q, fdm_nhf *d, int ii, int jj)
    {
        const double wmin = q->A545*q->A544;
        auto ok = [&](int a, int b)
        {
            const int n = lij(q,a,b);
            return q->flagslice4[n]>0 && q->wet[n]==1 && d->WL(a,b)>wmin;
        };

        // breaking (G 23): a cell with breaking viscosity
        if(flagbreak && q->flagslice4[lij(q,ii,jj)]>0 && d->vb(ii,jj)>0.0)
        return true;

        // shoreline (G 22): a fluid cell with a cell of the other wet state within nshore cells
        if(nshore>0 && q->flagslice4[lij(q,ii,jj)]>0)
        {
            const int w0 = q->wet[lij(q,ii,jj)];
            for(int a=-nshore; a<=nshore; ++a)
            for(int b=-nshore; b<=nshore; ++b)
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
            if(ok(ii+di[k],jj+dj[k]) && fabs(d->eta(ii+di[k],jj+dj[k])-ec)>tol_eta)
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

    // level 0: the native grid without the ring, then (multi-block) the blocks 1.. in lev[0]
    if(l==0)
    {
        for(int ii=0; ii<NX0; ++ii)
        for(int jj=0; jj<NY0; ++jj)
        if(!ring0(ii,jj) && test(p0,d0,ii,jj))
        tag_cell(lf,ii+O0i,jj+O0j,nbuf,M);
    }

    for(int n : lev[l])
    {
        nhflow_amr_patch *c = NP(n);
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EY; jj<EY+c->ny; ++jj)
        if(test(c->pp,c->d,ii,jj))
        tag_cell(lf,ii-EXT+c->I0,jj-EY+c->J0,nbuf,M);
    }
}

// a level-0 cell no patch may cover: dry or shallow (water level at most A 545 A 544, the deep
// criterion of the wetting and drying); the patches are fully wet
bool nhflow_amr::cell_unfit(int ii, int jj)
{
    const int n = lij(p0,ii,jj);
    return p0->flagslice4[n]>0 && (p0->wet[n]==0 || d0->WL(ii,jj)<=p0->A545*p0->A544);
}

// multi-block: a cell of a block 1..
bool nhflow_amr::cell_unfit_patch(reefamr_patch *c, int ii, int jj)
{
    lexer *q = c->pp;
    const int n = lij(q,ii,jj);
    return q->flagslice4[n]>0 && (q->wet[n]==0 || NP(c)->d->WL(ii,jj)<=q->A545*q->A544);
}

// --------------------------------------------------------------------- wetting and drying
// coarse cell (a,b) of lexer q usable in a stencil: mode 0 fluid, 1 wet, 2 wet and deep
bool nhflow_amr::wet_at(lexer *q, int a, int b, int mode) const
{
    const int n = lij(q,a,b);
    if(q->flagslice4[n]<=0)
    return false;
    if(mode>=1 && q->wet[n]!=1)
    return false;
    if(mode>=2 && q->deep[n]!=1)
    return false;
    return true;
}

// a dry cell as nhflow_fsf_f::wetdry leaves it: water level A 544, the surface on the bed, no
// velocity and no momentum
void nhflow_amr::dry_cell(lexer *pp, fdm_nhf *d, slice &WL, double *UH, double *VH, double *WH, int ii, int jj)
{
    d->eta(ii,jj) = pp->A544 - d->depth(ii,jj);
    WL(ii,jj) = pp->A544;
    for(int kk=0; kk<pp->knoz; ++kk)
    {
        const int n = cidx(pp,ii,jj,kk);
        d->U[n] = d->V[n] = d->W[n] = 0.0;
        UH[n] = VH[n] = WH[n] = 0.0;
    }
}

// flags of the interior of a fresh patch from its water level (the threshold of the wetting and
// drying); the patch's own wetdry then sets deep and applies its rules
void nhflow_amr::patch_flags(nhflow_amr_patch &c)
{
    lexer *pp = c.pp;
    fdm_nhf *d = c.d;
    for(int ii=EXT; ii<EXT+c.nx; ++ii)
    for(int jj=EY; jj<EY+c.ny; ++jj)
    {
        const int n = lij(pp,ii,jj);
        if(pp->flagslice4[n]<=0)
        continue;
        pp->wet[n] = pp->deep[n] = (d->WL(ii,jj)>pp->A544+1.0e-6) ? 1 : 0;
    }
}

// deep as at the end of nhflow_fsf_f::wetdry (no dry cell within 2 cells in x, 3 in y, the
// diagonal neighbours, and a water level above A 545 A 544), over the computed cells of q
void nhflow_amr::deep_rule(lexer *q, slice &WL)
{
    auto w = [&](int a, int b) { return q->wet[lij(q,a,b)]; };

    for(int ii=0; ii<q->knox; ++ii)
    for(int jj=0; jj<q->knoy; ++jj)
    {
        const int n = lij(q,ii,jj);
        if(q->flagslice4[n]<=0)
        continue;

        int dp = q->wet[n];
        if(w(ii+1,jj)==0 || w(ii+2,jj)==0 || w(ii-1,jj)==0 || w(ii-2,jj)==0)
        dp = 0;
        if(q->j_dir==1)
        {
            if(w(ii,jj+1)==0 || w(ii,jj+2)==0 || w(ii,jj+3)==0 || w(ii,jj-1)==0 || w(ii,jj-2)==0 || w(ii,jj-3)==0)
            dp = 0;
            if(w(ii+1,jj+1)==0 || w(ii+1,jj-1)==0 || w(ii-1,jj+1)==0 || w(ii-1,jj-1)==0)
            dp = 0;
        }
        if(WL(ii,jj)<=q->A545*q->A544)
        dp = 0;
        q->deep[n] = dp;
    }
}

// wet and deep of the covered coarse cells after the restriction of stage s (s<0: end of the
// step): wet only if all four children are wet - a coarse cell over the fine shoreline would
// take a surface WL - depth of the coarse cell above the water (the dry children hold A 544),
// rewet its coarse neighbours and drive them with that false slope; the fluxes across the patch
// box come from the fine faces anyway.  Then deep on every grid that holds covered cells, as its
// own wetdry would set it, so that the coarse rewetting rule, the flux zeroing and the Poisson
// rows next to a patch see the front of the patch
void nhflow_amr::restrict_flags(ghostcell *pgc, int s)
{
    if(!shore)
    return;

    for(int l=wtop(); l>=wlo+1; --l)
    block_up(l,1,7300+l,
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
                 q->wet[lij(q,B.ic,B.jc)] = (int)v[0];
             });

    // patches that hold covered cells (no partition exchange), then level 0 with its halo
    {
    comms_off guard(pgc);
    for(int l=MAX(wlo,1); l<wtop(); ++l)
    for(int id : lev[l])
    {
        nhflow_amr_patch *c = NP(id);
        lexer *pp = c->pp;
        slice &WL = (s<0) ? c->d->WL : *stage_out(id,s).WL;
        pgc->gcsl_start4Vint(pp,pp->wet,50);
        deep_rule(pp,WL);
        pgc->gcsl_start4Vint(pp,pp->deep,50);
    }
    }

    if(wlo>0)
    return;

    slice &WL0 = (s<0) ? d0->WL : *stage_out(-1,s).WL;
    pgc->gcsl_start4Vint(p0,p0->wet,50);
    deep_rule(p0,WL0);
    pgc->gcsl_start4Vint(p0,p0->deep,50);
}

// --------------------------------------------------------------------- interpolation
// interpolation from coarse cell (ic,jc) to the centre of its child (ox,oy), a quarter coarse
// cell away in x and y (as fnpf_amr): bicubic where all 16 cells are fluid, otherwise
// biquadratic on the 3x3 cells, slopes and curvatures switched off next to solids
void nhflow_amr::pweights(lexer *q, int ic, int jc, int ox, int oy, double *w, int mode)
{
    const int md = shore ? mode : 0;
    auto fl = [&](int a, int b) { return wet_at(q,a,b,md); };

    for(int m=0; m<25; ++m)
    w[m] = 0.0;

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
        const double cc = 0.25*x*y;
        W(1,1) += cc; W(-1,1) -= cc; W(1,-1) -= cc; W(-1,-1) += cc;
    }
}

double nhflow_amr::pq(slice &f, lexer *q, int ic, int jc, int ox, int oy, int mode)
{
    double w[25];
    pweights(q,ic,jc,ox,oy,w,mode);
    return pqw(f,q,ic,jc,w);
}

// with the weights w of pweights
double nhflow_amr::pqw(slice &f, lexer *q, int ic, int jc, const double *w)
{
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

// layer kk of a cell-layout array (cl true) or node kk of an F-layout array
double nhflow_amr::pq3(const double *f, lexer *q, int ic, int jc, int ox, int oy, int kk, bool cl, int mode)
{
    double w[25];
    pweights(q,ic,jc,ox,oy,w,mode);
    return pq3w(f,q,ic,jc,kk,cl,w);
}

double nhflow_amr::pq3w(const double *f, lexer *q, int ic, int jc, int kk, bool cl, const double *w)
{
    double r = 0.0;
    for(int di=-2; di<=2; ++di)
    for(int dj=-2; dj<=2; ++dj)
    {
        const double c = w[(di+2)*5+(dj+2)];
        if(c!=0.0)
        r += c*f[cl ? cidx(q,ic+di,jc+dj,kk) : fidx(q,ic+di,jc+dj,kk)];
    }
    return r;
}

// linear with limited slopes: the mean of the four children is the coarse value (the bed, so that
// the restricted water level of still water is still water)
double nhflow_amr::plin(slice &f, lexer *q, int ic, int jc, int ox, int oy)
{
    auto fl = [&](int a, int b) { return q->flagslice4[lij(q,a,b)]>0; };
    const double c = f(ic,jc);
    double sx=0.0, sy=0.0;
    if(fl(ic+1,jc) && fl(ic-1,jc))
    sx = mmod(f(ic+1,jc)-c,c-f(ic-1,jc));
    if(fl(ic,jc+1) && fl(ic,jc-1))
    sy = mmod(f(ic,jc+1)-c,c-f(ic,jc-1));
    return c + 0.25*ox*sx + 0.25*oy*sy;
}

// column of a child of coarse cell (ic,jc) of grid g: nodes 0..knf, the sigma nodes of grid g
// (knf = knoz of g) or, with vertical refinement (knf = 2 knoz of g), the nested nodes: the coarse
// nodes are the even nodes, the odd nodes linear between them
void nhflow_amr::pcol(int g, int ic, int jc, int ox, int oy, const double *src, int knf, double *v, int mode)
{
    lexer *q = glex(g);
    const int kc = q->knoz;
    const int fz = knf/kc;
    const int sI = q->jmax*q->kmaxF;
    const int sJ = q->kmaxF;

    double w[25];
    pweights(q,ic,jc,ox,oy,w,mode);

    const double *s[25];
    double ww[25];
    int nw=0;
    const int n0 = fidx(q,ic,jc,0);
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

    for(int kk=0; kk<=kc; ++kk)
    {
        double r = 0.0;
        for(int m=0; m<nw; ++m)
        r += ww[m]*s[m][kk];
        v[fz*kk] = r;
    }

    if(fz==2)
    for(int kk=0; kk<kc; ++kk)
    v[2*kk+1] = 0.5*(v[2*kk]+v[2*kk+2]);
}

// cell values of the nc layers of a column of grid lexer qc -> the fz*nc layers of the finer
// grid (fz 1 or 2): with vertical refinement linear in sigma with the minmod-limited slope of
// the neighbouring layers (one-sided in the bed and the surface layer); the two halves of a
// layer have the same thickness, so their mean is the coarse value (conservative)
void nhflow_amr::vcell(lexer *qc, const double *c, int nc, int fz, double *f)
{
    if(fz==1)
    {
        for(int k=0; k<nc; ++k)
        f[k] = c[k];
        return;
    }

    const double *ZP = qc->ZP + marge;
    const double *DZ = qc->DZN + marge;
    for(int k=0; k<nc; ++k)
    {
        double s = 0.0;
        if(nc>1)
        {
            if(k==0)
            s = (c[1]-c[0])/(ZP[1]-ZP[0]);
            else if(k==nc-1)
            s = (c[k]-c[k-1])/(ZP[k]-ZP[k-1]);
            else
            s = mmod((c[k+1]-c[k])/(ZP[k+1]-ZP[k]),(c[k]-c[k-1])/(ZP[k]-ZP[k-1]));
        }
        const double h = 0.25*DZ[k]*s;
        f[2*k] = c[k] - h;
        f[2*k+1] = c[k] + h;
    }
}

// --------------------------------------------------------------------- fill
// the cells around the level-l patches with the input of stage s (s<0: end of the step):
// eta, wet, deep, U, V, W per layer and P per node; water level and momentum follow from the
// depth of the patch (well balanced)
// the end of the step: the cells around the level-l patches filled from a patch of the same level (a
// sibling on this rank, or across a partition edge) take its eta, U, V, W and P of the end of the step,
// so that an output on the finest grid next to such a cell reads current values (reefamr::finest_point
// with fills_current).  The cells prolonged from a coarser grid keep theirs.  The next stage fills all
// of them again (fill_stage), with the same values for these cells
void nhflow_amr::fill_out(ghostcell *pgc)
{
    if(nblk>0)
    return;

    for(int l=1; l<=maxlev; ++l)
    {
        const int K = klev(l);
        const int nv = 1 + 3*K + K+1;

        fill_run(l,nv,7980+l,
                 [&](const reefamr_fill &f, double *v)
                 {
                     for(int m=0; m<nv; ++m)
                     v[m] = 0.0;
                     if(f.kind!=0)
                     return;
                     lexer *q = glex(f.g);
                     fdm_nhf *d = gfd(f.g);
                     v[0] = d->eta(f.si,f.sj);
                     for(int k=0; k<K; ++k)
                     {
                         const int n = cidx(q,f.si,f.sj,k);
                         v[1+k] = d->U[n];
                         v[1+K+k] = d->V[n];
                         v[1+2*K+k] = d->W[n];
                     }
                     for(int k=0; k<=K; ++k)
                     v[1+3*K+k] = d->P[fidx(q,f.si,f.sj,k)];
                 },
                 [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
                 {
                     // a cell of a sibling: on this rank (kind 0), or a level-l patch holds it (kind 2)
                     if(!(f.kind==0 || (f.kind==2 && covered(l,f.di-EXT+q->I0,f.dj-EY+q->J0))))
                     return;
                     lexer *pp = q->pp;
                     fdm_nhf *d = NP(q)->d;
                     d->eta(f.di,f.dj) = w[0];
                     for(int k=0; k<K; ++k)
                     {
                         const int n = cidx(pp,f.di,f.dj,k);
                         d->U[n] = w[1+k];
                         d->V[n] = w[1+K+k];
                         d->W[n] = w[1+2*K+k];
                     }
                     for(int k=0; k<=K; ++k)
                     d->P[fidx(pp,f.di,f.dj,k)] = w[1+3*K+k];
                 });
    }
}

void nhflow_amr::fill_stage(ghostcell *pgc, int l, int s)
{
    const int K = klev(l);
    const int Kc = klev(l-1);
    const bool brk = (p0->A550==1);
    const int nb = 3 + 3*K + K+1;                // position of vb (A 550)
    // sediment (S 10): the water level of the cell around the patch from the water level of its source, not
    // from its eta. After a bed update eta is that of the old bed until the first stage (the uniform grid
    // keeps WL and eta of the cells as they are): a copy takes the stage input WL of the source cell, a
    // coarser source WL - depth of the new bed (the stage after the bed update, s <= 0)
    const bool sedw = (p0->S10>0);
    const int nwl = nb + (brk ? 1 : 0);          // position of the water level and its kind
    const int nv = nwl + (sedw ? 2 : 0);
    vector<double> uc(3*Kc);

    fill_run(l,nv,7100+l,
             [&](const reefamr_fill &f, double *v)
             {
                 if(f.kind==0 || f.kind==1)
                 {
                     lexer *q = glex(f.g);
                     fdm_nhf *d = gfd(f.g);
                     if(f.kind==0)
                     {
                         v[0] = d->eta(f.si,f.sj);
                         v[1] = q->wet[lij(q,f.si,f.sj)];
                         v[2] = q->deep[lij(q,f.si,f.sj)];
                         for(int k=0; k<K; ++k)
                         {
                             const int n = cidx(q,f.si,f.sj,k);
                             v[3+k] = d->U[n];
                             v[3+K+k] = d->V[n];
                             v[3+2*K+k] = d->W[n];
                         }
                         for(int k=0; k<=K; ++k)
                         v[3+3*K+k] = d->P[fidx(q,f.si,f.sj,k)];
                     }
                     else
                     {
                         // G 30: surface and velocities from wet cells, pressure from wet and
                         // deep cells; the flags are the parent's (a fine face on the patch box
                         // carries mass only where the coarse face can)
                         // the weights of the wet cells once for the column
                         double w1[25];
                         pweights(q,f.si,f.sj,f.ox,f.oy,w1,1);
                         v[0] = pqw(d->eta,q,f.si,f.sj,w1);
                         v[1] = q->wet[lij(q,f.si,f.sj)];
                         v[2] = q->deep[lij(q,f.si,f.sj)];
                         // layer by layer on the coarser grid, then (G 6) into the halves
                         for(int k=0; k<Kc; ++k)
                         {
                             uc[k] = pq3w(d->U,q,f.si,f.sj,k,true,w1);
                             uc[Kc+k] = pq3w(d->V,q,f.si,f.sj,k,true,w1);
                             uc[2*Kc+k] = pq3w(d->W,q,f.si,f.sj,k,true,w1);
                         }
                         for(int m=0; m<3; ++m)
                         vcell(q,&uc[m*Kc],Kc,K/Kc,&v[3+m*K]);
                         pcol(f.g,f.si,f.sj,f.ox,f.oy,d->P,K,&v[3+3*K],2);
                     }
                     // A 550: the breaking viscosity of the source cell
                     if(brk)
                     v[nb] = d->vb(f.si,f.sj);
                     if(sedw)
                     {
                         if(f.kind==0)
                         {
                         v[nwl] = (*stage_in(f.g,s).WL)(f.si,f.sj);
                         v[nwl+1] = 0.0;
                         }
                         else
                         {
                         double w1[25];
                         pweights(q,f.si,f.sj,f.ox,f.oy,w1,1);
                         v[nwl] = (s<=0) ? pqw(d->WL,q,f.si,f.sj,w1) - pqw(d->depth,q,f.si,f.sj,w1) : v[0];
                         v[nwl+1] = 1.0;
                         }
                     }
                     return;
                 }
                 for(int m=0; m<nv; ++m)
                 v[m] = 0.0;
             },
             [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
             {
                 nhflow_amr_patch *c = NP(q);
                 lexer *pp = c->pp;
                 fdm_nhf *d = c->d;
                 stg S = stage_in(id,s);
                 const int ii=f.di, jj=f.dj;

                 d->eta(ii,jj) = w[0];
                 double wl = MAX(w[0] + d->depth(ii,jj), pp->A544);
                 if(sedw)
                 wl = MAX((w[nwl+1]>0.5) ? w[nwl] + d->depth(ii,jj) : w[nwl], pp->A544);
                 (*S.WL)(ii,jj) = wl;
                 if(s<=0)
                 d->WL(ii,jj) = wl;
                 pp->wet[lij(pp,ii,jj)] = (int)w[1];
                 pp->deep[lij(pp,ii,jj)] = (int)w[2];

                 for(int k=0; k<K; ++k)
                 {
                     const int n = cidx(pp,ii,jj,k);
                     d->U[n] = w[3+k];
                     d->V[n] = w[3+K+k];
                     d->W[n] = w[3+2*K+k];
                     S.UH[n] = wl*w[3+k];
                     S.VH[n] = wl*w[3+K+k];
                     S.WH[n] = wl*w[3+2*K+k];
                 }
                 for(int k=0; k<=K; ++k)
                 d->P[fidx(pp,ii,jj,k)] = w[3+3*K+k];

                 if(brk)
                 c->vbfill[lij(pp,ii,jj)] = w[nb];

                 // G 30: the flag of the source cell is kept through the patch's wetdry; a dry
                 // source cell gives a dry cell on the bed of the patch
                 if(shore)
                 c->wfix[lij(pp,ii,jj)] = (int)w[1];
                 if(shore && w[1]==0.0)
                 {
                     dry_cell(pp,d,*S.WL,S.UH,S.VH,S.WH,ii,jj);
                     if(s<=0)
                     d->WL(ii,jj) = pp->A544;
                 }
             });

    comms_off guard(pgc);

    for(int n : lev[l])
    bc_patch(pgc,*NP(n),s);
}

// wall ghost cells of a patch from its (filled) fluid cells, as at the end of a level-0 stage
void nhflow_amr::bc_patch(ghostcell *pgc, nhflow_amr_patch &c, int s)
{
    lexer *pp = c.pp;
    fdm_nhf *d = c.d;
    pscope ps(pgc,d,d0);
    int id=-1;
    for(int n=0; n<(int)P.size(); ++n)
    if(P[n]==&c)
    id=n;
    stg S = stage_in(id,s);

    pgc->gcsl_start4(pp,d->eta,gcval_eta);
    pgc->gcsl_start4(pp,*S.WL,gcval_eta);
    pgc->gcsl_start4Vint(pp,pp->wet,50);
    pgc->gcsl_start4Vint(pp,pp->deep,50);
    pgc->start4V(pp,d->U,10);
    pgc->start4V(pp,d->V,11);
    pgc->start4V(pp,d->W,12);
    pgc->start4V(pp,S.UH,14);
    pgc->start4V(pp,S.VH,15);
    pgc->start4V(pp,S.WH,16);
    pgc->start7P(pp,d->P,540);
}

// bed of the level-l patches around their interior: a patch of the same level or the coarser
// level (linear, mean preserving)
void nhflow_amr::fill_bed(ghostcell *pgc, int l)
{
    fill_run(l,1,7050+l,
             [&](const reefamr_fill &f, double *v)
             {
                 if(f.kind==0)
                 v[0] = gfd(f.g)->bed(f.si,f.sj);
                 else if(f.kind==1)
                 v[0] = plin(gfd(f.g)->bed,glex(f.g),f.si,f.sj,f.ox,f.oy);
                 else
                 v[0] = 0.0;
             },
             [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
             {
                 NP(q)->d->bed(f.di,f.dj) = w[0];
             });
}

// interior state of the fresh level-l patches from their parents: eta, U, V, W per layer and P
// (bicubic), water level and momentum with the depth of the patch.  The parent of every 2x2 block
// (on its rank) gives the interpolated columns of the four children and its own water level, wet
// flag and momentum; the patch does the rest.
void nhflow_amr::prolong_patch(ghostcell *pgc, int l)
{
    const int K = klev(l);
    const int Kc = klev(l-1);
    const int fz = K/Kc;
    const int cs = 3+4*K;               // per child: eta, detadt, U, V, W (K each), P (K+1)
    const int t0 = 4*cs;                // parent: WL, wet, UH VH WH per layer
    const int nv = t0+2+3*Kc;

    block_down_if(l,nv,7600+l,[](reefamr_patch *q) { return q->fresh; },
               [&](const reefamr_block &B, int key, double *v)
               {
                   lexer *q = glex(B.g);
                   fdm_nhf *dc = gfd(B.g);
                   const int ic = B.ic, jc = B.jc;
                   vector<double> uc(3*Kc);

                   for(int a=0; a<2; ++a)
                   for(int b=0; b<=sy; ++b)
                   {
                       const int ox = a==0?-1:1, oy = (sy==0) ? 0 : (b==0?-1:1);
                       double *w = &v[(2*a+b)*cs];

                       double w1[25];
                       pweights(q,ic,jc,ox,oy,w1,1);
                       w[0] = pqw(dc->eta,q,ic,jc,w1);
                       w[1] = pqw(dc->detadt,q,ic,jc,w1);

                       for(int kc=0; kc<Kc; ++kc)
                       {
                           uc[kc] = pq3w(dc->U,q,ic,jc,kc,true,w1);
                           uc[Kc+kc] = pq3w(dc->V,q,ic,jc,kc,true,w1);
                           uc[2*Kc+kc] = pq3w(dc->W,q,ic,jc,kc,true,w1);
                       }
                       for(int m=0; m<3; ++m)
                       vcell(q,&uc[m*Kc],Kc,fz,&w[2+m*K]);

                       pcol(B.g,ic,jc,ox,oy,dc->P,K,&w[2+3*K],2);
                   }

                   v[t0] = dc->WL(ic,jc);
                   v[t0+1] = q->wet[lij(q,ic,jc)];
                   for(int kc=0; kc<Kc; ++kc)
                   {
                       const int nc = cidx(q,ic,jc,kc);
                       v[t0+2+3*kc] = dc->UH[nc];
                       v[t0+3+3*kc] = dc->VH[nc];
                       v[t0+4+3*kc] = dc->WH[nc];
                   }
               },
               [&](reefamr_patch *qq, int id, int k, const double *v)
               {
                   if(!qq->fresh)
                   return;

                   nhflow_amr_patch &c = *NP(qq);
                   lexer *pp = c.pp;
                   fdm_nhf *d = c.d;
                   const int nby = c.ny>>sy;
                   const int bi = k/nby, bj = k%nby;

                   for(int a=0; a<2; ++a)
                   for(int b=0; b<=sy; ++b)
                   {
                       const int ii = EXT+2*bi+a, jj = EY+ychild(bj)+b;
                       const double *w = &v[(2*a+b)*cs];

                       d->eta(ii,jj) = w[0];
                       d->detadt(ii,jj) = w[1];
                       const double wl = MAX(d->eta(ii,jj) + d->depth(ii,jj), pp->A544);
                       d->WL(ii,jj) = wl;

                       for(int kk=0; kk<K; ++kk)
                       {
                           const int n = cidx(pp,ii,jj,kk);
                           d->U[n] = w[2+kk];
                           d->V[n] = w[2+K+kk];
                           d->W[n] = w[2+2*K+kk];
                           d->UH[n] = wl*d->U[n];
                           d->VH[n] = wl*d->V[n];
                           d->WH[n] = wl*d->W[n];
                       }

                       for(int kk=0; kk<=K; ++kk)
                       d->P[fidx(pp,ii,jj,kk)] = w[2+3*K+kk];
                   }

                   // conservative: the 2x2 block keeps the water level and, layer by layer, the momentum
                   // of its coarse cell (the restriction is the block mean, with G 6 over the 2x2x2
                   // cells of a coarse layer), a constant shift of the interpolated shape; the velocities
                   // follow
                   const int i0 = EXT+2*bi, j0 = EY+ychild(bj);
                   double cw[4];                       // weights of the children (a,b): cw[a+2b]
                   cw4(pp,i0,j0,cw);
                   double wm = 0.0;
                   for(int a=0; a<2; ++a)
                   for(int b=0; b<=sy; ++b)
                   wm += cw[a+2*b]*d->WL(i0+a,j0+b);
                   const double dwl = v[t0] - wm;
                   const bool pwet = ((int)v[t0+1]==1);

                   // G 30, well balanced at the shoreline: a dry parent gives dry children on the bed of
                   // the patch; a block where a child falls dry (the interpolated surface, or the shift,
                   // below the bed + A 544) keeps the surface of the wet parent without the shift (eta
                   // flat at rest; the water volume of the block changes by the clipping) and its dry
                   // children lose their momentum
                   if(shore)
                   {
                       const double wdry = pp->A544 + 1.0e-6;
                       bool drych = !pwet;
                       for(int a=0; a<2 && !drych; ++a)
                       for(int b=0; b<=sy && !drych; ++b)
                       if(d->eta(i0+a,j0+b) + d->depth(i0+a,j0+b)<=wdry || d->WL(i0+a,j0+b) + dwl<=wdry)
                       drych = true;

                       if(drych)
                       {
                           const bool pdry = !pwet;

                           // at t = 0 the shoreline of the patch is set on its own bed with the initial
                           // surface of NHFLOW (the still water level and the F 72 boxes, as nhflow_f::ini
                           // and nhflow_fsf_f::ini on level 0): a fine cell below the water is wet even
                           // where its coarse parent is dry
                           if(p0->count==0)
                           for(int a=0; a<2; ++a)
                           for(int b=0; b<=sy; ++b)
                           {
                               const int ii = i0+a, jj = j0+b;
                               double e = 0.0;
                               for(int qn=0; qn<p0->F72; ++qn)
                               if(cxc(pp,ii,jj)>=p0->F72_xs[qn] && cxc(pp,ii,jj)<p0->F72_xe[qn]
                               && cyc(pp,ii,jj)>=p0->F72_ys[qn] && cyc(pp,ii,jj)<p0->F72_ye[qn])
                               e = p0->F72_h[qn] - p0->F60;
                               if(pdry || d->eta(ii,jj) + d->depth(ii,jj)<=wdry)
                               {
                                   d->eta(ii,jj) = e;
                                   d->WL(ii,jj) = MAX(e + d->depth(ii,jj),pp->A544);
                                   for(int kk=0; kk<K; ++kk)
                                   {
                                       const int n = cidx(pp,ii,jj,kk);
                                       d->U[n] = d->V[n] = d->W[n] = 0.0;
                                       d->UH[n] = d->VH[n] = d->WH[n] = 0.0;
                                   }
                               }
                           }

                           for(int a=0; a<2; ++a)
                           for(int b=0; b<=sy; ++b)
                           if((pdry && p0->count>0) || d->eta(i0+a,j0+b) + d->depth(i0+a,j0+b)<=wdry)
                           {
                               dry_cell(pp,d,d->WL,d->UH,d->VH,d->WH,i0+a,j0+b);
                               d->detadt(i0+a,j0+b) = 0.0;
                           }
                           return;
                       }
                   }
                   for(int a=0; a<2; ++a)
                   for(int b=0; b<=sy; ++b)
                   {
                       d->WL(i0+a,j0+b) += dwl;
                       d->eta(i0+a,j0+b) += dwl;
                   }

                   for(int kc=0; kc<Kc; ++kc)
                   {
                       double um=0.0, vm=0.0, hm=0.0;
                       for(int kk=fz*kc; kk<fz*kc+fz; ++kk)
                       for(int a=0; a<2; ++a)
                       for(int b=0; b<=sy; ++b)
                       {
                           const double wb = cw[a+2*b]/double(fz);
                           const int n = cidx(pp,i0+a,j0+b,kk);
                           um += wb*d->UH[n];
                           vm += wb*d->VH[n];
                           hm += wb*d->WH[n];
                       }
                       for(int kk=fz*kc; kk<fz*kc+fz; ++kk)
                       for(int a=0; a<2; ++a)
                       for(int b=0; b<=sy; ++b)
                       {
                           const int n = cidx(pp,i0+a,j0+b,kk);
                           const double wl = d->WL(i0+a,j0+b);
                           d->UH[n] += v[t0+2+3*kc] - um;
                           d->VH[n] += v[t0+3+3*kc] - vm;
                           d->WH[n] += v[t0+4+3*kc] - hm;
                           const double wlvl = wl>pp->A544 ? wl : 1.0e20;
                           d->U[n] = d->UH[n]/wlvl;
                           d->V[n] = d->VH[n]/wlvl;
                           d->W[n] = d->WH[n]/wlvl;
                       }
                   }
               });
}

// initial water level boxes (F 72) on the interior of the fresh level-l patches at t = 0.  The
// bicubic prolongation of a box edge overshoots (12 % for the hump of 5f), so the blocks whose
// parent stencil (parent cell +-2) or whose own cells see more than one box state are set on the
// patch grid itself: the box height in a box, else the surface of the nearest parent-level cell
// outside the boxes (piecewise constant).  All other blocks keep the prolonged state, and so do U,
// V, W, P and detadt everywhere (initial waves or currents are not touched); UH, VH, WH follow the
// new water level.  The parent gives the boxes of its 5x5 cells and the surface of its 3x3 cells.
void nhflow_amr::ini_boxes(int l)
{
    // last box that holds the point (the later box wins, as on level 0), -1 none
    auto boxof = [&](double x, double y)
    {
        int r=-1;
        for(int qn=0; qn<p0->F72; ++qn)
        if(x>=p0->F72_xs[qn] && x<p0->F72_xe[qn] && y>=p0->F72_ys[qn] && y<p0->F72_ye[qn])
        r=qn;
        return r;
    };

    block_down_if(l,34,7610+l,[](reefamr_patch *q) { return q->fresh; },
               [&](const reefamr_block &B, int key, double *v)
               {
                   lexer *q = glex(B.g);
                   fdm_nhf *dq = gfd(B.g);
                   for(int a=-2; a<=2; ++a)
                   for(int b=-2; b<=2; ++b)
                   v[(a+2)*5+(b+2)] = boxof(cxc(q,B.ic+a,B.jc+b),cyc(q,B.ic+a,B.jc+b));
                   for(int a=-1; a<=1; ++a)
                   for(int b=-1; b<=1; ++b)
                   v[25+(a+1)*3+(b+1)] = dq->eta(B.ic+a,B.jc+b);
               },
               [&](reefamr_patch *qq, int id, int k, const double *v)
               {
                   if(!qq->fresh)
                   return;

                   nhflow_amr_patch &c = *NP(qq);
                   lexer *pp = c.pp;
                   fdm_nhf *d = c.d;
                   const int K = pp->knoz;
                   const int nby = c.ny>>sy;
                   const int bi = k/nby, bj = k%nby;
                   auto cbox = [&](int a, int b) { return (int)v[(a+2)*5+(b+2)]; };
                   auto ceta = [&](int a, int b) { return v[25+(a+1)*3+(b+1)]; };

                   const int b0 = cbox(0,0);
                   bool aff = false;
                   for(int a=-2; a<=2; ++a)
                   for(int b=-2; b<=2; ++b)
                   if(cbox(a,b)!=b0)
                   aff = true;

                   for(int a=0; a<2; ++a)
                   for(int b=0; b<=sy; ++b)
                   if(boxof(cxc(pp,EXT+2*bi+a,EY+ychild(bj)+b),cyc(pp,EXT+2*bi+a,EY+ychild(bj)+b))!=b0)
                   aff = true;

                   if(!aff)
                   return;

                   for(int a=0; a<2; ++a)
                   for(int b=0; b<=sy; ++b)
                   {
                       const int ii = EXT+2*bi+a, jj = EY+ychild(bj)+b;
                       const int bf = boxof(cxc(pp,ii,jj),cyc(pp,ii,jj));
                       const int da = a ? 1 : -1, db = b ? 1 : -1;

                       double e = ceta(0,0);
                       if(bf>=0)
                       e = p0->F72_h[bf] - p0->F60;
                       else
                       {
                           // the parent cell, else its neighbour on the side of this child
                           const int ca[4] = {0,da,0,da}, cb[4] = {0,0,db,db};
                           for(int n=0; n<4; ++n)
                           if(cbox(ca[n],cb[n])<0)
                           {
                               e = ceta(ca[n],cb[n]);
                               break;
                           }
                       }

                       d->eta(ii,jj) = e;
                       d->WL(ii,jj) = MAX(e + d->depth(ii,jj), pp->A544);
                       for(int kk=0; kk<K; ++kk)
                       {
                           const int n = cidx(pp,ii,jj,kk);
                           d->UH[n] = d->WL(ii,jj)*d->U[n];
                           d->VH[n] = d->WL(ii,jj)*d->V[n];
                           d->WH[n] = d->WL(ii,jj)*d->W[n];
                       }
                   }
               });
}

// --------------------------------------------------------------------- restriction
// water level and surface of the covered coarse cells after the continuity part of stage s
void nhflow_amr::restrict_surface(int s)
{
    for(int l=wtop(); l>=wlo+1; --l)
    block_up(l,2,7310+l,
             [&](reefamr_patch *q, int id, int k, double *v)
             {
                 nhflow_amr_patch *c = NP(q);
                 stg F = stage_out(id,s);
                 fdm_nhf *df = c->d;
                 const int nby = c->ny>>sy;
                 const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                 slice &WLf = *F.WL;
                 if(c->pp->geo_curv==1)
                 {
                 // curvilinear grid: per area (the volume of the children)
                 double w[4];
                 cw4(c->pp,i0,j0,w);
                 v[0] = (sy==0 ? w[0]*WLf(i0,j0)+w[1]*WLf(i0+1,j0) : w[0]*WLf(i0,j0) + w[1]*WLf(i0+1,j0) + w[2]*WLf(i0,j0+1) + w[3]*WLf(i0+1,j0+1));
                 v[1] = (sy==0 ? w[0]*df->detadt(i0,j0)+w[1]*df->detadt(i0+1,j0) : w[0]*df->detadt(i0,j0) + w[1]*df->detadt(i0+1,j0) + w[2]*df->detadt(i0,j0+1) + w[3]*df->detadt(i0+1,j0+1));
                 return;
                 }
                 v[0] = (sy==0 ? 0.5*(WLf(i0,j0)+WLf(i0+1,j0)) : 0.25*(WLf(i0,j0)+WLf(i0+1,j0)+WLf(i0,j0+1)+WLf(i0+1,j0+1)));
                 v[1] = (sy==0 ? 0.5*(df->detadt(i0,j0)+df->detadt(i0+1,j0)) : 0.25*(df->detadt(i0,j0)+df->detadt(i0+1,j0)+df->detadt(i0,j0+1)+df->detadt(i0+1,j0+1)));
             },
             [&](const reefamr_block &B, int key, const double *v)
             {
                 stg C = stage_out(B.g,s);
                 fdm_nhf *dc = gfd(B.g);
                 (*C.WL)(B.ic,B.jc) = v[0];
                 dc->eta(B.ic,B.jc) = v[0] - dc->depth(B.ic,B.jc);
                 dc->detadt(B.ic,B.jc) = v[1];
             });
}

// UH, VH, WH (and with P the pressure) of the covered coarse cells, U = UH/WL
void nhflow_amr::restrict_momentum(int s, bool withP)
{
    for(int l=wtop(); l>=wlo+1; --l)
    {
        const int K = klev(l-1);
        block_up(l,3*K,7320+l,
                 [&](reefamr_patch *q, int id, int k, double *v)
                 {
                     nhflow_amr_patch *c = NP(q);
                     lexer *pp = c->pp;
                     stg F = stage_out(id,s);
                     const int nby = c->ny>>sy;
                     const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                     const int fz = pp->knoz/K;

                     for(int kk=0; kk<K; ++kk)
                     {
                         // the 2x2 children of the layer, with G 6 in both fine layers (equal thickness)
                         auto avg = [&](const double *f)
                         {
                             double r = 0.0;
                             // curvilinear grid: per area of the children
                             if(pp->geo_curv==1)
                             {
                                 double w[4];
                                 cw4(pp,i0,j0,w);
                                 for(int kf=fz*kk; kf<fz*kk+fz; ++kf)
                                 r += (sy==0 ? w[0]*f[cidx(pp,i0,j0,kf)]+w[1]*f[cidx(pp,i0+1,j0,kf)] : w[0]*f[cidx(pp,i0,j0,kf)]+w[1]*f[cidx(pp,i0+1,j0,kf)]+w[2]*f[cidx(pp,i0,j0+1,kf)]+w[3]*f[cidx(pp,i0+1,j0+1,kf)]);
                                 return r/double(fz);
                             }
                             if(sy==0)      // 2D grid: the 2x1 children
                             {
                                 for(int kf=fz*kk; kf<fz*kk+fz; ++kf)
                                 r += f[cidx(pp,i0,j0,kf)]+f[cidx(pp,i0+1,j0,kf)];
                                 return 0.5*r/double(fz);
                             }
                             for(int kf=fz*kk; kf<fz*kk+fz; ++kf)
                             r += f[cidx(pp,i0,j0,kf)]+f[cidx(pp,i0+1,j0,kf)]+f[cidx(pp,i0,j0+1,kf)]+f[cidx(pp,i0+1,j0+1,kf)];
                             return 0.25*r/double(fz);
                         };
                         v[3*kk] = avg(F.UH);
                         v[3*kk+1] = avg(F.VH);
                         v[3*kk+2] = avg(F.WH);
                     }
                 },
                 [&](const reefamr_block &B, int key, const double *v)
                 {
                     lexer *q = glex(B.g);
                     stg C = stage_out(B.g,s);
                     fdm_nhf *dc = gfd(B.g);
                     const double wl = (*C.WL)(B.ic,B.jc);
                     const double wlvl = fabs(wl)>p0->A544 ? wl : 1.0e20;
                     for(int kk=0; kk<K; ++kk)
                     {
                         const int n = cidx(q,B.ic,B.jc,kk);
                         C.UH[n] = v[3*kk];
                         C.VH[n] = v[3*kk+1];
                         C.WH[n] = v[3*kk+2];
                         dc->U[n] = C.UH[n]/wlvl;
                         dc->V[n] = C.VH[n]/wlvl;
                         dc->W[n] = C.WH[n]/wlvl;
                     }
                 });
    }

    if(withP)
    restrict_col([&](int g) -> double* { return gfd(g)->P; });
}

// level-0 ghost cells and partition halo after the restriction
void nhflow_amr::halo0(ghostcell *pgc, int s)
{
    stg C = stage_out(-1,s);
    pgc->gcsl_start4(p0,*C.WL,gcval_eta);
    pgc->gcsl_start4(p0,d0->eta,gcval_eta);
    pgc->start4V(p0,C.UH,14);
    pgc->start4V(p0,C.VH,15);
    pgc->start4V(p0,C.WH,16);
    pgc->start4V(p0,d0->U,10);
    pgc->start4V(p0,d0->V,11);
    pgc->start4V(p0,d0->W,12);
    pgc->start7P(p0,d0->P,540);
}

// --------------------------------------------------------------------- regridding
void nhflow_amr::regrid_prepare(ghostcell *pgc)
{
    if(nblk>0)
    fill_level0(pgc,-1);
    for(int l=1; l<=maxlev; ++l)
    fill_stage(pgc,l,-1);
}

void nhflow_amr::regrid_static(ghostcell *pgc)
{
}

// state of the new patches, coarse to fine: bed and depth, eta, water level, velocities,
// momentum and pressure from the parent level, then the sigma grid and the derived fields
void nhflow_amr::regrid_state(ghostcell *pgc, vector<reefamr_patch*> &oldP)
{
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

    for(int l=1; l<=maxlev; ++l)
    {
        // bed: linear from the parent (on its rank); curvilinear grid: shifted so that the mean of the
        // children per area is the parent's bed (the lake at rest stays at rest: the restriction of
        // the water level is the mean per area)
        block_down_if(l,5,7620+l,[](reefamr_patch *q) { return q->fresh; },
                   [&](const reefamr_block &B, int key, double *v)
                   {
                       for(int a=0; a<2; ++a)
                       for(int b=0; b<=sy; ++b)
                       v[2*a+b] = plin(gfd(B.g)->bed,glex(B.g),B.ic,B.jc,a==0?-1:1,(sy==0) ? 0 : (b==0?-1:1));
                       v[4] = gfd(B.g)->bed(B.ic,B.jc);
                   },
                   [&](reefamr_patch *q, int id, int k, const double *v)
                   {
                       if(!q->fresh)
                       return;
                       const int nby = q->ny>>sy;
                       const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                       double sh = 0.0;
                       if(q->pp->geo_curv==1)
                       {
                           double w[4];
                           cw4(q->pp,i0,j0,w);
                           sh = v[4] - (w[0]*v[0] + w[1]*v[2] + w[2]*v[1] + w[3]*v[3]);
                       }
                       for(int a=0; a<2; ++a)
                       for(int b=0; b<=sy; ++b)
                       NP(q)->d->bed(i0+a,j0+b) = (q->pp->geo_curv==1) ? v[2*a+b] + sh : v[2*a+b];
                   });
        fill_bed(pgc,l);

        for(int id : lev[l])
        if(P[id]->fresh)
        {
            nhflow_amr_patch *c = NP(id);
            lexer *pp = c->pp;
            fdm_nhf *d = c->d;
            pscope ps(pgc,d,d0);

            pgc->gcsl_start4(pp,d->bed,50);
            for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
            for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
            {
                d->depth(ii,jj) = pp->wd - d->bed(ii,jj);
                pp->bed[lij(pp,ii,jj)] = d->bed(ii,jj);
            }
            pgc->gcsl_start4(pp,d->depth,50);
        }

        prolong_patch(pgc,l);

        // G 30: the flags from the prolonged water level
        if(shore)
        for(int id : lev[l])
        if(P[id]->fresh)
        patch_flags(*NP(id));

        // where an old patch of the same level was (on any rank), its state is taken over (G 30:
        // with its flags, which carry the history of the wetting and drying)
        {
            const int K = klev(l);
            // S 10: and the bed (the last value)
            const int nv0 = 3 + 2 + 6*K + K+1 + 1;
            // X 15 8: and the state of the hull (nhflow_hullbody::amr_get: sides with their hysteresis, thin body,
            // Rhie-Chow face corrections), so that a rebuilt patch continues the hull of the old one
            const int nh = hb ? (K + (K+1) + (K+2) + 2*K) : 0;
            const int nv = nv0 + nh;
            old_run(l,nv,7630+l,oldP,
                    [&](reefamr_patch *o, int io, int jo, double *v)
                    {
                        fdm_nhf *od = NP(o)->d;
                        lexer *op = o->pp;
                        v[0] = op->wet[lij(op,io,jo)];
                        v[1] = op->deep[lij(op,io,jo)];
                        v[2] = od->eta(io,jo);
                        v[3] = od->WL(io,jo);
                        v[4] = od->detadt(io,jo);
                        for(int kk=0; kk<K; ++kk)
                        {
                            const int m = cidx(op,io,jo,kk);
                            v[5+6*kk] = od->U[m];
                            v[6+6*kk] = od->V[m];
                            v[7+6*kk] = od->W[m];
                            v[8+6*kk] = od->UH[m];
                            v[9+6*kk] = od->VH[m];
                            v[10+6*kk] = od->WH[m];
                        }
                        for(int kk=0; kk<=K; ++kk)
                        v[5+6*K+kk] = od->P[fidx(op,io,jo,kk)];
                        v[5+6*K+K+1] = od->bed(io,jo);
                        if(nh>0)
                        {
                        if(od->hullbody!=nullptr)
                        od->hullbody->amr_get(op,od,io,jo,v+nv0);
                        else
                        for(int n=0; n<nh; ++n)
                        v[nv0+n] = 0.0;
                        }
                    },
                    [&](reefamr_patch *c, int id, int ii, int jj, const double *v)
                    {
                        fdm_nhf *d = NP(c)->d;
                        lexer *pp = c->pp;
                        if(shore)
                        {
                            pp->wet[lij(pp,ii,jj)] = (int)v[0];
                            pp->deep[lij(pp,ii,jj)] = (int)v[1];
                        }
                        d->eta(ii,jj) = v[2];
                        d->WL(ii,jj) = v[3];
                        d->detadt(ii,jj) = v[4];
                        for(int kk=0; kk<K; ++kk)
                        {
                            const int n = cidx(pp,ii,jj,kk);
                            d->U[n] = v[5+6*kk];
                            d->V[n] = v[6+6*kk];
                            d->W[n] = v[7+6*kk];
                            d->UH[n] = v[8+6*kk];
                            d->VH[n] = v[9+6*kk];
                            d->WH[n] = v[10+6*kk];
                        }
                        for(int kk=0; kk<=K; ++kk)
                        d->P[fidx(pp,ii,jj,kk)] = v[5+6*K+kk];
                        // S 10: the bed of the old patch (it evolved), the depth with it
                        if(pp->S10==1)
                        {
                        d->bed(ii,jj) = v[5+6*K+K+1];
                        d->depth(ii,jj) = pp->wd - d->bed(ii,jj);
                        pp->bed[lij(pp,ii,jj)] = d->bed(ii,jj);
                        }
                        if(nh>0 && d->hullbody!=nullptr)
                        d->hullbody->amr_set(pp,d,ii,jj,v+nv0);
                    });
        }

        // at t = 0 the initial water level boxes (F 72) on the patch grid itself, as
        // nhflow_fsf_ini on level 0 (the interpolated coarse box would overshoot at its edges);
        // level 0 takes their mean in regrid_finish
        if(p0->count==0 && p0->F72>0)
        {
            ini_boxes(l);
            if(shore)
            for(int id : lev[l])
            if(P[id]->fresh)
            patch_flags(*NP(id));
        }

        fill_stage(pgc,l,-1);

        for(int id : lev[l])
        if(P[id]->fresh)
        {
            nhflow_amr_patch *c = NP(id);
            lexer *pp = c->pp;
            fdm_nhf *d = c->d;
            pscope ps(pgc,d,d0);

            for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
            for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
            {
                d->eta_n(ii,jj) = d->eta(ii,jj);
                pp->WL[lij(pp,ii,jj)] = d->WL(ii,jj);
            }

            // sigma grid and derived fields as driver_ini_nhflow
            c->pmom->sigma_ini(pp,d,pgc,d->eta);
            c->pmom->inidisc(pp,d,pgc,c->pfsf);
            c->pfsf->kinematic_fsf(pp,d,d->U,d->V,d->W,d->eta);

            // floating bodies: the hull on the new sigma grid
            body_patch(pgc,*c);

            // fixed solids: the solids on the new sigma grid
            solid_patch(pgc,*c);
        }

        // S 10: the sediment of the fresh patches starts from their bed
        for(int id : lev[l])
        if(P[id]->fresh)
        sed_patch_ini(pgc,*NP(id));

        // A 560: k, epsilon/omega and the eddy viscosity of the fresh patches
        turb_regrid(pgc,l,oldP);
    }
}

// level 0 consistent with the patches
void nhflow_amr::regrid_finish(ghostcell *pgc, int old_total)
{
    for(int l=maxlev; l>=1; --l)
    block_up(l,1,7330+l,
             [&](reefamr_patch *q, int id, int k, double *v)
             {
                 fdm_nhf *df = NP(q)->d;
                 const int nby = q->ny>>sy;
                 const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                 // curvilinear grid: per area of the children
                 if(q->pp->geo_curv==1)
                 {
                 double w[4];
                 cw4(q->pp,i0,j0,w);
                 v[0] = (sy==0 ? w[0]*df->WL(i0,j0)+w[1]*df->WL(i0+1,j0) : w[0]*df->WL(i0,j0) + w[1]*df->WL(i0+1,j0) + w[2]*df->WL(i0,j0+1) + w[3]*df->WL(i0+1,j0+1));
                 return;
                 }
                 v[0] = (sy==0 ? 0.5*(df->WL(i0,j0)+df->WL(i0+1,j0)) : 0.25*(df->WL(i0,j0)+df->WL(i0+1,j0)+df->WL(i0,j0+1)+df->WL(i0+1,j0+1)));
             },
             [&](const reefamr_block &B, int key, const double *v)
             {
                 fdm_nhf *dc = gfd(B.g);
                 dc->WL(B.ic,B.jc) = v[0];
                 dc->eta(B.ic,B.jc) = v[0] - dc->depth(B.ic,B.jc);
             });

    // momentum, velocities and pressure: the end-of-step arrays are the last stage's output
    restrict_momentum(mom0->stages()-1,true);

    if(patches_total>0 || old_total>0)
    {
    restrict_flags(pgc,-1);
    halo0(pgc,mom0->stages()-1);
    }
}

// --------------------------------------------------------------------- flux matching
// face fluxes of grid id (0: level 0, n: patch n-1), ipol 1-3 momentum, 4 continuity
void nhflow_amr::flux_hook(lexer *p, fdm_nhf *d, int id, int ipol, double *Fx, double *Fy)
{
    const int K = p->knoz;

    // sharp thin bodies (X 15 8 on every grid): their wall and local fluxes first, as the hook they set on a
    // single grid (driver_ini_nhflow), then the face fluxes of the patch box are recorded and matched
    if(d->thinbody!=nullptr)
    d->thinbody->flux_hook(p,d,id,ipol,Fx,Fy);


    // a patch records its box faces
    if(id>0)
    {
        nhflow_amr_patch &c = *NP(id-1);
        const int il = EXT-1, ih = EXT+c.nx-1;
        const int jl = EY-1, jh = EY+c.ny-1;

        for(int side=0; side<4; ++side)
        {
            const int ns = side<2 ? c.ny : c.nx;
            c.rec[ipol][side].resize((size_t)ns*K);
            if(ipol==4)
            c.rec[0][side].resize(ns);
        }

        for(int r=0; r<c.ny; ++r)
        for(int k=0; k<K; ++k)
        {
            c.rec[ipol][0][r*K+k] = Fx[cidx(p,il,EY+r,k)];
            c.rec[ipol][1][r*K+k] = Fx[cidx(p,ih,EY+r,k)];
        }
        for(int r=0; r<c.nx; ++r)
        for(int k=0; k<K; ++k)
        {
            c.rec[ipol][2][r*K+k] = Fy[cidx(p,EXT+r,jl,k)];
            c.rec[ipol][3][r*K+k] = Fy[cidx(p,EXT+r,jh,k)];
        }

        if(ipol==4)
        {
            for(int r=0; r<c.ny; ++r)
            {
                c.rec[0][0][r] = d->dfx(il,EY+r);
                c.rec[0][1][r] = d->dfx(ih,EY+r);
            }
            for(int r=0; r<c.nx; ++r)
            {
                c.rec[0][2][r] = d->dfy(EXT+r,jl);
                c.rec[0][3][r] = d->dfy(EXT+r,jh);
            }
        }
    }

    // G 7 1: the fine register of the patch; the coarse faces keep their own flux (the fine one is
    // not there yet, refluxing after the fine steps) and add it to the coarse register, with the
    // fine face depth as in the synchronous coupling
    if(sub==1)
    {
        const double wdt = rkw(hstage)*p->dt;

        if(id>0)
        {
            nhflow_amr_patch &c = *NP(id-1);
            for(int side=0; side<4; ++side)
            {
                vector<double> &R = c.rec[ipol][side];
                vector<double> &G = c.freg[ipol][side];
                for(size_t n=0; n<R.size() && n<G.size(); ++n)
                G[n] += wdt*R[n];
            }
        }

        if(id>=(int)match.size())
        return;

        for(size_t n=0; n<match[id].size(); ++n)
        {
            const reefamr_match &m = match[id][n];
            if(ipol==4)
            {
                nhflow_amr_patch &c = *NP(m.child);
                const double df = 0.5*(c.rec[0][m.side][m.r]+c.rec[0][m.side][m.r+fr1(m.side)]);
                if(m.dir==0)
                d->dfx(m.fi,m.fj) = df;
                else
                d->dfy(m.fi,m.fj) = df;
            }
            double *cr = &cregL[id][(4*n+ipol-1)*K];
            for(int k=0; k<K; ++k)
            cr[k] += wdt*((m.dir==0) ? Fx[cidx(p,m.fi,m.fj,k)] : Fy[cidx(p,m.fi,m.fj,k)]);
        }

        for(size_t n=0; n<rmatch[id].size(); ++n)
        {
            const reefamr_match &m = rmatch[id][n];
            if(ipol==4 && id<(int)rval.size())
            {
                const double df = rval[id][n*NF+4*K];
                if(m.dir==0)
                d->dfx(m.fi,m.fj) = df;
                else
                d->dfy(m.fi,m.fj) = df;
            }
            double *cr = &cregR[id][(4*n+ipol-1)*K];
            for(int k=0; k<K; ++k)
            cr[k] += wdt*((m.dir==0) ? Fx[cidx(p,m.fi,m.fj,k)] : Fy[cidx(p,m.fi,m.fj,k)]);
        }
        return;
    }

    if(id>=(int)match.size())
    return;

    // coarse faces next to the patches of this rank: mean of the two fine faces, layer by layer
    // (G 6: of the two fine faces in both fine layers of the coarse layer; the fluxes are per
    // unit sigma and the halves have the same thickness)
    for(auto &m : match[id])
    {
        nhflow_amr_patch &c = *NP(m.child);
        const vector<double> &R = c.rec[ipol][m.side];
        const int Kf = c.pp->knoz;
        const int fz = Kf/K;

        for(int k=0; k<K; ++k)
        {
            double val;
            if(fz==1)
            val = 0.5*(R[m.r*Kf+k]+R[(m.r+fr1(m.side))*Kf+k]);
            else
            val = 0.25*(R[m.r*Kf+2*k]+R[m.r*Kf+2*k+1]+R[(m.r+fr1(m.side))*Kf+2*k]+R[(m.r+fr1(m.side))*Kf+2*k+1]);
            // multi-block: the coarse face in the neighbour block, its normal against the fine faces
            if(m.sgn<0)
            val = -val;
            if(m.dir==0)
            Fx[cidx(p,m.fi,m.fj,k)] = val;
            else
            Fy[cidx(p,m.fi,m.fj,k)] = val;
        }

        if(ipol==4)
        {
            const double df = 0.5*(c.rec[0][m.side][m.r]+c.rec[0][m.side][m.r+fr1(m.side)]);
            if(m.dir==0)
            d->dfx(m.fi,m.fj) = df;
            else
            d->dfy(m.fi,m.fj) = df;
        }
    }

    // coarse faces next to patches of other ranks
    if(id<(int)rval.size())
    for(size_t n=0; n<rmatch[id].size(); ++n)
    {
        const reefamr_match &m = rmatch[id][n];
        const double *v = &rval[id][n*NF];

        for(int k=0; k<K; ++k)
        {
            const double val = (m.sgn<0) ? -v[(ipol-1)*K+k] : v[(ipol-1)*K+k];
            if(m.dir==0)
            Fx[cidx(p,m.fi,m.fj,k)] = val;
            else
            Fy[cidx(p,m.fi,m.fj,k)] = val;
        }

        if(ipol==4)
        {
            if(m.dir==0)
            d->dfx(m.fi,m.fj) = v[4*K];
            else
            d->dfy(m.fi,m.fj) = v[4*K];
        }
    }
}

// fine face values of the level-l patches on the partition edges, to the rank of the coarse cell
void nhflow_amr::exchange_fluxes(int l)
{
    // layers of the coarse side and of the patches of level l
    const int K = klev(l-1);
    const int Kf = klev(l);

    rval.resize(rmatch.size());
    for(size_t g=0; g<rmatch.size(); ++g)
    rval[g].resize(rmatch[g].size()*NF);

    // G 31 1 exchanges after phase_F as well: a fresh patch has no momentum records yet (zero
    // until its phase_M)
    for(int id : lev[l])
    for(int side=0; side<4; ++side)
    {
        nhflow_amr_patch *c = NP(id);
        const size_t ns = side<2 ? c->ny : c->nx;
        for(int ip=1; ip<=4; ++ip)
        if(c->rec[ip][side].size()<ns*Kf)
        c->rec[ip][side].resize(ns*Kf,0.0);
        if(c->rec[0][side].size()<ns)
        c->rec[0][side].resize(ns,0.0);
    }

    face_run(l,NF,7200+l,
             [&](reefamr_patch *q, int side, int r, double *sb)
             {
                 nhflow_amr_patch *c = NP(q);
                 for(int ip=1; ip<=4; ++ip)
                 {
                     const vector<double> &R = c->rec[ip][side];
                     for(int k=0; k<K; ++k)
                     {
                         if(Kf==K)
                         sb[(ip-1)*K+k] = 0.5*(R[r*K+k]+R[(r+fr1(side))*K+k]);
                         else
                         sb[(ip-1)*K+k] = 0.25*(R[r*Kf+2*k]+R[r*Kf+2*k+1]+R[(r+fr1(side))*Kf+2*k]+R[(r+fr1(side))*Kf+2*k+1]);
                     }
                 }
                 sb[4*K] = 0.5*(c->rec[0][side][r]+c->rec[0][side][r+fr1(side)]);
                 for(int q=4*K+1; q<NF; ++q)
                 sb[q] = 0.0;
             },
             [&](reefamr_match &m, const double *v)
             {
                 for(size_t g=0; g<rmatch.size(); ++g)
                 if(!rmatch[g].empty() && &m>=&rmatch[g][0] && &m<=&rmatch[g].back())
                 {
                     const size_t n = &m-&rmatch[g][0];
                     for(int q=0; q<NF; ++q)
                     rval[g][n*NF+q] = v[q];
                 }
             });
}

// --------------------------------------------------------------------- time stepping
void nhflow_amr::step(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_momentum_func *mom, nhflow_stage_obj &S0)
{
    // A 534 1: hydrostatic semi-implicit free surface, one eta system over the blocks (nhflow_amr_si.cpp)
    if(p->A534==1)
    {
        si_step(p,d,pgc,mom,S0);
        return;
    }

    const int ns = mom->stages();

    if(!active())
    {
        mom->step_begin(p,d,pgc,S0);
        for(int s=0; s<ns; ++s)
        {
        mom->phase_F(p,d,pgc,S0,s);
        mom->phase_M(p,d,pgc,S0,s);
        mom->phase_P(p,d,pgc,S0,s);
        mom->phase_E(p,d,pgc,S0,s);
        }

        // adaptive flags: the first patches appear when the flags do
        regrid_step(p,pgc);
        return;
    }

    // G 7 1: subcycling
    if(sub==1)
    {
        sub_step(p,d,pgc,mom,S0);
        regrid_step(p,pgc);
        fill_out(pgc);
        // multi-block: the ring and the cells around the blocks from the end of the step
        if(nblk>0)
        fill_level0(pgc,-1);
        return;
    }

    {
    comms_off guard(pgc);
    for(auto q : P)
    {
        nhflow_amr_patch *c = NP(q);
        lexer *pp = c->pp;
        pp->dt = p->dt;
        pp->dt_old = p->dt_old;
        pp->simtime = p->simtime;
        pp->count = p->count;
        pp->amr_follow(*p);
        pscope ps(pgc,c->d,d0);
        c->pmom->step_begin(pp,c->d,pgc,c->S);
    }
    }

    // A 560: the turbulence of the blocks and the patches (level 0 has made its step before the momentum)
    turb_step(p,pgc);

    mom->step_begin(p,d,pgc,S0);
    S0p = &S0;

    for(int s=0; s<ns; ++s)
    {
        cur_stage = s;
        double t0 = MPI_Wtime();

        // stage input around the blocks and the patches, coarse to fine
        if(nblk>0)
        fill_level0(pgc,s);
        for(int l=1; l<=maxlev; ++l)
        fill_stage(pgc,l,s);
        tm[0] += MPI_Wtime()-t0;

        // continuity and momentum fluxes, finest first: the coarser grids take the recorded fine fluxes
        if(!cdiff)
        {
        for(int l=maxlev; l>=1; --l)
        {
            t0 = MPI_Wtime();
            {
            comms_off guard(pgc);
            for(int n : lev[l])
            {
                nhflow_amr_patch *c = NP(n);
                pscope ps(pgc,c->d,d0);
                c->pmom->phase_F(c->pp,c->d,pgc,c->S,s);
                c->pmom->phase_M(c->pp,c->d,pgc,c->S,s);
            }
            }
            tm[1] += MPI_Wtime()-t0;

            t0 = MPI_Wtime();
            exchange_fluxes(l);
            tm[2] += MPI_Wtime()-t0;
        }

        if(brk0)
        {
        // multi-block breaking: the continuity of level 0, its breaking, then its momentum
        blocks_phase(pgc,s,0);
        mom->phase_F(p,d,pgc,S0,s);
        t0 = MPI_Wtime();
        breaking_level0(pgc,-1);
        tm[1] += MPI_Wtime()-t0;
        blocks_phase(pgc,s,1);
        mom->phase_M(p,d,pgc,S0,s);
        }
        else
        {
        blocks_phase(pgc,s,0);
        blocks_phase(pgc,s,1);
        mom->phase_F(p,d,pgc,S0,s);
        mom->phase_M(p,d,pgc,S0,s);
        }
        }
        else
        {
        // G 31 1: continuity of all grids, the implicit diffusion of all grids in one solve,
        // then their momentum (the momentum fluxes of the finer grids still before the coarser)
        for(int l=maxlev; l>=1; --l)
        {
            t0 = MPI_Wtime();
            {
            comms_off guard(pgc);
            for(int n : lev[l])
            {
                nhflow_amr_patch *c = NP(n);
                pscope ps(pgc,c->d,d0);
                c->pmom->phase_F(c->pp,c->d,pgc,c->S,s);
            }
            }
            tm[1] += MPI_Wtime()-t0;

            t0 = MPI_Wtime();
            exchange_fluxes(l);
            tm[2] += MPI_Wtime()-t0;
        }
        blocks_phase(pgc,s,0);
        mom->phase_F(p,d,pgc,S0,s);

        t0 = MPI_Wtime();
        if(brk_split)
        breaking_patches(pgc);
        tm[1] += MPI_Wtime()-t0;

        t0 = MPI_Wtime();
        diff_solve(p,pgc,s);
        tdiff += MPI_Wtime()-t0;

        for(int l=maxlev; l>=1; --l)
        {
            t0 = MPI_Wtime();
            {
            comms_off guard(pgc);
            for(int n : lev[l])
            {
                nhflow_amr_patch *c = NP(n);
                pscope ps(pgc,c->d,d0);
                nhflow_stage_obj S = c->S;
                S.pdiff = dkeep;
                c->pmom->phase_M(c->pp,c->d,pgc,S,s);
            }
            }
            tm[1] += MPI_Wtime()-t0;

            t0 = MPI_Wtime();
            exchange_fluxes(l);
            tm[2] += MPI_Wtime()-t0;
        }

        {
        comms_off guard(pgc);
        for(int n : lev[0])
        {
            nhflow_amr_patch *c = NP(n);
            pscope ps(pgc,c->d,d0);
            nhflow_stage_obj S = c->S;
            S.pdiff = dkeep;
            c->pmom->phase_M(c->pp,c->d,pgc,S,s);
        }
        }

        nhflow_stage_obj S = S0;
        S.pdiff = dkeep;
        mom->phase_M(p,d,pgc,S,s);
        }

        t0 = MPI_Wtime();
        {
        comms_off guard(pgc);
        restrict_surface(s);
        restrict_momentum(s,false);
        }
        restrict_flags(pgc,s);
        halo0(pgc,s);
        tm[3] += MPI_Wtime()-t0;

        // pressure projection on all grids; level 0 first: it advances the floating bodies,
        // whose forcing the patches take
        mom->phase_P1(p,d,pgc,S0,s);
        {
        comms_off guard(pgc);
        for(auto q : P)
        {
        pscope ps(pgc,NP(q)->d,d0);
        NP(q)->pmom->phase_P1(q->pp,NP(q)->d,pgc,NP(q)->S,s);
        }
        }

        press_solve(p,pgc,s);

        // level 0 last: the loads on the floating bodies sample the patches
        {
        comms_off guard(pgc);
        for(auto q : P)
        {
        pscope ps(pgc,NP(q)->d,d0);
        NP(q)->pmom->phase_P2(q->pp,NP(q)->d,pgc,NP(q)->S,s);
        }
        }
        mom->phase_P2(p,d,pgc,S0,s);

        t0 = MPI_Wtime();
        {
        comms_off guard(pgc);
        restrict_momentum(s,true);
        }
        halo0(pgc,s);
        tm[3] += MPI_Wtime()-t0;

        // relaxation zones and ghost cells
        {
        comms_off guard(pgc);
        for(auto q : P)
        {
        pscope ps(pgc,NP(q)->d,d0);
        NP(q)->pmom->phase_E(q->pp,NP(q)->d,pgc,NP(q)->S,s);
        }
        }
        mom->phase_E(p,d,pgc,S0,s);
    }
    cur_stage = -1;

    regrid_step(p,pgc);

    // the cells around the patches filled from a patch of the same level: the end of the step (the output)
    {
        const double t0 = MPI_Wtime();
        fill_out(pgc);
        tm[0] += MPI_Wtime()-t0;
    }

    // multi-block: the ring and the cells around the blocks from the end of the step (the time step,
    // the output; after a regrid the restricted values)
    if(nblk>0)
    {
        const double t0 = MPI_Wtime();
        fill_level0(pgc,-1);
        tm[0] += MPI_Wtime()-t0;
    }
}

// the zone around the body moves, the flags move: new patches every G 2 steps, from the state
// at the end of the step
void nhflow_amr::regrid_step(lexer *p, ghostcell *pgc)
{
    // G 7 1: G 2 counts steps of the finest level (a level-0 step is 2^G1 of them)
    const int rint = (sub==1) ? MAX(regrid_int>>maxlev,1) : regrid_int;

    // (replace_due: G 41 0 with a static refinement, the placement once more with the measured model)
    if(maxlev>0 && ((regrid_int>0 && p->count%rint==0) || replace_due(p)))
    {
        const double t0 = MPI_Wtime();
        regrid(p,pgc,false);
        tm[7] += MPI_Wtime()-t0;
    }

    // load statistics of the ranks (G 43)
    load_step(p);
}

// --------------------------------------------------------------------- floating bodies
void nhflow_amr::zone_bodies(vector<sixdof_obj*> &obj)
{
    if(b6!=nullptr)
    for(int nb=0; nb<b6->objects(); ++nb)
    obj.push_back(b6->object(nb));
}

// the finest local grid whose interior holds (x,y), -1: level 0
int nhflow_amr::finest_at(double x, double y)
{
    // multi-block grid: the pieces of the blocks 1.. are level-0 grids of their own (a point in the ring of the
    // native grid lies in one of them)
    int g=-1, l=(nblk>0) ? -1 : 0;
    int haveX=0;
    bool inX=false;
    double X=0.0, Y=0.0;
    for(int n=0; n<(int)P.size(); ++n)
    {
        reefamr_patch *c = P[n];
        lexer *pp = c->pp;
        if(c->lev<=l)
        continue;
        // curvilinear grid (G 50): the logical range of the patch, as finest_rank (the same patch for the
        // ownership of a hull triangle and its loads)
        if(pp->geo_curv==1 && p0->geo_mb==0)
        {
            if(!haveX)
            {
                haveX = 1;
                inX = p0->curv_logical(x,y,X,Y);
            }
            const double r = double(1<<c->lev);
            if(inX && X*r>=double(c->I0) && X*r<double(c->I0+c->nx) && Y*r>=double(c->J0) && Y*r<double(c->J0+c->ny))
            {
                g=n;
                l=c->lev;
            }
            continue;
        }
        if(pp->geo_curv==1)
        {
            int a,b;
            if(pp->geo_cell(x,y,a,b) && a>=EXT && a<EXT+c->nx && b>=EY && b<EY+c->ny)
            {
                g=n;
                l=c->lev;
            }
            continue;
        }
        if(x>=pp->XN[EXT+marge] && x<pp->XN[EXT+c->nx+marge] && (sy==0 || (y>=pp->YN[EY+marge] && y<pp->YN[EY+c->ny+marge])))
        {
            g=n;
            l=c->lev;
        }
    }
    return g;
}

// the hull on the sigma grid of a patch: level set FB, solid flags
void nhflow_amr::body_patch(ghostcell *pgc, nhflow_amr_patch &c)
{
    if(b6==nullptr)
    return;

    pscope ps(pgc,c.d,d0);
    static_cast<nhflow_amr_6dof*>(c.p6dof)->body(c.pp,c.d,pgc);
}

// --------------------------------------------------------------------- fixed solids
// the refinement boxes around the fixed solids: the x-y extent of the surface triangles of every
// solid (the level-0 geometry, the same on all ranks) grown by G 12, refined to the finest level
void nhflow_amr::solid_zone(lexer *p, ghostcell *pgc)
{
    nhflow_forcing *f0 = mom0->pmfrc;
    if(maxlev<1 || f0==nullptr)
    return;

    if(p->G12==0)
    {
        if(p->mpirank==0)
        cout<<"NHFLOW AMR: fixed solids without G 12: no refinement around them (G 10 boxes)"<<endl;
        return;
    }

    const double r = MAX(p->G12_r,0.0);
    for(int n=0; n<f0->solid_entities(); ++n)
    {
        double b[4];
        f0->entity_box(n,b);
        b[0] = pgc->globalmin(b[0]);
        b[1] = pgc->globalmax(b[1]);
        b[2] = pgc->globalmin(b[2]);
        b[3] = pgc->globalmax(b[3]);
        if(b[0]>b[1])
        continue;

        // 2D grid: the strip in y
        if(p->j_dir==0)
        {
            b[2] = -1.0e20;
            b[3] = 1.0e20;
        }
        par.rbox.push_back(b[0]-r);
        par.rbox.push_back(b[1]+r);
        par.rbox.push_back(b[2]-r);
        par.rbox.push_back(b[3]+r);

        if(p->mpirank==0)
        cout<<"NHFLOW AMR: solid "<<n+1<<" refined in x "<<b[0]-r<<" "<<b[1]+r<<(p->j_dir==1 ? ", y " : "")
            <<(p->j_dir==1 ? std::to_string(b[2]-r)+" "+std::to_string(b[3]+r) : std::string(""))<<endl;
    }
}

// the solids on the sigma grid of a fresh patch: level set SOLID, solid flags (the forcing of the
// stages casts them again on the moving sigma grid)
void nhflow_amr::solid_patch(ghostcell *pgc, nhflow_amr_patch &c)
{
    if(!solids || maxlev<1)
    return;

    pscope ps(pgc,c.d,d0);
    c.pdf->forcing_ini(c.pp,c.d,pgc);
}

// the force objects of patch c: one per force box (P 81) its interior reaches into, without output
// files, as long as the patch lives
void nhflow_amr::force_objects(nhflow_amr_patch *c)
{
    lexer *pp = c->pp;
    const int m = marge;
    if(p0->P81<=0)
    return;

    const double x0 = pp->XN[EXT+m], x1 = pp->XN[EXT+c->nx+m];
    const double y0 = (sy==1) ? pp->YN[EY+m] : -1.0e20, y1 = (sy==1) ? pp->YN[EY+c->ny+m] : 1.0e20;

    c->forces.assign(p0->P81,nullptr);
    for(int ID=0; ID<p0->P81; ++ID)
    {
        if(x1<p0->P81_xs[ID] || x0>p0->P81_xe[ID] || y1<p0->P81_ys[ID] || y0>p0->P81_ye[ID])
        continue;

        c->forces[ID] = new nhflow_force(pp,c->d,pgc0,ID,false);
        c->forces[ID]->amr_attach(this,c);
    }
}

// a piece of the surface at (x,y) of grid g (nullptr: level 0, else the patch): in the interior of
// the grid (level 0: the cells of this rank) and no patch of the next finer level there
bool nhflow_amr::force_leaf(const void *g, double x, double y) const
{
    nhflow_amr *self = const_cast<nhflow_amr*>(this);
    const int m = marge;

    lexer *q;
    int l, i, j, I, J;
    if(g==nullptr)
    {
        q = p0;
        l = 0;
        if(x<q->XN[m] || x>=q->XN[q->knox+m] || (sy==1 && (y<q->YN[m] || y>=q->YN[q->knoy+m])))
        return true;
        i = q->posc_i(x);
        j = (sy==1) ? q->posc_j(y) : 0;
        I = i+O0i;
        J = j+O0j;
    }
    else
    {
        const reefamr_patch *c = static_cast<const reefamr_patch*>(g);
        q = c->pp;
        l = c->lev;
        if(x<q->XN[EXT+m] || x>=q->XN[EXT+c->nx+m])
        return false;
        if(sy==1 && (y<q->YN[EY+m] || y>=q->YN[EY+c->ny+m]))
        return false;
        i = q->posc_i(x);
        j = (sy==1) ? q->posc_j(y) : EY;
        I = i-EXT+c->I0;
        J = j-EY+c->J0;
    }

    if(l>=maxlev || patches_total==0)
    return true;

    // the child of the cell that holds the point (the levels halve the cells of their parents)
    const int Ic = 2*I + ((x-q->XN[i+m]) >= 0.5*(q->XN[i+1+m]-q->XN[i+m]) ? 1 : 0);
    const int Jc = (sy==1) ? 2*J + ((y-q->YN[j+m]) >= 0.5*(q->YN[j+1+m]-q->YN[j+m]) ? 1 : 0) : J;

    return !self->covered(l+1,Ic,Jc);
}

// the part of the patches of this rank in force box ID
void nhflow_amr::solid_forces(int ID, double *F) const
{
    nhflow_amr *self = const_cast<nhflow_amr*>(this);
    if(!solids || maxlev<1)
    return;

    for(auto q : self->P)
    {
        nhflow_amr_patch *c = NP(q);
        if(ID>=(int)c->forces.size() || c->forces[ID]==nullptr)
        continue;

        pscope ps(pgc0,c->d,d0);
        c->forces[ID]->patch_sum(c->pp,c->d,pgc0,F);
    }
}

// loads on the floating bodies from the current state of all grids (outside the time step)
void nhflow_amr::body_loads(lexer *p, ghostcell *pgc)
{
    for(int nb=0; nb<b6->objects(); ++nb)
    b6->object(nb)->hydrodynamic_forces_nhflow(p,d0,pgc,d0->WL,false);
}

// --------------------------------------------------------------------- time step
double nhflow_amr::dt_local_max(int m)
{
    double r=0.0;
    for(auto q : P)
    {
        nhflow_amr_patch *c = NP(q);
        lexer *pp = c->pp;
        fdm_nhf *d = c->d;
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EY; jj<EY+c->ny; ++jj)
        {
            if(pp->flagslice4[lij(pp,ii,jj)]<0)
            continue;

            if(m==0)
            r = MAX(r,d->WL(ii,jj));

            if(m>=1 && m<=3)
            for(int kk=0; kk<pp->knoz; ++kk)
            {
                const int n = cidx(pp,ii,jj,kk);
                const double v = m==1 ? d->U[n] : (m==2 ? d->V[n] : d->W[n]);
                r = MAX(r,fabs(v));
            }

            if(m==4)
            for(int kk=0; kk<=pp->knoz; ++kk)
            r = MAX(r,fabs(d->omegaF[fidx(pp,ii,jj,kk)]));

            // the horizontal speed of a cell (curvilinear grids)
            if(m==5)
            for(int kk=0; kk<pp->knoz; ++kk)
            {
                const int n = cidx(pp,ii,jj,kk);
                r = MAX(r,sqrt(d->U[n]*d->U[n] + d->V[n]*d->V[n]));
            }
        }
    }
    return r;
}

void nhflow_amr::dt_cell_size(int wetonly, double &dx, double &dz)
{
    for(auto q : P)
    {
        nhflow_amr_patch *c = NP(q);
        lexer *pp = c->pp;
        fdm_nhf *d = c->d;

        // G 7 1: level l takes 2^l steps per level-0 step, its cells limit the level-0 step as
        // cells 2^l times larger
        const double sc = (sub==1) ? double(1<<c->lev) : 1.0;

        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EY; jj<EY+c->ny; ++jj)
        {
            if(pp->flagslice4[lij(pp,ii,jj)]<0)
            continue;
            if(wetonly==1 && pp->wet[lij(pp,ii,jj)]!=1)
            continue;

            double h;
            // curvilinear grid: the widths of the cell across its faces (as nhflow_timestep)
            if(pp->geo_curv==1)
            h = MIN(nhflow_curv::width1(pp,ii,jj),nhflow_curv::width2(pp,ii,jj));
            else
            if(pp->j_dir==1 && pp->knoy>1)
            h = MIN(pp->DXN[ii+marge],pp->DYN[jj+marge]);
            else
            h = pp->DXN[ii+marge];
            dx = MIN(dx,sc*h);

            for(int kk=0; kk<pp->knoz; ++kk)
            dz = MIN(dz,sc*pp->DZN[kk+marge]*d->WL(ii,jj));
        }
    }
}
