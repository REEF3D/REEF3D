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
#include"bc_zone.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice4.h"
#include"nhflow_momentum_RK2.h"
#include"nhflow_momentum_RK3.h"
#include"nhflow_signal_speed.h"
#include"nhflow_reconstruct_hires.h"
#include"nhflow_reconstruct_weno.h"
#include"nhflow_HLL.h"
#include"nhflow_HLLC.h"
#include"nhflow_diff_void.h"
#include"nhflow_pjm.h"
#include"nhflow_pjm_corr.h"
#include"nhflow_pjm_hs.h"
#include"nhflow_komega_void.h"
#include"vrans_nhflow_v.h"
#include"nhflow_fsf_f.h"
#include"nhflow_forcing.h"
#include"nhflow_amr_6dof.h"
#include"6DOF_nhflow.h"
#include"6DOF_obj.h"
#include"sediment_void.h"
#include"ioflow_void.h"
#include"patchBC_void.h"
#include"reefmg_core.h"
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

    reefamr_param q;
    q.name = "NHFLOW AMR";
    q.maxlev = p->A270;
    q.regrid = 0;
    q.nbuf = MAX(p->A272,0);
    q.tile = MAX(p->A275,4);
    q.tile += q.tile%2;
    q.keep = MIN(MAX(p->A280,0),250);

    // cells computed beyond the patch box: the WENO5 reconstruction reaches three cells; the patch
    // arrays reach EXT + margin fine cells beyond the patch, which the coarser level has to cover
    q.ext = 3;
    q.nest = MAX(3,(q.ext+p->margin+1)/2);

    for(int k=0; k<p->A276; ++k)
    {
        q.rbox.push_back(p->A276_xs[k]); q.rbox.push_back(p->A276_xe[k]);
        q.rbox.push_back(p->A276_ys[k]); q.rbox.push_back(p->A276_ye[k]);
    }
    for(int k=0; k<p->A277; ++k)
    {
        q.fbox.push_back(p->A277_xs[k]); q.fbox.push_back(p->A277_xe[k]);
        q.fbox.push_back(p->A277_ys[k]); q.fbox.push_back(p->A277_ye[k]);
    }

    // no refinement in the relaxation zones of the wave generation (B 98 2) and the numerical
    // beach (B 99 1, 2), measured from the ends of the domain in x, and in the B 520 zones (bc_zone)
    bc_zone_set::from_legacy(p,pgc).norefine_boxes(p,q.fbox);
    q.ioband = 4;

    // refinement around the floating body: margin A 278 around the wetted hull, a rectangle
    // aligned with x and y (moored or oscillating bodies), built at t = 0
    q.zones = (b6!=nullptr && p->A278>0);
    q.zr = p->A278_r;
    q.zalign = true;

    // optional wake wedge (A 279 L a, as SFLOW): the zone is oriented along the direction of
    // motion (the yaw angle at rest) and a wedge of half angle a reaches back from the bow by L
    // at most, where the bow has been; without A 279 the rectangle aligned with x and y
    if(p->A279_L>0.0)
    {
    q.zalign = false;
    q.zL = p->A279_L;
    q.za = p->A279_a;
    }

    // solution-adaptive flags (as sflow_amr): a cell of level l gets level l+1 where the surface
    // jumps by more than A 273 to a neighbour cell or its second difference along x or y exceeds
    // A 282 (both in m, so a finer level flags only steeper or shorter waves)
    tol_eta = MAX(p->A273,0.0);
    tol_curv = MAX(p->A282,0.0);
    adaptive = (tol_eta>0.0 || tol_curv>0.0);

    // the zone follows the body, the patches follow the flags: new patches every A 271 steps
    // (A 271 0: static), the layout is kept as long as it covers the flagged tiles with at most
    // 50 % excess (as fnpf_amr)
    if(q.zones || adaptive)
    {
    q.regrid = MAX(p->A271,0);
    q.lazy = 1.5;
    }

    // patches are fully wet: no patch within 4 level-0 cells of a dry or shallow cell, checked at
    // every regrid (cell_unfit; a static box over dry cells loses those tiles as well)
    q.dryband = 4;

    // the zone holds the hull: its triangles are sized for the finest level (sixdof_obj::amr_hfac,
    // before the 6DOF initialisation builds them), as on the uniform fine grid
    if(q.zones)
    for(int nb=0; nb<b6->objects(); ++nb)
    b6->object(nb)->amr_hfac = 1.0/double(1<<MAX(p->A270,0));

    // vertical refinement (A 281 1): every level doubles the sigma layers of its parent, the
    // coarse nodes are nodes of the fine grid
    vr = (p->A281==1) ? 2 : 1;
    q.vref.assign(q.maxlev+1,vr);

    configure(q);

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

    for(auto v : kv0)
    delete [] v;
    delete mg0;
}

nhflow_amr::pscope::pscope(ghostcell *gg, fdm_nhf *dp, fdm_nhf *d00) : g(gg), dl0(d00)
{
    old = g->set_comms(false);
    g->fdm_nhf_update(dp);
}

nhflow_amr::pscope::~pscope()
{
    g->fdm_nhf_update(dl0);
    g->set_comms(old);
}

// --------------------------------------------------------------------- grids
int nhflow_amr::fidx(lexer *q, int ii, int jj, int kk) const
{
    return (ii-q->imin)*q->jmax*q->kmaxF + (jj-q->jmin)*q->kmaxF + kk - q->kmin;
}

int nhflow_amr::cidx(lexer *q, int ii, int jj, int kk) const
{
    return (ii-q->imin)*q->jmax*q->kmax + (jj-q->jmin)*q->kmax + kk - q->kmin;
}

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

// --------------------------------------------------------------------- setup
void nhflow_amr::ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    int ok=1;
    if(mom0==nullptr) ok=0;
    if(p->A510!=2 && p->A510!=3) ok=0;
    if(p->A511!=1 && p->A511!=2) ok=0;
    if(p->A520<0 || p->A520>2) ok=0;
    if(p->A512!=0 || p->A560!=0 || p->A550!=0) ok=0;
    if(p->B200!=0 || p->S10!=0 || p->X330!=0 || p->A599!=0) ok=0;
    if(p->X10!=0 && (b6==nullptr || p->X60!=1 || p->X16!=0 || p->X320!=0 || p->A516==2 || p->A516==4)) ok=0;
    if(p->A581>0 || p->A583>0 || p->A584>0 || p->A585>0 || p->A586>0 || p->A587>0 || p->A588>0 || p->A589>0 || p->A590>0) ok=0;
    if(p->A580==1) ok=0;
    if(p->E10>0 || p->L10>0 || p->Z20>0) ok=0;
    if(p->j_dir!=1) ok=0;

    if(ok==0)
    {
        if(p->mpirank==0)
        cout<<"NHFLOW AMR (A 270): only for A 510 2/3, A 511 1/2, A 520 0/1/2, A 512 0, A 560 0, A 550 0, B 200 0, S 10 0, "
            <<"X 10 0/1/2 (X 60 1, X 16 0, A 516 0/1/3), no solids, membranes, nets, DEM, particles or rods, and 3D grids "
            <<"-- refinement switched off"<<endl;
        maxlev=0;
        return;
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

    setup(p,pgc);

    mkdir("./REEF3D_NHFLOW_AMR",0777);

    if(p->mpirank==0)
    {
        logout.open("./REEF3D_NHFLOW_AMR/REEF3D_NHFLOW_AMR_log.dat");
        logout<<"# count \t simtime \t dt \t patches \t cells \t pressure iterations \t residual \t water volume \t relative change \t layout"<<endl;
    }

    for(int it=0; it<maxlev; ++it)
    regrid(p,pgc,true);

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
            return sixdof_obj_nhflow::nhflow_grid{glex(g),NP(g)->d,(cur_stage<0) ? &NP(g)->d->WL : stage_out(g,cur_stage).WL};
        };

        cur_stage = -1;
        body_loads(p,pgc);
    }

    m0 = mass(p,d,pgc);

    if(p->mpirank==0)
    {
    cout<<"NHFLOW AMR: "<<maxlev<<" level(s), "<<patches_total<<" patch(es), "<<cells_total<<" refined columns, dt "<<p->dt;
    if(vr==2)
    cout<<", sigma layers doubled on every level (A 281)";
    if(adaptive)
    {
    cout<<", flags:";
    if(tol_eta>0.0)
    cout<<" surface jump "<<tol_eta<<" m (A 273)";
    if(tol_curv>0.0)
    cout<<" second difference "<<tol_curv<<" m (A 282)";
    }
    if(regrid_int>0)
    cout<<", regrid every "<<regrid_int<<" steps (A 271)"<<(par.zones ? ", the zone follows the body" : "");
    cout<<endl;
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

    build_bc2D(pgc,*c);
    build_lexer3D(*c);

    lexer *pp = c->pp;

    // vectors of the Poisson rows (nhflow_poisson fills them over FBASELOOP and LOOP)
    pp->veclength = pp->imax*pp->jmax*pp->kmaxF;

    c->d = new fdm_nhf(pp);
    fdm_nhf *d = c->d;
    pscope ps(pgc,d,d0);

    c->pBC = pBCv;
    c->pflow = pflowv;
    // the floating bodies on the patch grid (all calls do nothing without bodies)
    c->p6dof = new nhflow_amr_6dof(pp,b6,(b6!=nullptr) ? b6->object(0)->nhflow_dsm()*pp->DXM/glex(-1)->DXM : 0.0);
    c->psed = psedv;

    c->pss = new nhflow_signal_speed(pp);

    if(pp->A511==1)
    c->pconv = new nhflow_HLL(pp,pgc,c->pBC);
    if(pp->A511==2)
    c->pconv = new nhflow_HLLC(pp,pgc,c->pBC);

    c->pdiff = new nhflow_diff_void(pp);

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

    c->pturb = new nhflow_komega_func_void(pp,d,pgc);
    c->pvrans = new vrans_nhflow_v(pp,d,pgc);
    c->pfsf = new nhflow_fsf_f(pp,d,pgc,c->pflow,c->pBC);
    c->pdf = new nhflow_forcing(pp,d,pgc);

    if(pp->A510==2)
    c->pmom = new nhflow_momentum_RK2(pp,d,pgc,c->p6dof,c->pvrans,c->pdf,c->psed);
    if(pp->A510==3)
    c->pmom = new nhflow_momentum_RK3(pp,d,pgc,c->p6dof,c->pvrans,c->pdf);

    c->S = {c->pflow,c->pss,c->precon,c->pconv,c->pdiff,c->ppress,nullptr,nullptr,nullptr,c->pfsf,c->pturb,c->pvrans};

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
}

void nhflow_amr::patch_delete(reefamr_patch *q)
{
    nhflow_amr_patch *c = NP(q);

    for(auto v : c->kv)
    delete [] v;
    delete c->mg;

    delete c->pmom;
    delete c->p6dof;
    delete c->pdf;
    delete c->pfsf;
    delete c->pvrans;
    delete c->pturb;
    delete c->ppress;
    delete c->precon;
    delete c->pdiff;
    delete c->pconv;
    delete c->pss;

    // fdm_nhf has no destructor (on level 0 it lives for the whole run): its 3D arrays are
    // freed here, otherwise every regrid that replaces a patch leaks them
    {
        fdm_nhf *d = c->d;
        lexer *pp = c->pp;
        const int n7 = pp->imax*pp->jmax*(pp->kmax+2);
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
}

// cells of level l that need level l+1 (the boxes A 276 and the body zone are marked by the
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

    if(l==0)
    {
        for(int ii=0; ii<NX0; ++ii)
        for(int jj=0; jj<NY0; ++jj)
        if(test(p0,d0,ii,jj))
        tag_cell(lf,ii+O0i,jj+O0j,nbuf,M);
        return;
    }

    for(int n : lev[l])
    {
        nhflow_amr_patch *c = NP(n);
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
        if(test(c->pp,c->d,ii,jj))
        tag_cell(lf,ii-EXT+c->I0,jj-EXT+c->J0,nbuf,M);
    }
}

// a level-0 cell no patch may cover: dry or shallow (water level at most A 545 A 544, the deep
// criterion of the wetting and drying); the patches are fully wet
bool nhflow_amr::cell_unfit(int ii, int jj)
{
    const int n = lij(p0,ii,jj);
    return p0->flagslice4[n]>0 && (p0->wet[n]==0 || d0->WL(ii,jj)<=p0->A545*p0->A544);
}

// --------------------------------------------------------------------- interpolation
// interpolation from coarse cell (ic,jc) to the centre of its child (ox,oy), a quarter coarse
// cell away in x and y (as fnpf_amr): bicubic where all 16 cells are fluid, otherwise
// biquadratic on the 3x3 cells, slopes and curvatures switched off next to solids
void nhflow_amr::pweights(lexer *q, int ic, int jc, int ox, int oy, double *w)
{
    auto fl = [&](int a, int b) { return q->flagslice4[lij(q,a,b)]>0; };

    for(int m=0; m<25; ++m)
    w[m] = 0.0;

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

double nhflow_amr::pq(slice &f, lexer *q, int ic, int jc, int ox, int oy)
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

// layer kk of a cell-layout array (cl true) or node kk of an F-layout array
double nhflow_amr::pq3(const double *f, lexer *q, int ic, int jc, int ox, int oy, int kk, bool cl)
{
    double w[25];
    pweights(q,ic,jc,ox,oy,w);

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
void nhflow_amr::pcol(int g, int ic, int jc, int ox, int oy, const double *src, int knf, double *v)
{
    lexer *q = glex(g);
    const int kc = q->knoz;
    const int fz = knf/kc;
    const int sI = q->jmax*q->kmaxF;
    const int sJ = q->kmaxF;

    double w[25];
    pweights(q,ic,jc,ox,oy,w);

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
void nhflow_amr::fill_stage(ghostcell *pgc, int l, int s)
{
    const int K = klev(l);
    const int Kc = klev(l-1);
    const int nv = 3 + 3*K + K+1;
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
                         v[0] = pq(d->eta,q,f.si,f.sj,f.ox,f.oy);
                         v[1] = q->wet[lij(q,f.si,f.sj)];
                         v[2] = q->deep[lij(q,f.si,f.sj)];
                         // layer by layer on the coarser grid, then (A 281) into the halves
                         for(int k=0; k<Kc; ++k)
                         {
                             uc[k] = pq3(d->U,q,f.si,f.sj,f.ox,f.oy,k,true);
                             uc[Kc+k] = pq3(d->V,q,f.si,f.sj,f.ox,f.oy,k,true);
                             uc[2*Kc+k] = pq3(d->W,q,f.si,f.sj,f.ox,f.oy,k,true);
                         }
                         for(int m=0; m<3; ++m)
                         vcell(q,&uc[m*Kc],Kc,K/Kc,&v[3+m*K]);
                         pcol(f.g,f.si,f.sj,f.ox,f.oy,d->P,K,&v[3+3*K]);
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
                 const double wl = MAX(w[0] + d->depth(ii,jj), pp->A544);
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

// interior state of a fresh patch from its parent: eta, U, V, W per layer and P (bicubic),
// water level and momentum with the depth of the patch
void nhflow_amr::prolong_patch(ghostcell *pgc, nhflow_amr_patch &c)
{
    lexer *pp = c.pp;
    fdm_nhf *d = c.d;
    const int K = pp->knoz;
    const int nby = c.ny/2;
    vector<double> uc, uf(3*K), v(K+1);

    for(int bi=0; bi<c.nx/2; ++bi)
    for(int bj=0; bj<nby; ++bj)
    {
        const int k = bi*nby+bj;
        const int g = c.rgrid[k];
        if(g<-1)
        continue;

        lexer *q = glex(g);
        fdm_nhf *dc = gfd(g);
        const int ic = c.ric[k], jc = c.rjc[k];
        const int Kc = q->knoz;
        const int fz = K/Kc;
        uc.resize(3*Kc);

        for(int a=0; a<2; ++a)
        for(int b=0; b<2; ++b)
        {
            const int ox = a==0?-1:1, oy = b==0?-1:1;
            const int ii = EXT+2*bi+a, jj = EXT+2*bj+b;

            d->eta(ii,jj) = pq(dc->eta,q,ic,jc,ox,oy);
            d->detadt(ii,jj) = pq(dc->detadt,q,ic,jc,ox,oy);
            const double wl = MAX(d->eta(ii,jj) + d->depth(ii,jj), pp->A544);
            d->WL(ii,jj) = wl;

            for(int kc=0; kc<Kc; ++kc)
            {
                uc[kc] = pq3(dc->U,q,ic,jc,ox,oy,kc,true);
                uc[Kc+kc] = pq3(dc->V,q,ic,jc,ox,oy,kc,true);
                uc[2*Kc+kc] = pq3(dc->W,q,ic,jc,ox,oy,kc,true);
            }
            for(int m=0; m<3; ++m)
            vcell(q,&uc[m*Kc],Kc,fz,&uf[m*K]);

            for(int kk=0; kk<K; ++kk)
            {
                const int n = cidx(pp,ii,jj,kk);
                d->U[n] = uf[kk];
                d->V[n] = uf[K+kk];
                d->W[n] = uf[2*K+kk];
                d->UH[n] = wl*d->U[n];
                d->VH[n] = wl*d->V[n];
                d->WH[n] = wl*d->W[n];
            }

            pcol(g,ic,jc,ox,oy,dc->P,K,&v[0]);
            for(int kk=0; kk<=K; ++kk)
            d->P[fidx(pp,ii,jj,kk)] = v[kk];
        }

        // conservative: the 2x2 block keeps the water level and, layer by layer, the momentum
        // of its coarse cell (the restriction is the block mean, with A 281 over the 2x2x2
        // cells of a coarse layer), a constant shift of the interpolated shape; the velocities
        // follow
        const int i0 = EXT+2*bi, j0 = EXT+2*bj;
        double wm = 0.0;
        for(int a=0; a<2; ++a)
        for(int b=0; b<2; ++b)
        wm += 0.25*d->WL(i0+a,j0+b);
        const double dwl = dc->WL(ic,jc) - wm;
        for(int a=0; a<2; ++a)
        for(int b=0; b<2; ++b)
        {
            d->WL(i0+a,j0+b) += dwl;
            d->eta(i0+a,j0+b) += dwl;
        }

        const double wb = 0.25/double(fz);
        for(int kc=0; kc<Kc; ++kc)
        {
            const int nc = cidx(q,ic,jc,kc);
            double um=0.0, vm=0.0, hm=0.0;
            for(int kk=fz*kc; kk<fz*kc+fz; ++kk)
            for(int a=0; a<2; ++a)
            for(int b=0; b<2; ++b)
            {
                const int n = cidx(pp,i0+a,j0+b,kk);
                um += wb*d->UH[n];
                vm += wb*d->VH[n];
                hm += wb*d->WH[n];
            }
            for(int kk=fz*kc; kk<fz*kc+fz; ++kk)
            for(int a=0; a<2; ++a)
            for(int b=0; b<2; ++b)
            {
                const int n = cidx(pp,i0+a,j0+b,kk);
                const double wl = d->WL(i0+a,j0+b);
                d->UH[n] += dc->UH[nc] - um;
                d->VH[n] += dc->VH[nc] - vm;
                d->WH[n] += dc->WH[nc] - hm;
                d->U[n] = d->UH[n]/wl;
                d->V[n] = d->VH[n]/wl;
                d->W[n] = d->WH[n]/wl;
            }
        }
    }
}

// initial water level boxes (F 72) on the interior of a fresh patch at t = 0.  The bicubic
// prolongation of a box edge overshoots (12 % for the hump of 5f), so the blocks whose parent
// stencil (parent cell +-2) or whose own cells see more than one box state are set on the patch
// grid itself: the box height in a box, else the surface of the nearest parent-level cell
// outside the boxes (piecewise constant).  All other blocks keep the prolonged state, and so
// do U, V, W, P and detadt everywhere (initial waves or currents are not touched); UH, VH, WH
// follow the new water level.
void nhflow_amr::ini_boxes(nhflow_amr_patch &c)
{
    lexer *pp = c.pp;
    fdm_nhf *d = c.d;
    const int K = pp->knoz;
    const int nby = c.ny/2;

    // last box that holds the point (the later box wins, as on level 0), -1 none
    auto boxof = [&](double x, double y)
    {
        int r=-1;
        for(int qn=0; qn<p0->F72; ++qn)
        if(x>=p0->F72_xs[qn] && x<p0->F72_xe[qn] && y>=p0->F72_ys[qn] && y<p0->F72_ye[qn])
        r=qn;
        return r;
    };

    for(int bi=0; bi<c.nx/2; ++bi)
    for(int bj=0; bj<nby; ++bj)
    {
        const int k = bi*nby+bj;
        const int g = c.rgrid[k];
        if(g<-1)
        continue;

        lexer *q = glex(g);
        fdm_nhf *dq = gfd(g);
        const int ic = c.ric[k], jc = c.rjc[k];
        auto cbox = [&](int a, int b) { return boxof(q->XP[ic+a+marge],q->YP[jc+b+marge]); };

        const int b0 = cbox(0,0);
        bool aff = false;
        for(int a=-2; a<=2; ++a)
        for(int b=-2; b<=2; ++b)
        if(cbox(a,b)!=b0)
        aff = true;

        for(int a=0; a<2; ++a)
        for(int b=0; b<2; ++b)
        if(boxof(pp->XP[EXT+2*bi+a+marge],pp->YP[EXT+2*bj+b+marge])!=b0)
        aff = true;

        if(!aff)
        continue;

        for(int a=0; a<2; ++a)
        for(int b=0; b<2; ++b)
        {
            const int ii = EXT+2*bi+a, jj = EXT+2*bj+b;
            const int bf = boxof(pp->XP[ii+marge],pp->YP[jj+marge]);
            const int da = a ? 1 : -1, db = b ? 1 : -1;

            double e = dq->eta(ic,jc);
            if(bf>=0)
            e = p0->F72_h[bf] - p0->F60;
            else
            {
                // the parent cell, else its neighbour on the side of this child
                const int ca[4] = {0,da,0,da}, cb[4] = {0,0,db,db};
                for(int n=0; n<4; ++n)
                if(cbox(ca[n],cb[n])<0)
                {
                    e = dq->eta(ic+ca[n],jc+cb[n]);
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
    }
}

// cells of the fresh patch c that an old patch of the same level held (the zone moved with the
// body): fn(old patch, ii, jj, io, jo) with the lexer indices on c and on the old patch
template<class F>
void nhflow_amr::from_old(nhflow_amr_patch &c, vector<reefamr_patch*> &oldP, F fn)
{
    for(auto q : oldP)
    {
        if(q==&c || q->lev!=c.lev)
        continue;

        // kept and removed patches both hold the state of the end of the step
        const int I0 = MAX(c.I0,q->I0), I1 = MIN(c.I1,q->I1);
        const int J0 = MAX(c.J0,q->J0), J1 = MIN(c.J1,q->J1);
        if(I0>I1 || J0>J1)
        continue;

        nhflow_amr_patch *o = NP(q);
        for(int I=I0; I<=I1; ++I)
        for(int J=J0; J<=J1; ++J)
        fn(*o, I-c.I0+EXT, J-c.J0+EXT, I-q->I0+EXT, J-q->J0+EXT);
    }
}

// --------------------------------------------------------------------- restriction
// water level and surface of the covered coarse cells after the continuity part of stage s
void nhflow_amr::restrict_surface(int s)
{
    for(int l=maxlev; l>=1; --l)
    for(int id : lev[l])
    {
        nhflow_amr_patch *c = NP(id);
        stg F = stage_out(id,s);
        fdm_nhf *df = c->d;
        const int nby = c->ny/2;

        for(int bi=0; bi<c->nx/2; ++bi)
        for(int bj=0; bj<nby; ++bj)
        {
            const int k = bi*nby+bj;
            const int g = c->rgrid[k];
            if(g<-1)
            continue;

            stg C = stage_out(g,s);
            fdm_nhf *dc = gfd(g);
            const int ic = c->ric[k], jc = c->rjc[k];
            const int i0 = EXT+2*bi, j0 = EXT+2*bj;
            slice &WLf = *F.WL;

            const double wl = 0.25*(WLf(i0,j0)+WLf(i0+1,j0)+WLf(i0,j0+1)+WLf(i0+1,j0+1));
            (*C.WL)(ic,jc) = wl;
            dc->eta(ic,jc) = wl - dc->depth(ic,jc);
            dc->detadt(ic,jc) = 0.25*(df->detadt(i0,j0)+df->detadt(i0+1,j0)+df->detadt(i0,j0+1)+df->detadt(i0+1,j0+1));
        }
    }
}

// UH, VH, WH (and with P the pressure) of the covered coarse cells, U = UH/WL
void nhflow_amr::restrict_momentum(int s, bool withP)
{
    for(int l=maxlev; l>=1; --l)
    for(int id : lev[l])
    {
        nhflow_amr_patch *c = NP(id);
        lexer *pp = c->pp;
        stg F = stage_out(id,s);
        const int nby = c->ny/2;

        for(int bi=0; bi<c->nx/2; ++bi)
        for(int bj=0; bj<nby; ++bj)
        {
            const int k = bi*nby+bj;
            const int g = c->rgrid[k];
            if(g<-1)
            continue;

            lexer *q = glex(g);
            stg C = stage_out(g,s);
            fdm_nhf *dc = gfd(g);
            const int ic = c->ric[k], jc = c->rjc[k];
            const int i0 = EXT+2*bi, j0 = EXT+2*bj;
            const double wl = (*C.WL)(ic,jc);
            const double wlvl = fabs(wl)>p0->A544 ? wl : 1.0e20;
            const int K = q->knoz;
            const int fz = pp->knoz/K;

            for(int kk=0; kk<K; ++kk)
            {
                // the 2x2 children of the layer, with A 281 in both fine layers (equal thickness)
                auto avg = [&](const double *f)
                {
                    double r = 0.0;
                    for(int kf=fz*kk; kf<fz*kk+fz; ++kf)
                    r += f[cidx(pp,i0,j0,kf)]+f[cidx(pp,i0+1,j0,kf)]+f[cidx(pp,i0,j0+1,kf)]+f[cidx(pp,i0+1,j0+1,kf)];
                    return 0.25*r/double(fz);
                };
                const int n = cidx(q,ic,jc,kk);
                C.UH[n] = avg(F.UH);
                C.VH[n] = avg(F.VH);
                C.WH[n] = avg(F.WH);
                dc->U[n] = C.UH[n]/wlvl;
                dc->V[n] = C.VH[n]/wlvl;
                dc->W[n] = C.WH[n]/wlvl;
            }
        }
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
        // bed
        for(int id : lev[l])
        if(P[id]->fresh)
        {
            nhflow_amr_patch *c = NP(id);
            const int nby = c->ny/2;
            for(int bi=0; bi<c->nx/2; ++bi)
            for(int bj=0; bj<nby; ++bj)
            {
                const int k = bi*nby+bj;
                const int g = c->rgrid[k];
                if(g<-1)
                continue;
                for(int a=0; a<2; ++a)
                for(int b=0; b<2; ++b)
                c->d->bed(EXT+2*bi+a,EXT+2*bj+b) = plin(gfd(g)->bed,glex(g),c->ric[k],c->rjc[k],a==0?-1:1,b==0?-1:1);
            }
        }
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

            prolong_patch(pgc,*c);

            // where an old patch of the same level was, its state is taken over
            const int K = pp->knoz;
            from_old(*c,oldP,[&](nhflow_amr_patch &o, int ii, int jj, int io, int jo)
            {
                fdm_nhf *od = o.d;
                lexer *op = o.pp;
                d->eta(ii,jj) = od->eta(io,jo);
                d->WL(ii,jj) = od->WL(io,jo);
                d->detadt(ii,jj) = od->detadt(io,jo);
                for(int kk=0; kk<K; ++kk)
                {
                    const int n = cidx(pp,ii,jj,kk), m = cidx(op,io,jo,kk);
                    d->U[n] = od->U[m];
                    d->V[n] = od->V[m];
                    d->W[n] = od->W[m];
                    d->UH[n] = od->UH[m];
                    d->VH[n] = od->VH[m];
                    d->WH[n] = od->WH[m];
                }
                for(int kk=0; kk<=K; ++kk)
                d->P[fidx(pp,ii,jj,kk)] = od->P[fidx(op,io,jo,kk)];
            });

            // at t = 0 the initial water level boxes (F 72) on the patch grid itself, as
            // nhflow_fsf_ini on level 0 (the interpolated coarse box would overshoot at its edges);
            // level 0 takes their mean in regrid_finish
            if(p0->count==0 && p0->F72>0)
            ini_boxes(*c);
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
        }
    }
}

// level 0 consistent with the patches
void nhflow_amr::regrid_finish(ghostcell *pgc, int old_total)
{
    for(int l=maxlev; l>=1; --l)
    for(int id : lev[l])
    {
        nhflow_amr_patch *c = NP(id);
        fdm_nhf *df = c->d;
        const int nby = c->ny/2;
        for(int bi=0; bi<c->nx/2; ++bi)
        for(int bj=0; bj<nby; ++bj)
        {
            const int k = bi*nby+bj;
            const int g = c->rgrid[k];
            if(g<-1)
            continue;
            fdm_nhf *dc = gfd(g);
            const int ic = c->ric[k], jc = c->rjc[k];
            const int i0 = EXT+2*bi, j0 = EXT+2*bj;
            const double wl = 0.25*(df->WL(i0,j0)+df->WL(i0+1,j0)+df->WL(i0,j0+1)+df->WL(i0+1,j0+1));
            dc->WL(ic,jc) = wl;
            dc->eta(ic,jc) = wl - dc->depth(ic,jc);
        }
    }

    // momentum, velocities and pressure: the end-of-step arrays are the last stage's output
    restrict_momentum(mom0->stages()-1,true);

    if(patches_total>0 || old_total>0)
    halo0(pgc,mom0->stages()-1);
}

// --------------------------------------------------------------------- flux matching
// face fluxes of grid id (0: level 0, n: patch n-1), ipol 1-3 momentum, 4 continuity
void nhflow_amr::flux_hook(lexer *p, fdm_nhf *d, int id, int ipol, double *Fx, double *Fy)
{
    const int K = p->knoz;

    // a patch records its box faces
    if(id>0)
    {
        nhflow_amr_patch &c = *NP(id-1);
        const int il = EXT-1, ih = EXT+c.nx-1;
        const int jl = EXT-1, jh = EXT+c.ny-1;

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
            c.rec[ipol][0][r*K+k] = Fx[cidx(p,il,EXT+r,k)];
            c.rec[ipol][1][r*K+k] = Fx[cidx(p,ih,EXT+r,k)];
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
                c.rec[0][0][r] = d->dfx(il,EXT+r);
                c.rec[0][1][r] = d->dfx(ih,EXT+r);
            }
            for(int r=0; r<c.nx; ++r)
            {
                c.rec[0][2][r] = d->dfy(EXT+r,jl);
                c.rec[0][3][r] = d->dfy(EXT+r,jh);
            }
        }
    }

    if(id>=(int)match.size())
    return;

    // coarse faces next to the patches of this rank: mean of the two fine faces, layer by layer
    // (A 281: of the two fine faces in both fine layers of the coarse layer; the fluxes are per
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
            val = 0.5*(R[m.r*Kf+k]+R[(m.r+1)*Kf+k]);
            else
            val = 0.25*(R[m.r*Kf+2*k]+R[m.r*Kf+2*k+1]+R[(m.r+1)*Kf+2*k]+R[(m.r+1)*Kf+2*k+1]);
            if(m.dir==0)
            Fx[cidx(p,m.fi,m.fj,k)] = val;
            else
            Fy[cidx(p,m.fi,m.fj,k)] = val;
        }

        if(ipol==4)
        {
            const double df = 0.5*(c.rec[0][m.side][m.r]+c.rec[0][m.side][m.r+1]);
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
            const double val = v[(ipol-1)*K+k];
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
                         sb[(ip-1)*K+k] = 0.5*(R[r*K+k]+R[(r+1)*K+k]);
                         else
                         sb[(ip-1)*K+k] = 0.25*(R[r*Kf+2*k]+R[r*Kf+2*k+1]+R[(r+1)*Kf+2*k]+R[(r+1)*Kf+2*k+1]);
                     }
                 }
                 sb[4*K] = 0.5*(c->rec[0][side][r]+c->rec[0][side][r+1]);
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
        pscope ps(pgc,c->d,d0);
        c->pmom->step_begin(pp,c->d,pgc,c->S);
    }
    }

    mom->step_begin(p,d,pgc,S0);
    S0p = &S0;

    for(int s=0; s<ns; ++s)
    {
        cur_stage = s;
        double t0 = MPI_Wtime();

        // stage input around the patches, coarse to fine
        for(int l=1; l<=maxlev; ++l)
        fill_stage(pgc,l,s);
        tm[0] += MPI_Wtime()-t0;

        // continuity and momentum fluxes, finest first: the coarser grids take the recorded fine fluxes
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

        mom->phase_F(p,d,pgc,S0,s);
        mom->phase_M(p,d,pgc,S0,s);

        t0 = MPI_Wtime();
        {
        comms_off guard(pgc);
        restrict_surface(s);
        restrict_momentum(s,false);
        }
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
}

// the zone around the body moves, the flags move: new patches every A 271 steps, from the state
// at the end of the step
void nhflow_amr::regrid_step(lexer *p, ghostcell *pgc)
{
    if(maxlev>0 && regrid_int>0 && p->count%regrid_int==0)
    {
        const double t0 = MPI_Wtime();
        regrid(p,pgc,false);
        tm[7] += MPI_Wtime()-t0;
    }
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
    int g=-1, l=0;
    for(int n=0; n<(int)P.size(); ++n)
    {
        reefamr_patch *c = P[n];
        lexer *pp = c->pp;
        if(c->lev<=l)
        continue;
        if(x>=pp->XN[EXT+marge] && x<pp->XN[EXT+c->nx+marge] && y>=pp->YN[EXT+marge] && y<pp->YN[EXT+c->ny+marge])
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
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
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
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
        {
            if(pp->flagslice4[lij(pp,ii,jj)]<0)
            continue;
            if(wetonly==1 && pp->wet[lij(pp,ii,jj)]!=1)
            continue;

            double h;
            if(pp->j_dir==1 && pp->knoy>1)
            h = MIN(pp->DXN[ii+marge],pp->DYN[jj+marge]);
            else
            h = pp->DXN[ii+marge];
            dx = MIN(dx,h);

            for(int kk=0; kk<pp->knoz; ++kk)
            dz = MIN(dz,pp->DZN[kk+marge]*d->WL(ii,jj));
        }
    }
}
