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

#include"nhflow_momentum_SI.h"
#include"seastate_nhflow.h"
#include"si_eta2D.h"
#include"bc_zone.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"reefckpt.h"
#include"nhflow_reconstruct.h"
#include"nhflow_curv.h"
#include"sflow_curv.h"
#include"nhflow_diffusion.h"
#include"nhflow_fsf.h"
#include"nhflow_turbulence.h"
#include"vrans_nhflow.h"
#include"solver.h"
#include"solver2D.h"
#include"reefmg2D.h"
#ifdef REEF3D_USE_HYPRE
#include"hypre_struct2D.h"
#endif
#include"wind_f.h"
#include"wind_v.h"
#include"met_forcing.h"
#include"nhflow_tracer_SI.h"
#include<cstring>
#include<sys/stat.h>
#include<sys/types.h>

nhflow_momentum_SI::nhflow_momentum_SI(lexer *p, fdm_nhf *d, ghostcell *pgc, sixdof *pp6dof, vrans_nhflow *ppvrans,
                                       nhflow_forcing *ppnhfdf, bool block)
                                    : nhflow_momentum_func(p,d,pgc),
                                      qx(p),qxn(p),qxE(p),Hx(p),alx(p),gmx(p),ax(p),kx(p),gx0(p),cx0(p),
                                      qy(p),qyn(p),qyE(p),Hy(p),aly(p),gmy(p),ay(p),ky(p),gy0(p),cy0(p),
                                      eta0(p),WL0(p),etaI(p),gmc(p),rhsW(p),bdiag(p),tS(p),tE(p),
                                      WLk(p),etaB(p),WLb(p),wsrc(p)
{
    check(p);

    const int N = p->imax*p->jmax*(p->kmax+2);

    double **arr[] = {&Qx,&Qy,&Qxn,&Qyn,&Tx,&Ty,&QAx,&QAy,&Rx,&Ry,&UH0,&VH0,&WH0,&UHa,&VHa,&WHa,&UHb,&VHb,&WHb,&UHp,&VHp,&WHp,
                      &UHD,&VHD,&WHD,&dUH,&dVH};

    for(double **a : arr)
    p->Darray(*a,N);

    sigma_ini(p,d,pgc,d->eta);

    p6dof = pp6dof;
    pnhfdf = ppnhfdf;
    pmfrc = ppnhfdf;
    pm6dof = pp6dof;
    pvrans = ppvrans;

    if(p->A570==0)
    pwind = new wind_v(p);

    if(p->A570>0)
    pwind = new wind_f(p);

    // atmospheric forcing (W 200) and Coriolis force (W 211)
    if(met_forcing::active(p))
    {
    pmet = new met_forcing(p,pgc);
    d->met = pmet;
    }

    psi = new si_eta2D(p);

    // tracers (A 600)
    if(p->A600>0)
    ptrc = new nhflow_tracer_SI(p,d,this);

    // the 2D solver of the eta system: created with the first solve of one grid (new_solver2D); the
    // multi-block grid solves the eta system of all blocks together (nhflow_amr_si.cpp)
    psolv2D = nullptr;
    psi_iter = 0;
    psi_time = 0.0;
    patch = block;

    theta = p->A535;
    g = fabs(p->W22);
    rk = (p->A510==2)?2:3;
    qini = false;
    wd_iter = 0;
    wd_cells = 0;
    Qin = Qout = vrelax = 0.0;
    gcval_eta = (p->F50>=1 && p->F50<=4) ? 50+p->F50 : 50;

    // vertical diffusion in the column operator: as nhflow_idiff (A 512 2, or A 560 > 0)
    visc = !(p->A512==0 && p->A560==0) || (p->A547_c>0.0 && p->A560==0) || (p->A547_nu0>0.0 && p->A560==0);

    wcu.assign(p->knoz,0.0);
    wcd.assign(p->knoz,0.0);
    wcl.assign(p->knoz,0.0);
    wnu.assign(p->knoz,0.0);
    wa.assign(p->knoz,0.0);
    wb.assign(p->knoz,0.0);
    wr.assign(p->knoz,0.0);
    wx.assign(p->knoz,0.0);

    bc_method(p);

    if(p->mpirank==0 && !patch)
    {
    cout<<"NHFLOW hydrostatic semi-implicit free surface (A 534 1): theta "<<theta<<", SSP-RK"<<rk<<" predictor";

    if(p->A537==1)
    cout<<", wetting and drying iteration";

    if(p->A538>0.0)
    cout<<", mass flux dissipation "<<p->A538;

    cout<<", thin-water velocities below "<<p->A546<<" m (A 546)";

    cout<<endl;

    mkdir("./REEF3D_NHFLOW_Log",0777);
    if(reefckpt_open(silog,"./REEF3D_NHFLOW_Log/REEF3D_NHFLOW_SI.dat",p->I45))   // I 45: continued
    {
        silog<<"# semi-implicit NHFLOW (A 534 1): volume and open boundary discharges"<<endl;
        silog<<"# count  simtime  dt  volume  Qin  Qout  relax_volume  eta_iterations  wetdry_iterations  dry_cells"<<endl;
    }
    }
}

nhflow_momentum_SI::~nhflow_momentum_SI()
{
    double *arr[] = {Qx,Qy,Qxn,Qyn,Tx,Ty,QAx,QAy,Rx,Ry,UH0,VH0,WH0,UHa,VHa,WHa,UHb,VHb,WHb,UHp,VHp,WHp,UHD,VHD,WHD,dUH,dVH};

    for(double *a : arr)
    delete [] a;

    delete psi;
    delete ptrc;

    if(psolv2D!=nullptr)
    delete psolv2D;
    delete pwind;
    delete pmet;
}

void nhflow_momentum_SI::check(lexer *p)
{
    // the options the semi-implicit path does not have yet
    int stop=0;

    auto msg = [&](const char *s)
    {
        if(p->mpirank==0)
        cout<<"!!! NHFLOW A 534 1 (semi-implicit): "<<s<<" not available yet !!!"<<endl;
        stop=1;
    };

    if(p->A520!=0)
    msg("the non-hydrostatic pressure (A 520 > 0) is; set A 520 0,");

    if(p->A510!=2 && p->A510!=3)
    msg("a predictor other than SSP-RK2/3 (A 510 2/3) is");

    if(p->G1>0)
    msg("mesh refinement (G 1) is");



    if(p->A519==2)
    msg("wall friction at every wall face (A 519 2) is");

    if(p->B200>0)
    msg("VRANS (B 200) is");

    if(p->X10>0 || p->X13>0)
    msg("6DOF, floating bodies and direct forcing (X 10, X 13) are");

    if(p->X330>0)
    msg("membranes (X 330) are");

    if(p->Z30>0)
    msg("FEM coupling (Z 30) is");

    // SEASTATE coupling (A 750 1): the wave force in the predictor (the same force in every layer, as the explicit
    // NHFLOW), the Stokes mass source (A 751 2) in the eta system; surfbeat (A 770 1) needs the host's step for the
    // wave groups and the long-wave boundary
    if(p->A750>0 && p->A770==1)
    msg("SEASTATE surfbeat (A 770 1) is");

    if(p->B411>0 || p->B412>0 || p->B413>0 || p->B414>0 || p->B415>0 || p->B416>0 || p->B417>0 || p->B418>0 || p->B421>0 || p->B422>0)
    msg("patch boundary conditions (B 411-422) are");

    if(stop==1)
    exit(1);

    if(p->A535<0.5 && p->mpirank==0)
    cout<<"NHFLOW A 535 "<<p->A535<<": theta < 0.5 is unstable for large time steps"<<endl;

    if(p->A512==1 && p->mpirank==0)
    cout<<"NHFLOW A 534 1: explicit diffusion (A 512 1) limits the time step, use A 512 2 (implicit)"<<endl;
}

void nhflow_momentum_SI::bc_method(lexer *p)
{
    // iowave open edge method per side (B 520 id method prio, B 521 id edge ...): 3 Riemann,
    // 4 Flather, 5 clamped level, 6 clamped discharge; B 521 edge 1 x-, 2 x+, 3 y-, 4 y+
    for(int q=0;q<5;++q)
    method[q]=0;

    const int edge_cs[5] = {0,1,4,3,2};

    for(int q=0;q<p->B521;++q)
    for(int qn=0;qn<p->B520;++qn)
    if(p->B520_id[qn]==p->B521_id[q] && p->B520_method[qn]>=3 && p->B520_method[qn]<=6
       && bc_zone_set::zone_block(p,p->B520_id[qn])==p->amr_block)
    {
        const int e = p->B521_edge[q];

        if(e<1 || e>4)
        continue;

        const int cs = edge_cs[e];

        if(method[cs]!=0 && method[cs]!=p->B520_method[qn])
        {
            if(p->mpirank==0)
            cout<<"!!! NHFLOW A 534 1 (semi-implicit): different B 520 methods on one edge (B 521 edge "<<e<<") are not available yet !!!"<<endl;
            exit(1);
        }

        method[cs] = p->B520_method[qn];
    }
}

// ---------------------------------------------------------------------------
// the stage interface (nhflow_amr) is refused in check()
// ---------------------------------------------------------------------------
void nhflow_momentum_SI::step_begin(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&) {}
void nhflow_momentum_SI::phase_F(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) {}
void nhflow_momentum_SI::phase_M(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) {}
void nhflow_momentum_SI::phase_D(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int,int) {}
void nhflow_momentum_SI::phase_P1(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) {}
void nhflow_momentum_SI::phase_P2(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) {}
void nhflow_momentum_SI::phase_E(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) {}
slice& nhflow_momentum_SI::stage_WL(fdm_nhf *d, int) { return d->WL; }
double* nhflow_momentum_SI::stage_UH(fdm_nhf *d, int, int m) { return m==0?d->UH:(m==1?d->VH:d->WH); }

// ---------------------------------------------------------------------------
// one time step
// ---------------------------------------------------------------------------
void nhflow_momentum_SI::start(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *ppflow, nhflow_signal_speed *pss,
                               nhflow_reconstruct *precon, nhflow_convection *pconvec, nhflow_diffusion *pnhfdiff,
                               nhflow_pressure *ppress, solver *ppoissonsolv, solver *ppsolv, nhflow *pnhf, nhflow_fsf *pfsf,
                               nhflow_turbulence *pnhfturb, vrans_nhflow *ppvrans)
{
    nhflow_stage_obj S = {ppflow,pss,precon,pconvec,pnhfdiff,ppress,ppoissonsolv,ppsolv,pnhf,pfsf,pnhfturb,ppvrans};

    // multi-block grid (G 50 1): the step of all blocks with one composite eta system (nhflow_amr_si.cpp)
    if(prun!=nullptr)
    {
    prun->step(p,d,pgc,this,S);
    return;
    }

    // one grid: the parts of the step in sequence (with the multi-block grid the stage runner fills the
    // cells around the blocks between them)
    si_restore(p,d,pgc);
    si_begin(p,d,pgc,S);

    for(int s=0; s<rk; ++s)
    {
    si_stage_vel(p,d,pgc,S,s);
    si_stage_rhs(p,d,pgc,S,s);
    }

    si_predictor_end(p,d,pgc);
    si_assemble(p,d,pgc);

    // initial guess eta^n; Newton iteration on WL = max(eta + d, 0) (Casulli 2009)
    if(psolv2D==nullptr)
    psolv2D = new_solver2D(p,pgc);

    wd_iter = psi->solve_wetdry(p,pgc,psolv2D,ax,ay,rhsW,d->depth,bdiag,etaI,p->A537==1 ? 20 : 1,wd_cells);
    psi_iter = psi->iter;
    psi_time = psi->time;

    si_end(p,d,pgc,S);

    // tracers (A 600) with the transport of the step
    if(ptrc!=nullptr)
    tracer_step(p,d,pgc);
}

// tracers of one grid (the multi-block grid: nhflow_amr::si_tracer): the substeps with the partition
// halo between them, the log and the gauges
void nhflow_momentum_SI::tracer_step(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const int nt = ptrc->nt;

    double r = ptrc->prepare(p,d);
    r = pgc->globalmax(r);
    const int ns = nhflow_tracer_SI::substeps(p,r);

    // the fluid mask of the halo (static, once)
    if(!trc_mask)
    {
    pgc->gcparaxV(p,ptrc->mask,14);
    pgc->gcparacoxV(p,ptrc->mask,14);
    trc_mask = true;
    }

    for(int s=0; s<ns; ++s)
    {
        for(int n=0; n<nt; ++n)
        {
        pgc->gcparaxV(p,nhflow_tracer_SI::field(d,n),14);
        pgc->gcparacoxV(p,nhflow_tracer_SI::field(d,n),14);
        }

        ptrc->flushed(p,d,s,ns);

        for(int n=0; n<nt; ++n)
        {
        pgc->gcparaxV(p,ptrc->Cx[n],14);
        pgc->gcparacoxV(p,ptrc->Cx[n],14);
        }

        ptrc->substep(p,d,s,ns);
    }

    for(int n=0; n<nt; ++n)
    {
    pgc->gcparaxV(p,nhflow_tracer_SI::field(d,n),14);
    pgc->gcparacoxV(p,nhflow_tracer_SI::field(d,n),14);
    }

    ptrc->ghosts(p,d);

    nhflow_tracer_SI::stats st;
    nhflow_tracer_SI::stats_zero(st);
    ptrc->sums(p,d,[](int,int){ return true; },st);
    nhflow_tracer_SI::reduce(p,st);
    ptrc->write_log(p,st,ns,r);

    if(p->P51>0)
    {
        static vector<double> gdist;
        static vector<int> gcell;
        static bool found = false;

        if(!found)
        {
        ptrc->gauge_find(p,[](int,int){ return true; },gdist,gcell);
        found = true;
        }

        vector<double> val, dist = gdist;
        ptrc->gauge_values(p,d,gcell,val);
        nhflow_tracer_SI::gauge_combine(p,dist,val);
        nhflow_tracer_SI::gauge_write(p,val);
    }
}

void nhflow_momentum_SI::tracer_ckpt_fields(const std::function<void(const char*, char, double*)> &fd)
{
    if(ptrc!=nullptr)
    ptrc->ckpt_fields(fd);
}

solver2D* nhflow_momentum_SI::new_solver2D(lexer *p, ghostcell *pgc)
{
    // REEFMG on the slice (N 10 1), otherwise hypre
#ifdef REEF3D_USE_HYPRE
    if(p->N10!=1)
    return new hypre_struct2D(p,pgc);
#endif
    return new reefmg2D(p,pgc);
}

void nhflow_momentum_SI::si_restore(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // the water column of the end of the previous step: the main loop passes WL and eta through
    // ghost cell routines (depth_update) which may reach the fluid behind a thin solid; the
    // semi-implicit step keeps its volume
    if(qini)
    {
        SLICELOOP4
        {
        d->WL(i,j) = WLk(i,j);
        d->eta(i,j) = WLk(i,j) - d->depth(i,j);
        }

    exchange4(p,pgc,d->WL);
    exchange4(p,pgc,d->eta);
    }

}

void nhflow_momentum_SI::si_begin(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S)
{
    pflow = S.pflow;
    pturb = S.pturb;
    psolv = S.psolv;
    nhflow_fsf *pfsf = S.pfsf;

    const int N = p->imax*p->jmax*(p->kmax+2);

    // ghost cells of the open boundaries (ioflow discharge, iowave open edges), as step_begin of
    // the explicit integrators
    pflow->discharge_nhflow(p,d,pgc);
    pflow->inflow_nhflow(p,d,pgc,d->U,d->V,d->W,d->UH,d->VH,d->WH,d->WL);

    // state at t^n (all cells incl. ghost cells)
    memcpy(UH0,d->UH,N*sizeof(double));
    memcpy(VH0,d->VH,N*sizeof(double));
    memcpy(WH0,d->WH,N*sizeof(double));

    if(p->j_dir==0)
    for(int q=0;q<N;++q)
    VH0[q] = 0.0;

    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    {
    eta0(i,j) = d->eta(i,j);
    WL0(i,j) = d->WL(i,j);
    }

    SLICELOOP4
    d->eta_n(i,j) = d->eta(i,j);

    exchange4(p,pgc,eta0);
    exchange4(p,pgc,WL0);

    // sigma grid of t^n (metrics of the predictor) and the kinematic boundary values
    sigma_update(p,d,pgc,d->WL);

    // A 547: the algebraic eddy viscosity of t^n (without a turbulence model)
    if(p->A560==0 && (p->A547_c>0.0 || p->A547_nu0>0.0))
    eddy_algebraic(p,d);
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);

    // 1. face depths
    face_depth(p,d,pgc);

    // first step: layer face discharges from the cell values, the layer transport of the
    // "previous step" with them
    if(!qini)
    {
        for(int q=0;q<N;++q)
        Qx[q] = Qy[q] = 0.0;

        double nx,ny;

        SLICELOOP1
        {
        qx(i,j) = 0.0;
        fnorm(p,1,i,j,nx,ny);

            if(Hx(i,j)>0.0)
            for(k=0;k<p->knoz;++k)
            {
            Qx[IJK] = 0.5*(d->UH[IJK] + d->UH[Ip1JK])*nx + 0.5*(d->VH[IJK] + d->VH[Ip1JK])*ny*p->y_dir;
            qx(i,j) += p->DZN[KP]*Qx[IJK];
            }
        }

        if(p->j_dir==1)
        SLICELOOP2
        {
        qy(i,j) = 0.0;
        fnorm(p,2,i,j,nx,ny);

            if(Hy(i,j)>0.0)
            for(k=0;k<p->knoz;++k)
            {
            Qy[IJK] = 0.5*(d->UH[IJK] + d->UH[IJp1K])*nx + 0.5*(d->VH[IJK] + d->VH[IJp1K])*ny;
            qy(i,j) += p->DZN[KP]*Qy[IJK];
            }
        }

    exchange1(p,pgc,qx);
    exchange2(p,pgc,qy);
    exchangeF(p,pgc,Qx);
    exchangeF(p,pgc,Qy);

    // face depths upwind of the new discharges
    face_depth(p,d,pgc);
    }

    // implicit bed friction (A 519 1)
    friction(p,d);

    // 2. open boundaries (the first step also sets their face discharges)
    bc_setup(p,d,pgc);

    if(!qini)
    {
    memcpy(Tx,Qx,N*sizeof(double));
    memcpy(Ty,Qy,N*sizeof(double));
    }

    qini = true;

    // Q^n, q^n
    memcpy(Qxn,Qx,N*sizeof(double));
    memcpy(Qyn,Qy,N*sizeof(double));

    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    {
    qxn(i,j) = qx(i,j);
    qyn(i,j) = qy(i,j);
    }

    // the inflow ghost cells of the stage arrays from those of t^n
    S.pflow->rkinflow_nhflow(p,d,pgc,d->U,d->V,d->W,UHa,VHa,WHa,WL0);
    S.pflow->rkinflow_nhflow(p,d,pgc,d->U,d->V,d->W,UHb,VHb,WHb,WL0);
}

void nhflow_momentum_SI::si_end(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S)
{
    // 7. layer face discharges
    face_update(p,d,pgc);

    // 8. WL, eta, wetting and drying
    eta_update(p,d,pgc,S.pfsf);

    // 9.-10. cell momentum, relaxation zones, velocities, omega
    cell_update(p,d,pgc);
    relax(p,d,pgc);
    finish(p,d,pgc);

    SLICELOOP4
    WLk(i,j) = d->WL(i,j);

    p->poissoniter = psi_iter;
    p->poissontime = psi_time;

    // the multi-block grid writes the log of all blocks (nhflow_amr::si_step)
    if(patch || prun!=nullptr)
    return;

    log(p,d,pgc);

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
    cout<<"SI eta: iter "<<psi_iter<<"  time "<<setprecision(3)<<psi_time;

    if(p->A537==1)
    cout<<"  wetdry iterations "<<wd_iter<<"  dry cells "<<wd_cells;

    cout<<endl;
    }
}

void nhflow_momentum_SI::face_depth(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // total water depth of the face: still water depth of the face plus eta upwind of the column
    // discharge (Stelling and Duinmeijer 2003), max(eta) where it is zero; a face is open between
    // two fluid cells when this depth exceeds the wetting criterion A 544. Next to a dry cell the
    // still water depth of the face is the smaller one (the higher bed) and eta that of the wet cell:
    // the face opens when the water of the wet cell stands above the bed of the dry cell (flooding);
    // the film of a dry cell does not flow (its eta is its bed level, a lake at rest stays at rest).
    // Between two dry cells the face is closed.
    // A 548: between two wet cells the still water depth of the face is the mean of the two (0, default), the
    // smaller one (the higher bed, Stelling and Duinmeijer 2003) everywhere (1), or only where one of the two
    // cells holds a film thinner than 10 A 544 (2). The mean gives a face over a bed step the depth of the deep
    // side, so that a film on the high side drains and floods by turns (flip-flop of shoreline cells).
    double df,eu,h;

    SLICELOOP1
    {
    Hx(i,j) = 0.0;

        if(p->flagslice4[IJ]>0 && p->flagslice4[Ip1J]>0)
        {
        const int wL = p->wet[IJ], wR = p->wet[Ip1J];

        if(wL==1 && wR==1)
        {
        df = (p->A548==1 || (p->A548==2 && MIN(d->WL(i,j),d->WL(i+1,j))<10.0*p->A544)) ? MIN(d->depth(i,j),d->depth(i+1,j)) : 0.5*(d->depth(i,j) + d->depth(i+1,j));

        if(qini && qx(i,j)>1.0e-12)
        eu = d->eta(i,j);

        else
        if(qini && qx(i,j)<-1.0e-12)
        eu = d->eta(i+1,j);

        else
        eu = MAX(d->eta(i,j),d->eta(i+1,j));

        h = eu + df;
        }

        else
        if(wL==1 || wR==1)
        h = (wL==1 ? d->eta(i,j) : d->eta(i+1,j)) + MIN(d->depth(i,j),d->depth(i+1,j));

        else
        h = 0.0;

        Hx(i,j) = h>p->A544 ? h : 0.0;
        }
    }

    if(p->j_dir==1)
    SLICELOOP2
    {
    Hy(i,j) = 0.0;

        if(p->flagslice4[IJ]>0 && p->flagslice4[IJp1]>0)
        {
        const int wL = p->wet[IJ], wR = p->wet[IJp1];

        if(wL==1 && wR==1)
        {
        df = (p->A548==1 || (p->A548==2 && MIN(d->WL(i,j),d->WL(i,j+1))<10.0*p->A544)) ? MIN(d->depth(i,j),d->depth(i,j+1)) : 0.5*(d->depth(i,j) + d->depth(i,j+1));

        if(qini && qy(i,j)>1.0e-12)
        eu = d->eta(i,j);

        else
        if(qini && qy(i,j)<-1.0e-12)
        eu = d->eta(i,j+1);

        else
        eu = MAX(d->eta(i,j),d->eta(i,j+1));

        h = eu + df;
        }

        else
        if(wL==1 || wR==1)
        h = (wL==1 ? d->eta(i,j) : d->eta(i,j+1)) + MIN(d->depth(i,j),d->depth(i,j+1));

        else
        h = 0.0;

        Hy(i,j) = h>p->A544 ? h : 0.0;
        }
    }

    exchange1(p,pgc,Hx);
    exchange2(p,pgc,Hy);
}

double nhflow_momentum_SI::bed_gamma(lexer *p, double ut, double dn, double ks)
{
    // implicit form of the log-law bed friction of the bottom layer (nhflow_bcmom, A 519 1):
    // d(UH)/dt = -gamma UH,  gamma = |u| / (u+^2 dn),  u+ = 1/kappa max(1, ln(30 z0/ks)),
    // z0 = max(dn/2, ks/30)
    const double kappa=0.4;

    if(ks<=0.0)
    ks=0.0001;

    double z0 = 0.5*dn;

    if(30.0*z0<ks)
    z0=ks/30.0;

    const double uplus = (1.0/kappa)*MAX(1.0,std::log(30.0*(z0/ks)));

    return dn>1.0e-12 ? ut/(uplus*uplus*dn) : 0.0;
}

void nhflow_momentum_SI::friction(lexer *p, fdm_nhf *d)
{
    // bed friction rate gamma of the bottom layer (A 519 1), implicit in the column operator; the
    // velocity of a face: normal from the face discharge, tangential the mean of the two cells
    double un,ut,nx,ny;

    SLICELOOP1
    gmx(i,j) = 0.0;

    SLICELOOP2
    gmy(i,j) = 0.0;

    SLICELOOP4
    gmc(i,j) = 0.0;

    if(p->A519!=1)
    return;

    k=0;
    const double ds0 = p->DZN[KP];

    SLICELOOP1
    if(Hx(i,j)>0.0)
    {
    fnorm(p,1,i,j,nx,ny);
    un = Qx[IJK]/Hx(i,j);
    ut = (-ny*0.5*(d->U[IJK] + d->U[Ip1JK]) + nx*0.5*(d->V[IJK] + d->V[Ip1JK]))*p->y_dir;

    gmx(i,j) = bed_gamma(p,sqrt(un*un + ut*ut),ds0*Hx(i,j),0.5*(d->ks(i,j) + d->ks(i+1,j)));
    }

    if(p->j_dir==1)
    SLICELOOP2
    if(Hy(i,j)>0.0)
    {
    fnorm(p,2,i,j,nx,ny);
    un = Qy[IJK]/Hy(i,j);
    ut = -ny*0.5*(d->U[IJK] + d->U[IJp1K]) + nx*0.5*(d->V[IJK] + d->V[IJp1K]);

    gmy(i,j) = bed_gamma(p,sqrt(un*un + ut*ut),ds0*Hy(i,j),0.5*(d->ks(i,j) + d->ks(i,j+1)));
    }

    SLICELOOP4
    WETDRY
    if(d->WL(i,j)>p->A544)
    {
    un = d->U[IJK];
    ut = d->V[IJK]*p->y_dir;

    gmc(i,j) = bed_gamma(p,sqrt(un*un + ut*ut),ds0*d->WL(i,j),d->ks(i,j));
    }
}

void nhflow_momentum_SI::eddy_algebraic(lexer *p, fdm_nhf *d)
{
    // A 547 c nu0: nu_t = c kappa (u*_b + u*_s) z' (1 - z'/h) + nu0 at the centre of every layer, z' = sigma h
    // above the bed (the parabolic profiles of the bed and of the surface boundary layer, Davies 1985);
    // u*_b = |u_b|/u+ from the log law of the bottom layer (as bed_gamma: u+ = 1/kappa max(1, ln(30 z0/ks)),
    // z0 = max(dn/2, ks/30)), so that the eddy viscosity belongs to the bed stress of the column operator;
    // u*_s = sqrt(|tau_s|/rho_w) of the met forcing (W 200, wind stress of the step before), 0 without it.
    // All columns incl. the ghost cells (each column from its own state), dry columns nu0.
    const double kappa = 0.4;
    const int ii=i, jj=j, kk=k;

    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    {
        const double h = d->WL(i,j);
        double ustar = 0.0;

        if(h>p->A544 && p->wet[IJ]==1)
        {
            k=0;
            const double dn = p->DZN[KP]*h;
            double ks = d->ks(i,j);
            if(ks<=0.0)
            ks = 0.0001;
            double z0 = 0.5*dn;
            if(30.0*z0<ks)
            z0 = ks/30.0;
            const double uplus = (1.0/kappa)*MAX(1.0,std::log(30.0*(z0/ks)));
            const double ub = sqrt(d->U[IJK]*d->U[IJK] + d->V[IJK]*d->V[IJK]*p->y_dir);
            ustar = ub/uplus;

            if(pmet!=nullptr && pmet->wind)
            ustar += sqrt(sqrt(pmet->Tx(i,j)*pmet->Tx(i,j) + pmet->Ty(i,j)*pmet->Ty(i,j)));
        }

        double zs = 0.0;
        for(k=0; k<p->knoz; ++k)
        {
            const double sc = zs + 0.5*p->DZN[KP];
            zs += p->DZN[KP];
            d->EV[IJK] = p->A547_c*kappa*ustar*sc*(1.0-sc)*MAX(h,0.0) + p->A547_nu0;
        }
    }

    i=ii; j=jj; k=kk;
}

void nhflow_momentum_SI::col_nu(lexer *p, fdm_nhf *d, int i1, int j1, int i2, int j2, double *nu)
{
    // layer viscosities of the column (i1,j1), or the mean of two columns (a face)
    const int ii=i, jj=j, kk=k;

    for(k=0;k<p->knoz;++k)
    {
        if(!visc)
        nu[k] = 0.0;

        else
        {
        i=i1; j=j1;
        const double n1 = d->VISC[IJK] + d->EV[IJK];
        i=i2; j=j2;
        const double n2 = d->VISC[IJK] + d->EV[IJK];
        nu[k] = 0.5*(n1 + n2);
        }
    }

    i=ii; j=jj; k=kk;
}

void nhflow_momentum_SI::col_coef(lexer *p, fdm_nhf *d, double h, int i1, int j1, int i2, int j2, const double *nu,
                                  double *cu, double *cl, double *cd)
{
    // vertical diffusion of the layer variable (UH per unit sigma), as the vertical rows of nhflow_idiff:
    // dt nu_{k+1/2} sigz^2/(dsigP_k dsig_k), nu at the level the mean of the two layers, no flux at the
    // bed and the surface (the bed shear is gamma);
    // vertical advection (sigma velocity omega of the step before, the mean of the two columns at a
    // face), upwind and implicit: the flux omega u at k+1/2 over dsig_k h, zero at the bed and the
    // surface; conservative (the columns of the matrix weighted with dsig sum to one), so omega does
    // not limit the time step
    const int ii=i, jj=j, kk=k;
    const double hh = MAX(h,p->A544);
    const double s2 = 1.0/(hh*hh);
    const int K = p->knoz;

    auto om = [&](int kf) -> double
    {
        if(kf<=0 || kf>=K)
        return 0.0;

        k=kf;
        i=i1; j=j1;
        const double o1 = d->omegaF[FIJK];
        i=i2; j=j2;
        const double o2 = d->omegaF[FIJK];
        return 0.5*(o1 + o2);
    };

    for(int q=0;q<K;++q)
    {
    k=q;
    cu[q] = q<K-1 ? p->dt*0.5*(nu[q]+nu[q+1])*s2/(p->DZP[KP]*p->DZN[KP]) : 0.0;
    cl[q] = q>0 ? p->dt*0.5*(nu[q]+nu[q-1])*s2/(p->DZP[KM1]*p->DZN[KP]) : 0.0;

    const double a = p->dt/(p->DZN[KP]*hh);
    const double ot = om(q+1);
    const double ob = om(q);

    cu[q] += a*MAX(-ot,0.0);
    cl[q] += a*MAX(ob,0.0);
    cd[q] = a*(ot - ob);
    }

    i=ii; j=jj; k=kk;
}

void nhflow_momentum_SI::col_solve(lexer *p, const double *cu, const double *cl, const double *cd, double dtgam, const double *b, double *x)
{
    // A x = b, A = I - dt L_v + dt V_omega + dt gamma e_0 e_0^T: tridiagonal, upper -cu, lower -cl,
    // diagonal 1 + cu + cl + cd (an M-matrix, column dominant with the weights dsig): Thomas algorithm
    const int K = p->knoz;
    double *cp = wa.data();

    double diag = 1.0 + cu[0] + cl[0] + cd[0] + dtgam;
    cp[0] = -cu[0]/diag;
    x[0] = b[0]/diag;

    for(int q=1;q<K;++q)
    {
    diag = 1.0 + cu[q] + cl[q] + cd[q] + cl[q]*cp[q-1];
    cp[q] = -cu[q]/diag;
    x[q] = (b[q] + cl[q]*x[q-1])/diag;
    }

    for(int q=K-2;q>=0;--q)
    x[q] -= cp[q]*x[q+1];
}

void nhflow_momentum_SI::bc_setup(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // open boundary faces of the domain edges, from the ghost cells set by ioflow / iowave
    //  discharge (inflow B 60, wave generation B 98 >= 3, clamped discharge edge B 520 6):
    //     the layer discharges of the ghost state, Q_k = u_n,k h_g
    //  water level (outflow with F 50 2/4, clamped level edge B 520 5): Dirichlet eta_b
    //  Riemann edge (B 520 3), active absorption (B 99 >= 3), free outflow: the incoming
    //     invariant of the ghost column R = u_n,g + 2 sqrt(g h_g), linearised column discharge
    //     q_in = h (R - 2c) + (R - 3c)(eta^{n+1} - eta^n), implicit in eta^{n+1}
    //  Flather edge (B 520 4): q_in = h (u_n,g + sqrt(g/h) eta^n) - sqrt(g h) eta^{n+1}
    //  radiating faces distribute q over the layers with the profile of the upwind side
    //  the ghost state is that of the start of the step (one step behind the background)
    const double dt = p->dt;
    const int K = p->knoz;
    double eu,df;

    bcf.clear();

    // WL, eta with the boundary values: water level outflow (ioflow), open edges (iowave); the
    // routines write eta (ioflow also the outflow cells), which is restored
    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    {
    WLb(i,j) = d->WL(i,j);
    tE(i,j) = d->eta(i,j);
    }

    pflow->fsfinflow_nhflow(p,d,pgc,WLb);

    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    {
    etaB(i,j) = d->eta(i,j);
    d->eta(i,j) = tE(i,j);
    }

    const int inflow_open  = (p->B98>=3 || p->B60>=1) ? 1 : 0;
    const int outflow_open = (p->B99>=3 || p->B60>=1) ? 1 : 0;
    const int open_side[5] = {0,p->open_xm,p->open_yp,p->open_ym,p->open_xp};

    GCSL4LOOP
    {
    bc_face f;

    f.i  = i = p->gcbsl4[n][0];
    f.j  = j = p->gcbsl4[n][1];
    f.cs = p->gcbsl4[n][3];
    const int bc = p->gcbsl4[n][4];

        if(p->j_dir==0 && (f.cs==2 || f.cs==3))
        continue;

        if(f.cs==1) {f.fi=i-1; f.fj=j; f.gi=i-1; f.gj=j; f.s= 1.0;}
        if(f.cs==4) {f.fi=i;   f.fj=j; f.gi=i+1; f.gj=j; f.s=-1.0;}
        if(f.cs==3) {f.fi=i; f.fj=j-1; f.gi=i; f.gj=j-1; f.s= 1.0;}
        if(f.cs==2) {f.fi=i; f.fj=j;   f.gi=i; f.gj=j+1; f.s=-1.0;}

        const bool xf = (f.cs==1 || f.cs==4);
        const int dir = xf ? 1 : 2;

        // geometry of the face: normal distance of the centres, length, cell width across it, the unit
        // normal towards +i / +j (curvilinear grid: from the geometry layer)
        f.dxp = fdn(p,dir,f.fi,f.fj);
        f.ds = flen(p,dir,f.fi,f.fj);
        f.dxn = carea(p,i,j)/MAX(f.ds,1.0e-20);
        fnorm(p,dir,f.fi,f.fj,f.nx,f.ny);

        if(p->geo_curv!=1)
        f.dxn = xf ? p->DXN[IP] : p->DYN[JP];

        slice &Hf  = xf ? (slice&)Hx  : (slice&)Hy;
        slice &alf = xf ? (slice&)alx : (slice&)aly;
        slice &gmf = xf ? (slice&)gmx : (slice&)gmy;
        slice &af  = xf ? (slice&)ax  : (slice&)ay;
        slice &kf  = xf ? (slice&)kx  : (slice&)ky;
        slice &qEf = xf ? (slice&)qxE : (slice&)qyE;
        slice &gf  = xf ? (slice&)gx0 : (slice&)gy0;
        slice &qf  = xf ? (slice&)qx  : (slice&)qy;
        double *Qf = xf ? Qx : Qy;

        Hf(f.fi,f.fj) = 0.0;
        alf(f.fi,f.fj) = 1.0;
        gmf(f.fi,f.fj) = 0.0;
        af(f.fi,f.fj) = 0.0;
        kf(f.fi,f.fj) = 0.0;
        qEf(f.fi,f.fj) = 0.0;
        gf(f.fi,f.fj) = 0.0;

        auto Qface = [&](int kk) -> double& { return Qf[(f.fi-p->imin)*p->jmax*p->kmax + (f.fj-p->jmin)*p->kmax + kk-p->kmin]; };
        auto cidx = [&](int ii, int jj, int kk) { return (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin; };

        // boundary type
        f.type = 0;

        if(open_side[f.cs]==1 && method[f.cs]>0)
        {
        if(method[f.cs]==3) f.type = bc_riemann;
        if(method[f.cs]==4) f.type = bc_flather;
        if(method[f.cs]==5) f.type = bc_level;
        if(method[f.cs]==6) f.type = bc_flux;
        }

        else
        if((bc==1 || bc==6) && inflow_open==1)
        f.type = bc_flux;

        else
        if((bc==2 || bc==7 || bc==8) && outflow_open==1)
        f.type = ((p->F50==2 || p->F50==4) && p->B99<3) ? bc_level : bc_riemann;

        // wall: the face stays closed
        if(f.type==0)
        {
        qf(f.fi,f.fj) = 0.0;

        for(k=0;k<K;++k)
        Qface(k) = 0.0;

        continue;
        }

        // interior and ghost state of t^n
        const double hi = WL0(i,j);
        const double etai = eta0(i,j);
        const double hg = WLb(f.gi,f.gj);

        f.Qg.assign(K,0.0);
        f.wo.assign(K,1.0);
        f.ug.assign(K,0.0);
        f.vg.assign(K,0.0);
        f.wg.assign(K,0.0);

        double ung=0.0, uni=0.0;

        for(k=0;k<K;++k)
        {
        const int g0 = cidx(f.gi,f.gj,k);
        const int c0 = cidx(i,j,k);

        f.ug[k] = d->U[g0];
        f.vg[k] = d->V[g0]*p->y_dir;
        f.wg[k] = d->W[g0];

        ung += p->DZN[KP]*f.s*(d->U[g0]*f.nx + d->V[g0]*f.ny);
        uni += p->DZN[KP]*f.s*(d->U[c0]*f.nx + d->V[c0]*f.ny);
        }

        f.qin = 0.0;
        f.A = 0.0;
        f.B = 0.0;

        // a dry boundary cell: only an inflowing discharge
        if((p->wet[IJ]==0 || hi<=p->A544) && f.type!=bc_flux)
        {
        qf(f.fi,f.fj) = 0.0;

        for(k=0;k<K;++k)
        Qface(k) = 0.0;

        continue;
        }

        // layer profiles of the radiating faces: inflow of the ghost column, outflow of the interior
        auto profile = [&](vector<double> &w, bool ghost, double ubar)
        {
            bool ok = fabs(ubar)>1.0e-6;

            for(k=0;k<K && ok;++k)
            {
            const int c = ghost ? cidx(f.gi,f.gj,k) : cidx(i,j,k);
            w[k] = f.s*(d->U[c]*f.nx + d->V[c]*f.ny)/ubar;

            if(w[k]<0.0 || w[k]>3.0)
            ok = false;
            }

            if(!ok)
            for(k=0;k<K;++k)
            w[k] = 1.0;
        };

        if(f.type==bc_flux)
        {
        const double h = hg>p->A544 ? hg : hi;

            for(k=0;k<K;++k)
            {
            const int g0 = cidx(f.gi,f.gj,k);
            f.Qg[k] = f.s*(d->U[g0]*f.nx + d->V[g0]*f.ny)*h;
            f.qin += p->DZN[KP]*f.Qg[k];
            }

            if(!qini)
            {
            for(k=0;k<K;++k)
            Qface(k) = f.s*f.Qg[k];

            qf(f.fi,f.fj) = f.s*f.qin;
            }
        }

        if(f.type==bc_riemann || f.type==bc_flather)
        {
        profile(f.Qg,true,ung);
        profile(f.wo,false,uni);

            if(f.type==bc_riemann)
            {
            const double ci = sqrt(g*hi);
            const double R = ung + 2.0*sqrt(g*MAX(hg,1.0e-6));

            f.B = MAX(3.0*ci - R, 0.0);
            f.A = hi*(R - 2.0*ci) + f.B*etai;
            }

            if(f.type==bc_flather)
            {
            f.B = sqrt(g*hi);
            f.A = hi*(ung + sqrt(g/hi)*etai);
            }

            if(!qini)
            {
            qf(f.fi,f.fj) = 0.0;

                for(k=0;k<K;++k)
                {
                Qface(k) = hi*(d->U[cidx(i,j,k)]*f.nx + d->V[cidx(i,j,k)]*f.ny);
                qf(f.fi,f.fj) += p->DZN[KP]*Qface(k);
                }
            }
        }

        if(f.type==bc_level)
        {
        const double etab = etaB(f.gi,f.gj);

        df = 0.5*(d->depth(i,j) + d->depth(f.gi,f.gj));

            if(!qini)
            {
            qf(f.fi,f.fj) = 0.0;

                for(k=0;k<K;++k)
                {
                Qface(k) = UH0[cidx(i,j,k)]*f.nx + VH0[cidx(i,j,k)]*f.ny;
                qf(f.fi,f.fj) += p->DZN[KP]*Qface(k);
                }
            }

        // eta upwind of the face discharge (inward: the boundary value)
        const double qin_n = f.s*qf(f.fi,f.fj);

        if(qin_n>1.0e-12)
        eu = etab;

        else
        if(qin_n<-1.0e-12)
        eu = etai;

        else
        eu = MAX(etab,etai);

        const double h = eu + df;

            if(h<=p->A544)
            {
            qf(f.fi,f.fj) = 0.0;

            for(k=0;k<K;++k)
            Qface(k) = 0.0;

            continue;
            }

        Hf(f.fi,f.fj) = h;

            if(p->A519==1)
            {
            k=0;
            const double ds0 = p->DZN[KP];
            const double un = Qface(0)/h;
            const int c0 = cidx(i,j,0);
            const double ut = (-f.ny*d->U[c0] + f.nx*d->V[c0])*p->y_dir;

            gmf(f.fi,f.fj) = bed_gamma(p,sqrt(un*un + ut*ut),ds0*h,d->ks(i,j));
            }

        // the boundary value for the face gradients of t^n and t^{n+1} (held over the step)
        etaI(f.gi,f.gj) = etab;

        // face gradient in the face direction: (right - left)/dxp
        gf(f.fi,f.fj) = f.s*(etai - etab)/f.dxp;
        }

    bcf.push_back(f);
    }
}

void nhflow_momentum_SI::stage_arrays(int s, double *&UHs, double *&VHs, double *&WHs, double *&UHo, double *&VHo, double *&WHo, double &a)
{
    // SSP-RK of the explicit terms on the sigma grid of t^n (A 510 2: stages 0, 1; 3: 0, 1, 2), the
    // layer transport and omega of the previous step frozen
    double *in[3][3]  = {{UH0,VH0,WH0},{UHa,VHa,WHa},{UHb,VHb,WHb}};
    double *out[3][3] = {{UHa,VHa,WHa},{UHb,VHb,WHb},{UHp,VHp,WHp}};
    const int so = (s==rk-1) ? 2 : s;

    UHs = in[s][0]; VHs = in[s][1]; WHs = in[s][2];
    UHo = out[so][0]; VHo = out[so][1]; WHo = out[so][2];

    if(rk==2)
    a = s==0 ? 0.0 : 0.5;

    else
    a = s==0 ? 0.0 : (s==1 ? 0.75 : 1.0/3.0);
}

void nhflow_momentum_SI::si_stage_vel(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S, int s)
{
    double *UHs,*VHs,*WHs,*UHo,*VHo,*WHo,a;
    stage_arrays(s,UHs,VHs,WHs,UHo,VHo,WHo,a);

    // ghost cells of the stage input, velocities of the stage on the water depth of t^n
    pgc->start4V(p,UHs,gcval_uh);
    pgc->start4V(p,VHs,gcval_vh);
    pgc->start4V(p,WHs,gcval_wh);

    // below the depth A 546 desingularised (Kurganov and Petrova 2007), as sflow_momentum_SI::rk_stage:
    // u = sqrt(2) h (uh)/sqrt(h^4 + max(h^4, h_v^4)); a wetting cell gets the momentum of its deep
    // neighbours with the inflow while its depth is that of t^n
    const double hv4 = p->A546*p->A546*p->A546*p->A546;

    LOOP
    {
        if(p->wet[IJ]==1 && WL0(i,j)>p->A544)
        {
            if(WL0(i,j)>=p->A546)
            {
            d->U[IJK] = UHs[IJK]/WL0(i,j);
            d->V[IJK] = VHs[IJK]/WL0(i,j)*p->y_dir;
            d->W[IJK] = WHs[IJK]/WL0(i,j);
            }

            else
            {
            const double h = WL0(i,j);
            const double r = sqrt(2.0)*h/sqrt(h*h*h*h + hv4);
            d->U[IJK] = r*UHs[IJK];
            d->V[IJK] = r*VHs[IJK]*p->y_dir;
            d->W[IJK] = r*WHs[IJK];
            }
        }

        else
        {
        d->U[IJK] = 0.0;
        d->V[IJK] = 0.0;
        d->W[IJK] = 0.0;
        }
    }

    pgc->start4V(p,d->U,gcval_u);
    pgc->start4V(p,d->V,gcval_v);
    pgc->start4V(p,d->W,gcval_w);

}

void nhflow_momentum_SI::si_stage_rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_stage_obj &S, int s)
{
    double *UHs,*VHs,*WHs,*UHo,*VHo,*WHo,a;
    stage_arrays(s,UHs,VHs,WHs,UHo,VHo,WHo,a);
    const double alpha = 1.0-a;

    // sources (as the explicit stage, without the pressure terms and the bed friction)
    pturb->isource(p,d);
    pflow->isource_nhflow(p,d,pgc,pvrans,WL0);
    pwind->wind_forcing_nhf_x(p,d,pgc,d->U,d->V,d->F,WL0,d->eta);
    if(pmet!=nullptr)
    pmet->source_nhflow(p,pgc,d->F,0,WL0,UHs,VHs);
    if(d->wave!=nullptr)
    d->wave->u_source(p,d);
    irhs(p,d,pgc);

    pturb->jsource(p,d);
    pflow->jsource_nhflow(p,d,pgc,pvrans,WL0);
    pwind->wind_forcing_nhf_y(p,d,pgc,d->U,d->V,d->G,WL0,d->eta);
    if(pmet!=nullptr)
    pmet->source_nhflow(p,pgc,d->G,1,WL0,UHs,VHs);
    if(d->wave!=nullptr)
    d->wave->v_source(p,d);
    jrhs(p,d,pgc);

    pturb->ksource(p,d);
    pflow->ksource_nhflow(p,d,pgc,pvrans,WL0);
    krhs(p,d,pgc);

    // advection with the layer transport and omega of the previous step
    advection(p,d,pgc,S.precon);

    // horizontal diffusion along the sigma layers (explicit); the vertical diffusion is implicit in
    // the column operator
    if(visc)
    hdiffusion(p,d);

    LOOP
    {
        if(p->wet[IJ]==1)
        {
        UHo[IJK] = a*UH0[IJK] + alpha*(UHs[IJK] + p->dt*d->F[IJK]);
        VHo[IJK] = p->j_dir==1 ? a*VH0[IJK] + alpha*(VHs[IJK] + p->dt*d->G[IJK]) : 0.0;
        WHo[IJK] = a*WH0[IJK] + alpha*(WHs[IJK] + p->dt*d->H[IJK]);
        }

        else
        {
        UHo[IJK] = 0.0;
        VHo[IJK] = 0.0;
        WHo[IJK] = 0.0;
        }
    }
}

void nhflow_momentum_SI::advection(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_reconstruct *precon)
{
    // conservative horizontal advection of UH, VH, WH with the layer transport T of the previous
    // step, velocities upwind of the flux (reconstructed, A 514); no gravity wave speed in the
    // dissipation: the Courant number of the flow alone limits the time step; the vertical
    // advection (omega) is implicit in the column operator
    const int N = p->imax*p->jmax*(p->kmax+2);

    precon->reconstruct_3D_x(p, pgc, d, d->U, d->Us, d->Un);
    precon->reconstruct_3D_x(p, pgc, d, d->W, d->Ws, d->Wn);

    if(p->j_dir==1)
    {
    precon->reconstruct_3D_y(p, pgc, d, d->U, d->Ue, d->Uw);
    precon->reconstruct_3D_x(p, pgc, d, d->V, d->Vs, d->Vn);
    precon->reconstruct_3D_y(p, pgc, d, d->V, d->Ve, d->Vw);
    precon->reconstruct_3D_y(p, pgc, d, d->W, d->We, d->Ww);
    }

    // comp 0 U into F, 1 V into G, 2 W into H
    auto component = [&](int comp, double *R, double *Ls, double *Ln, double *Le, double *Lw)
    {
        for(int q=0;q<N;++q)
        {
        d->Fx[q] = 0.0;
        d->Fy[q] = 0.0;
        d->Fz[q] = 0.0;
        }

        SLICELOOP1
        if(Hx(i,j)>0.0)
        for(k=0;k<p->knoz;++k)
        d->Fx[IJK] = Tx[IJK]>=0.0 ? Tx[IJK]*Ls[IJK] : Tx[IJK]*Ln[IJK];

        if(p->j_dir==1)
        SLICELOOP2
        if(Hy(i,j)>0.0)
        for(k=0;k<p->knoz;++k)
        d->Fy[IJK] = Ty[IJK]>=0.0 ? Ty[IJK]*Le[IJK] : Ty[IJK]*Lw[IJK];

        exchangeF(p,pgc,d->Fx);
        exchangeF(p,pgc,d->Fy);

        // open boundary faces: first-order upwind, inflow with the ghost velocities
        for(const bc_face &f : bcf)
        {
            const bool xf = (f.cs==1 || f.cs==4);
            double *T = xf ? Tx : Ty;
            double *Fb = xf ? d->Fx : d->Fy;

            for(k=0;k<p->knoz;++k)
            {
            i=f.fi;
            j=f.fj;
            const double t = T[IJK];
            const bool in = f.s*t>0.0;

            i=f.i;
            j=f.j;
            const double uc = comp==0 ? d->U[IJK] : (comp==1 ? d->V[IJK] : d->W[IJK]);
            const double ug = comp==0 ? f.ug[k] : (comp==1 ? f.vg[k] : f.wg[k]);

            i=f.fi;
            j=f.fj;
            Fb[IJK] = t*(in ? ug : uc);
            }
        }

        // the vertical advection is implicit in the column operator (col_coef)
        LOOP
        WETDRY
        R[IJK] -= layer_div(p,d->Fx,d->Fy);
    };

    component(0,d->F,d->Us,d->Un,d->Ue,d->Uw);

    if(p->j_dir==1)
    component(1,d->G,d->Vs,d->Vn,d->Ve,d->Vw);

    component(2,d->H,d->Ws,d->Wn,d->We,d->Ww);
}

void nhflow_momentum_SI::hdiffusion(lexer *p, fdm_nhf *d)
{
    // WL div_h( nu_h grad_h u ) per layer, u = U, V, W of the stage, nu_h = VISC + EV (EVH with the
    // HLES split A 561), the mean of the two cells at a face; compact Laplacian (curvilinear grid: the
    // orthogonal part, nhflow_curv::lapcoef); no flux through walls and to dry cells (slip walls, as
    // the domain walls of nhflow_idiff)
    double *EVh = (p->A561>0) ? d->EVH : d->EV;

    auto nu = [&](int q) { return d->VISC[q] + EVh[q]; };

    LOOP
    WETDRY
    {
    const double wl = WL0(i,j);
    double fu=0.0, fv=0.0, fw=0.0, c;
    double cE,cW,cN,cS;

        if(p->geo_curv==1)
        nhflow_curv::lapcoef(p,i,j,cE,cW,cN,cS);

        else
        {
        cE = 1.0/(p->DXP[IP]*p->DXN[IP]);
        cW = 1.0/(p->DXP[IM1]*p->DXN[IP]);
        cN = p->j_dir==1 ? 1.0/(p->DYP[JP]*p->DYN[JP]) : 0.0;
        cS = p->j_dir==1 ? 1.0/(p->DYP[JM1]*p->DYN[JP]) : 0.0;
        }

        if(p->flagslice4[Ip1J]>0 && p->wet[Ip1J]==1)
        {
        c = 0.5*(nu(IJK) + nu(Ip1JK))*cE;
        fu += c*(d->U[Ip1JK] - d->U[IJK]);
        fv += c*(d->V[Ip1JK] - d->V[IJK]);
        fw += c*(d->W[Ip1JK] - d->W[IJK]);
        }

        if(p->flagslice4[Im1J]>0 && p->wet[Im1J]==1)
        {
        c = 0.5*(nu(IJK) + nu(Im1JK))*cW;
        fu += c*(d->U[Im1JK] - d->U[IJK]);
        fv += c*(d->V[Im1JK] - d->V[IJK]);
        fw += c*(d->W[Im1JK] - d->W[IJK]);
        }

        if(p->j_dir==1)
        {
            if(p->flagslice4[IJp1]>0 && p->wet[IJp1]==1)
            {
            c = 0.5*(nu(IJK) + nu(IJp1K))*cN;
            fu += c*(d->U[IJp1K] - d->U[IJK]);
            fv += c*(d->V[IJp1K] - d->V[IJK]);
            fw += c*(d->W[IJp1K] - d->W[IJK]);
            }

            if(p->flagslice4[IJm1]>0 && p->wet[IJm1]==1)
            {
            c = 0.5*(nu(IJK) + nu(IJm1K))*cS;
            fu += c*(d->U[IJm1K] - d->U[IJK]);
            fv += c*(d->V[IJm1K] - d->V[IJK]);
            fw += c*(d->W[IJm1K] - d->W[IJK]);
            }
        }

    d->F[IJK] += wl*fu;
    d->G[IJK] += wl*fv*p->y_dir;
    d->H[IJK] += wl*fw;
    }
}

double nhflow_momentum_SI::dissipation(lexer *p, double q, double h, double dxp)
{
    // optional mass flux dissipation (A 538): K = A538 (|u| + c) dx/2, the Rusanov dissipation
    // of the continuity equation, switched on from Froude number 0.5 to 1; implicit (in the eta
    // system), for hydraulic jumps and transcritical flow
    if(p->A538<=0.0 || h<=p->A544)
    return 0.0;

    const double u = fabs(q)/h;
    const double c = sqrt(g*h);
    const double fr = u/c;
    const double phi = MAX(MIN((fr-0.5)/0.5,1.0),0.0);

    return p->A538*0.5*(u + c)*dxp*phi;
}

void nhflow_momentum_SI::si_predictor_end(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    clearrhs(p,d,pgc);

    // the cell predictor to the faces and to the cell update (UHD, VHD)
    LOOP
    {
    UHD[IJK] = p->wet[IJ]==1 ? UHp[IJK] : 0.0;
    VHD[IJK] = (p->wet[IJ]==1 && p->j_dir==1) ? VHp[IJK] : 0.0;
    }

    exchangeC(p,pgc,UHD);
    exchangeC(p,pgc,VHD);

}

void nhflow_momentum_SI::si_assemble(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const double dt = p->dt;
    double K,qh;

    // column of the face predictor (momentum interpolation: the mean of the cell predictors,
    // which carry no gravity term of this step, along the face normal), explicit part of the column
    // discharge, eta-system coefficients a_f = theta dt (g theta dt alpha_f H_f + K_f)/dxp;
    // face gradient of eta^n with its non-orthogonal part (curvilinear grid, explicit in both parts)
    double nx,ny;

    SLICELOOP1
    {
        if(Hx(i,j)>0.0)
        {
        fnorm(p,1,i,j,nx,ny);
        const double dn = fdn(p,1,i,j);
        cx0(i,j) = fcorr(p,d->eta,1,i,j);
        gx0(i,j) = (d->eta(i+1,j) - d->eta(i,j))/dn + cx0(i,j);

        // column: Q~ = A^-1 (I - dt L_v) Q^, r = A^-1 1
        col_nu(p,d,i,j,i+1,j,wnu.data());
        col_coef(p,d,Hx(i,j),i,j,i+1,j,wnu.data(),wcu.data(),wcl.data(),wcd.data());

        for(k=0;k<p->knoz;++k)
        wb[k] = 0.5*(UHD[IJK] + UHD[Ip1JK])*nx + 0.5*(VHD[IJK] + VHD[Ip1JK])*ny;

        col_solve(p,wcu.data(),wcl.data(),wcd.data(),dt*gmx(i,j),wb.data(),wr.data());

        qh = 0.0;
        for(k=0;k<p->knoz;++k)
        {
        QAx[IJK] = wr[k];
        qh += p->DZN[KP]*wr[k];
        wb[k] = 1.0;
        }

        col_solve(p,wcu.data(),wcl.data(),wcd.data(),dt*gmx(i,j),wb.data(),wr.data());

        alx(i,j) = 0.0;
        for(k=0;k<p->knoz;++k)
        {
        Rx[IJK] = wr[k];
        alx(i,j) += p->DZN[KP]*wr[k];
        }

        qxE(i,j) = qh - dt*g*Hx(i,j)*alx(i,j)*((1.0-theta)*gx0(i,j) + theta*cx0(i,j));

        K = dissipation(p,qxn(i,j),Hx(i,j),dn);
        kx(i,j) = K;

        ax(i,j) = theta*dt*(g*theta*dt*alx(i,j)*Hx(i,j) + K)/dn;
        }

        else
        {
        gx0(i,j) = 0.0;
        cx0(i,j) = 0.0;
        qxE(i,j) = 0.0;
        kx(i,j) = 0.0;
        ax(i,j) = 0.0;
        }
    }

    if(p->j_dir==1)
    SLICELOOP2
    {
        if(Hy(i,j)>0.0)
        {
        fnorm(p,2,i,j,nx,ny);
        const double dn = fdn(p,2,i,j);
        cy0(i,j) = fcorr(p,d->eta,2,i,j);
        gy0(i,j) = (d->eta(i,j+1) - d->eta(i,j))/dn + cy0(i,j);

        col_nu(p,d,i,j,i,j+1,wnu.data());
        col_coef(p,d,Hy(i,j),i,j,i,j+1,wnu.data(),wcu.data(),wcl.data(),wcd.data());

        for(k=0;k<p->knoz;++k)
        wb[k] = 0.5*(UHD[IJK] + UHD[IJp1K])*nx + 0.5*(VHD[IJK] + VHD[IJp1K])*ny;

        col_solve(p,wcu.data(),wcl.data(),wcd.data(),dt*gmy(i,j),wb.data(),wr.data());

        qh = 0.0;
        for(k=0;k<p->knoz;++k)
        {
        QAy[IJK] = wr[k];
        qh += p->DZN[KP]*wr[k];
        wb[k] = 1.0;
        }

        col_solve(p,wcu.data(),wcl.data(),wcd.data(),dt*gmy(i,j),wb.data(),wr.data());

        aly(i,j) = 0.0;
        for(k=0;k<p->knoz;++k)
        {
        Ry[IJK] = wr[k];
        aly(i,j) += p->DZN[KP]*wr[k];
        }

        qyE(i,j) = qh - dt*g*Hy(i,j)*aly(i,j)*((1.0-theta)*gy0(i,j) + theta*cy0(i,j));

        K = dissipation(p,qyn(i,j),Hy(i,j),dn);
        ky(i,j) = K;

        ay(i,j) = theta*dt*(g*theta*dt*aly(i,j)*Hy(i,j) + K)/dn;
        }

        else
        {
        gy0(i,j) = 0.0;
        cy0(i,j) = 0.0;
        qyE(i,j) = 0.0;
        ky(i,j) = 0.0;
        ay(i,j) = 0.0;
        }
    }

    exchange1(p,pgc,gx0);
    exchange1(p,pgc,qxE);
    exchange1(p,pgc,ax);
    exchange1(p,pgc,kx);
    exchange1(p,pgc,cx0);
    exchange2(p,pgc,cy0);
    exchange2(p,pgc,gy0);
    exchange2(p,pgc,qyE);
    exchange2(p,pgc,ay);
    exchange2(p,pgc,ky);

    // open boundary faces
    SLICELOOP4
    bdiag(i,j) = 0.0;

    for(const bc_face &f : bcf)
    {
    i=f.i;
    j=f.j;

        const bool xf = (f.cs==1 || f.cs==4);

        slice &Hf  = xf ? (slice&)Hx  : (slice&)Hy;
        slice &alf = xf ? (slice&)alx : (slice&)aly;
        slice &gmf = xf ? (slice&)gmx : (slice&)gmy;
        slice &af  = xf ? (slice&)ax  : (slice&)ay;
        slice &qEf = xf ? (slice&)qxE : (slice&)qyE;
        slice &gf  = xf ? (slice&)gx0 : (slice&)gy0;
        double *QA = xf ? QAx : QAy;
        double *R  = xf ? Rx : Ry;

        if(f.type==bc_level && Hf(f.fi,f.fj)>0.0)
        {
        // column of the boundary face with the interior predictor and viscosity
        col_nu(p,d,i,j,i,j,wnu.data());
        col_coef(p,d,Hf(f.fi,f.fj),i,j,i,j,wnu.data(),wcu.data(),wcl.data(),wcd.data());

        for(k=0;k<p->knoz;++k)
        wb[k] = UHD[IJK]*f.nx + VHD[IJK]*f.ny;

        col_solve(p,wcu.data(),wcl.data(),wcd.data(),dt*gmf(f.fi,f.fj),wb.data(),wr.data());

        const int fi=f.fi, fj=f.fj;
        auto fidx = [&](int kk) { return (fi-p->imin)*p->jmax*p->kmax + (fj-p->jmin)*p->kmax + kk-p->kmin; };

        qh = 0.0;
        for(k=0;k<p->knoz;++k)
        {
        QA[fidx(k)] = wr[k];
        qh += p->DZN[KP]*wr[k];
        wb[k] = 1.0;
        }

        col_solve(p,wcu.data(),wcl.data(),wcd.data(),dt*gmf(f.fi,f.fj),wb.data(),wr.data());

        alf(f.fi,f.fj) = 0.0;
        for(k=0;k<p->knoz;++k)
        {
        R[fidx(k)] = wr[k];
        alf(f.fi,f.fj) += p->DZN[KP]*wr[k];
        }

        qEf(f.fi,f.fj) = qh - dt*g*Hf(f.fi,f.fj)*alf(f.fi,f.fj)*(1.0-theta)*gf(f.fi,f.fj);
        af(f.fi,f.fj) = theta*dt*g*theta*dt*alf(f.fi,f.fj)*Hf(f.fi,f.fj)/f.dxp;
        }

        if(f.type==bc_flux)
        qEf(f.fi,f.fj) = f.s*f.qin;

        // radiating: q_in = A - B eta^{n+1}; A in the right-hand side, B on the diagonal
        if(f.type==bc_riemann || f.type==bc_flather)
        {
        qEf(f.fi,f.fj) = f.s*f.A;
        bdiag(i,j) += theta*dt*f.B/f.dxn;
        }
    }

    // a face closed in this step (dry, H_f = 0) carries no old discharge either: the explicit
    // part (1-theta) q^n would drain a cell without coupling in the eta system below zero
    SLICELOOP1
    if(Hx(i,j)<=0.0)
    {
    qxn(i,j) = 0.0;

    for(k=0;k<p->knoz;++k)
    Qxn[IJK] = 0.0;
    }

    if(p->j_dir==1)
    SLICELOOP2
    if(Hy(i,j)<=0.0)
    {
    qyn(i,j) = 0.0;

    for(k=0;k<p->knoz;++k)
    Qyn[IJK] = 0.0;
    }

    exchange1(p,pgc,qxn);
    exchange2(p,pgc,qyn);
    exchangeF(p,pgc,Qxn);
    exchangeF(p,pgc,Qyn);

    for(const bc_face &f : bcf)
    if(f.type==bc_level && (f.cs==1 || f.cs==4 ? Hx(f.fi,f.fj) : Hy(f.fi,f.fj))<=0.0)
    {
        const bool xf = (f.cs==1 || f.cs==4);
        double *Qn = xf ? Qxn : Qyn;

        if(xf)
        qxn(f.fi,f.fj) = 0.0;

        else
        qyn(f.fi,f.fj) = 0.0;

        i=f.fi;
        j=f.fj;
        for(k=0;k<p->knoz;++k)
        Qn[IJK] = 0.0;
    }

    // right-hand side of the water depth: WL^n - dt div( theta q^E + (1-theta) q^n )
    SLICELOOP4
    rhsW(i,j) = WL0(i,j) - dt*(theta*div(p,qxE,qyE) + (1.0-theta)*div(p,qxn,qyn));

    // + dt K: the Stokes mass source of the SEASTATE vortex force (A 751 2), K = -div(M) of the wet cells, held
    // over the step (also in eta_update)
    mass_source = (d->wave!=nullptr && p->A751==2);

    if(mass_source)
    {
        SLICELOOP4
        wsrc(i,j) = 0.0;

    d->wave->mass_source(p,d,wsrc);

        SLICELOOP4
        rhsW(i,j) += dt*wsrc(i,j);
    }

    // initial guess eta^n
    SLICELOOP4
    etaI(i,j) = eta0(i,j);
}

void nhflow_momentum_SI::face_update(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // Q_k^{n+1} = r_k ( Q^_k - dt g H_f (theta G^{n+1} + (1-theta) G^n) ) - K_f G^{n+1}
    const double dt = p->dt;
    const int N = p->imax*p->jmax*(p->kmax+2);
    double G1,G;

    SLICELOOP1
    {
    qx(i,j) = 0.0;

        if(Hx(i,j)>0.0)
        {
        G1 = (etaI(i+1,j) - etaI(i,j))/fdn(p,1,i,j) + cx0(i,j);
        G = dt*g*Hx(i,j)*(theta*G1 + (1.0-theta)*gx0(i,j));

            for(k=0;k<p->knoz;++k)
            {
            Qx[IJK] = QAx[IJK] - G*Rx[IJK] - kx(i,j)*G1;
            qx(i,j) += p->DZN[KP]*Qx[IJK];
            }
        }

        else
        for(k=0;k<p->knoz;++k)
        Qx[IJK] = 0.0;
    }

    if(p->j_dir==1)
    SLICELOOP2
    {
    qy(i,j) = 0.0;

        if(Hy(i,j)>0.0)
        {
        G1 = (etaI(i,j+1) - etaI(i,j))/fdn(p,2,i,j) + cy0(i,j);
        G = dt*g*Hy(i,j)*(theta*G1 + (1.0-theta)*gy0(i,j));

            for(k=0;k<p->knoz;++k)
            {
            Qy[IJK] = QAy[IJK] - G*Ry[IJK] - ky(i,j)*G1;
            qy(i,j) += p->DZN[KP]*Qy[IJK];
            }
        }

        else
        for(k=0;k<p->knoz;++k)
        Qy[IJK] = 0.0;
    }

    // open boundary faces
    for(const bc_face &f : bcf)
    {
        const bool xf = (f.cs==1 || f.cs==4);

        slice &Hf  = xf ? (slice&)Hx  : (slice&)Hy;
        slice &gf  = xf ? (slice&)gx0 : (slice&)gy0;
        slice &qf  = xf ? (slice&)qx  : (slice&)qy;
        double *Qf = xf ? Qx : Qy;
        double *QA = xf ? QAx : QAy;
        double *R  = xf ? Rx : Ry;

        double *lay = wb.data();
        const int K = p->knoz;

        if(f.type==bc_level)
        {
            if(Hf(f.fi,f.fj)>0.0)
            {
            G1 = f.s*(etaI(f.i,f.j) - etaI(f.gi,f.gj))/f.dxp;
            G = dt*g*Hf(f.fi,f.fj)*(theta*G1 + (1.0-theta)*gf(f.fi,f.fj));

            i=f.fi;
            j=f.fj;
            for(k=0;k<K;++k)
            lay[k] = QA[IJK] - G*R[IJK];
            }

            else
            for(k=0;k<K;++k)
            lay[k] = 0.0;
        }

        if(f.type==bc_flux)
        for(k=0;k<K;++k)
        lay[k] = f.s*f.Qg[k];

        if(f.type==bc_riemann || f.type==bc_flather)
        {
        const double q = f.s*(f.A - f.B*etaI(f.i,f.j));
        const bool in = f.s*q>0.0;

        for(k=0;k<K;++k)
        lay[k] = q*(in ? f.Qg[k] : f.wo[k]);
        }

    i=f.fi;
    j=f.fj;
    qf(i,j) = 0.0;

        for(k=0;k<K;++k)
        {
        Qf[IJK] = lay[k];
        qf(i,j) += p->DZN[KP]*lay[k];
        }
    }

    exchange1(p,pgc,qx);
    exchange2(p,pgc,qy);
    exchangeF(p,pgc,Qx);
    exchangeF(p,pgc,Qy);

    // layer transport of the step
    for(int q=0;q<N;++q)
    {
    Tx[q] = theta*Qx[q] + (1.0-theta)*Qxn[q];
    Ty[q] = theta*Qy[q] + (1.0-theta)*Qyn[q];
    }
}

double nhflow_momentum_SI::layer_div(lexer *p, double *fx, double *fy)
{
    if(p->geo_curv==1)
    return nhflow_curv::divh(p,fx,fy,i,j,k);

    double s = (fx[IJK] - fx[Im1JK])/p->DXN[IP];

    if(p->j_dir==1)
    s += (fy[IJK] - fy[IJm1K])/p->DYN[JP];

    return s;
}

void nhflow_momentum_SI::eta_update(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_fsf *pfsf)
{
    starttime=pgc->timer();

    // WL from the divergence of the layer transport: the volume is conserved exactly, independent
    // of the tolerance of the eta solve, and omega (layer continuity with the same transport)
    // vanishes at the surface; with the wetting and drying iteration WL >= 0 up to round-off
    SLICELOOP4
    {
    double s=0.0;

    for(k=0;k<p->knoz;++k)
    s += p->DZN[KP]*layer_div(p,Tx,Ty);

    double w = WL0(i,j) - p->dt*s;

    if(mass_source)
    w += p->dt*wsrc(i,j);

    d->WL(i,j) = MAX(w, 0.0);
    d->eta(i,j) = d->WL(i,j) - d->depth(i,j);
    d->detadt(i,j) = (d->WL(i,j) - WL0(i,j))/p->dt;
    }

    // discharge over the open boundaries (theta weighted, inward positive)
    Qin = Qout = 0.0;

    for(const bc_face &f : bcf)
    {
        const bool xf = (f.cs==1 || f.cs==4);
        double *T = xf ? Tx : Ty;
        double q=0.0;

        i=f.fi;
        j=f.fj;
        for(k=0;k<p->knoz;++k)
        q += p->DZN[KP]*T[IJK];

        q *= f.s*f.ds;

        if(q>0.0)
        Qin += q;

        else
        Qout += q;
    }

    Qin = pgc->globalsum(Qin);
    Qout = pgc->globalsum(Qout);

    // ghost cells; the interior values are kept
    SLICELOOP4
    {
    tE(i,j) = d->eta(i,j);
    tS(i,j) = d->WL(i,j);
    }

    pgc->gcsl_start4(p,d->eta,gcval_eta);
    pgc->gcsl_start4(p,d->WL,gcval_eta);
    pgc->gcsl_start4(p,d->detadt,1);

    SLICELOOP4
    {
    d->eta(i,j) = tE(i,j);
    d->WL(i,j) = tS(i,j);
    }

    level_ghosts(p,d,pgc);

    // wetting & drying flags (A 540) on a copy of the water depth: the explicit wetdry lifts dry
    // cells to the film A 544; the semi-implicit step keeps its volume
    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    {
    WLb(i,j) = d->WL(i,j);
    tE(i,j) = d->eta(i,j);
    }

    pfsf->wetdry(p,d,pgc,d->UH,d->VH,d->WH,WLb);

    // a dry cell that took water through a flooding face of the step (water depth above the film A 544)
    // is wet: the explicit rule (A 540 1: wet when a wet neighbour stands higher) keeps such a cell dry
    // once its level stands above the neighbour's, and the flooding face of the next steps fills it
    // further (a one-way valve at the shoreline under onshore forcing, e.g. wave setup)
    {
    int nw = 0;

    SLICELOOP4
    if(p->wet[IJ]==0 && p->flagslice4[IJ]>0 && d->WL(i,j)>2.0*p->A544
       && (p->wetfix==nullptr || p->wetfix[IJ]<0))
    {
    p->wet[IJ] = 1;
    ++nw;
    }

    if(pgc->globalisum(nw)>0)
    pgc->gcsl_start4Vint(p,p->wet,50);
    }

    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    d->eta(i,j) = tE(i,j);

    p->fsftime+=pgc->timer()-starttime;
}

void nhflow_momentum_SI::level_ghosts(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // ghost values of the water level and radiating open boundaries (the boundary state of the step)
    for(const bc_face &f : bcf)
    if(f.type!=bc_flux)
    {
        const int di = f.gi-f.i, dj = f.gj-f.j;

        for(int q=1;q<=3;++q)
        {
        d->eta(f.i+q*di,f.j+q*dj) = etaB(f.gi,f.gj);
        d->WL(f.i+q*di,f.j+q*dj) = WLb(f.gi,f.gj);
        }
    }
}

void nhflow_momentum_SI::cell_update(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // cell momentum: the gravity term of the cell is the mean of that of its open faces, H_f times the
    // face gradient (theta weighted), the same in every layer; with H_f instead of the depth of the
    // cell, a thin face gives the cell no more momentum than its own discharge gets. A face towards a
    // cell that is dry at the end of the step (eta^{n+1} + d <= A 544: the eta of the Newton iteration
    // is not a water level there, it may lie far below the bed) does not count: its gradient carries
    // the drainage of the last film of that cell, not a pressure gradient of the wet cell
    const double dt = p->dt;
    double we,ww,wsum,gu,gv;

    auto wn = [&](int a, int b) { return etaI(a,b) + d->depth(a,b) > p->A544; };

    SLICELOOP4
    {
        if(p->wet[IJ]==1)
        {
        we = (Hx(i,j)>0.0 && wn(i+1,j)) ? 1.0 : 0.0;
        ww = (Hx(i-1,j)>0.0 && wn(i-1,j)) ? 1.0 : 0.0;
        wsum = MAX(we+ww,1.0);

        gu = dt*g*(ww*Hx(i-1,j)*(theta*(etaI(i,j) - etaI(i-1,j))/p->DXP[IM1] + (1.0-theta)*gx0(i-1,j))
                 + we*Hx(i,j)*(theta*(etaI(i+1,j) - etaI(i,j))/p->DXP[IP] + (1.0-theta)*gx0(i,j)))/wsum;
        gv = 0.0;

            if(p->j_dir==1)
            {
            we = (Hy(i,j)>0.0 && wn(i,j+1)) ? 1.0 : 0.0;
            ww = (Hy(i,j-1)>0.0 && wn(i,j-1)) ? 1.0 : 0.0;
            wsum = MAX(we+ww,1.0);

            gv = dt*g*(ww*Hy(i,j-1)*(theta*(etaI(i,j) - etaI(i,j-1))/p->DYP[JM1] + (1.0-theta)*gy0(i,j-1))
                     + we*Hy(i,j)*(theta*(etaI(i,j+1) - etaI(i,j))/p->DYP[JP] + (1.0-theta)*gy0(i,j)))/wsum;
            }

            // curvilinear grid: per index direction the mean normal and the mean of H_f times the normal
            // gradient of the open faces, the cell vector from the two (on Cartesian nodes the mean above);
            // with open faces in one direction only its component along that normal (no tangential
            // part: two nearly parallel faces of a skewed cell do not fix it)
            if(p->geo_curv==1)
            {
            double nx, ny;
            double ax_=0.0, ay_=0.0, ag=0.0, bx_=0.0, by_=0.0, bg=0.0;
            int na=0, nb=0;

            auto face = [&](int dir, int fi, int fj, double H, double G0, double G1)
            {
                // the other cell of the face
                const int oi = dir==1 ? (fi==i ? i+1 : i-1) : i;
                const int oj = dir==2 ? (fj==j ? j+1 : j-1) : j;

                if(H<=0.0 || !wn(oi,oj))
                return;

                fnorm(p,dir,fi,fj,nx,ny);
                const double G = H*(theta*G1 + (1.0-theta)*G0);

                if(dir==1)
                {
                ax_ += nx; ay_ += ny; ag += G; ++na;
                }

                else
                {
                bx_ += nx; by_ += ny; bg += G; ++nb;
                }
            };

            face(1,i,j,Hx(i,j),gx0(i,j),(etaI(i+1,j) - etaI(i,j))/fdn(p,1,i,j) + cx0(i,j));
            face(1,i-1,j,Hx(i-1,j),gx0(i-1,j),(etaI(i,j) - etaI(i-1,j))/fdn(p,1,i-1,j) + cx0(i-1,j));
            face(2,i,j,Hy(i,j),gy0(i,j),(etaI(i,j+1) - etaI(i,j))/fdn(p,2,i,j) + cy0(i,j));
            face(2,i,j-1,Hy(i,j-1),gy0(i,j-1),(etaI(i,j) - etaI(i,j-1))/fdn(p,2,i,j-1) + cy0(i,j-1));

            double gx=0.0, gy=0.0;

            // unit mean normals, mean gradients
            auto unit = [&](double &x, double &y, double &gg, int n)
            {
                const double l = sqrt(x*x + y*y);
                x = l>0.0 ? x/l : 0.0;
                y = l>0.0 ? y/l : 0.0;
                gg = n>0 ? gg/double(n) : 0.0;
            };

            unit(ax_,ay_,ag,na);
            unit(bx_,by_,bg,nb);

            const double det = ax_*by_ - ay_*bx_;

            if(na>0 && nb>0 && fabs(det)>1.0e-3)
            {
            gx = (ag*by_ - ay_*bg)/det;
            gy = (ax_*bg - ag*bx_)/det;
            }

            else
            if(na>0)
            {
            gx = ag*ax_;
            gy = ag*ay_;
            }

            else
            if(nb>0)
            {
            gx = bg*bx_;
            gy = bg*by_;
            }

            gu = dt*g*gx;
            gv = dt*g*gy;

            }

        // a cell that is dry at the end of the step gets no gravity term (its eta is not a level)
        if(!wn(i,j))
        gu = gv = 0.0;

        // column: UH = A^-1 (I - dt L_v) UH* - dt g [H G] r,  r = A^-1 1
        col_nu(p,d,i,j,i,j,wnu.data());
        col_coef(p,d,WL0(i,j),i,j,i,j,wnu.data(),wcu.data(),wcl.data(),wcd.data());

        for(k=0;k<p->knoz;++k)
        wb[k] = 1.0;

        col_solve(p,wcu.data(),wcl.data(),wcd.data(),dt*gmc(i,j),wb.data(),wr.data());

        for(k=0;k<p->knoz;++k)
        wb[k] = UHD[IJK];

        col_solve(p,wcu.data(),wcl.data(),wcd.data(),dt*gmc(i,j),wb.data(),wx.data());

            for(k=0;k<p->knoz;++k)
            {
            d->UH[IJK] = wx[k] - gu*wr[k];
            d->VH[IJK] = 0.0;
            }

            if(p->j_dir==1)
            {
            for(k=0;k<p->knoz;++k)
            wb[k] = VHD[IJK];

            col_solve(p,wcu.data(),wcl.data(),wcd.data(),dt*gmc(i,j),wb.data(),wx.data());

            for(k=0;k<p->knoz;++k)
            d->VH[IJK] = wx[k] - gv*wr[k];
            }

        // WH: the vertical advection alone
        for(k=0;k<p->knoz;++k)
        {
        wnu[k] = 0.0;
        wb[k] = WHp[IJK];
        }

        col_coef(p,d,WL0(i,j),i,j,i,j,wnu.data(),wcu.data(),wcl.data(),wcd.data());
        col_solve(p,wcu.data(),wcl.data(),wcd.data(),0.0,wb.data(),wx.data());

        for(k=0;k<p->knoz;++k)
        d->WH[IJK] = wx[k];
        }

        else
        for(k=0;k<p->knoz;++k)
        {
        d->UH[IJK] = 0.0;
        d->VH[IJK] = 0.0;
        d->WH[IJK] = 0.0;
        }
    }
}

void nhflow_momentum_SI::relax(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // relaxation zones (B 98 2, B 99 1/2, B 520 1/2), once per step: WL, then U, V, W with the
    // new water depth; the targets of the end of the step with B 95 (the zone weight acts once
    // per step instead of once per RK stage)
    vrelax = 0.0;

    if(p->B98!=2 && p->B99!=1 && p->B99!=2 && p->B90==0)
    return;

    const int N = p->imax*p->jmax*(p->kmax+2);

    if(p->B95!=0.0)
    pflow->wavegen_stage_nhflow(p,d,pgc,p->simtime+p->dt);

    memcpy(dUH,d->UH,N*sizeof(double));
    memcpy(dVH,d->VH,N*sizeof(double));

    SLICELOOP4
    tS(i,j) = d->WL(i,j);

    pflow->WL_relax(p,pgc,d->WL,d->depth);

    SLICELOOP4
    {
    d->WL(i,j) = MAX(d->WL(i,j),0.0);
    d->eta(i,j) = d->WL(i,j) - d->depth(i,j);
    }

    // U, V, W on the new water depth
    LOOP
    {
        if(p->wet[IJ]==1 && d->WL(i,j)>p->A544)
        {
        d->U[IJK] = d->UH[IJK]/d->WL(i,j);
        d->V[IJK] = d->VH[IJK]/d->WL(i,j);
        d->W[IJK] = d->WH[IJK]/d->WL(i,j);
        }
    }

    pflow->U_relax(p,pgc,d->U,d->UH);
    pflow->V_relax(p,pgc,d->V,d->VH);
    pflow->W_relax(p,pgc,d->W,d->WH);

    SLICELOOP4
    vrelax += (d->WL(i,j) - tS(i,j))*carea(p,i,j);

    // the layer face discharges follow the change of the cell discharges
    LOOP
    {
    dUH[IJK] = d->UH[IJK] - dUH[IJK];
    dVH[IJK] = (d->VH[IJK] - dVH[IJK])*p->y_dir;
    }

    exchangeC(p,pgc,dUH);
    exchangeC(p,pgc,dVH);

    double nx,ny;

    SLICELOOP1
    if(Hx(i,j)>0.0)
    {
    fnorm(p,1,i,j,nx,ny);
    qx(i,j) = 0.0;

        for(k=0;k<p->knoz;++k)
        {
        Qx[IJK] += 0.5*(dUH[IJK] + dUH[Ip1JK])*nx + 0.5*(dVH[IJK] + dVH[Ip1JK])*ny;
        qx(i,j) += p->DZN[KP]*Qx[IJK];
        }
    }

    if(p->j_dir==1)
    SLICELOOP2
    if(Hy(i,j)>0.0)
    {
    fnorm(p,2,i,j,nx,ny);
    qy(i,j) = 0.0;

        for(k=0;k<p->knoz;++k)
        {
        Qy[IJK] += 0.5*(dUH[IJK] + dUH[IJp1K])*nx + 0.5*(dVH[IJK] + dVH[IJp1K])*ny;
        qy(i,j) += p->DZN[KP]*Qy[IJK];
        }
    }

    exchange1(p,pgc,qx);
    exchange2(p,pgc,qy);
    exchangeF(p,pgc,Qx);
    exchangeF(p,pgc,Qy);

    exchange4(p,pgc,d->WL);
    exchange4(p,pgc,d->eta);

    vrelax = pgc->globalsum(vrelax);
}

void nhflow_momentum_SI::finish(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const int N = p->imax*p->jmax*(p->kmax+2);

    // velocities (Froude limiter A 531), ghost cells
    velcalc(p,d,pgc,d->UH,d->VH,d->WH,d->WL,1.0);

    pgc->start4V(p,d->UH,gcval_uh);
    pgc->start4V(p,d->VH,gcval_vh);
    pgc->start4V(p,d->WH,gcval_wh);

    // the layer transport of the step is the continuity flux: omega, and the advection of the next step
    memcpy(d->FEx,Tx,N*sizeof(double));
    memcpy(d->FEy,Ty,N*sizeof(double));

    omega_update(p,d,pgc,d->WL,d->U,d->V,d->W);

    // turbulence (A 560 > 0): metrics at the new water level for the k, eps/omega step
    if(p->A560>0)
    sigma_update(p,d,pgc,d->WL);
}

void nhflow_momentum_SI::log(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    double vol=0.0;

    SLICELOOP4
    vol += d->WL(i,j)*carea(p,i,j);

    vol = pgc->globalsum(vol);

    write_log(p,vol,Qin,Qout,vrelax,psi_iter,wd_iter,wd_cells);
}

void nhflow_momentum_SI::write_log(lexer *p, double vol, double qi, double qo, double vr, int it, int wdit, int dry)
{
    if(p->mpirank==0 && !patch)
    silog<<p->count<<" "<<setprecision(10)<<p->simtime+p->dt<<" "<<p->dt<<" "<<setprecision(16)<<vol<<" "
         <<qi<<" "<<qo<<" "<<vr<<" "<<it<<" "<<wdit<<" "<<dry<<"\n"<<flush;
}

double nhflow_momentum_SI::div(lexer *p, slice &fx, slice &fy)
{
    // divergence of column face discharges (per unit face length, along the face normals)
    if(p->geo_curv==1)
    return (fx(i,j)*flen(p,1,i,j) - fx(i-1,j)*flen(p,1,i-1,j) + fy(i,j)*flen(p,2,i,j) - fy(i,j-1)*flen(p,2,i,j-1))/carea(p,i,j);

    double s = (fx(i,j) - fx(i-1,j))/p->DXN[IP];

    if(p->j_dir==1)
    s += (fy(i,j) - fy(i,j-1))/p->DYN[JP];

    return s;
}

double nhflow_momentum_SI::flen(lexer *p, int dir, int ii, int jj)
{
    // length of face i+1/2 (dir 1) or j+1/2 (dir 2) of cell (ii,jj)
    if(p->geo_curv==1)
    return dir==1 ? p->LEN1[(ii-p->imin)*p->jmax + (jj-p->jmin)] : p->LEN2[(ii-p->imin)*p->jmax + (jj-p->jmin)];

    if(dir==1)
    return p->j_dir==1 ? p->DYN[jj+marge] : 1.0;

    return p->DXN[ii+marge];
}

double nhflow_momentum_SI::fdn(lexer *p, int dir, int ii, int jj)
{
    // distance of the two cell centres of the face along its normal
    if(p->geo_curv==1)
    {
    const double dn = sflow_curv::dnorm(p,dir,ii,jj);
    return dn>1.0e-12 ? dn : 1.0e-12;
    }

    return dir==1 ? p->DXP[ii+marge] : p->DYP[jj+marge];
}

double nhflow_momentum_SI::carea(lexer *p, int ii, int jj)
{
    if(p->geo_curv==1)
    return p->AREA2D[(ii-p->imin)*p->jmax + (jj-p->jmin)];

    return p->DXN[ii+marge]*(p->j_dir==1 ? p->DYN[jj+marge] : 1.0);
}

void nhflow_momentum_SI::fnorm(lexer *p, int dir, int ii, int jj, double &nx, double &ny)
{
    // unit normal of the face towards +i (dir 1) / +j (dir 2)
    if(p->geo_curv==1)
    {
    double len;
    nhflow_curv::normal(p,dir,ii,jj,nx,ny,len);
    return;
    }

    nx = dir==1 ? 1.0 : 0.0;
    ny = dir==1 ? 0.0 : 1.0;
}

void nhflow_momentum_SI::cgrad(lexer *p, slice &f, int ii, int jj, double &fx, double &fy)
{
    // cell gradient of eta from the centres of the neighbours (sflow_curv::grad5), from the wet
    // neighbours only: one-sided at a wall or a dry cell (the eta of a dry cell is its bed level),
    // without the ghost cells behind a wall (their values and centres are not the same on every block
    // of a multi-block grid); no wet neighbour along an index direction: no gradient along it
    auto lij = [&](int a, int b) { return (a-p->imin)*p->jmax + (b-p->jmin); };
    auto in  = [&](int a, int b) { return p->flagslice4[lij(a,b)]>0; };
    auto wt  = [&](int a, int b) { return in(a,b) && p->wet[lij(a,b)]==1; };
    auto xc  = [&](int a, int b) { return p->XC2D[lij(a,b)]; };
    auto yc  = [&](int a, int b) { return p->YC2D[lij(a,b)]; };

    // values: the wet neighbours; geometry of a direction without them: the fluid neighbours
    const int iE = wt(ii+1,jj) ? ii+1 : ii, iW = wt(ii-1,jj) ? ii-1 : ii;
    const int jN = wt(ii,jj+1) ? jj+1 : jj, jS = wt(ii,jj-1) ? jj-1 : jj;
    const int gE = in(ii+1,jj) ? ii+1 : ii, gW = in(ii-1,jj) ? ii-1 : ii;
    const int gN = in(ii,jj+1) ? jj+1 : jj, gS = in(ii,jj-1) ? jj-1 : jj;

    if(gE==gW || gN==gS)
    {
    sflow_curv::grad(p,f,ii,jj,fx,fy);
    return;
    }

    double xe,ye,fe,xn,yn,fn;

    if(iE!=iW)
    {
    xe = xc(iE,jj) - xc(iW,jj);
    ye = yc(iE,jj) - yc(iW,jj);
    fe = f(iE,jj) - f(iW,jj);
    }

    else
    {
    xe = xc(gE,jj) - xc(gW,jj);
    ye = yc(gE,jj) - yc(gW,jj);
    fe = 0.0;
    }

    if(jN!=jS)
    {
    xn = xc(ii,jN) - xc(ii,jS);
    yn = yc(ii,jN) - yc(ii,jS);
    fn = f(ii,jN) - f(ii,jS);
    }

    else
    {
    xn = xc(ii,gN) - xc(ii,gS);
    yn = yc(ii,gN) - yc(ii,gS);
    fn = 0.0;
    }

    const double det = xe*yn - ye*xn;

    if(fabs(det)>0.0)
    {
    fx = (fe*yn - ye*fn)/det;
    fy = (xe*fn - fe*xn)/det;
    }

    else
    fx = fy = 0.0;
}

double nhflow_momentum_SI::fcorr(lexer *p, slice &f, int dir, int ii, int jj)
{
    // curvilinear grid: the non-orthogonal part of the normal gradient of the face, k.grad(f)/L with
    // k = S - d |S|^2/(S.d) (sflow_curv::kvec) and the mean of the cell gradients; the orthogonal part
    // (f_R - f_L)/dn is implicit in the eta system, this part explicit (deferred correction)
    if(p->geo_curv!=1)
    return 0.0;

    if(p->A539!=1)
    return 0.0;

    double kx,ky,fx1,fy1,fx2,fy2;
    sflow_curv::kvec(p,dir,ii,jj,kx,ky);
    cgrad(p,f,ii,jj,fx1,fy1);

    if(dir==1)
    cgrad(p,f,ii+1,jj,fx2,fy2);

    else
    cgrad(p,f,ii,jj+1,fx2,fy2);

    const double l = flen(p,dir,ii,jj);

    return l>0.0 ? 0.5*(kx*(fx1+fx2) + ky*(fy1+fy2))/l : 0.0;
}

void nhflow_momentum_SI::exchange1(lexer *p, ghostcell *pgc, slice &f)
{
    // MPI only: the first face of a partition belongs to the neighbour, physical
    // boundary faces are set by bc_setup
    if(p->mpi_size>1)
    {
    pgc->gcslparax(p,f,1);
    pgc->gcslparacox(p,f,14);
    }
}

void nhflow_momentum_SI::exchange2(lexer *p, ghostcell *pgc, slice &f)
{
    if(p->mpi_size>1 && p->j_dir==1)
    {
    pgc->gcslparax(p,f,2);
    pgc->gcslparacox(p,f,14);
    }
}

void nhflow_momentum_SI::exchange4(lexer *p, ghostcell *pgc, slice &f)
{
    if(p->mpi_size>1)
    {
    pgc->gcslparax(p,f,4);
    pgc->gcslparacox(p,f,14);
    }
}

void nhflow_momentum_SI::exchangeF(lexer *p, ghostcell *pgc, double *f)
{
    // layer face arrays (x- and y-faces), MPI only
    if(p->mpi_size>1)
    {
    pgc->gcparaxV1(p,f,14);
    pgc->gcparacoxV1(p,f,14);
    }
}

void nhflow_momentum_SI::exchangeC(lexer *p, ghostcell *pgc, double *f)
{
    // layer cell arrays, MPI only
    if(p->mpi_size>1)
    pgc->start4V_par(p,f,14);
}
