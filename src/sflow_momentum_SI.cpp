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

#include"sflow_momentum_SI.h"
#include"si_eta2D.h"
#include"bc_zone.h"
#include"geo_mesh.h"
#include"lexer.h"
#include"fdm2D.h"
#include"met_forcing.h"
#include"seastate_sflow.h"
#include"hydro_sflow.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"reefckpt.h"
#include"sflow_fsf.h"
#include"sflow_reconstruct.h"
#include"sflow_diffusion.h"
#include"sflow_rheology.h"
#include"sflow_curv.h"
#include"sflow_amr.h"
#include"6DOF.h"
#include"sflow_forcing.h"
#include<sys/stat.h>
#include<sys/types.h>

sflow_momentum_SI::sflow_momentum_SI(lexer *p, fdm2D *b, ghostcell *pgc, sflow_HLL *pphll, sflow_signal_speed *ppss,
                                     sflow_reconstruct *pprecon, sflow_diffusion *ppdiff, sflow_pressure *pppress,
                                     solver2D *ppsolv, solver2D *pppoissonsolv, ioflow *ppflow, sflow_fsf *ppfsf,
                                     sflow_forcing *ppsfdf, sixdof *pp6dof, bool block)
                                    : sflow_momentum_func(p,b,pgc,pphll,ppss,pprecon,ppdiff,pppress,ppsolv,pppoissonsolv,ppflow,ppfsf,ppsfdf,pp6dof),
                                      qx(p),qxn(p),qxE(p),Hx(p),alx(p),ax(p),kx(p),gx0(p),
                                      qy(p),qyn(p),qyE(p),Hy(p),aly(p),ay(p),ky(p),gy0(p),
                                      fxU(p),fxV(p),fyU(p),fyV(p),
                                      eta0(p),WL0(p),UH0(p),VH0(p),UHa(p),VHa(p),UHb(p),VHb(p),UHp(p),VHp(p),
                                      WHz(p),dUH(p),dVH(p),etaI(p),alc(p),rhs(p),rhsW(p),bdiag(p),wsrc(p),hsrc(p),
                                      tU(p),tV(p),tUH(p),tVH(p),WLk(p),etaB(p),
                                      qxN(p),Hxo(p),qxo(p),qxa(p),nxR(p),nxL(p),nxB(p),nx0(p),kxEL(p),kxER(p),kxPL(p),kxPR(p),kx1(p),kx2(p),
                                      qyN(p),Hyo(p),qyo(p),qya(p),nyR(p),nyL(p),nyB(p),ny0(p),kyEL(p),kyER(p),kyPL(p),kyPR(p),ky1(p),ky2(p),
                                      phiN(p),WHp(p),Lc0(p),Unm1(p),Vnm1(p),WLh(p),WLo(p),AWo(p),WLq(p),
                                      nPW(p),nEb(p),nEo(p),nPo(p),nU(p),nV(p),nUs(p),nVs(p),nA(p),nsx(p),nsy(p),nE(p),nP(p),nL(p),nZ(p),nAE(p),nAP(p),nBE(p),nBP(p),nGx(p),nGy(p),
                                      nF1x(p),nF2x(p),nF1y(p),nF2y(p),nH3(p),nW2(p),nCs(p),nRm(p),nRc(p),nZE(p),nZP(p),sps(p),bcA(p),bcS(p),FSx(p),FSy(p),dNx(p),dNy(p),mfx(p),mfy(p),acU(p),acV(p)
{
    check(p);

    psi = new si_eta2D(p);

    theta = p->A227;
    g = fabs(p->W22);
    rk = (p->A210==2)?2:3;
    qini = false;

    // A 210 11 (inside A210 = 1): one stage, the BCG advection with its source term (second order in time also for the
    // nonlinear terms: the RK stages of A 210 12/13 see neither the gravity nor the non-hydrostatic pressure of the step)
    if(p->A210==1)
    {
    rk = 1;
    onestep = true;
    }

    wd_iter = 0;
    wd_cells = 0;
    Qin = Qout = 0.0;

    bc_method(p);

    nh_ini(p);

    // the time levels of the nonlinear terms (nh_levels), every predictor: the depths of the gravity and A 220 3 terms
    // at t^{n+1/2}, the advecting discharge at t^{n+1/2}; with the SSP-RK predictors of A 210 12/13 (before only with
    // the non-hydrostatic pressure) also the acceleration of the step before in the later stages (stage_vel). With t^n
    // and without that acceleration the nonlinear terms of A 210 12/13 were first order in time (hydrostatic and
    // non-hydrostatic, SFLOW_NH_SI validation 01)
    // (REEF3D::HYDRO, rain and point sources: A 210 12/13 as before, the time levels only with the non-hydrostatic
    // pressure: on the rain-fed sheet flow of a tilted plane the extrapolated levels gave smaller steps and a 10 times
    // larger volume residual, 1e-6 m3)
    tl12 = !onestep && !hydro_sflow::wetting(p);
    tlev = (nh || onestep || tl12);

    patch = block;

    // X 10 3: the ship pressure of the step is the mean of its start and its end (ship_level; a refined patch: its
    // sflow_amr_ship); X 10 2: the body moves once per step (sixdof_sflow::start_oneway)
    if(p->X10==2 || p->X10==3)
    p6dof->sflow_midpoint(true);

    if(p->mpirank==0 && !patch)
    {
    if(onestep)
    cout<<"SFLOW semi-implicit free surface (A 210 11): theta "<<theta<<", one-step predictor (BCG advection)";

    else
    cout<<"SFLOW semi-implicit free surface (A 210 "<<10+rk<<"): theta "<<theta<<", SSP-RK"<<rk<<" predictor";

    if(p->A229==1)
    cout<<", wetting and drying iteration";

    if(p->A231>0.0)
    cout<<", mass flux dissipation "<<p->A231;

    cout<<", thin-water velocities below "<<p->A232<<" m (A 232)";

    if(nh)
    cout<<", non-hydrostatic quadratic pressure (A 220 "<<p->A220<<") coupled implicitly, alpha "<<adisp;

    cout<<endl;

    mkdir("./REEF3D_SFLOW_Log",0777);
    if(reefckpt_open(silog,"./REEF3D_SFLOW_Log/REEF3D_SFLOW_SI.dat",p->I45))   // I 45: continued
    {
        silog<<"# semi-implicit SFLOW (A 210 11/12/13): volume and open boundary discharges"<<endl;
        silog<<"# count  simtime  dt  volume  Qin  Qout  relax_volume  eta_iterations  wetdry_iterations  dry_cells"<<endl;
    }
    }
}

sflow_momentum_SI::~sflow_momentum_SI()
{
    delete psi;
}

void sflow_momentum_SI::check(lexer *p)
{
    // the options the semi-implicit path does not have yet
    int stop=0;

    auto msg = [&](const char *s)
    {
        if(p->mpirank==0)
        cout<<"!!! SFLOW A 210 11/12/13 (semi-implicit): "<<s<<" not available yet !!!"<<endl;
        stop=1;
    };

    // A 220 defaults to 2: the quadratic non-hydrostatic pressure (A 220 2/3) is coupled into the eta system, on one
    // Cartesian grid, without the wetting and drying iteration
    if(p->A220!=0 && p->A220!=2 && p->A220!=3)
    msg("Boussinesq and the linear non-hydrostatic pressure (A 220 1, 4) are; set A 220 0, 2 or 3,");

    // curvilinear grid (G 50): one block without mesh refinement (the composite coupled system of sflow_amr_sinh.cpp
    // is Cartesian)
    if(p->A220==2 || p->A220==3)
    {
        if(p->geo_mb>0)
        msg("the non-hydrostatic pressure (A 220 2/3) on multi-block grids is; set A 220 0,");

        else
        if(p->geo_curv==1 && p->G1>0)
        msg("the non-hydrostatic pressure (A 220 2/3) on curvilinear grids with mesh refinement (G 1) is; set A 220 0,");

        if(p->geo_curv==1 && p->j_dir==0)
        msg("the non-hydrostatic pressure (A 220 2/3) on a curvilinear grid with one cell in y is");
    }

    // mesh refinement (G 1) and multi-block grids (G 50 1): one eta system over the leaf cells of all grids
    // (reefamr_si.cpp), one global time step
    if(p->G1>0 && p->G7==1)
    msg("subcycling of the mesh refinement (G 7 1) is");


    if(p->B411>0 || p->B412>0 || p->B413>0 || p->B414>0 || p->B415>0 || p->B416>0 || p->B417>0 || p->B418>0 || p->B421>0 || p->B422>0)
    msg("patch boundary conditions (B 411-422) are");

    // X 10 3: ship waves from the pressure of the moving body (sixdof_sflow, X 400), with the gravity term; X 10 2: the
    // direct forcing of the moving body on the cells after the step (cell_update), one grid, Cartesian
    if(p->X10>0 && p->X10!=2 && p->X10!=3)
    msg("6DOF floating bodies (X 10 1; SFLOW has the moving body X 10 2 and the ship pressure X 10 3) are");

    if(p->X10==2 && (p->G1>0 || p->G50>0 || p->geo_curv==1))
    msg("the direct forcing of a moving body (X 10 2) with mesh refinement, curvilinear and multi-block grids is");

    // X 10 3 with mesh refinement (G 1, Cartesian): the pressure of the step on the patches (sflow_amr::si_step)
    if(p->X10==3 && (p->G50>0 || p->geo_curv==1))
    msg("the ship pressure (X 10 3) on curvilinear and multi-block grids is");

    // SEASTATE coupling (A 750 1): the wave force in the predictor, the Stokes mass source (A 751 2) in the eta
    // system; surfbeat (A 770 1) needs the host's step for the wave groups and the long-wave boundary
    if(p->A750>0 && p->A770==1)
    msg("SEASTATE surfbeat (A 770 1) is");

    if(stop==1)
    exit(1);

    if(p->A227<0.5 && p->mpirank==0)
    cout<<"SFLOW A 227 "<<p->A227<<": theta < 0.5 is unstable for large time steps"<<endl;
}

void sflow_momentum_SI::bc_method(lexer *p)
{
    // iowave open edge method per side (B 520 id method prio, B 521 id edge ...): 3 Riemann,
    // 4 Flather, 5 clamped level, 6 clamped discharge; B 521 edge 1 x-, 2 x+, 3 y-, 4 y+
    for(int k=0;k<5;++k)
    method[k]=0;

    const int edge_cs[5] = {0,1,4,3,2};

    for(int k=0;k<p->B521;++k)
    for(int qn=0;qn<p->B520;++qn)
    if(p->B520_id[qn]==p->B521_id[k] && p->B520_method[qn]>=3 && p->B520_method[qn]<=6
       && bc_zone_set::zone_block(p,p->B520_id[qn])==p->amr_block)
    {
        // a zone of a block (B 522): its face after turning the block (geo_mesh::turn_blocks)
        int e = p->B521_edge[k];
        if(p->amr_block>0 && p->gridgeo!=nullptr && e>=1 && e<=4)
        e = p->gridgeo->face_turned(p->amr_block,e-1) + 1;

        if(e<1 || e>4)
        continue;

        const int cs = edge_cs[e];

        if(method[cs]!=0 && method[cs]!=p->B520_method[qn])
        {
            if(p->mpirank==0)
            cout<<"!!! SFLOW A 210 11/12/13 (semi-implicit): different B 520 methods on one edge (B 521 edge "<<e<<") are not available yet !!!"<<endl;
            exit(1);
        }

        method[cs] = p->B520_method[qn];
    }
}

void sflow_momentum_SI::start(lexer *p, fdm2D* b, ghostcell* pgc)
{
    // multi-block grid (G 50 1) and mesh refinement (G 1): the step of all grids with one eta system (sflow_amr::si_step)
    if(prun!=nullptr)
    {
    prun->si_step(p,b,pgc);
    return;
    }

    step_single(p,b,pgc);
}

void sflow_momentum_SI::step_single(lexer *p, fdm2D* b, ghostcell* pgc)
{
    si_restore(p,b,pgc);

    // X 10 3: the body moves to the end of the step, its pressure of the start is kept (sixdof_sflow::sflow_midpoint);
    // the mean of the two goes with the gravity term (ship_level, si_assemble, cell_update)
    if(p->X10==3)
    ship_begin(p,b,pgc,true);

    si_begin(p,b,pgc);

    // time levels of the nonlinear terms (every predictor)
    if(tlev)
    nh_levels(p,b,pgc);

    // 3. predictor of the explicit terms
    for(int s=0; s<rk; ++s)
    {
    si_stage_vel(p,b,pgc,s);
    si_stage_rhs(p,b,pgc,s);
    }

    // 4.-7. face predictor, eta system, face discharges
    si_predictor_end(p,b,pgc);
    si_assemble(p,b,pgc);
    si_rhs(p,b,pgc);

    // non-hydrostatic pressure (A 220 2/3): eta and P from one coupled system
    if(nh)
    {
    nh_predictor(p,b,pgc);
    nh_solve(p,b,pgc);
    wd_iter = nh_wd_iter;
    wd_cells = nh_wd_dry;
    }

    else
    wd_iter = psi->solve_wetdry(p,pgc,ppoissonsolv,ax,ay,rhsW,b->depth,bdiag,etaI,p->A229==1 ? 20 : 1,wd_cells);

    si_end_faces(p,b,pgc);
    si_end_cells(p,b,pgc);
}

void sflow_momentum_SI::ship_begin(lexer *p, fdm2D* b, ghostcell* pgc, bool move)
{
    // X 10 3: the body of level 0 moves to the end of the step (move), the pressure of the step on this grid (a refined
    // patch: from the body at the start and the end of the step, sflow_amr::si_step and sflow_amr_ship)
    if(move)
    psfdf->forcing(p,b,pgc,p6dof,p->count,1.0,b->UH,b->VH,b->WH,b->WL,true);

    ship_level(p,b,pgc);
}

void sflow_momentum_SI::si_restore(lexer *p, fdm2D* b, ghostcell* pgc)
{
    // the water column of the end of the previous step: sflow_eta::depth_update (main loop,
    // after the momentum step) sets eta = WL - depth for a changed bed, fixes eta of the
    // outflow cells (F 50 2/3) and passes eta through gcsl_start4, whose ghost cells reach the
    // fluid behind a solid one or two cells thick (gcsl_neumann, three cells deep); the
    // semi-implicit step keeps its volume and treats the water level as a boundary value
    if(qini)
    {
        SLICELOOP4
        {
        b->WL(i,j) = WLk(i,j);
        b->eta(i,j) = WLk(i,j) - b->depth(i,j);
        }

    exchange4(p,pgc,b->WL);
    exchange4(p,pgc,b->eta);
    }
}

void sflow_momentum_SI::si_begin(lexer *p, fdm2D* b, ghostcell* pgc)
{
    // ghost cells of the open boundaries (ioflow discharge, iowave open edges), as the
    // explicit integrators. B 98 4: the wave generation (iowave::wavegen2D) also sets eta of the first cell inside;
    // here the wave stays in the ghost cells (the target of the Flather condition) and the cell keeps its level
    std::vector<double> eta_in;

    if(p->B98==4 && p->B90>0)
    for(n=0;n<p->gcslin_count;n++)
    eta_in.push_back(b->eta(p->gcslin[n][0],p->gcslin[n][1]));

    inflow(p,b,pgc,pflow);

    if(p->B98==4 && p->B90>0)
    for(n=0;n<p->gcslin_count;n++)
    b->eta(p->gcslin[n][0],p->gcslin[n][1]) = eta_in[n];

    // state at t^n (one grid of several: all cells of the arrays, the cells around the grid filled)
    if(amrgrid)
    {
        for(i=p->imin; i<p->imin+p->imax; ++i)
        for(j=p->jmin; j<p->jmin+p->jmax; ++j)
        {
        eta0(i,j) = b->eta(i,j);
        WL0(i,j) = b->WL(i,j);
        UH0(i,j) = b->UH(i,j);
        VH0(i,j) = b->VH(i,j)*p->y_dir;
        WHz(i,j) = 0.0;
        }
    }

    else
    SLICELOOP4
    {
    eta0(i,j) = b->eta(i,j);
    WL0(i,j) = b->WL(i,j);
    UH0(i,j) = b->UH(i,j);
    VH0(i,j) = b->VH(i,j)*p->y_dir;
    WHz(i,j) = 0.0;
    }

    // 1. face depths and implicit bed friction
    face_depth(p,b,pgc);

    // first step: face discharges from the cell values (curvilinear grid: along the face normal)
    if(!qini)
    {
        double nx,ny;

        SLICELOOP1
        {
            if(p->geo_curv==1)
            {
            fnorm(p,1,i,j,nx,ny);
            qx(i,j) = Hx(i,j)>0.0 ? 0.5*(b->UH(i,j) + b->UH(i+1,j))*nx + 0.5*(b->VH(i,j) + b->VH(i+1,j))*ny : 0.0;
            }

            else
            qx(i,j) = Hx(i,j)>0.0 ? 0.5*(b->UH(i,j) + b->UH(i+1,j)) : 0.0;
        }

        if(p->j_dir==1)
        SLICELOOP2
        {
            if(p->geo_curv==1)
            {
            fnorm(p,2,i,j,nx,ny);
            qy(i,j) = Hy(i,j)>0.0 ? 0.5*(b->UH(i,j) + b->UH(i,j+1))*nx + 0.5*(b->VH(i,j) + b->VH(i,j+1))*ny : 0.0;
            }

            else
            qy(i,j) = Hy(i,j)>0.0 ? 0.5*(b->VH(i,j) + b->VH(i,j+1)) : 0.0;
        }

    exchange1(p,pgc,qx);
    exchange2(p,pgc,qy);

    // face depths upwind of the new discharges
    face_depth(p,b,pgc);
    }

    friction(p,b);

    // 2. open boundaries (the first step also sets their face discharges)
    bc_setup(p,b,pgc);

    // 3. the face depths of the gravity term (multilayer flux)
    mlflux_faces(p,b,pgc);

    qini = true;

    SLICELOOP1
    qxn(i,j) = qx(i,j);

    SLICELOOP2
    qyn(i,j) = qy(i,j);

    exchange1(p,pgc,qxn);
    exchange2(p,pgc,qyn);

    // boundary faces lie outside SLICELOOP1/2
    for(const bc_face &f : bcf)
    {
        if(f.cs==1 || f.cs==4)
        qxn(f.fi,f.fj) = qx(f.fi,f.fj);

        else
        qyn(f.fi,f.fj) = qy(f.fi,f.fj);
    }

    // an open boundary face closed for the step (dry cell, no depth at the face) carries no old discharge: it kept
    // the discharge of the last step it was open, which drained the cell and left uncounted
    for(const bc_face &f : bco)
    {
        if(f.cs==1 || f.cs==4)
        qxn(f.fi,f.fj) = 0.0;

        else
        qyn(f.fi,f.fj) = 0.0;
    }
}

void sflow_momentum_SI::stage_io_si(int s, slice *&UHs, slice *&VHs, slice *&UHo, slice *&VHo)
{
    // SSP-RK of the explicit terms (A 210 2: stages 0, 1; 3: 0, 1, 2), depth H^n and face discharges q^n frozen
    slice *in[3][2]  = {{&UH0,&VH0},{&UHa,&VHa},{&UHb,&VHb}};
    slice *out[3][2] = {{&UHa,&VHa},{&UHb,&VHb},{&UHp,&VHp}};
    const int so = (s==rk-1) ? 2 : s;

    UHs = in[s][0]; VHs = in[s][1];
    UHo = out[so][0]; VHo = out[so][1];
}

void sflow_momentum_SI::si_stage_vel(lexer *p, fdm2D* b, ghostcell* pgc, int s)
{
    slice *UHs,*VHs,*UHo,*VHo;
    stage_io_si(s,UHs,VHs,UHo,VHo);

    // the time of the stage: SSP-RK2 0, 1; SSP-RK3 0, 1, 1/2. Not with REEF3D::HYDRO (rain and point sources, the
    // depth of the higher bed): on the wetting banks of a steep valley with rain the acceleration of the films drove
    // the stages to the N 61 stop (sflow_2d_sed_steep_valley_s67; tl12)
    stc = (!tl12 || s==0) ? 0.0 : (s==1 ? 1.0 : 0.5);

    stage_vel(p,b,pgc,*UHs,*VHs);
    stc = 0.0;
}

void sflow_momentum_SI::si_stage_rhs(lexer *p, fdm2D* b, ghostcell* pgc, int s)
{
    slice *UHs,*VHs,*UHo,*VHo;
    stage_io_si(s,UHs,VHs,UHo,VHo);

    double a;
    if(rk==2)
    a = s==0 ? 0.0 : 0.5;
    else
    a = s==0 ? 0.0 : (s==1 ? 0.75 : 1.0/3.0);

    stage_rhs(p,b,pgc,*UHs,*VHs,*UHo,*VHo,a);
}

void sflow_momentum_SI::si_predictor_end(lexer *p, fdm2D* b, ghostcell* pgc)
{
    // predictor of the cells, to the faces
    exchange4(p,pgc,UHp);
    exchange4(p,pgc,VHp);
}

void sflow_momentum_SI::si_end_faces(lexer *p, fdm2D* b, ghostcell* pgc)
{
    face_update(p,b,pgc);
}

void sflow_momentum_SI::water_column(lexer *p, fdm2D* b)
{
    // REEF3D::HYDRO: evaporation, infiltration and abstraction after the step (hydro_sflow::sinks) are part
    // of the water column the next step starts from
    // (the losses only: with the mesh refinement the fills and the restriction after the step change b->WL of the
    // cells around a grid and of the covered cells, which the next step must not take)
    if(b->hydro!=nullptr)
    SLICELOOP4
    WLk(i,j) -= b->hydro->L(i,j);
}

void sflow_momentum_SI::si_end_cells(lexer *p, fdm2D* b, ghostcell* pgc)
{
    // 8. WL, eta, wetting and drying
    eta_update(p,b,pgc);

    // non-hydrostatic pressure: breaking persistence (A 248 1) with eta^{n+1} and eta^n, as the explicit stages
    if(nh && p->A248==1)
    pfsf->breaking_persist(p,b,pgc,b->eta,eta0,1.0);

    // 9.-10. cell momentum, relaxation zones, velocities, ghost cells
    cell_update(p,b,pgc);

    SLICELOOP4
    {
    b->eta_n(i,j) = b->eta(i,j);
    WLk(i,j) = b->WL(i,j);
    }

    pgc->gcsl_start4(p,b->eta_n,gcval_eta);

    p->poissoniter = psi->iter;
    p->poissontime = psi->time;

    // a multi-block grid and the mesh refinement write the log of all grids (sflow_amr::si_step)
    if(patch || prun!=nullptr)
    return;

    log(p,b,pgc);

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
    cout<<"SI eta: iter "<<psi->iter<<"  time "<<setprecision(3)<<psi->time;

    if(p->A229==1)
    cout<<"  wetdry iterations "<<wd_iter<<"  dry cells "<<wd_cells;

    cout<<endl;
    }
}

void sflow_momentum_SI::face_depth(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // total water depth of the face: still water depth of the face plus eta upwind of the
    // face discharge (Stelling and Duinmeijer 2003), max(eta) where it is zero; a face is
    // open between two fluid cells when this depth exceeds the wetting criterion A 244.
    // Next to a dry cell the still water depth of the face is the smaller one (the higher
    // bed): water below the bed of the dry cell does not flow into it.
    // REEF3D::HYDRO rain and point sources (R 10, R 40): the higher bed on every face. Sheet flow on a slope
    // is thinner than the bed step between two cells (1 cm on a slope of 0.02 with 2 m cells against 4 cm),
    // and the mean of the still water depths adds half the step to the face depth: the friction of the face
    // drops and the sheet drains several times too fast.
    double df,eu,h;
    const bool hbed = hydro_sflow::wetting(p);

    SLICELOOP1
    {
    Hx(i,j) = 0.0;

        if(p->flagslice4[IJ]>0 && p->flagslice4[Ip1J]>0)
        {
        if(p->wet[IJ]==1 && p->wet[Ip1J]==1 && !hbed)
        df = 0.5*(b->depth(i,j) + b->depth(i+1,j));

        else
        df = MIN(b->depth(i,j),b->depth(i+1,j));

        int su;
        if(qini && qx(i,j)>1.0e-12)
        su = 0;

        else
        if(qini && qx(i,j)<-1.0e-12)
        su = 1;

        else
        su = (b->eta(i+1,j)>b->eta(i,j)) ? 1 : 0;

        eu = (qini && fabs(qx(i,j))>1.0e-12) ? b->eta(i+su,j) : MAX(b->eta(i,j),b->eta(i+1,j));

        h = eu + df;

        Hx(i,j) = h>p->A244 ? h : 0.0;

        // a source cell with only the film of the wetting criterion gives no water (dry_source)
        if(!hbed && dry_source(p,b,i+su,j))
        Hx(i,j) = 0.0;
        }
    }

    if(p->j_dir==1)
    SLICELOOP2
    {
    Hy(i,j) = 0.0;

        if(p->flagslice4[IJ]>0 && p->flagslice4[IJp1]>0)
        {
        if(p->wet[IJ]==1 && p->wet[IJp1]==1 && !hbed)
        df = 0.5*(b->depth(i,j) + b->depth(i,j+1));

        else
        df = MIN(b->depth(i,j),b->depth(i,j+1));

        int su;
        if(qini && qy(i,j)>1.0e-12)
        su = 0;

        else
        if(qini && qy(i,j)<-1.0e-12)
        su = 1;

        else
        su = (b->eta(i,j+1)>b->eta(i,j)) ? 1 : 0;

        eu = (qini && fabs(qy(i,j))>1.0e-12) ? b->eta(i,j+su) : MAX(b->eta(i,j),b->eta(i,j+1));

        h = eu + df;

        Hy(i,j) = h>p->A244 ? h : 0.0;

        if(!hbed && dry_source(p,b,i,j+su))
        Hy(i,j) = 0.0;
        }
    }

    exchange1(p,pgc,Hx);
    exchange2(p,pgc,Hy);
}

// The source cell of a face (upwind of its discharge, without flow the cell with the higher surface) holds only the
// film of the wetting criterion: the face stays closed. Between two dry cells on a slope the face depth
// max(eta) + min(depth) equals A 244 and opened by round-off, and the film (1 mm with A 244 1e-3) drained downhill
// with the bed slope as surface gradient; with the mesh refinement the films of the prolonged and restricted cells
// differ from A 244 by 1e-3 relative and opened many faces: sheets of a few mm ran down the dry canal banks behind a
// ship, up to 0.2 m above the water, at 0.5 m/s, until an eta solve set every cell dry (NOK canal, A 210 12, G 1 1).
// Water still flows out of a wet cell into a dry one (wetting) and out of a wet cell (drying).
bool sflow_momentum_SI::dry_source(lexer *p, fdm2D *b, int ii, int jj)
{
    return b->eta(ii,jj) + b->depth(ii,jj) <= 1.01*p->A244;
}

void sflow_momentum_SI::mlflux_faces(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // the face depth of the gravity term (multilayer flux, Popinet 2020; as NHFLOW v12), after face_depth, the friction
    // and the open boundaries: the mean depth of the two cells where both are wet and deep (WL above 10 A 244). Then
    // g H_f (eta_R - eta_L) is a difference of g h^2/2 over a flat bed and the momentum is conserved across a bore; with
    // the upwind face depth of face_depth the bore of a dam break runs ahead and the plateau behind it sinks. Shoreline
    // faces keep the depth of face_depth (wetting and drying), and so do the open boundary faces (bc_setup: the ghost
    // cells are not fluid cells) and REEF3D::HYDRO (the depth of the higher bed). mfx, mfy mark the faces (cell_update)
    const double hd = 10.0*p->A244;
    const bool hbed = hydro_sflow::wetting(p);

    // a fluid cell of the grid or of a neighbouring rank
    auto deep = [&](int a, int c)
    {
        const int q = (a-p->imin)*p->jmax + c-p->jmin;
        return p->flagslice4[q]>0 && p->wet[q]==1 && WL0(a,c)>hd;
    };

    SLICELOOP1
    {
    mfx(i,j) = 0.0;

        if(!hbed && Hx(i,j)>0.0 && deep(i,j) && deep(i+1,j))
        {
        Hx(i,j) = 0.5*(WL0(i,j) + WL0(i+1,j));
        mfx(i,j) = 1.0;
        }
    }

    if(p->j_dir==1)
    SLICELOOP2
    {
    mfy(i,j) = 0.0;

        if(!hbed && Hy(i,j)>0.0 && deep(i,j) && deep(i,j+1))
        {
        Hy(i,j) = 0.5*(WL0(i,j) + WL0(i,j+1));
        mfy(i,j) = 1.0;
        }
    }

    exchange1(p,pgc,Hx);
    exchange1(p,pgc,mfx);
    exchange2(p,pgc,Hy);
    exchange2(p,pgc,mfy);
}

void sflow_momentum_SI::face_weights(lexer *p, fdm2D *b, int dir, int ii, int jj, double &wl, double &wr)
{
    // weights of the two cell predictors in the face predictor of face ii+1/2 (dir 1) or jj+1/2 (dir 2), REEF3D::HYDRO:
    // w = the water of the cell above the higher bed of the face over the water depth of the cell (0..1), at t^n
    // (eta of t^n, its ghost cells exchanged); the cell with the higher bed has w = 1, a deep channel cell below the bed
    // of the bank w = 0, equal beds 1 and 1; the weights are normalised (wl + wr = 1: the steady sheet flow keeps the
    // friction balance), and the mean where neither cell holds water
    const int i2 = dir==1 ? ii+1 : ii;
    const int j2 = dir==1 ? jj : jj+1;
    const double df = MIN(b->depth(ii,jj),b->depth(i2,j2));

    auto w = [&](int i1, int j1)
    {
        const double h = b->eta(i1,j1) + b->depth(i1,j1);

        if(h<=p->A244)
        return 0.0;

        const double hc = b->eta(i1,j1) + df;
        return hc<=0.0 ? 0.0 : (hc>=h ? 1.0 : hc/h);
    };

    wl = w(ii,jj);
    wr = w(i2,j2);

    if(wl+wr>0.0)
    {
    const double s = wl+wr;
    wl /= s;
    wr /= s;
    }

    else
    {
    wl = 0.5;
    wr = 0.5;
    }
}

void sflow_momentum_SI::friction(lexer *p, fdm2D *b)
{
    // implicit Manning friction (A 218 1): d(UH)/dt = -gamma UH,  gamma = c_f |u|/H,
    // c_f = g n^2/H^(1/3), n = ks^(1/6)/20 (as sflow_rough_manning);  alpha = 1/(1 + dt gamma)
    double manning,cf,u,v,h;

    SLICELOOP1
    alx(i,j) = 1.0;

    SLICELOOP2
    aly(i,j) = 1.0;

    SLICELOOP4
    alc(i,j) = 1.0;

    if(p->A218!=1)
    return;

    SLICELOOP1
    if(Hx(i,j)>0.0)
    {
    h = Hx(i,j);
    manning = pow(0.5*(b->ks(i,j) + b->ks(i+1,j)),1.0/6.0)/20.0;
    cf = manning*manning*g/pow(h,1.0/3.0);
    u = qx(i,j)/h;
    v = 0.5*(b->V(i,j) + b->V(i+1,j))*p->y_dir;

    // curvilinear grid: the normal discharge and the tangential part of the mean cell velocity
    if(p->geo_curv==1)
    {
    double nx,ny;
    fnorm(p,1,i,j,nx,ny);
    v = -0.5*(b->U(i,j) + b->U(i+1,j))*ny + 0.5*(b->V(i,j) + b->V(i+1,j))*nx;
    }

    alx(i,j) = 1.0/(1.0 + p->dt*cf*sqrt(u*u + v*v)/h);
    }

    if(p->j_dir==1)
    SLICELOOP2
    if(Hy(i,j)>0.0)
    {
    h = Hy(i,j);
    manning = pow(0.5*(b->ks(i,j) + b->ks(i,j+1)),1.0/6.0)/20.0;
    cf = manning*manning*g/pow(h,1.0/3.0);
    u = 0.5*(b->U(i,j) + b->U(i,j+1));
    v = qy(i,j)/h;

    if(p->geo_curv==1)
    {
    double nx,ny;
    fnorm(p,2,i,j,nx,ny);
    u = -0.5*(b->U(i,j) + b->U(i,j+1))*ny + 0.5*(b->V(i,j) + b->V(i,j+1))*nx;
    }

    aly(i,j) = 1.0/(1.0 + p->dt*cf*sqrt(u*u + v*v)/h);
    }

    SLICELOOP4
    WETDRY
    if(b->WL(i,j)>p->A244)
    {
    h = b->WL(i,j);
    manning = pow(b->ks(i,j),1.0/6.0)/20.0;
    cf = manning*manning*g/pow(h,1.0/3.0);
    u = b->U(i,j);
    v = b->V(i,j)*p->y_dir;

    alc(i,j) = 1.0/(1.0 + p->dt*cf*sqrt(u*u + v*v)/h);
    }
}

void sflow_momentum_SI::bc_setup(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // open boundary faces of the domain edges, from the ghost cells set by ioflow / iowave
    //  discharge (inflow B 60, wave generation B 98 >= 3, clamped discharge edge B 520 6):
    //     the inward discharge of the ghost state, q_in = u_n,g h_g
    //  water level (outflow with F 50 2/3, clamped level edge B 520 5): Dirichlet eta_b
    //  Riemann edge (B 520 3), active absorption (B 99 >= 3), free outflow: the incoming
    //     invariant of the ghost state R = u_n,g + 2 sqrt(g h_g), linearised discharge
    //     q_in = h (R - 2c) + (R - 3c)(eta^{n+1} - eta^n), implicit in eta^{n+1}
    //  Flather edge (B 520 4): q_in = h (u_n,g + sqrt(g/h) eta^n) - sqrt(g h) eta^{n+1}
    //  the ghost state is that of the start of the step (one step behind the background)
    const double dt = p->dt;
    double manning,cf,u,h,eu,df;

    bcf.clear();
    bco.clear();

    // eta with the boundary values: water level outflow (ioflow), open edges (iowave)
    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    etaB(i,j) = b->eta(i,j);

    pflow->waterlevel2D(p,b,pgc,etaB);

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

        // face of the entry: reset (a boundary face is not in SLICELOOP1/2)
        if(f.cs==1) {f.fi=i-1; f.fj=j; f.gi=i-1; f.gj=j; f.s= 1.0; f.dxp=p->DXP[IM1]; f.dxn=p->DXN[IP]; f.ds=p->DYN[JP];}
        if(f.cs==4) {f.fi=i;   f.fj=j; f.gi=i+1; f.gj=j; f.s=-1.0; f.dxp=p->DXP[IP];  f.dxn=p->DXN[IP]; f.ds=p->DYN[JP];}
        if(f.cs==3) {f.fi=i; f.fj=j-1; f.gi=i; f.gj=j-1; f.s= 1.0; f.dxp=p->DYP[JM1]; f.dxn=p->DYN[JP]; f.ds=p->DXN[IP];}
        if(f.cs==2) {f.fi=i; f.fj=j;   f.gi=i; f.gj=j+1; f.s=-1.0; f.dxp=p->DYP[JP];  f.dxn=p->DYN[JP]; f.ds=p->DXN[IP];}

        const bool xf = (f.cs==1 || f.cs==4);

        // curvilinear grid: the normal distance of the centres, the face length, the cell width across it, the
        // unit normal towards +i / +j
        double fnx = xf ? 1.0 : 0.0, fny = xf ? 0.0 : 1.0;
        if(p->geo_curv==1)
        {
        const int dir = xf ? 1 : 2;
        f.dxp = fdn(p,dir,f.fi,f.fj);
        f.ds = flen(p,dir,f.fi,f.fj);
        f.dxn = carea(p,i,j)/MAX(f.ds,1.0e-20);
        fnorm(p,dir,f.fi,f.fj,fnx,fny);
        }
        f.nx = fnx;
        f.ny = fny;

        // the velocity normal to the face of a cell
        auto unorm = [&](int a, int c) { return p->geo_curv==1 ? b->U(a,c)*fnx + b->V(a,c)*fny : (xf ? b->U(a,c) : b->V(a,c)); };

        slice &Hf  = xf ? (slice&)Hx  : (slice&)Hy;
        slice &alf = xf ? (slice&)alx : (slice&)aly;
        slice &af  = xf ? (slice&)ax  : (slice&)ay;
        slice &kf  = xf ? (slice&)kx  : (slice&)ky;
        slice &qEf = xf ? (slice&)qxE : (slice&)qyE;
        slice &gf  = xf ? (slice&)gx0 : (slice&)gy0;
        slice &qf  = xf ? (slice&)qx  : (slice&)qy;

        Hf(f.fi,f.fj) = 0.0;
        alf(f.fi,f.fj) = 1.0;
        af(f.fi,f.fj) = 0.0;
        kf(f.fi,f.fj) = 0.0;
        qEf(f.fi,f.fj) = 0.0;
        gf(f.fi,f.fj) = 0.0;

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
        f.type = (p->B98==4 && p->B90>0) ? bc_flather : bc_flux;

        else
        if((bc==2 || bc==7 || bc==8) && outflow_open==1)
        f.type = ((p->F50==2 || p->F50==3) && p->B99<3) ? bc_level : bc_riemann;

        // B 98 4: generating-absorbing inflow (Flather with the wave of the ghost cell, sflow_momentum_SI::si_begin)
        const bool wgen = !(open_side[f.cs]==1 && method[f.cs]>0) && (bc==1 || bc==6) && inflow_open==1 && p->B98==4 && p->B90>0;

        // wall: the face stays closed
        if(f.type==0)
        {
        qf(f.fi,f.fj) = 0.0;
        continue;
        }

        // interior and ghost state of t^n
        const double hi = WL0(i,j);
        const double etai = eta0(i,j);
        const double ung = f.s*unorm(f.gi,f.gj);
        const double hg = MAX(b->eta(f.gi,f.gj) + b->depth(f.gi,f.gj), 0.0);
        const double uni = f.s*unorm(i,j);

        f.qin = 0.0;
        f.A = 0.0;
        f.B = 0.0;

        // a dry boundary cell: only an inflowing discharge. REEF3D::HYDRO (R 10, R 40): dry by the water column
        // of the step; the outflow cell takes the level of the outflow in sflow_eta::depth_update, and with a level at
        // or below the bed (a free outfall of rain-on-grid) its wet flag is that of the film, not of its water
        if((hydro_sflow::wetting(p) ? hi<=p->A244 : (p->wet[IJ]==0 || hi<=p->A244)) && f.type!=bc_flux)
        {
        qf(f.fi,f.fj) = 0.0;
        bco.push_back(f);
        continue;
        }

        if(f.type==bc_flux)
        {
        f.qin = ung*(hg>p->A244 ? hg : hi);

            if(!qini)
            qf(f.fi,f.fj) = f.s*f.qin;
        }

        if(f.type==bc_riemann)
        {
        const double ci = sqrt(g*hi);
        const double R = ung + 2.0*sqrt(g*MAX(hg,1.0e-6));

        f.B = MAX(3.0*ci - R, 0.0);
        f.A = hi*(R - 2.0*ci) + f.B*etai;

            if(!qini)
            qf(f.fi,f.fj) = f.s*hi*uni;
        }

        if(f.type==bc_flather && wgen)
        {
        // q_in = q_w + c (eta_w - eta^{n+1}): the wave (eta_w, u_w depth averaged) comes in, the reflected wave leaves
        const double etaw = b->eta(f.gi,f.gj);
        f.B = sqrt(g*hi);
        f.A = ung*(hg>p->A244 ? hg : hi) + f.B*etaw;

            if(!qini)
            qf(f.fi,f.fj) = f.s*hi*uni;
        }

        else
        if(f.type==bc_flather)
        {
        f.B = sqrt(g*hi);
        f.A = hi*(ung + sqrt(g/hi)*etai);

            if(!qini)
            qf(f.fi,f.fj) = f.s*hi*uni;
        }

        if(f.type==bc_level)
        {
        const double etab = etaB(f.gi,f.gj);

        df = 0.5*(b->depth(i,j) + b->depth(f.gi,f.gj));

        if(!qini)
        qf(f.fi,f.fj) = p->geo_curv==1 ? UH0(i,j)*fnx + VH0(i,j)*fny : (xf ? UH0(i,j) : VH0(i,j));

        // eta upwind of the face discharge (inward: the boundary value)
        const double qin_n = f.s*qf(f.fi,f.fj);

        if(qin_n>1.0e-12)
        eu = etab;

        else
        if(qin_n<-1.0e-12)
        eu = etai;

        else
        eu = MAX(etab,etai);

        h = eu + df;

            if(h<=p->A244)
            {
            qf(f.fi,f.fj) = 0.0;
            bco.push_back(f);
            continue;
            }

        Hf(f.fi,f.fj) = h;

            if(p->A218==1)
            {
            manning = pow(0.5*(b->ks(i,j) + b->ks(f.gi,f.gj)),1.0/6.0)/20.0;
            cf = manning*manning*g/pow(h,1.0/3.0);
            u = qf(f.fi,f.fj)/h;
            alf(f.fi,f.fj) = 1.0/(1.0 + dt*cf*fabs(u)/h);
            }

        // the boundary value for the face gradients of t^n and t^{n+1} (held over the step)
        etaI(f.gi,f.gj) = etab;

        // face gradient in the face direction: (right - left)/dxp
        gf(f.fi,f.fj) = f.s*(etai - etab)/f.dxp;
        }

    bcf.push_back(f);
    }
}

void sflow_momentum_SI::stage_vel(lexer *p, fdm2D *b, ghostcell *pgc, slice &UHs, slice &VHs)
{
    // stage velocities on the depth of t^n; below the depth A 232 desingularised (Kurganov and
    // Petrova 2007): u = sqrt(2) h (uh)/sqrt(h^4 + max(h^4, h_v^4)), continuous at h_v, u -> 0 with h.
    // A cell that is just wetting (h of a few mm) receives the momentum of its deep neighbours with
    // the inflow, but its depth is that of t^n: (uh)/h gave hundreds of m/s within the stages, carried
    // into the face predictor of the next faces (a blow-up of 57 m in the 1 km Michael hindcast)
    const double hv4 = p->A232*p->A232*p->A232*p->A232;

    SLICELOOP4
    {
        if(p->wet[IJ]==1 && WL0(i,j)>p->A244)
        {
            // a later stage of A 210 12/13 (stc > 0): the momentum with the acceleration of the implicit part of the
            // step before up to the time of the stage, on the depth of that time; thin water (the depth of the stage
            // or of t^n below A 232) as before
            if(stc>0.0 && WL0(i,j)>=p->A232 && (stc>0.75 ? WLq(i,j) : WLh(i,j))>=p->A232)
            {
            const double hs = stc>0.75 ? WLq(i,j) : WLh(i,j);
            b->U(i,j) = (UHs(i,j) + stc*p->dt*acU(i,j))/hs;
            b->V(i,j) = (VHs(i,j) + stc*p->dt*acV(i,j))/hs*p->y_dir;
            }

            else
            if(WL0(i,j)>=p->A232)
            {
            b->U(i,j) = UHs(i,j)/WL0(i,j);
            b->V(i,j) = VHs(i,j)/WL0(i,j)*p->y_dir;
            }

            else
            {
            const double h = WL0(i,j);
            const double r = sqrt(2.0)*h/sqrt(h*h*h*h + hv4);
            b->U(i,j) = r*UHs(i,j);
            b->V(i,j) = r*VHs(i,j)*p->y_dir;
            }
        }

        else
        {
        b->U(i,j) = 0.0;
        b->V(i,j) = 0.0;
        }

    b->W(i,j) = 0.0;
    }

    // MPI and boundary ghost cells of U, V and UH, VH
    // the interior is kept: the ghost cells of a wall reach three cells into a solid
    // (sflow_momentum_func::mirror, ::ghostcells), into the fluid behind a thin solid
    SLICELOOP4
    {
    tU(i,j) = b->U(i,j);
    tV(i,j) = b->V(i,j);
    tUH(i,j) = UHs(i,j);
    tVH(i,j) = VHs(i,j);
    }

    ghostcells(p,b,pgc,UHs,VHs,WHz,b->WL);

    SLICELOOP4
    {
    b->U(i,j) = tU(i,j);
    b->V(i,j) = tV(i,j);
    UHs(i,j) = tUH(i,j);
    VHs(i,j) = tVH(i,j);
    b->WL(i,j) = WL0(i,j);
    }
}

void sflow_momentum_SI::stage_rhs(lexer *p, fdm2D *b, ghostcell *pgc, slice &UHs, slice &VHs, slice &UHo, slice &VHo, double a)
{
    const double alpha = 1.0-a;

    clearrhs(p,b);

    // advection with the face discharges of t^n
    advection(p,b,pgc);

    // the sources of the stage (FSx, FSy) for the wall sides of cell_update (SSP-RK: those of the last stage)
    SLICELOOP4
    {
    FSx(i,j) = b->F(i,j);
    FSy(i,j) = b->G(i,j);
    }

    // sources and diffusion (as the explicit stage, without gravity and bed friction; the ship pressure X 10 3 with
    // the gravity term, si_assemble and cell_update)
    prheo->u_source(p,b,b->U,b->V,WL0);
    if(pmet!=nullptr)
    pmet->source_sflow(p,pgc,b->F,0,WL0,UHs,VHs);
    if(b->wave!=nullptr)
    b->wave->u_source(p,b);
    irhs(p,b);

    SLICELOOP4
    FSx(i,j) = b->F(i,j) - FSx(i,j);

	pdiff->diff_u(p,b,pgc,psolv,UHDIFF,UHs,b->U,b->V,WL0,alpha);

    if(p->j_dir==1)
    {
    prheo->v_source(p,b,b->U,b->V,WL0);
    if(pmet!=nullptr)
    pmet->source_sflow(p,pgc,b->G,1,WL0,UHs,VHs);
    if(b->wave!=nullptr)
    b->wave->v_source(p,b);
    jrhs(p,b);

        SLICELOOP4
        FSy(i,j) = b->G(i,j) - FSy(i,j);
    }

    else
    SLICELOOP4
    FSy(i,j) = 0.0;
    pdiff->diff_v(p,b,pgc,psolv,VHDIFF,VHs,b->U,b->V,WL0,alpha);

    SLICELOOP4
    {
        if(p->wet[IJ]==1)
        {
        UHo(i,j) = a*UH0(i,j) + alpha*(UHDIFF(i,j) + p->dt*b->F(i,j));
        VHo(i,j) = (a*VH0(i,j) + alpha*(VHDIFF(i,j) + p->dt*b->G(i,j)))*p->y_dir;
        }

        else
        {
        UHo(i,j) = 0.0;
        VHo(i,j) = 0.0;
        }
    }
}

void sflow_momentum_SI::advection(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // conservative advection of UH, VH with the face discharges of t^n (the time levels of nh_levels: extrapolated to
    // t^{n+1/2}), velocities upwind of the face discharge (reconstructed, A 211; A 210 11: the BCG face values); no
    // gravity wave speed in the dissipation: the Courant number of the flow alone limits the time step
    slice &qa = tlev ? (slice&)qxa : (slice&)qxn;
    slice &qb = tlev ? (slice&)qya : (slice&)qyn;

    if(onestep)
    advection_bcg(p,b,pgc,qa,qb);

    else
    {
        precon->reconstruct_x(p, pgc, b, b->U, b->Us, b->Un);

        if(p->j_dir==1)
        {
        precon->reconstruct_x(p, pgc, b, b->V, b->Vs, b->Vn);
        precon->reconstruct_y(p, pgc, b, b->U, b->Ue, b->Uw);
        precon->reconstruct_y(p, pgc, b, b->V, b->Ve, b->Vw);
        }

        SLICELOOP1
        {
        fxU(i,j) = qa(i,j)>=0.0 ? qa(i,j)*b->Us(i,j) : qa(i,j)*b->Un(i,j);
        fxV(i,j) = 0.0;

        if(p->j_dir==1)
        fxV(i,j) = qa(i,j)>=0.0 ? qa(i,j)*b->Vs(i,j) : qa(i,j)*b->Vn(i,j);
        }

        if(p->j_dir==1)
        SLICELOOP2
        {
        fyU(i,j) = qb(i,j)>=0.0 ? qb(i,j)*b->Ue(i,j) : qb(i,j)*b->Uw(i,j);
        fyV(i,j) = qb(i,j)>=0.0 ? qb(i,j)*b->Ve(i,j) : qb(i,j)*b->Vw(i,j);
        }
    }

    exchange1(p,pgc,fxU);
    exchange1(p,pgc,fxV);
    exchange2(p,pgc,fyU);
    exchange2(p,pgc,fyV);

    // open boundary faces: first-order upwind, inflow with the ghost velocities
    for(const bc_face &f : bcf)
    {
    i=f.i;
    j=f.j;

        const bool xf = (f.cs==1 || f.cs==4);
        const double q = xf ? qxn(f.fi,f.fj) : qyn(f.fi,f.fj);
        const bool in = f.s*q>0.0;
        const double uu = in ? b->U(f.gi,f.gj) : b->U(i,j);
        const double vv = (in ? b->V(f.gi,f.gj) : b->V(i,j))*p->y_dir;

        if(xf)
        {
        fxU(f.fi,f.fj) = q*uu;
        fxV(f.fi,f.fj) = q*vv;
        }

        else
        {
        fyU(f.fi,f.fj) = q*uu;
        fyV(f.fi,f.fj) = q*vv;
        }
    }

    // mesh refinement: the coarse faces next to a finer patch take the fine fluxes, a patch records its box faces
    if(fhook)
    fhook(p,fxU,fxV,fyU,fyV);

    if(p->geo_curv==1)
    {
        SLICELOOP4
        WETDRY
        {
        b->F(i,j) -= div(p,fxU,fyU);
        b->G(i,j) -= div(p,fxV,fyV);
        }
        return;
    }

    SLICELOOP4
    WETDRY
    {
    b->F(i,j) -= (fxU(i,j) - fxU(i-1,j))/p->DXN[IP];

        if(p->j_dir==1)
        {
        b->F(i,j) -= (fyU(i,j) - fyU(i,j-1))/p->DYN[JP];

        b->G(i,j) -= (fxV(i,j) - fxV(i-1,j))/p->DXN[IP]
                   + (fyV(i,j) - fyV(i,j-1))/p->DYN[JP];
        }
    }
}

void sflow_momentum_SI::wall_sides(lexer *p)
{
    // A 210 11: the sides of the cells with an open boundary face (iowave, ioflow; also those closed for the step), as
    // bits: 1 x-, 2 x+, 4 y-, 8 y+ (bc_face cs: 1 x-, 4 x+, 3 y-, 2 y+); wall() tells the other sides towards a solid
    // or the outside of the domain
    SLICELOOP4
    bcS(i,j) = 0.0;

    auto mark = [&](const bc_face &f)
    {
        const int bit = f.cs==1 ? 1 : (f.cs==4 ? 2 : (f.cs==3 ? 4 : 8));
        bcS(f.i,f.j) = double(int(bcS(f.i,f.j)) | bit);
    };

    for(const bc_face &f : bcf)
    mark(f);

    for(const bc_face &f : bco)
    mark(f);
}

double sflow_momentum_SI::wall(lexer *p, int ii, int jj, int side)
{
    // 1.0 if the side of cell (ii,jj) (1 x-, 2 x+, 3 y-, 4 y+) faces a solid or the outside of the domain and no open
    // boundary (wall_sides of the step)
    // A 210 11, and A 210 12/13 with the second-order time levels (tl12): the wall sides; 0 otherwise (as before)
    if(!(onestep || tl12))
    return 0.0;

    const int ai = side==1 ? ii-1 : (side==2 ? ii+1 : ii);
    const int aj = side==3 ? jj-1 : (side==4 ? jj+1 : jj);
    const int bit = side==1 ? 1 : (side==2 ? 2 : (side==3 ? 4 : 8));

    if(p->flagslice4[(ai-p->imin)*p->jmax + (aj-p->jmin)]>=0)
    return 0.0;

    return (int(bcS(ii,jj)) & bit) ? 0.0 : 1.0;
}

void sflow_momentum_SI::advection_bcg(lexer *p, fdm2D *b, ghostcell *pgc, slice &qa, slice &qb)
{
    // A 210 11: the fluxes q_f u_f of U and V with the face values of the BCG scheme (Bell, Colella & Glaz 1989), as
    // NHFLOW A 510 11 and the multilayer advection of Popinet (2020): the upwind cell extrapolated to the face and to
    // t^{n+1/2}, with the source term of the upwind cell,
    //   u_f = u_up + a (1 - c) g dx/2 - dt/2 v_up g_t + dt/2 S_up,    c = |u_f| dt/dx, a = sign(q_f),
    // q_f the discharge of the step at t^{n+1/2} (nh_levels), g the slope of the upwind cell normal to the face, g_t the
    // one along it (2D), both the generalised minmod (theta 1.3), v_up the velocity along the face in the
    // upwind cell (from its two faces). S = -g grad(eta^n) (with the ship pressure level of the step, X 10 3; a wall
    // side with the level mirrored) plus the non-hydrostatic increment of the cell momentum of the step before (dNx,
    // dNy of nh_cell) over H dt: without S the advected velocity misses the acceleration of the step and the nonlinear
    // terms are first order in time (standing wave a/h 0.094, A 210 12: order 1.2 at T/1280). Slopes next to solids,
    // dry cells and open boundaries: 0. Curvilinear grid (G 50): the cell velocities U, V are Cartesian, q_f the
    // normal discharge; the slopes along the index directions (per index), dx the index step 1, c with the centre
    // distance and the normal velocity of the face, g_t over the centre distance along the face, S by the chain rule
    const double dt = p->dt;
    const bool yd = (p->j_dir==1);
    const bool cv = (p->geo_curv==1);

    // the cells the slopes may use: wet fluid cells (MPI: from the neighbouring ranks, flagslice4 marks only the first
    // ghost column of a rank boundary)
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    {
    const int q = (ii-p->imin)*p->jmax + (jj-p->jmin);
    bcA(ii,jj) = (p->flagslice4[q]>0 && p->wet[q]==1 && WL0(ii,jj)>p->A244) ? 1.0 : 0.0;
    }

    exchange4(p,pgc,bcA);

    // the open sides of the cells (wall)
    wall_sides(p);

    // eta^n of the source term and the non-hydrostatic increment of the step before in the cells of the neighbouring
    // ranks (si_begin copies the interior)
    exchange4(p,pgc,eta0);

    if(nh)
    {
    exchange4(p,pgc,dNx);
    exchange4(p,pgc,dNy);
    }

    auto ok = [&](int a, int c) { return bcA(a,c)>0.5; };
    auto lev = [&](int a, int c) { return eta0(a,c) + (ship ? sps(a,c) : 0.0); };

    // the non-hydrostatic part of the source term of the cell for U (dir 1) or V (dir 2)
    auto snh = [&](int dir, int a, int c) -> double
    {
        if(!nh || WL0(a,c)<=p->A244)
        return 0.0;

        return (dir==1 ? dNx(a,c) : dNy(a,c))/(WL0(a,c)*MAX(p->dt_old,1.0e-20));
    };

    // the level of the neighbour of cell (a,c) on side sd (1 x-, 2 x+, 3 y-, 4 y+) for the gradient: its own if the
    // slopes may use it, mirrored at a wall, extrapolated through the cell from the other side towards a dry cell
    // (one-sided), the level of the cell if neither side is open
    auto nlev = [&](int a, int c, int sd) -> double
    {
        const int an = sd==1 ? a-1 : (sd==2 ? a+1 : a), cn = sd==3 ? c-1 : (sd==4 ? c+1 : c);
        const int ao = sd==1 ? a+1 : (sd==2 ? a-1 : a), co = sd==3 ? c+1 : (sd==4 ? c-1 : c);
        const int so = sd==1 ? 2 : (sd==2 ? 1 : (sd==3 ? 4 : 3));

        if(ok(an,cn))
        return lev(an,cn);

        if(wall(p,a,c,sd)>0.5)
        return lev(a,c);

        if(ok(ao,co) || wall(p,a,c,so)>0.5)
        return 2.0*lev(a,c) - (ok(ao,co) ? lev(ao,co) : lev(a,c));

        return lev(a,c);
    };

    // the source term of the cell for U (dir 1) or V (dir 2)
    auto src = [&](int dir, int a, int c) -> double
    {
        const double s = snh(dir,a,c);

        if(cv)
        {
        double gx,gy;
        sflow_curv::grad5(p,a,c,nlev(a,c,2),nlev(a,c,1),nlev(a,c,4),nlev(a,c,3),gx,gy);
        return s - g*(dir==1 ? gx : gy);
        }

        // a wall side counts with the mirrored level (as cell_update)
        const double w1 = wall(p,a,c,dir==1 ? 1 : 3), w2 = wall(p,a,c,dir==1 ? 2 : 4);

        if(dir==1)
        {
        const bool oe = ok(a+1,c), ow = ok(a-1,c), me = oe || w2>0.5, mw = ow || w1>0.5;
        const double le = oe ? lev(a+1,c) : lev(a,c), lw = ow ? lev(a-1,c) : lev(a,c);
        if(me && mw) return s - g*(le - lw)/(p->DXP[a-1+marge] + p->DXP[a+marge]);
        if(oe) return s - g*(le - lev(a,c))/p->DXP[a+marge];
        if(ow) return s - g*(lev(a,c) - lw)/p->DXP[a-1+marge];
        return s;
        }

        const bool on = ok(a,c+1), os = ok(a,c-1), mn = on || w2>0.5, ms = os || w1>0.5;
        const double ln = on ? lev(a,c+1) : lev(a,c), ls = os ? lev(a,c-1) : lev(a,c);
        if(mn && ms) return s - g*(ln - ls)/(p->DYP[c-1+marge] + p->DYP[c+marge]);
        if(on) return s - g*(ln - lev(a,c))/p->DYP[c+marge];
        if(os) return s - g*(lev(a,c) - ls)/p->DYP[c-1+marge];
        return s;
    };

    // the normal velocity of the face for its Courant number: the depth-weighted mean of the cell velocities (as
    // the face velocity of Popinet (2020); curvilinear: their component along the face normal)
    auto uface = [&](int dir, int a, int c)
    {
        const int a2 = dir==1 ? a+1 : a, c2 = dir==2 ? c+1 : c;
        const double hl = ok(a,c) ? WL0(a,c) : 0.0, hr = ok(a2,c2) ? WL0(a2,c2) : 0.0;

        if((hl+hr)<=0.0)
        return 0.0;

        if(cv)
        {
        double nx,ny;
        fnorm(p,dir,a,c,nx,ny);
        return (hl*(b->U(a,c)*nx + b->V(a,c)*ny) + hr*(b->U(a2,c2)*nx + b->V(a2,c2)*ny))/(hl + hr);
        }

        slice &u = dir==1 ? (slice&)b->U : (slice&)b->V;
        return (hl*u(a,c) + hr*u(a2,c2))/(hl + hr);
    };

    // the limited slope of a cell from the values and centre distances on its two sides (one-sided: 0)
    auto slope = [&](double fm, double f0, double fp, double dm, double dp, bool om, bool op)
    {
        if(!om || !op)
        return 0.0;

        const double sl = (f0 - fm)/dm, sr = (fp - f0)/dp, sc = (fp - fm)/(dm + dp);

        if(sl*sr<=0.0)
        return 0.0;

        const double th = 1.3;
        const double m = MIN(MIN(th*fabs(sl),th*fabs(sr)),fabs(sc));
        return sl>0.0 ? m : -m;
    };

    // x-faces (i+1/2)
    SLICELOOP1
    {
    fxU(i,j) = 0.0;
    fxV(i,j) = 0.0;

        if(Hx(i,j)>0.0)
        {
        const double t = qa(i,j);
        const int iu = t>=0.0 ? i : i+1;
        double fu = b->U(iu,j);
        double fv = (yd || cv) ? b->V(iu,j) : 0.0;

            if(ok(iu,j))
            {
            const double a = t>=0.0 ? 1.0 : -1.0;
            const double dx = cv ? 1.0 : p->DXN[iu+marge];
            const double c = MIN(fabs(uface(1,i,j))*dt/(cv ? fdn(p,1,i,j) : dx),1.0);
            const double dm = cv ? 1.0 : p->DXP[iu-1+marge], dp = cv ? 1.0 : p->DXP[iu+marge];
            const bool om = ok(iu-1,j), op = ok(iu+1,j);

            fu += a*(1.0 - c)*slope(b->U(iu-1,j),b->U(iu,j),b->U(iu+1,j),dm,dp,om,op)*0.5*dx + 0.5*dt*src(1,iu,j);

                if(yd || cv)
                fv += a*(1.0 - c)*slope(b->V(iu-1,j),b->V(iu,j),b->V(iu+1,j),dm,dp,om,op)*0.5*dx + 0.5*dt*src(2,iu,j);

                // transverse term: the velocity along the face in the upwind cell from its two y-faces
                if(yd)
                {
                const double hs = Hy(iu,j-1) + Hy(iu,j);

                    if(hs>p->A244)
                    {
                    const double v = (qb(iu,j-1) + qb(iu,j))/hs;
                    const double em = cv ? 1.0 : p->DYP[j-1+marge], ep = cv ? 1.0 : p->DYP[j+marge];
                    const double ds = cv ? 0.5*(fdn(p,2,iu,j-1) + fdn(p,2,iu,j)) : 1.0;
                    const bool on = ok(iu,j+1), os = ok(iu,j-1);
                    fu -= 0.5*dt*v*slope(b->U(iu,j-1),b->U(iu,j),b->U(iu,j+1),em,ep,os,on)/ds;
                    fv -= 0.5*dt*v*slope(b->V(iu,j-1),b->V(iu,j),b->V(iu,j+1),em,ep,os,on)/ds;
                    }
                }
            }

        fxU(i,j) = t*fu;
        fxV(i,j) = t*fv;
        }
    }

    // y-faces (j+1/2)
    if(yd)
    SLICELOOP2
    {
    fyU(i,j) = 0.0;
    fyV(i,j) = 0.0;

        if(Hy(i,j)>0.0)
        {
        const double t = qb(i,j);
        const int ju = t>=0.0 ? j : j+1;
        double fu = b->U(i,ju);
        double fv = b->V(i,ju);

            if(ok(i,ju))
            {
            const double a = t>=0.0 ? 1.0 : -1.0;
            const double dy = cv ? 1.0 : p->DYN[ju+marge];
            const double c = MIN(fabs(uface(2,i,j))*dt/(cv ? fdn(p,2,i,j) : dy),1.0);
            const double dm = cv ? 1.0 : p->DYP[ju-1+marge], dp = cv ? 1.0 : p->DYP[ju+marge];
            const bool om = ok(i,ju-1), op = ok(i,ju+1);

            fu += a*(1.0 - c)*slope(b->U(i,ju-1),b->U(i,ju),b->U(i,ju+1),dm,dp,om,op)*0.5*dy + 0.5*dt*src(1,i,ju);
            fv += a*(1.0 - c)*slope(b->V(i,ju-1),b->V(i,ju),b->V(i,ju+1),dm,dp,om,op)*0.5*dy + 0.5*dt*src(2,i,ju);

            // transverse term: the velocity along the face in the upwind cell from its two x-faces
            const double hs = Hx(i-1,ju) + Hx(i,ju);

                if(hs>p->A244)
                {
                const double u = (qa(i-1,ju) + qa(i,ju))/hs;
                const double em = cv ? 1.0 : p->DXP[i-1+marge], ep = cv ? 1.0 : p->DXP[i+marge];
                const double ds = cv ? 0.5*(fdn(p,1,i-1,ju) + fdn(p,1,i,ju)) : 1.0;
                const bool oe = ok(i+1,ju), ow = ok(i-1,ju);
                fu -= 0.5*dt*u*slope(b->U(i-1,ju),b->U(i,ju),b->U(i+1,ju),em,ep,ow,oe)/ds;
                fv -= 0.5*dt*u*slope(b->V(i-1,ju),b->V(i,ju),b->V(i+1,ju),em,ep,ow,oe)/ds;
                }
            }

        fyU(i,j) = t*fu;
        fyV(i,j) = t*fv;
        }
    }
}

double sflow_momentum_SI::dissipation(lexer *p, double q, double h, double dxp)
{
    // optional mass flux dissipation (A 231): K = A231 (|u| + c) dx/2, the Rusanov
    // dissipation of the continuity equation, switched on from Froude number 0.5 to 1;
    // implicit (in the eta system), for hydraulic jumps and transcritical flow
    if(p->A231<=0.0 || h<=p->A244)
    return 0.0;

    const double u = fabs(q)/h;
    const double c = sqrt(g*h);
    const double fr = u/c;
    const double phi = MAX(MIN((fr-0.5)/0.5,1.0),0.0);

    return p->A231*0.5*(u + c)*dxp*phi;
}

void sflow_momentum_SI::si_assemble(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const double dt = p->dt;
    double K;

    // face predictor (momentum interpolation: the mean of the cell predictors, which carry
    // no gravity term of this step), explicit part of the face discharge, eta-system
    // coefficients a_f = theta dt (g theta dt alpha_f H_f + K_f)/dxp
    // REEF3D::HYDRO rain and point sources (R 10, R 40, R 100 2): the face depth is that of the higher bed
    // (face_depth), and the mean is weighted by the water of each cell above that bed (face_weights): a deep
    // channel cell next to a bank whose bed lies above the channel level gave half its momentum to the thin
    // face of the bank, pushed water up the bank against the gravity of the thin face and set off surges
    const bool cv = (p->geo_curv==1);
    const bool hbed = hydro_sflow::wetting(p);
    double nx,ny,wl,wr;

    SLICELOOP1
    {
        if(Hx(i,j)>0.0 && cv)
        {
        const double dn = fdn(p,1,i,j);
        fnorm(p,1,i,j,nx,ny);
        gx0(i,j) = (b->eta(i+1,j) - b->eta(i,j))/dn;

        if(hbed)
        {
        face_weights(p,b,1,i,j,wl,wr);
        qxE(i,j) = alx(i,j)*(wl*(UHp(i,j)*nx + VHp(i,j)*ny) + wr*(UHp(i+1,j)*nx + VHp(i+1,j)*ny) - dt*g*Hx(i,j)*(1.0-theta)*gx0(i,j));
        }

        else
        qxE(i,j) = alx(i,j)*(0.5*(UHp(i,j) + UHp(i+1,j))*nx + 0.5*(VHp(i,j) + VHp(i+1,j))*ny - dt*g*Hx(i,j)*(1.0-theta)*gx0(i,j));

        K = dissipation(p,qxn(i,j),Hx(i,j),dn);
        kx(i,j) = K;

        ax(i,j) = theta*dt*(g*theta*dt*alx(i,j)*Hx(i,j) + K)/dn;
        }

        else
        if(Hx(i,j)>0.0)
        {
        gx0(i,j) = (b->eta(i+1,j) - b->eta(i,j))/p->DXP[IP];

        if(hbed)
        {
        face_weights(p,b,1,i,j,wl,wr);
        qxE(i,j) = alx(i,j)*(wl*UHp(i,j) + wr*UHp(i+1,j) - dt*g*Hx(i,j)*(1.0-theta)*gx0(i,j));
        }

        else
        qxE(i,j) = alx(i,j)*(0.5*(UHp(i,j) + UHp(i+1,j)) - dt*g*Hx(i,j)*(1.0-theta)*gx0(i,j));

        K = dissipation(p,qxn(i,j),Hx(i,j),p->DXP[IP]);
        kx(i,j) = K;

        ax(i,j) = theta*dt*(g*theta*dt*alx(i,j)*Hx(i,j) + K)/p->DXP[IP];
        }

        else
        {
        gx0(i,j) = 0.0;
        qxE(i,j) = 0.0;
        kx(i,j) = 0.0;
        ax(i,j) = 0.0;
        }
    }

    if(p->j_dir==1)
    SLICELOOP2
    {
        if(Hy(i,j)>0.0 && cv)
        {
        const double dn = fdn(p,2,i,j);
        fnorm(p,2,i,j,nx,ny);
        gy0(i,j) = (b->eta(i,j+1) - b->eta(i,j))/dn;

        if(hbed)
        {
        face_weights(p,b,2,i,j,wl,wr);
        qyE(i,j) = aly(i,j)*(wl*(UHp(i,j)*nx + VHp(i,j)*ny) + wr*(UHp(i,j+1)*nx + VHp(i,j+1)*ny) - dt*g*Hy(i,j)*(1.0-theta)*gy0(i,j));
        }

        else
        qyE(i,j) = aly(i,j)*(0.5*(UHp(i,j) + UHp(i,j+1))*nx + 0.5*(VHp(i,j) + VHp(i,j+1))*ny - dt*g*Hy(i,j)*(1.0-theta)*gy0(i,j));

        K = dissipation(p,qyn(i,j),Hy(i,j),dn);
        ky(i,j) = K;

        ay(i,j) = theta*dt*(g*theta*dt*aly(i,j)*Hy(i,j) + K)/dn;
        }

        else
        if(Hy(i,j)>0.0)
        {
        gy0(i,j) = (b->eta(i,j+1) - b->eta(i,j))/p->DYP[JP];

        if(hbed)
        {
        face_weights(p,b,2,i,j,wl,wr);
        qyE(i,j) = aly(i,j)*(wl*VHp(i,j) + wr*VHp(i,j+1) - dt*g*Hy(i,j)*(1.0-theta)*gy0(i,j));
        }

        else
        qyE(i,j) = aly(i,j)*(0.5*(VHp(i,j) + VHp(i,j+1)) - dt*g*Hy(i,j)*(1.0-theta)*gy0(i,j));

        K = dissipation(p,qyn(i,j),Hy(i,j),p->DYP[JP]);
        ky(i,j) = K;

        ay(i,j) = theta*dt*(g*theta*dt*aly(i,j)*Hy(i,j) + K)/p->DYP[JP];
        }

        else
        {
        gy0(i,j) = 0.0;
        qyE(i,j) = 0.0;
        ky(i,j) = 0.0;
        ay(i,j) = 0.0;
        }
    }

    // X 10 3: the ship pressure as the surface level it holds, -p/(rho g), its compact face gradient with the gravity
    // term (the mean of the pressure of the start and of the end of the step)
    if(ship)
    {
        SLICELOOP1
        if(Hx(i,j)>0.0)
        qxE(i,j) -= alx(i,j)*dt*g*Hx(i,j)*(sps(i+1,j) - sps(i,j))/p->DXP[IP];

        if(p->j_dir==1)
        SLICELOOP2
        if(Hy(i,j)>0.0)
        qyE(i,j) -= aly(i,j)*dt*g*Hy(i,j)*(sps(i,j+1) - sps(i,j))/p->DYP[JP];
    }

    exchange1(p,pgc,gx0);
    exchange1(p,pgc,qxE);
    exchange1(p,pgc,ax);
    exchange1(p,pgc,kx);
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
        slice &af  = xf ? (slice&)ax  : (slice&)ay;
        slice &qEf = xf ? (slice&)qxE : (slice&)qyE;
        slice &gf  = xf ? (slice&)gx0 : (slice&)gy0;
        slice &Up  = xf ? (slice&)UHp : (slice&)VHp;

        if(f.type==bc_level && Hf(f.fi,f.fj)>0.0)
        {
        const double up = p->geo_curv==1 ? UHp(i,j)*f.nx + VHp(i,j)*f.ny : Up(i,j);
        qEf(f.fi,f.fj) = alf(f.fi,f.fj)*(up - dt*g*Hf(f.fi,f.fj)*(1.0-theta)*gf(f.fi,f.fj));
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
    qxn(i,j) = 0.0;
    
    if(p->j_dir==1)
    SLICELOOP2
    if(Hy(i,j)<=0.0)
    qyn(i,j) = 0.0;
    
    exchange1(p,pgc,qxn);
    exchange2(p,pgc,qyn);
    
    for(const bc_face &f : bcf)
    if(f.type==bc_level && (f.cs==1 || f.cs==4 ? Hx(f.fi,f.fj) : Hy(f.fi,f.fj))<=0.0)
    {
        if(f.cs==1 || f.cs==4)
        qxn(f.fi,f.fj) = 0.0;
        
        else
        qyn(f.fi,f.fj) = 0.0;
    }
    
}

void sflow_momentum_SI::si_rhs(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const double dt = p->dt;

    // right-hand side of the water depth: WL^n - dt div( theta q^E + (1-theta) q^n )
    SLICELOOP4
    rhsW(i,j) = WL0(i,j) - dt*(theta*div(p,qxE,qyE) + (1.0-theta)*div(p,qxn,qyn));

    // + dt K: the Stokes mass source of the SEASTATE vortex force (A 751 2), K = -div(M) of the wet cells, held
    // over the step (also in eta_update)
    mass_source = (b->wave!=nullptr && p->A751==2);

    if(mass_source)
    {
        SLICELOOP4
        wsrc(i,j) = 0.0;

    b->wave->mass_source(p,b,wsrc);

        SLICELOOP4
        rhsW(i,j) += dt*wsrc(i,j);
    }

    // + dt S: REEF3D::HYDRO rain and point sources on wet and dry cells, held over the step (also in eta_update)
    if(b->hydro!=nullptr)
    {
        SLICELOOP4
        hsrc(i,j) = 0.0;

    b->hydro->mass_source(p,b,hsrc);

        SLICELOOP4
        rhsW(i,j) += dt*hsrc(i,j);
    }

    // initial guess eta^n; wet pattern: all wet (the linear system of P1), then the
    // Newton iteration on WL = max(eta + d, 0) (Casulli 2009), converging from above
    SLICELOOP4
    etaI(i,j) = eta0(i,j);
}

void sflow_momentum_SI::face_update(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // q^{n+1} = q^E - (alpha_f dt g H_f theta + K_f) grad_f eta^{n+1} = q^E - a_f/(theta dt) (eta_r - eta_l)
    const double tdt = theta*p->dt;

    SLICELOOP1
    qx(i,j) = Hx(i,j)>0.0 ? qxE(i,j) - ax(i,j)/tdt*(etaI(i+1,j) - etaI(i,j)) : 0.0;

    if(p->j_dir==1)
    SLICELOOP2
    qy(i,j) = Hy(i,j)>0.0 ? qyE(i,j) - ay(i,j)/tdt*(etaI(i,j+1) - etaI(i,j)) : 0.0;

    // non-hydrostatic part of the face discharges (A 220 2/3, interior faces)
    if(nh)
    {
        SLICELOOP1
        if(Hx(i,j)>0.0)
        qx(i,j) += qxN(i,j);

        if(p->j_dir==1)
        SLICELOOP2
        if(Hy(i,j)>0.0)
        qy(i,j) += qyN(i,j);
    }

    exchange1(p,pgc,qx);
    exchange2(p,pgc,qy);

    // open boundary faces
    for(const bc_face &f : bcf)
    {
    i=f.i;
    j=f.j;

        const bool xf = (f.cs==1 || f.cs==4);

        slice &Hf = xf ? (slice&)Hx  : (slice&)Hy;
        slice &af = xf ? (slice&)ax  : (slice&)ay;
        slice &qEf = xf ? (slice&)qxE : (slice&)qyE;
        slice &qf = xf ? (slice&)qx  : (slice&)qy;

        if(f.type==bc_level)
        qf(f.fi,f.fj) = Hf(f.fi,f.fj)>0.0 ? qEf(f.fi,f.fj) - af(f.fi,f.fj)/tdt*f.s*(etaI(i,j) - etaI(f.gi,f.gj)) : 0.0;

        if(f.type==bc_flux)
        qf(f.fi,f.fj) = f.s*f.qin;

        if(f.type==bc_riemann || f.type==bc_flather)
        qf(f.fi,f.fj) = f.s*(f.A - f.B*etaI(i,j));
    }
}

void sflow_momentum_SI::eta_update(lexer *p, fdm2D *b, ghostcell *pgc)
{
    starttime=pgc->timer();

    // WL from the divergence of the face discharges: the volume is conserved exactly,
    // independent of the tolerance of the eta solve; with the wetting and drying iteration
    // WL >= 0 up to round-off
    SLICELOOP4
    {
    double w = WL0(i,j) - p->dt*(theta*div(p,qx,qy) + (1.0-theta)*div(p,qxn,qyn));

    if(mass_source)
    w += p->dt*wsrc(i,j);

    if(b->hydro!=nullptr)
    w += p->dt*hsrc(i,j);

    b->WL(i,j) = MAX(w, 0.0);
    b->eta(i,j) = b->WL(i,j) - b->depth(i,j);
    }
    
    // depth-averaged suspended load (S 120 > 0): the face discharges of this continuity step,
    // theta q + (1-theta) q^n, the mass flux of the suspended load (as FEacc of the explicit scheme)
    if(b->FEacc_x!=nullptr && p->geo_curv==0)
    {
        for(i=-1; i<p->knox; ++i)
        for(j=0; j<p->knoy; ++j)
        (*b->FEacc_x)(i,j) = theta*qx(i,j) + (1.0-theta)*qxn(i,j);
        
        for(i=0; i<p->knox; ++i)
        for(j=-1; j<p->knoy; ++j)
        (*b->FEacc_y)(i,j) = (p->j_dir==1) ? theta*qy(i,j) + (1.0-theta)*qyn(i,j) : 0.0;
        
    b->FEacc_step = p->count;
    }

    // discharge over the open boundaries (theta weighted, inward positive)
    Qin = Qout = 0.0;

    for(const bc_face &f : bcf)
    {
        const bool xf = (f.cs==1 || f.cs==4);
        const double q = f.s*(theta*(xf ? qx(f.fi,f.fj) : qy(f.fi,f.fj)) + (1.0-theta)*(xf ? qxn(f.fi,f.fj) : qyn(f.fi,f.fj)))*f.ds;

        if(q>0.0)
        Qin += q;

        else
        Qout += q;

        // REEF3D::HYDRO coupled (R 100 2): the outflow volume of the boundary cell, returned to the catchment model
        if(q<0.0 && b->hydro!=nullptr && b->hydro->cpl)
        b->hydro->Bq(f.i,f.j) -= q*p->dt;
    }

    // a block or a patch of the mesh refinement: its own sum (sflow_amr::si_log sums the grids over the ranks; a global
    // sum per patch hung the run where the ranks held different numbers of patches)
    if(!patch)
    {
    Qin = pgc->globalsum(Qin);
    Qout = pgc->globalsum(Qout);
    }

    // ghost cells; the interior values are kept: gcsl_start4 extends a cell next to a solid
    // three cells deep (gcsl_neumann), into the fluid behind a solid one or two cells thick
    SLICELOOP4
    {
    tU(i,j) = b->eta(i,j);
    tV(i,j) = b->WL(i,j);
    }

    pgc->gcsl_start4(p,b->eta,gcval_eta);
    pgc->gcsl_start4(p,b->WL,gcval_eta);

    // ghost values of the open water level boundaries
    for(const bc_face &f : bcf)
    if(f.type==bc_level)
    {
        const int di = f.gi-f.i, dj = f.gj-f.j;

        for(int q=1;q<=3;++q)
        b->eta(f.i+q*di,f.j+q*dj) = etaB(f.gi,f.gj);
    }

    // wetting & drying (flags, dry cells without momentum); the interior eta and WL are
    // kept, wetdry passes them through gcsl_start4 as well
    pfsf->wetdry(p,b,pgc,b->WL);

    SLICELOOP4
    {
    b->eta(i,j) = tU(i,j);
    b->WL(i,j) = tV(i,j);
    }

    p->lsmtime+=pgc->timer()-starttime;
}

void sflow_momentum_SI::cell_update(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // cell momentum with the mean of the open face gradients of the cell. A face towards a cell
    // that is dry at the end of the step (eta^{n+1} + d <= A 244) does not count, as in
    // nhflow_momentum_SI: the eta of the Newton iteration is not a water level there, its row of the
    // eta system holds only the coefficient of the thin face (H_f just above A 244), and it can lie
    // metres below the bed; its gradient gave the wet neighbour a spurious jet of tens of m/s, the
    // flicker of the gauges near shallow banks (1 km Michael hindcast, validation/09_michael)
    const double dt = p->dt;
    double we,ww,gn,gn1,wsum;

    auto wn = [&](int a, int c) { return etaI(a,c) + b->depth(a,c) > p->A244; };

    wall_sides(p);

    // the depth of the gravity term (the non-hydrostatic step and A 210 11: extrapolated to t^{n+1/2}, nh_levels)
    slice &Hg = tlev ? (slice&)WLh : (slice&)WL0;

    SLICELOOP4
    {
        if(p->wet[IJ]==1)
        {
        we = (Hx(i,j)>0.0 && wn(i+1,j)) ? 1.0 : 0.0;
        ww = (Hx(i-1,j)>0.0 && wn(i-1,j)) ? 1.0 : 0.0;
        wsum = MAX(we+ww,1.0);

        // a wall face (every predictor; A 210 12/13 since the second-order time levels, before A 210 11 only) counts in
        // the mean with the term that keeps the face discharge zero there: g H d_n eta
        // = the sources of the predictor normal to the wall (wind, Coriolis, atmospheric pressure; FSx, FSy of
        // stage_rhs), zero in free waves (the mirrored level). With zero also under wind a wall cell of a closed basin
        // got half the gradient that balances the stress, and the wind drove a flow into the walls (NHFLOW A 510 11,
        // nhflow_2d_si_wind_setup_a547: the wall cells 4 times the depth-mean velocity of the interior). Only with an
        // open face in the direction (as NHFLOW): a cell between two walls has no gradient
        double fw = 0.0;

        if(we+ww>0.0)
        {
        const double wx = wall(p,i,j,1) + wall(p,i,j,2);
        wsum = MAX(we + ww + wx, 1.0);
        fw = wx*FSx(i,j);
        }

        // a face between two wet deep cells (mfx, mlflux_faces) with its own depth H_f, the mean of the two cells: the
        // mean of g H_f d_x eta of the faces conserves the momentum across a bore (a difference of g h^2/2)
        const double hw = mfx(i-1,j)>0.5 ? Hx(i-1,j)/Hg(i,j) : 1.0, he = mfx(i,j)>0.5 ? Hx(i,j)/Hg(i,j) : 1.0;

        gn  = (ww*hw*gx0(i-1,j) + we*he*gx0(i,j))/wsum;
        gn1 = (ww*hw*(etaI(i,j) - etaI(i-1,j))/p->DXP[IM1] + we*he*(etaI(i+1,j) - etaI(i,j))/p->DXP[IP])/wsum;

        if(ship)
        gn1 += (ww*hw*(sps(i,j) - sps(i-1,j))/p->DXP[IM1] + we*he*(sps(i+1,j) - sps(i,j))/p->DXP[IP])/(wsum*theta);

        b->UH(i,j) = alc(i,j)*(UHp(i,j) - dt*g*Hg(i,j)*(theta*gn1 + (1.0-theta)*gn) - dt*fw/wsum);

        b->VH(i,j) = 0.0;

            if(p->j_dir==1)
            {
            we = (Hy(i,j)>0.0 && wn(i,j+1)) ? 1.0 : 0.0;
            ww = (Hy(i,j-1)>0.0 && wn(i,j-1)) ? 1.0 : 0.0;
            wsum = MAX(we+ww,1.0);

            fw = 0.0;

            if(we+ww>0.0)
            {
            const double wy = wall(p,i,j,3) + wall(p,i,j,4);
            wsum = MAX(we + ww + wy, 1.0);
            fw = wy*FSy(i,j);
            }

            const double hs = mfy(i,j-1)>0.5 ? Hy(i,j-1)/Hg(i,j) : 1.0, hn = mfy(i,j)>0.5 ? Hy(i,j)/Hg(i,j) : 1.0;

            gn  = (ww*hs*gy0(i,j-1) + we*hn*gy0(i,j))/wsum;
            gn1 = (ww*hs*(etaI(i,j) - etaI(i,j-1))/p->DYP[JM1] + we*hn*(etaI(i,j+1) - etaI(i,j))/p->DYP[JP])/wsum;

            if(ship)
            gn1 += (ww*hs*(sps(i,j) - sps(i,j-1))/p->DYP[JM1] + we*hn*(sps(i,j+1) - sps(i,j))/p->DYP[JP])/(wsum*theta);

            b->VH(i,j) = alc(i,j)*(VHp(i,j) - dt*g*Hg(i,j)*(theta*gn1 + (1.0-theta)*gn) - dt*fw/wsum);
            }

            // curvilinear grid: per index direction the mean normal and the mean normal gradient of the open faces,
            // the cell vector from the two (as nhflow_momentum_SI); with open faces in one direction only its component
            // along that normal
            if(p->geo_curv==1)
            {
            double nx,ny;
            double ax_=0.0, ay_=0.0, ag=0.0, bx_=0.0, by_=0.0, bg=0.0;
            int na=0, nb=0;

            // A 210 11: the index directions with an open face (the wall faces carry the sources only there)
            const bool o1 = (Hx(i,j)>0.0 && wn(i+1,j)) || (Hx(i-1,j)>0.0 && wn(i-1,j));
            const bool o2 = (Hy(i,j)>0.0 && wn(i,j+1)) || (Hy(i,j-1)>0.0 && wn(i,j-1));

            auto face = [&](int dir, int fi, int fj, int oi, int oj, double H, double G0, double G1)
            {
                // a wall face counts with its normal and the gradient that keeps the discharge through it
                // zero, g H d_n eta = the sources normal to the face (as the Cartesian grid above; zero in free waves,
                // the mirrored level)
                const int sd = dir==1 ? (fi==i ? 2 : 1) : (fj==j ? 4 : 3);
                const bool wf = (H<=0.0 || !wn(oi,oj)) && wall(p,i,j,sd)>0.5;

                if((H<=0.0 || !wn(oi,oj)) && !wf)
                return;

                fnorm(p,dir,fi,fj,nx,ny);
                const bool wo = wf && (dir==1 ? o1 : o2) && Hg(i,j)>0.0;
                // a face between two wet deep cells with its own depth (mlflux_faces), as the Cartesian grid
                const double mf = dir==1 ? mfx(fi,fj) : mfy(fi,fj);
                const double hr = (!wf && mf>0.5) ? H/Hg(i,j) : 1.0;
                const double G = wf ? (wo ? (FSx(i,j)*nx + FSy(i,j)*ny)/(g*Hg(i,j)) : 0.0) : hr*(theta*G1 + (1.0-theta)*G0);

                if(dir==1)
                {
                ax_ += nx; ay_ += ny; ag += G; ++na;
                }

                else
                {
                bx_ += nx; by_ += ny; bg += G; ++nb;
                }
            };

            face(1,i,j,i+1,j,Hx(i,j),gx0(i,j),(etaI(i+1,j) - etaI(i,j))/fdn(p,1,i,j));
            face(1,i-1,j,i-1,j,Hx(i-1,j),gx0(i-1,j),(etaI(i,j) - etaI(i-1,j))/fdn(p,1,i-1,j));
            face(2,i,j,i,j+1,Hy(i,j),gy0(i,j),(etaI(i,j+1) - etaI(i,j))/fdn(p,2,i,j));
            face(2,i,j-1,i,j-1,Hy(i,j-1),gy0(i,j-1),(etaI(i,j) - etaI(i,j-1))/fdn(p,2,i,j-1));

            auto unit = [&](double &x, double &y, double &gg, int n)
            {
                const double l = sqrt(x*x + y*y);
                x = l>0.0 ? x/l : 0.0;
                y = l>0.0 ? y/l : 0.0;
                gg = n>0 ? gg/double(n) : 0.0;
            };

            unit(ax_,ay_,ag,na);
            unit(bx_,by_,bg,nb);

            double gx=0.0, gy=0.0;
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

            b->UH(i,j) = alc(i,j)*(UHp(i,j) - dt*g*Hg(i,j)*gx);
            b->VH(i,j) = alc(i,j)*(VHp(i,j) - dt*g*Hg(i,j)*gy);
            }
        }

        else
        {
        b->UH(i,j) = 0.0;
        b->VH(i,j) = 0.0;
        }

    b->WH(i,j) = 0.0;
    }

    // non-hydrostatic pressure (A 220 2/3): UH with the mean of the face terms, WH
    if(nh)
    nh_cell(p,b,pgc);

    // A 210 12/13: the acceleration of the implicit part of the step (gravity, non-hydrostatic pressure, bed friction)
    // for the later predictor stages of the next step (stage_vel)
    // in wet deep cells (WL of t^n and t^{n+1} above 10 A 244, not next to the shoreline); 0 in thin water, where
    // the momentum of the step is that of wetting, drying and the friction of a film
    if(tl12)
    SLICELOOP4
    {
    const bool dp = p->wet[IJ]==1 && p->deep[IJ]==1 && WL0(i,j)>10.0*p->A244 && etaI(i,j) + b->depth(i,j)>10.0*p->A244;
    acU(i,j) = dp ? (b->UH(i,j) - UHp(i,j))/dt : 0.0;
    acV(i,j) = dp ? (b->VH(i,j) - VHp(i,j))/dt*p->y_dir : 0.0;
    }

    // X 10 2: the body moves over the step, the direct forcing sets the momentum of its cells to the body velocity
    // (sflow_forcing with the velocities of the new momentum, as the stages of the explicit step); the face discharges
    // of the next step follow from the cells (si_assemble)
    if(p->X10==2)
    {
        SLICELOOP4
        {
        tU(i,j) = b->U(i,j);
        tV(i,j) = b->V(i,j);
        const bool w = (p->wet[IJ]==1 && b->WL(i,j)>p->A244);
        b->U(i,j) = w ? b->UH(i,j)/b->WL(i,j) : 0.0;
        b->V(i,j) = (w && p->j_dir==1) ? b->VH(i,j)/b->WL(i,j) : 0.0;
        }

    psfdf->forcing(p,b,pgc,p6dof,0,1.0,b->UH,b->VH,b->WH,b->WL,true);

        SLICELOOP4
        {
        b->U(i,j) = tU(i,j);
        b->V(i,j) = tV(i,j);
        }
    }

    // velocities of the interior (Froude limiter A 531)
    velcalc(p,b,pgc,b->UH,b->VH,b->WH,b->WL,0);

    // relaxation zones (once per step)
    relax(p,b,pgc);

    // ghost cells and face diagnostics; the interior is kept, as in rk_stage
    SLICELOOP4
    {
    tU(i,j) = b->U(i,j);
    tV(i,j) = b->V(i,j);
    tUH(i,j) = b->UH(i,j);
    tVH(i,j) = b->VH(i,j);
    rhs(i,j) = b->WL(i,j);
    }

    velcalc(p,b,pgc,b->UH,b->VH,b->WH,b->WL,2);

    SLICELOOP4
    {
    b->U(i,j) = tU(i,j);
    b->V(i,j) = tV(i,j);
    b->UH(i,j) = tUH(i,j);
    b->VH(i,j) = tVH(i,j);
    b->WL(i,j) = rhs(i,j);
    }

    face_velocities(p,b,pgc,b->WL);
}

void sflow_momentum_SI::relax(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // relaxation zones (B 98 2, B 99 1/2, B 520 1/2), once per step: eta, then U, V with
    // the new water depth; the targets of the end of the step with B 95
    // (the zone weight acts once per step instead of once per RK stage)
    vrelax = 0.0;

    if(p->B98!=2 && p->B99!=1 && p->B99!=2 && p->B90==0)
    return;

    if(p->B95!=0.0 && p->B98==2)
    pflow->wavegen_2D_stage(p,b,pgc,p->simtime+p->dt);

    SLICELOOP4
    {
    tUH(i,j) = b->UH(i,j);
    tVH(i,j) = b->VH(i,j);
    rhs(i,j) = b->WL(i,j);
    }

    pflow->eta_relax(p,pgc,b->eta);

    SLICELOOP4
    {
    b->WL(i,j) = MAX(b->eta(i,j) + b->depth(i,j), 0.0);
    b->eta(i,j) = b->WL(i,j) - b->depth(i,j);
    }

    pflow->um_relax(p,pgc,b->U,b->UH,b->WL);
    pflow->vm_relax(p,pgc,b->V,b->VH,b->WL);

    // non-hydrostatic pressure (A 220 2/3): W and P as the explicit stages
    if(nh)
    {
    pflow->wm_relax(p,pgc,b->W,b->WH,b->WL);
    pflow->pm_relax(p,pgc,b->press);
    }

    SLICELOOP4
    vrelax += (b->WL(i,j) - rhs(i,j))*carea(p,i,j);

    // the face discharges follow the change of the cell discharges
    SLICELOOP4
    {
    dUH(i,j) = b->UH(i,j) - tUH(i,j);
    dVH(i,j) = (b->VH(i,j) - tVH(i,j))*p->y_dir;
    }

    exchange4(p,pgc,dUH);
    exchange4(p,pgc,dVH);

    if(p->geo_curv==1)
    {
        double nx,ny;

        SLICELOOP1
        if(Hx(i,j)>0.0)
        {
        fnorm(p,1,i,j,nx,ny);
        qx(i,j) += 0.5*(dUH(i,j) + dUH(i+1,j))*nx + 0.5*(dVH(i,j) + dVH(i+1,j))*ny;
        }

        SLICELOOP2
        if(Hy(i,j)>0.0)
        {
        fnorm(p,2,i,j,nx,ny);
        qy(i,j) += 0.5*(dUH(i,j) + dUH(i,j+1))*nx + 0.5*(dVH(i,j) + dVH(i,j+1))*ny;
        }
    }

    else
    {
    SLICELOOP1
    if(Hx(i,j)>0.0)
    qx(i,j) += 0.5*(dUH(i,j) + dUH(i+1,j));

    if(p->j_dir==1)
    SLICELOOP2
    if(Hy(i,j)>0.0)
    qy(i,j) += 0.5*(dVH(i,j) + dVH(i,j+1));
    }

    exchange1(p,pgc,qx);
    exchange2(p,pgc,qy);

    if(!patch)
    vrelax = pgc->globalsum(vrelax);
}

void sflow_momentum_SI::log(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double vol=0.0;

    SLICELOOP4
    vol += b->WL(i,j)*carea(p,i,j);

    vol = pgc->globalsum(vol);

    if(p->mpirank==0)
    silog<<p->count<<" "<<setprecision(10)<<p->simtime+p->dt<<" "<<p->dt<<" "<<setprecision(16)<<vol<<" "
         <<Qin<<" "<<Qout<<" "<<vrelax<<" "<<psi->iter<<" "<<wd_iter<<" "<<wd_cells<<"\n";
}

double sflow_momentum_SI::div(lexer *p, slice &fx, slice &fy)
{
    // divergence of face discharges (per unit face length, along the face normals)
    if(p->geo_curv==1)
    return (fx(i,j)*flen(p,1,i,j) - fx(i-1,j)*flen(p,1,i-1,j) + fy(i,j)*flen(p,2,i,j) - fy(i,j-1)*flen(p,2,i,j-1))/carea(p,i,j);

    double d = (fx(i,j) - fx(i-1,j))/p->DXN[IP];

    if(p->j_dir==1)
    d += (fy(i,j) - fy(i,j-1))/p->DYN[JP];

    return d;
}

void sflow_momentum_SI::exchange1(lexer *p, ghostcell *pgc, slice &f)
{
    // MPI only: the first face of a partition belongs to the neighbour, physical
    // boundary faces are set by bc_setup
    if(p->mpi_size>1)
    {
    pgc->gcslparax(p,f,1);
    pgc->gcslparacox(p,f,14);
    }
}

void sflow_momentum_SI::exchange2(lexer *p, ghostcell *pgc, slice &f)
{
    if(p->mpi_size>1 && p->j_dir==1)
    {
    pgc->gcslparax(p,f,2);
    pgc->gcslparacox(p,f,14);
    }
}

void sflow_momentum_SI::exchange4(lexer *p, ghostcell *pgc, slice &f)
{
    if(p->mpi_size>1)
    {
    pgc->gcslparax(p,f,4);
    pgc->gcslparacox(p,f,14);
    }
}

double sflow_momentum_SI::flen(lexer *p, int dir, int ii, int jj)
{
    // length of face i+1/2 (dir 1) or j+1/2 (dir 2) of cell (ii,jj)
    if(p->geo_curv==1)
    return dir==1 ? p->LEN1[(ii-p->imin)*p->jmax + (jj-p->jmin)] : p->LEN2[(ii-p->imin)*p->jmax + (jj-p->jmin)];

    if(dir==1)
    return p->j_dir==1 ? p->DYN[jj+marge] : 1.0;

    return p->DXN[ii+marge];
}

double sflow_momentum_SI::fdn(lexer *p, int dir, int ii, int jj)
{
    // distance of the two cell centres of the face along its normal
    if(p->geo_curv==1)
    {
    const double dn = sflow_curv::dnorm(p,dir,ii,jj);
    return dn>1.0e-12 ? dn : 1.0e-12;
    }

    return dir==1 ? p->DXP[ii+marge] : p->DYP[jj+marge];
}

double sflow_momentum_SI::carea(lexer *p, int ii, int jj)
{
    if(p->geo_curv==1)
    return p->AREA2D[(ii-p->imin)*p->jmax + (jj-p->jmin)];

    return p->DXN[ii+marge]*p->DYN[jj+marge];
}

void sflow_momentum_SI::fnorm(lexer *p, int dir, int ii, int jj, double &nx, double &ny)
{
    // unit normal of the face towards +i (dir 1) / +j (dir 2)
    if(p->geo_curv==1)
    {
    double len;
    sflow_curv::normal(p,dir,ii,jj,nx,ny,len);
    return;
    }

    nx = dir==1 ? 1.0 : 0.0;
    ny = dir==1 ? 0.0 : 1.0;
}

void sflow_momentum_SI::ship_level(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // X 10 3: the ship pressure of the step (sixdof_sflow::sflow_press_step, the mean of its start and end) as the
    // surface level it holds, -p/(rho g); the stored pressure is minus the pressure on the water
    ship = p6dof->sflow_press_step(p,sps);

    if(ship)
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    sps(ii,jj) = -sps(ii,jj)/(p->W1*g);

    exchange4(p,pgc,sps);
}
