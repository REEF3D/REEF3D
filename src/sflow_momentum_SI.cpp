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
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"reefckpt.h"
#include"sflow_fsf.h"
#include"sflow_reconstruct.h"
#include"sflow_diffusion.h"
#include"sflow_rheology.h"
#include<sys/stat.h>
#include<sys/types.h>

sflow_momentum_SI::sflow_momentum_SI(lexer *p, fdm2D *b, ghostcell *pgc, sflow_HLL *pphll, sflow_signal_speed *ppss,
                                     sflow_reconstruct *pprecon, sflow_diffusion *ppdiff, sflow_pressure *pppress,
                                     solver2D *ppsolv, solver2D *pppoissonsolv, ioflow *ppflow, sflow_fsf *ppfsf,
                                     sflow_forcing *ppsfdf, sixdof *pp6dof)
                                    : sflow_momentum_func(p,b,pgc,pphll,ppss,pprecon,ppdiff,pppress,ppsolv,pppoissonsolv,ppflow,ppfsf,ppsfdf,pp6dof),
                                      qx(p),qxn(p),qxE(p),Hx(p),alx(p),ax(p),kx(p),gx0(p),
                                      qy(p),qyn(p),qyE(p),Hy(p),aly(p),ay(p),ky(p),gy0(p),
                                      fxU(p),fxV(p),fyU(p),fyV(p),
                                      eta0(p),WL0(p),UH0(p),VH0(p),UHa(p),VHa(p),UHb(p),VHb(p),UHp(p),VHp(p),
                                      WHz(p),dUH(p),dVH(p),etaI(p),alc(p),rhs(p),rhsW(p),bdiag(p),
                                      tU(p),tV(p),tUH(p),tVH(p),WLk(p),etaB(p)
{
    check(p);

    psi = new si_eta2D(p);

    theta = p->A227;
    g = fabs(p->W22);
    rk = (p->A210==2)?2:3;
    qini = false;
    wd_iter = 0;
    wd_cells = 0;
    Qin = Qout = 0.0;

    bc_method(p);

    if(p->mpirank==0)
    {
    cout<<"SFLOW semi-implicit free surface (A 226 1): theta "<<theta<<", SSP-RK"<<rk<<" predictor";

    if(p->A229==1)
    cout<<", wetting and drying iteration";

    if(p->A231>0.0)
    cout<<", mass flux dissipation "<<p->A231;

    cout<<endl;

    mkdir("./REEF3D_SFLOW_Log",0777);
    if(reefckpt_open(silog,"./REEF3D_SFLOW_Log/REEF3D_SFLOW_SI.dat",p->I45))   // I 45: continued
    {
        silog<<"# semi-implicit SFLOW (A 226 1): volume and open boundary discharges"<<endl;
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
        cout<<"!!! SFLOW A 226 1 (semi-implicit): "<<s<<" not available yet !!!"<<endl;
        stop=1;
    };

    // A 220 defaults to 2: the semi-implicit path needs A 220 0 in ctrl.txt
    if(p->A220!=0)
    msg("non-hydrostatic pressure and Boussinesq (A 220 > 0, default 2) are; set A 220 0,");

    if(p->G1>0 || p->geo_mb>0)
    msg("mesh refinement and multi-block grids (G 1, G 50) are");

    if(p->geo_curv==1)
    msg("the curvilinear grid is");

    if(p->B411>0 || p->B412>0 || p->B413>0 || p->B414>0 || p->B415>0 || p->B416>0 || p->B417>0 || p->B418>0 || p->B421>0 || p->B422>0)
    msg("patch boundary conditions (B 411-422) are");

    if(p->X10>0)
    msg("6DOF and direct forcing (X 10) are");

    if(p->A750>0)
    msg("SEASTATE coupling (A 750) is");

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
        const int e = p->B521_edge[k];

        if(e<1 || e>4)
        continue;

        const int cs = edge_cs[e];

        if(method[cs]!=0 && method[cs]!=p->B520_method[qn])
        {
            if(p->mpirank==0)
            cout<<"!!! SFLOW A 226 1 (semi-implicit): different B 520 methods on one edge (B 521 edge "<<e<<") are not available yet !!!"<<endl;
            exit(1);
        }

        method[cs] = p->B520_method[qn];
    }
}

void sflow_momentum_SI::start(lexer *p, fdm2D* b, ghostcell* pgc)
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

    // ghost cells of the open boundaries (ioflow discharge, iowave open edges), as the
    // explicit integrators
    inflow(p,b,pgc,pflow);

    // state at t^n
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

    // first step: face discharges from the cell values
    if(!qini)
    {
        SLICELOOP1
        qx(i,j) = Hx(i,j)>0.0 ? 0.5*(b->UH(i,j) + b->UH(i+1,j)) : 0.0;

        if(p->j_dir==1)
        SLICELOOP2
        qy(i,j) = Hy(i,j)>0.0 ? 0.5*(b->VH(i,j) + b->VH(i,j+1)) : 0.0;

    exchange1(p,pgc,qx);
    exchange2(p,pgc,qy);

    // face depths upwind of the new discharges
    face_depth(p,b,pgc);
    }

    friction(p,b);

    // 2. open boundaries (the first step also sets their face discharges)
    bc_setup(p,b,pgc);

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

    // 3. predictor of the explicit terms
    predictor(p,b,pgc);

    // 4.-7. face predictor, eta system, face discharges
    eta_solve(p,b,pgc);
    face_update(p,b,pgc);

    // 8. WL, eta, wetting and drying
    eta_update(p,b,pgc);

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
    double df,eu,h;

    SLICELOOP1
    {
    Hx(i,j) = 0.0;

        if(p->flagslice4[IJ]>0 && p->flagslice4[Ip1J]>0)
        {
        if(p->wet[IJ]==1 && p->wet[Ip1J]==1)
        df = 0.5*(b->depth(i,j) + b->depth(i+1,j));

        else
        df = MIN(b->depth(i,j),b->depth(i+1,j));

        if(qini && qx(i,j)>1.0e-12)
        eu = b->eta(i,j);

        else
        if(qini && qx(i,j)<-1.0e-12)
        eu = b->eta(i+1,j);

        else
        eu = MAX(b->eta(i,j),b->eta(i+1,j));

        h = eu + df;

        Hx(i,j) = h>p->A244 ? h : 0.0;
        }
    }

    if(p->j_dir==1)
    SLICELOOP2
    {
    Hy(i,j) = 0.0;

        if(p->flagslice4[IJ]>0 && p->flagslice4[IJp1]>0)
        {
        if(p->wet[IJ]==1 && p->wet[IJp1]==1)
        df = 0.5*(b->depth(i,j) + b->depth(i,j+1));

        else
        df = MIN(b->depth(i,j),b->depth(i,j+1));

        if(qini && qy(i,j)>1.0e-12)
        eu = b->eta(i,j);

        else
        if(qini && qy(i,j)<-1.0e-12)
        eu = b->eta(i,j+1);

        else
        eu = MAX(b->eta(i,j),b->eta(i,j+1));

        h = eu + df;

        Hy(i,j) = h>p->A244 ? h : 0.0;
        }
    }

    exchange1(p,pgc,Hx);
    exchange2(p,pgc,Hy);
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
        f.type = bc_flux;

        else
        if((bc==2 || bc==7 || bc==8) && outflow_open==1)
        f.type = ((p->F50==2 || p->F50==3) && p->B99<3) ? bc_level : bc_riemann;

        // wall: the face stays closed
        if(f.type==0)
        {
        qf(f.fi,f.fj) = 0.0;
        continue;
        }

        // interior and ghost state of t^n
        const double hi = WL0(i,j);
        const double etai = eta0(i,j);
        const double ung = f.s*(xf ? b->U(f.gi,f.gj) : b->V(f.gi,f.gj));
        const double hg = MAX(b->eta(f.gi,f.gj) + b->depth(f.gi,f.gj), 0.0);
        const double uni = f.s*(xf ? b->U(i,j) : b->V(i,j));

        f.qin = 0.0;
        f.A = 0.0;
        f.B = 0.0;

        // a dry boundary cell: only an inflowing discharge
        if((p->wet[IJ]==0 || hi<=p->A244) && f.type!=bc_flux)
        {
        qf(f.fi,f.fj) = 0.0;
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
        qf(f.fi,f.fj) = xf ? UH0(i,j) : VH0(i,j);

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

void sflow_momentum_SI::predictor(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // SSP-RK of the explicit terms, depth H^n and face discharges q^n frozen
    if(rk==2)
    {
    rk_stage(p,b,pgc,UH0,VH0,UHa,VHa,0.0);
    rk_stage(p,b,pgc,UHa,VHa,UHp,VHp,0.5);
    }

    if(rk==3)
    {
    rk_stage(p,b,pgc,UH0,VH0,UHa,VHa,0.0);
    rk_stage(p,b,pgc,UHa,VHa,UHb,VHb,0.75);
    rk_stage(p,b,pgc,UHb,VHb,UHp,VHp,1.0/3.0);
    }
}

void sflow_momentum_SI::rk_stage(lexer *p, fdm2D *b, ghostcell *pgc, slice &UHs, slice &VHs, slice &UHo, slice &VHo, double a)
{
    const double alpha = 1.0-a;

    // stage velocities on the depth of t^n
    SLICELOOP4
    {
        if(p->wet[IJ]==1 && WL0(i,j)>p->A244)
        {
        b->U(i,j) = UHs(i,j)/WL0(i,j);
        b->V(i,j) = VHs(i,j)/WL0(i,j)*p->y_dir;
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

    clearrhs(p,b);

    // advection with the face discharges of t^n
    advection(p,b,pgc);

    // sources and diffusion (as the explicit stage, without gravity and bed friction)
    prheo->u_source(p,b,b->U,b->V,WL0);
    irhs(p,b);
	pdiff->diff_u(p,b,pgc,psolv,UHDIFF,UHs,b->U,b->V,WL0,alpha);

    if(p->j_dir==1)
    {
    prheo->v_source(p,b,b->U,b->V,WL0);
    jrhs(p,b);
    }
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
    // conservative advection of UH, VH with the face discharges of t^n, velocities
    // upwind of the face discharge (reconstructed, A 211); no gravity wave speed in the
    // dissipation: the Courant number of the flow alone limits the time step
    precon->reconstruct_x(p, pgc, b, b->U, b->Us, b->Un);

    if(p->j_dir==1)
    {
    precon->reconstruct_x(p, pgc, b, b->V, b->Vs, b->Vn);
    precon->reconstruct_y(p, pgc, b, b->U, b->Ue, b->Uw);
    precon->reconstruct_y(p, pgc, b, b->V, b->Ve, b->Vw);
    }

    SLICELOOP1
    {
    fxU(i,j) = qxn(i,j)>=0.0 ? qxn(i,j)*b->Us(i,j) : qxn(i,j)*b->Un(i,j);
    fxV(i,j) = 0.0;

    if(p->j_dir==1)
    fxV(i,j) = qxn(i,j)>=0.0 ? qxn(i,j)*b->Vs(i,j) : qxn(i,j)*b->Vn(i,j);
    }

    if(p->j_dir==1)
    SLICELOOP2
    {
    fyU(i,j) = qyn(i,j)>=0.0 ? qyn(i,j)*b->Ue(i,j) : qyn(i,j)*b->Uw(i,j);
    fyV(i,j) = qyn(i,j)>=0.0 ? qyn(i,j)*b->Ve(i,j) : qyn(i,j)*b->Vw(i,j);
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

void sflow_momentum_SI::eta_solve(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const double dt = p->dt;
    double K;

    // predictor of the cells, to the faces
    exchange4(p,pgc,UHp);
    exchange4(p,pgc,VHp);

    // face predictor (momentum interpolation: the mean of the cell predictors, which carry
    // no gravity term of this step), explicit part of the face discharge, eta-system
    // coefficients a_f = theta dt (g theta dt alpha_f H_f + K_f)/dxp
    SLICELOOP1
    {
        if(Hx(i,j)>0.0)
        {
        gx0(i,j) = (b->eta(i+1,j) - b->eta(i,j))/p->DXP[IP];

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
        if(Hy(i,j)>0.0)
        {
        gy0(i,j) = (b->eta(i,j+1) - b->eta(i,j))/p->DYP[JP];

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
        qEf(f.fi,f.fj) = alf(f.fi,f.fj)*(Up(i,j) - dt*g*Hf(f.fi,f.fj)*(1.0-theta)*gf(f.fi,f.fj));
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
    
    // right-hand side of the water depth: WL^n - dt div( theta q^E + (1-theta) q^n )
    SLICELOOP4
    rhsW(i,j) = WL0(i,j) - dt*(theta*div(p,qxE,qyE) + (1.0-theta)*div(p,qxn,qyn));

    // initial guess eta^n; wet pattern: all wet (the linear system of P1), then the
    // Newton iteration on WL = max(eta + d, 0) (Casulli 2009), converging from above
    SLICELOOP4
    etaI(i,j) = eta0(i,j);

    wd_iter = psi->solve_wetdry(p,pgc,ppoissonsolv,ax,ay,rhsW,b->depth,bdiag,etaI,p->A229==1 ? 20 : 1,wd_cells);
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
    b->WL(i,j) = MAX(WL0(i,j) - p->dt*(theta*div(p,qx,qy) + (1.0-theta)*div(p,qxn,qyn)), 0.0);
    b->eta(i,j) = b->WL(i,j) - b->depth(i,j);
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
    }

    Qin = pgc->globalsum(Qin);
    Qout = pgc->globalsum(Qout);

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
    // cell momentum with the mean of the open face gradients of the cell
    const double dt = p->dt;
    double we,ww,gn,gn1,wsum;

    SLICELOOP4
    {
        if(p->wet[IJ]==1)
        {
        we = Hx(i,j)>0.0 ? 1.0 : 0.0;
        ww = Hx(i-1,j)>0.0 ? 1.0 : 0.0;
        wsum = MAX(we+ww,1.0);

        gn  = (ww*gx0(i-1,j) + we*gx0(i,j))/wsum;
        gn1 = (ww*(etaI(i,j) - etaI(i-1,j))/p->DXP[IM1] + we*(etaI(i+1,j) - etaI(i,j))/p->DXP[IP])/wsum;

        b->UH(i,j) = alc(i,j)*(UHp(i,j) - dt*g*WL0(i,j)*(theta*gn1 + (1.0-theta)*gn));

        b->VH(i,j) = 0.0;

            if(p->j_dir==1)
            {
            we = Hy(i,j)>0.0 ? 1.0 : 0.0;
            ww = Hy(i,j-1)>0.0 ? 1.0 : 0.0;
            wsum = MAX(we+ww,1.0);

            gn  = (ww*gy0(i,j-1) + we*gy0(i,j))/wsum;
            gn1 = (ww*(etaI(i,j) - etaI(i,j-1))/p->DYP[JM1] + we*(etaI(i,j+1) - etaI(i,j))/p->DYP[JP])/wsum;

            b->VH(i,j) = alc(i,j)*(VHp(i,j) - dt*g*WL0(i,j)*(theta*gn1 + (1.0-theta)*gn));
            }
        }

        else
        {
        b->UH(i,j) = 0.0;
        b->VH(i,j) = 0.0;
        }

    b->WH(i,j) = 0.0;
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

    SLICELOOP4
    vrelax += (b->WL(i,j) - rhs(i,j))*p->DXN[IP]*p->DYN[JP];

    // the face discharges follow the change of the cell discharges
    SLICELOOP4
    {
    dUH(i,j) = b->UH(i,j) - tUH(i,j);
    dVH(i,j) = (b->VH(i,j) - tVH(i,j))*p->y_dir;
    }

    exchange4(p,pgc,dUH);
    exchange4(p,pgc,dVH);

    SLICELOOP1
    if(Hx(i,j)>0.0)
    qx(i,j) += 0.5*(dUH(i,j) + dUH(i+1,j));

    if(p->j_dir==1)
    SLICELOOP2
    if(Hy(i,j)>0.0)
    qy(i,j) += 0.5*(dVH(i,j) + dVH(i,j+1));

    exchange1(p,pgc,qx);
    exchange2(p,pgc,qy);

    vrelax = pgc->globalsum(vrelax);
}

void sflow_momentum_SI::log(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double vol=0.0;

    SLICELOOP4
    vol += b->WL(i,j)*p->DXN[IP]*p->DYN[JP];

    vol = pgc->globalsum(vol);

    if(p->mpirank==0)
    silog<<p->count<<" "<<setprecision(10)<<p->simtime+p->dt<<" "<<p->dt<<" "<<setprecision(16)<<vol<<" "
         <<Qin<<" "<<Qout<<" "<<vrelax<<" "<<psi->iter<<" "<<wd_iter<<" "<<wd_cells<<"\n";
}

double sflow_momentum_SI::div(lexer *p, slice &fx, slice &fy)
{
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
