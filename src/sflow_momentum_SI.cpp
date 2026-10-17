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
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"sflow_fsf.h"
#include"sflow_reconstruct.h"
#include"sflow_diffusion.h"
#include"sflow_rheology.h"

sflow_momentum_SI::sflow_momentum_SI(lexer *p, fdm2D *b, ghostcell *pgc, sflow_HLL *pphll, sflow_signal_speed *ppss,
                                     sflow_reconstruct *pprecon, sflow_diffusion *ppdiff, sflow_pressure *pppress,
                                     solver2D *ppsolv, solver2D *pppoissonsolv, ioflow *ppflow, sflow_fsf *ppfsf,
                                     sflow_forcing *ppsfdf, sixdof *pp6dof)
                                    : sflow_momentum_func(p,b,pgc,pphll,ppss,pprecon,ppdiff,pppress,ppsolv,pppoissonsolv,ppflow,ppfsf,ppsfdf,pp6dof),
                                      qx(p),qxn(p),qxE(p),Hx(p),alx(p),ax(p),gx0(p),
                                      qy(p),qyn(p),qyE(p),Hy(p),aly(p),ay(p),gy0(p),
                                      fxU(p),fxV(p),fyU(p),fyV(p),
                                      eta0(p),WL0(p),UH0(p),VH0(p),UHa(p),VHa(p),UHb(p),VHb(p),UHp(p),VHp(p),
                                      WHz(p),dUH(p),dVH(p),etaI(p),alc(p),rhs(p),tU(p),tV(p),tUH(p),tVH(p),WLk(p)
{
    check(p);

    psi = new si_eta2D(p);

    theta = p->A227;
    g = fabs(p->W22);
    rk = (p->A210==2)?2:3;
    qini = false;

    if(p->mpirank==0)
    cout<<"SFLOW semi-implicit free surface (A 226 1): theta "<<theta<<", SSP-RK"<<rk<<" predictor"<<endl;
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

    if(p->B60>0 || p->B90>0 || p->B98>=2 || p->B99>=1)
    msg("open boundaries (B 60, B 90, B 98, B 99: discharge, waves, relaxation, water level) are");

    if(p->X10>0)
    msg("6DOF and direct forcing (X 10) are");

    if(p->A750>0)
    msg("SEASTATE coupling (A 750) is");

    if(stop==1)
    exit(1);

    if(p->A227<0.5 && p->mpirank==0)
    cout<<"SFLOW A 227 "<<p->A227<<": theta < 0.5 is unstable for large time steps"<<endl;
}

void sflow_momentum_SI::start(lexer *p, fdm2D* b, ghostcell* pgc)
{
    // the water column of the end of the previous step: sflow_eta::depth_update (main loop,
    // after the momentum step) sets eta = WL - depth for a changed bed and passes eta through
    // gcsl_start4, whose solid bed forcing (ghostcell::solid_forcing_bed) also changes fluid
    // cells behind a solid one cell thick; the volume of the semi-implicit step is kept
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

    // ghost cells, as the explicit integrators
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

    qini = true;

    // face depths upwind of the new discharges
    face_depth(p,b,pgc);
    }

    SLICELOOP1
    qxn(i,j) = qx(i,j);

    SLICELOOP2
    qyn(i,j) = qy(i,j);

    exchange1(p,pgc,qxn);
    exchange2(p,pgc,qyn);

    friction(p,b);

    // 2. predictor of the explicit terms
    predictor(p,b,pgc);

    // 3.-6. face predictor, eta system, face discharges
    eta_solve(p,b,pgc);
    face_update(p,b,pgc);

    // 7. eta, WL, wetting and drying
    eta_update(p,b,pgc);

    // 8. cell momentum, velocities, ghost cells, relaxation zones
    cell_update(p,b,pgc);

    SLICELOOP4
    {
    b->eta_n(i,j) = b->eta(i,j);
    WLk(i,j) = b->WL(i,j);
    }

    pgc->gcsl_start4(p,b->eta_n,gcval_eta);

    p->poissoniter = psi->iter;
    p->poissontime = psi->time;

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"SI eta: iter "<<psi->iter<<"  time "<<setprecision(3)<<psi->time<<endl;
}

void sflow_momentum_SI::face_depth(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // total water depth of the face: still water depth of the face plus eta upwind of the
    // face discharge (Stelling and Duinmeijer 2003), max(eta) where it is zero; a face is
    // open between two fluid cells when this depth exceeds the wetting criterion A 244
    double df,eu,h;

    SLICELOOP1
    {
    Hx(i,j) = 0.0;

        if(p->flagslice4[IJ]>0 && p->flagslice4[Ip1J]>0)
        {
        df = 0.5*(b->depth(i,j) + b->depth(i+1,j));

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
        df = 0.5*(b->depth(i,j) + b->depth(i,j+1));

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
    // (sflow_momentum_func::mirror), into the fluid behind a solid only one cell thick
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

void sflow_momentum_SI::eta_solve(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const double dt = p->dt;
    const double gt2 = g*theta*theta*dt*dt;

    // increments of the explicit terms, to the faces
    SLICELOOP4
    {
    dUH(i,j) = UHp(i,j) - UH0(i,j);
    dVH(i,j) = (VHp(i,j) - VH0(i,j))*p->y_dir;
    }

    exchange4(p,pgc,dUH);
    exchange4(p,pgc,dVH);

    // face predictor, explicit part of the face discharge, eta-system coefficients
    SLICELOOP1
    {
        if(Hx(i,j)>0.0)
        {
        gx0(i,j) = (b->eta(i+1,j) - b->eta(i,j))/p->DXP[IP];

        qxE(i,j) = alx(i,j)*(qxn(i,j) + 0.5*(dUH(i,j) + dUH(i+1,j)) - dt*g*Hx(i,j)*(1.0-theta)*gx0(i,j));

        ax(i,j) = gt2*alx(i,j)*Hx(i,j)/p->DXP[IP];
        }

        else
        {
        gx0(i,j) = 0.0;
        qxE(i,j) = 0.0;
        ax(i,j) = 0.0;
        }
    }

    if(p->j_dir==1)
    SLICELOOP2
    {
        if(Hy(i,j)>0.0)
        {
        gy0(i,j) = (b->eta(i,j+1) - b->eta(i,j))/p->DYP[JP];

        qyE(i,j) = aly(i,j)*(qyn(i,j) + 0.5*(dVH(i,j) + dVH(i,j+1)) - dt*g*Hy(i,j)*(1.0-theta)*gy0(i,j));

        ay(i,j) = gt2*aly(i,j)*Hy(i,j)/p->DYP[JP];
        }

        else
        {
        gy0(i,j) = 0.0;
        qyE(i,j) = 0.0;
        ay(i,j) = 0.0;
        }
    }

    exchange1(p,pgc,gx0);
    exchange1(p,pgc,qxE);
    exchange1(p,pgc,ax);
    exchange2(p,pgc,gy0);
    exchange2(p,pgc,qyE);
    exchange2(p,pgc,ay);

    // right-hand side, initial guess eta^n
    SLICELOOP4
    {
    rhs(i,j) = eta0(i,j) - dt*(theta*div(p,qxE,qyE) + (1.0-theta)*div(p,qxn,qyn));
    etaI(i,j) = eta0(i,j);
    }

    psi->solve(p,pgc,ppoissonsolv,ax,ay,rhs,etaI);

    exchange4(p,pgc,etaI);
}

void sflow_momentum_SI::face_update(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const double dt = p->dt;

    SLICELOOP1
    qx(i,j) = Hx(i,j)>0.0 ? qxE(i,j) - alx(i,j)*dt*g*Hx(i,j)*theta*(etaI(i+1,j) - etaI(i,j))/p->DXP[IP] : 0.0;

    if(p->j_dir==1)
    SLICELOOP2
    qy(i,j) = Hy(i,j)>0.0 ? qyE(i,j) - aly(i,j)*dt*g*Hy(i,j)*theta*(etaI(i,j+1) - etaI(i,j))/p->DYP[JP] : 0.0;

    exchange1(p,pgc,qx);
    exchange2(p,pgc,qy);
}

void sflow_momentum_SI::eta_update(lexer *p, fdm2D *b, ghostcell *pgc)
{
    starttime=pgc->timer();

    // eta from the divergence of the face discharges: the volume is conserved exactly,
    // independent of the tolerance of the eta solve
    SLICELOOP4
    b->eta(i,j) = eta0(i,j) - p->dt*(theta*div(p,qx,qy) + (1.0-theta)*div(p,qxn,qyn));

    // boundary conditions, as sflow_eta::update
    pflow->waterlevel2D(p,b,pgc,b->eta);
    pflow->eta_relax(p,pgc,b->eta);

    // ghost cells; the interior values are kept: gcsl_start4 extends the cell next to a
    // solid bed cell two cells into the solid (ghostcell::solid_forcing_bed), which reaches
    // the fluid cell behind a solid one cell thick and would change the volume
    SLICELOOP4
    rhs(i,j) = b->eta(i,j);

    pgc->gcsl_start4(p,b->eta,gcval_eta);

    SLICELOOP4
    {
    b->eta(i,j) = rhs(i,j);
    b->WL(i,j) = MAX(b->eta(i,j) + b->depth(i,j), 0.0);
    rhs(i,j) = b->WL(i,j);
    }

    pgc->gcsl_start4(p,b->WL,gcval_eta);

    SLICELOOP4
    b->WL(i,j) = rhs(i,j);

    // wetting & drying (flags, dry cells without momentum); the interior eta and WL are
    // kept, wetdry passes them through gcsl_start4 as well
    SLICELOOP4
    {
    tU(i,j) = b->eta(i,j);
    tV(i,j) = b->WL(i,j);
    }

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

    // ghost cells, relaxation zones (none with A 226 1 yet), face diagnostics; the interior is
    // kept, as in rk_stage
    SLICELOOP4
    {
    tU(i,j) = b->U(i,j);
    tV(i,j) = b->V(i,j);
    tUH(i,j) = b->UH(i,j);
    tVH(i,j) = b->VH(i,j);
    rhs(i,j) = b->WL(i,j);
    }

    stage_finish(p,b,pgc,b->UH,b->VH,b->WH,b->WL);

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
    // boundary faces stay closed
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
