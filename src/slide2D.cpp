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

#include"slide2D.h"
#include"sflow_level_file.h"
#include"lexer.h"
#include"fdm2D.h"
#include"fdm_nhf.h"
#include"matrix2D.h"
#include"vec2D.h"
#include"ghostcell.h"
#include"gcsl_yskip.h"
#include"sflow_bicgstab.h"
#include"reefckpt.h"
#include<cmath>
#include<iomanip>
#include<sys/stat.h>
#include<sys/types.h>

slide2D::slide2D(lexer *p, fdm_nhf *d, ghostcell *pgc) : slide2D(p,(fdm2D*)nullptr,pgc)
{
}

slide2D::slide2D(lexer *p, fdm2D *b, ghostcell *pgc) : h(p), qx(p), qy(p), lambda(p), zb(p),
                                                    h0(p), qx0(p), qy0(p), h1(p), qx1(p), qy1(p), Rh(p), Rx(p), Ry(p),
                                                    FhX(p), FxLX(p), FxRX(p), FyX(p), FhY(p), FxY(p), FyLY(p), FyRY(p),
                                                    hE(p), hW(p), zE(p), zW(p), hN(p), hS(p), zN(p), zS(p), lnew(p), gr(p), hstep(p), hf(p), dfil(p), wb(p), uwb(p), vwb(p), sgc(p), suc(p), svc(p), pnb(p), phis(p), tpsi(p), fsx(p), fsy(p)
{
    using namespace rheo;

    // the basal law of the slide: the keys of the SFLOW adapter of the core (W 120 - W 128)
    q.law = p->W120;
    q.tandelta = tan(deg2rad(p->W121));
    q.xi = p->W122;
    q.tz1 = tan(deg2rad(p->W123_z1));
    q.tz2 = tan(deg2rad(p->W123_z2));
    q.tz3 = tan(deg2rad(p->W123_z3));
    q.beta = p->W123_beta;
    q.L = p->W123_L;
    q.mus = tan(deg2rad(p->W124_ds));
    q.mu2 = tan(deg2rad(p->W124_d2));
    q.I0 = p->W124_I0;
    q.d = p->W124_d;
    q.phis = p->W124_phis;
    q.tauy = p->W125_tauy;
    q.K = p->W125_K;
    q.n = p->W125_n;

    g = fabs(p->W22)>0.0 ? fabs(p->W22) : 9.81;
    rhos = p->W154;
    rhow = p->W1;
    r = rhow/rhos;
    Cd = p->W155;
    dsub = MAX(p->W156,1.0e-12);
    cfl = p->W157;
    hmin = p->W158;
    nsub = 1;
    nsub_max = 100000;
    vol0 = 0.0;
    p_W162 = p->W162;

    // W 164: dilatancy (the parameters of the CFD three-phase slide: W 145, W 144 k, W 124 d)
    dil = p->W164>0.0;
    vtp_dil = dil;
    cgam = p->W164;
    etaf = p->W1*p->W2;
    K3 = p->W145_K3;
    phic = p->W145_phic;
    phi0 = p->W145_phi0;
    kperm = p->W144_k;
    dgrain = p->W124_d;

    if(dil && (phi0<=0.0 || (kperm<=0.0 && dgrain<=0.0)))
    {
        if(p->mpirank==0)
        cout<<"!!! slide layer W 164 (dilatancy): needs the initial solid fraction W 145 phi_0 > 0 and the permeability W 144 k > 0 or the grain size W 124 d > 0 !!!"<<endl;
        pgc->final(true);
    }

    if(r>=1.0)
    {
        if(p->mpirank==0)
        cout<<"!!! slide layer W 150: the bulk density of the slide W 154 ("<<rhos<<") must be larger than the water density W 1 ("<<rhow<<") !!!"<<endl;
        pgc->final(true);
    }

    if(p->mpirank==0)
    {
    const char *name[7] = {"off","Coulomb","Voellmy","Pouliquen-Forterre","mu(I)","Herschel-Bulkley","Herschel-Bulkley-Coulomb"};
    cout<<"slide layer (W 150): rho_s "<<rhos<<", rho_w "<<rhow<<", basal law "<<name[q.law>=0 && q.law<=6 ? q.law : 0];
    if(p->W126>0)
    cout<<", pore pressure "<<(p->W126==1 ? "constant" : "advected/dissipating");
    cout<<", C_d "<<Cd<<", CFL "<<cfl<<endl;
    if(dil)
    cout<<"slide layer dilatancy (W 164): shear rate "<<cgam<<" |u|/h, phi_0 "<<phi0<<", K3 "<<K3<<", phi_c "<<phic
        <<", k "<<(kperm>0.0 ? kperm : permeability(phi0))<<" m2"<<(kperm>0.0 ? "" : " (Carman-Kozeny at phi_0)")<<", eta_f "<<etaf<<" Pa s"<<endl;
    }
}

bool slide2D::vtp_dil = false;

// ---------------------------------------------------------------------------------------------
// initial slide
// ---------------------------------------------------------------------------------------------

// the water: SFLOW (fdm2D) or NHFLOW (fdm_nhf, the velocity of the bottom layer for the drag)
void slide2D::bind(lexer *p, fdm2D *b)
{
    wbed = &b->bed;
    weta = &b->eta;
    wWL = &b->WL;
    wU = &b->U;
    wV = &b->V;
    wdry = p->A244;
    nhf = false;
}

void slide2D::bind(lexer *p, fdm_nhf *d)
{
    wbed = &d->bed;
    weta = &d->eta;
    wWL = &d->WL;
    
    // the bottom layer of the water; at a shoreline the velocity of a thin column (UH/WL) can be large (the semi-implicit
    // NHFLOW at small time steps: 50 m/s at a cell of a few mm at the water line of case 32): at most |u| = 2 sqrt(g WL)
    k = 0;
    SLICELOOP4
    {
    double u = d->U[IJK];
    double v = p->j_dir==1 ? d->V[IJK] : 0.0;
    const double um = 2.0*sqrt(g*MAX(d->WL(i,j),0.0));
    const double ua = sqrt(u*u + v*v);
    
        if(ua>um)
        {
        u *= um/ua;
        v *= um/ua;
        }
    
    uwb(i,j) = u;
    vwb(i,j) = v;
    }
    
    wU = &uwb;
    wV = &vwb;
    wdry = p->A544;
    nhf = true;
}

void slide2D::ini(lexer *p, fdm2D *b, ghostcell *pgc)
{
    bind(p,b);
    ini_common(p,pgc);
}

void slide2D::ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    bind(p,d);
    ini_common(p,pgc);
}

void slide2D::start(lexer *p, fdm2D *b, ghostcell *pgc)
{
    bind(p,b);
    step(p,pgc);
}

// NHFLOW: the water column is kept, the surface moves with the bed (as sflow_eta::depth_update does for SFLOW)
void slide2D::start(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    bind(p,d);
    
    // W 163: the non-hydrostatic pressure of the bed node (explicit and semi-implicit NHFLOW: d->P = rho_w phi_0) of the
    // water step before
    if(p->W163>0 && p->A520>0)
    {
        k = 0;
        SLICELOOP4
        pnb(i,j) = p->wet[IJ]==1 ? d->P[FIJK]/p->W1 : 0.0;
        
        pgc->gcsl_start4(p,pnb,50);
    }
    
    step(p,pgc);
    
    SLICELOOP4
    {
    d->depth(i,j) = p->wd - d->bed(i,j);
    d->eta(i,j) = d->WL(i,j) - d->depth(i,j);
    }
    
    pgc->gcsl_start4(p,d->depth,50);
    pgc->gcsl_start4(p,d->eta,50);
}

void slide2D::ini_common(lexer *p, ghostcell *pgc)
{
    SLICELOOP4
    {
    zb(i,j) = (*wbed)(i,j);
    h(i,j) = 0.0;
    qx(i,j) = 0.0;
    qy(i,j) = 0.0;
    lambda(i,j) = p->W126>0 ? p->W127_lambda : 0.0;
    phis(i,j) = dil ? phi0 : 0.0;
    tpsi(i,j) = 0.0;
    }

    // W 161: the top of the initial slide from slide-level.dat (the format of sflow-level.dat, F 73), e.g. the
    // pre-event terrain over the sliding surface of a real slide
    level_grid lev;
    
    if(p->W161>0)
    lev.read(p,pgc,"slide-level.dat","W 161");
    
    SLICELOOP4
    {
        const double x = p->XP[IP];
        const double y = p->YP[JP];
        const bool oneD = (p->j_dir==0);
        double hs = 0.0;
        
        if(p->W161>0)
        {
        const double zt = lev.value(x,y);
        
            if(zt>-9000.0)
            hs = MAX(hs, zt - zb(i,j));
        }

        // W 151: the material between the bed and the level z_top in a box
        for(int n=0; n<p->W151; ++n)
        if(x>=p->W151_xs[n] && x<p->W151_xe[n] && (oneD || (y>=p->W151_ys[n] && y<p->W151_ye[n])))
        hs = MAX(hs, p->W151_z[n] - zb(i,j));

        // W 152: a layer of the thickness d over the bed in a box
        for(int n=0; n<p->W152; ++n)
        if(x>=p->W152_xs[n] && x<p->W152_xe[n] && (oneD || (y>=p->W152_ys[n] && y<p->W152_ye[n])))
        hs = MAX(hs, p->W152_d[n]);

        // W 153: a semi-ellipsoidal mound of the thickness T and the horizontal axes Lx, Ly
        for(int n=0; n<p->W153; ++n)
        {
        const double ax = (x - p->W153_xc[n])/(0.5*p->W153_Lx[n]);
        const double ay = oneD ? 0.0 : (y - p->W153_yc[n])/(0.5*p->W153_Ly[n]);
        const double rr = ax*ax + ay*ay;

        if(rr<1.0)
        hs = MAX(hs, p->W153_T[n]*sqrt(1.0 - rr));
        }

        h(i,j) = MAX(hs,0.0);

        // W 160: the initial velocity of the slide
        if(h(i,j)>hmin)
        {
        qx(i,j) = h(i,j)*p->W160_u;
        qy(i,j) = h(i,j)*p->W160_v*p->y_dir;
        }
    }

    pgc->gcsl_start4(p,zb,50);
    ghosts(p,pgc,h,qx,qy);
    pgc->gcsl_start4(p,lambda,50);

    // W 164: tan(psi) of the initial state at rest (J = 0)
    if(dil)
    {
        SLICELOOP4
        tpsi(i,j) = h(i,j)>hmin ? MIN(MAX(K3*(phis(i,j) - phic),-0.6),0.6) : 0.0;
        
        pgc->gcsl_start4(p,phis,50);
        pgc->gcsl_start4(p,tpsi,50);
    }

    // W 162: the diffusivity of the filtered transfer from the still water depth over the bed
    if(p->W162>0.0)
    {
        // the still water level: F 60 (F 62); p->wd is set later (sflow_f::ini_fsf)
        const double wl = p->F60>-1.0e19 ? p->F60 : (p->F62>-1.0e19 ? p->F62 : p->wd);
        
        SLICELOOP4
        {
        const double d = MAX(wl - zb(i,j), 0.0);
        dfil(i,j) = p->W162*d*d;
        hf(i,j) = h(i,j);
        }
        
        pgc->gcsl_start4(p,dfil,50);
        pfilt = new sflow_bicgstab(p,pgc);
        Mf = new matrix2D(p);
        xf = new vec2D(p);
        rf = new vec2D(p);
    }

    bed_update(p,pgc);

    vol0 = volume_slide(p,pgc);

    if(p->mpirank==0)
    {
    cout<<"slide layer: initial volume "<<setprecision(8)<<vol0<<" m^3"<<(p->j_dir==0 ? " (per m width)" : "")<<endl;

    mkdir(nhf ? "./REEF3D_NHFLOW_Log" : "./REEF3D_SFLOW_Log",0777);
    // a restart from a checkpoint (I 45) continues the log from its size at the checkpoint
    logfresh = reefckpt_open(logout, nhf ? "./REEF3D_NHFLOW_Log/REEF3D_NHFLOW_slide.dat" : "./REEF3D_SFLOW_Log/REEF3D_SFLOW_slide.dat", p->I45);
    if(logfresh)
    logout<<"# time  dt  substeps  V_slide  dV_slide/V0  V_water  h_max  |u|_max  x_min  x_max  y_min  y_max"<<(dil ? "  phi_min  phi_max  tanpsi_min  tanpsi_max" : "")
          <<"  (front: h > max(10 W 158, 1 mm))"<<endl;
    }

    logfresh = pgc->globalimax((p->mpirank==0 && logfresh) ? 1 : 0)==1;
    
    if(logfresh)
    log(p,pgc,0);
}

// ---------------------------------------------------------------------------------------------
// one water step
// ---------------------------------------------------------------------------------------------

void slide2D::step(lexer *p, ghostcell *pgc)
{
    // sub-steps: the slide CFL inside the water step
    const double dts_max = timestep(p,pgc);
    nsub = int(ceil(p->dt/dts_max - 1.0e-9));
    nsub = MAX(nsub,1);

    if(nsub>nsub_max)
    {
        if(p->mpirank==0)
        cout<<"slide layer: "<<nsub<<" sub-steps needed, limited to "<<nsub_max<<endl;
        nsub = nsub_max;
    }

    const double dts = p->dt/double(nsub);
    // the reduced gravity with the water state of the end of the water step, ghost cells included; the slide depth of
    // the start of the step (the water surface over the sub-steps, sources)
    for(i=-1; i<=p->knox; ++i)
    for(j=-1; j<=p->knoy; ++j)
    {
    gr(i,j) = gred(p,i,j);
    hstep(i,j) = h(i,j);
    }

    if(dil)
    for(i=-1; i<=p->knox; ++i)
    for(j=-1; j<=p->knoy; ++j)
    fsx(i,j) = fsy(i,j) = 0.0;

    for(int s=0; s<nsub; ++s)
    {
        // SSP-RK2, stage 1
        SLICELOOP4
        {
        h0(i,j) = h(i,j);
        qx0(i,j) = qx(i,j);
        qy0(i,j) = qy(i,j);
        }

        rhs(p,h,qx,qy);
        flux_sum(p,0.5*dts);

        SLICELOOP4
        {
        h1(i,j) = h0(i,j) + dts*Rh(i,j);
        qx1(i,j) = qx0(i,j) + dts*Rx(i,j);
        qy1(i,j) = (qy0(i,j) + dts*Ry(i,j))*p->y_dir;
        }

        positivity(p,h1,qx1,qy1);
        friction(p,h1,qx1,qy1,dts);
        ghosts(p,pgc,h1,qx1,qy1);

        // stage 2
        rhs(p,h1,qx1,qy1);
        flux_sum(p,0.5*dts);

        SLICELOOP4
        {
        h(i,j) = 0.5*h0(i,j) + 0.5*(h1(i,j) + dts*Rh(i,j));
        qx(i,j) = 0.5*qx0(i,j) + 0.5*(qx1(i,j) + dts*Rx(i,j));
        qy(i,j) = (0.5*qy0(i,j) + 0.5*(qy1(i,j) + dts*Ry(i,j)))*p->y_dir;
        }

        positivity(p,h,qx,qy);
        friction(p,h,qx,qy,0.5*dts);
        ghosts(p,pgc,h,qx,qy);
    }

    if(p->W126==2)
    pore_update(p,pgc);

    if(dil)
    dilatancy_update(p,pgc);

    // the bed of the water: z_b + h (sflow_eta::depth_update then moves the surface with the bed)
    bed_update(p,pgc);
    log(p,pgc,1);
}

double slide2D::dt_cap(lexer *p, ghostcell *pgc)
{
    return double(MAX(p->W159,1))*timestep(p,pgc);
}

double slide2D::timestep(lexer *p, ghostcell *pgc)
{
    double dtm = 1.0e20;

    SLICELOOP4
    if(h(i,j)>hmin)
    {
    const double c = sqrt(g*h(i,j));
    const double u = fabs(qx(i,j))/h(i,j);
    const double v = fabs(qy(i,j))/h(i,j);

    dtm = MIN(dtm, p->DXN[IP]/(u + c));

    if(p->j_dir==1)
    dtm = MIN(dtm, p->DYN[JP]/(v + c));
    }

    dtm = pgc->globalmin(dtm);

    return cfl*dtm;
}

// ---------------------------------------------------------------------------------------------
// ghost cells: zero gradient, walls at the solids and the domain boundaries (normal discharge mirrored)
// ---------------------------------------------------------------------------------------------

void slide2D::ghosts(lexer *p, ghostcell *pgc, slice &H, slice &QX, slice &QY)
{
    pgc->gcsl_start4(p,H,50);
    pgc->gcsl_start4(p,QX,50);
    pgc->gcsl_start4(p,QY,50);

    const bool ysk = gcsl_yskip(p) || p->j_dir==0;

    auto fl = [p](int ii, int jj){ return p->flagslice4[(ii-p->imin)*p->jmax + (jj-p->jmin)]; };

    for(int n=0; n<p->gcbsl4_count; ++n)
    {
        const int ii = p->gcbsl4[n][0];
        const int jj = p->gcbsl4[n][1];
        const int cs = p->gcbsl4[n][3];
        const int d  = p->gcbsl4[n][2];
        const int mg = (d>0 && d<3) ? d : 3;

        if(ysk && (cs==2 || cs==3))
        continue;

        const int di = (cs==1) ? -1 : (cs==4 ? 1 : 0);
        const int dj = (cs==3) ? -1 : (cs==2 ? 1 : 0);

        for(int qq=1; qq<=mg; ++qq)
        {
            // the mirror image of the ghost cell, the boundary cell when the mirror is not a fluid cell
            int mi = ii - di*(qq-1);
            int mj = jj - dj*(qq-1);

            if(fl(mi,mj)<=0)
            {
            mi = ii;
            mj = jj;
            }

            const int gi = ii + di*qq;
            const int gj = jj + dj*qq;

            H(gi,gj) = H(mi,mj);

            if(cs==1 || cs==4)
            {
            QX(gi,gj) = -QX(mi,mj);
            QY(gi,gj) = QY(mi,mj);
            }

            if(cs==2 || cs==3)
            {
            QX(gi,gj) = QX(mi,mj);
            QY(gi,gj) = -QY(mi,mj);
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// the bed of the water
// ---------------------------------------------------------------------------------------------

void slide2D::bed_update(lexer *p, ghostcell *pgc)
{
    // the bed of the step before: wb for the bed condition of the non-hydrostatic water
    SLICELOOP4
    wb(i,j) = (*wbed)(i,j);
    
    if(p->W162>0.0)
    {
        filter(p,pgc);
        
        SLICELOOP4
        (*wbed)(i,j) = zb(i,j) + hf(i,j);
    }
    
    else
    SLICELOOP4
    (*wbed)(i,j) = zb(i,j) + h(i,j);
    
    SLICELOOP4
    wb(i,j) = p->dt>0.0 && p->count>0 ? ((*wbed)(i,j) - wb(i,j))/p->dt : 0.0;

    pgc->gcsl_start4(p,(*wbed),50);
}

// W 162: h_f - div(D grad h_f) = h, D at the faces the smaller of the two cells (no filter across a shoreline),
// no flux through walls and the outer boundaries; conservative (the faces are symmetric)
void slide2D::filter(lexer *p, ghostcell *pgc)
{
    int n=0;
    
    SLICELOOP4
    {
        const double ce = p->flagslice4[Ip1J]>0 ? MIN(dfil(i,j),dfil(i+1,j))/(p->DXP[IP]*p->DXN[IP]) : 0.0;
        const double cw = p->flagslice4[Im1J]>0 ? MIN(dfil(i,j),dfil(i-1,j))/(p->DXP[IM1]*p->DXN[IP]) : 0.0;
        const double cn = p->y_dir==1 && p->flagslice4[IJp1]>0 ? MIN(dfil(i,j),dfil(i,j+1))/(p->DYP[JP]*p->DYN[JP]) : 0.0;
        const double cs = p->y_dir==1 && p->flagslice4[IJm1]>0 ? MIN(dfil(i,j),dfil(i,j-1))/(p->DYP[JM1]*p->DYN[JP]) : 0.0;
        
        Mf->p[n] = 1.0 + ce + cw + cn + cs;
        Mf->n[n] = -ce;
        Mf->s[n] = -cw;
        Mf->w[n] = -cn;
        Mf->e[n] = -cs;
        (*rf).V[n] = h(i,j);
        ++n;
    }
    
    // the filter is stiff in deep water (D/dx^2 up to 1e5): more iterations than the pressure (N 46), from the
    // filtered slide of the step before
    const int it = p->solveriter;
    const int n46 = p->N46;
    p->N46 = MAX(n46,20000);
    pfilt->start(p,pgc,hf,*Mf,*xf,*rf,4);
    
    if(p->solveriter>=p->N46 && p->mpirank==0)
    cout<<"slide layer W 162: the filter did not converge in "<<p->N46<<" iterations"<<endl;
    
    p->N46 = n46;
    p->solveriter = it;
    
    SLICELOOP4
    hf(i,j) = MAX(hf(i,j),0.0);
    
    pgc->gcsl_start4(p,hf,50);
}

// ---------------------------------------------------------------------------------------------
// friction and drag: backward-Euler steps on the stage discharges
// ---------------------------------------------------------------------------------------------

double slide2D::gn_cell(lexer *p)
{
    // gravity normal to the bed; W 128: cos^2(theta) of the slope of the bed under the slide
    double gg = g;

    if(p->W128==1)
    {
    const double bx = (zb(i+1,j) - zb(i-1,j))/(p->DXP[IM1] + p->DXP[IP]);
    const double by = p->j_dir==1 ? (zb(i,j+1) - zb(i,j-1))/(p->DYP[JM1] + p->DYP[JP]) : 0.0;
    gg /= 1.0 + bx*bx + by*by;
    }

    return gg;
}

double slide2D::sub_cell(lexer *p, int ii, int jj)
{
    // the submerged fraction: 0 subaerial, 1 under more than W 156 of water
    const double wl = (*wWL)(ii,jj);

    if(p->wet[(ii-p->imin)*p->jmax + (jj-p->jmin)]==0 || wl<=wdry)
    return 0.0;

    return MIN(wl/dsub,1.0);
}

double slide2D::gred(lexer *p, int ii, int jj)
{
    return g*(1.0 - r*sub_cell(p,ii,jj));
}

void slide2D::friction(lexer *p, slice &H, slice &QX, slice &QY, double dt)
{
    rheo::basal_state s;
    s.rho = rhos;
    s.g = g;

    SLICELOOP4
    {
        const double hh = H(i,j);

        if(hh<=hmin)
        {
        QX(i,j) = 0.0;
        QY(i,j) = 0.0;
        continue;
        }

        double u = QX(i,j)/hh;
        double v = p->j_dir==1 ? QY(i,j)/hh : 0.0;
        const double sb = sub_cell(p,i,j);

        // interface drag: relaxation to the water velocity, rho_w C_d |du| / (rho_s h)
        if(Cd>0.0 && sb>0.0)
        {
        const double du = (*wU)(i,j) - u;
        const double dv = p->j_dir==1 ? (*wV)(i,j) - v : 0.0;
        const double k = sb*r*Cd*sqrt(du*du + dv*dv)/hh;

        u = (u + dt*k*(*wU)(i,j))/(1.0 + dt*k);
        v = p->j_dir==1 ? (v + dt*k*(*wV)(i,j))/(1.0 + dt*k) : 0.0;
        }

        const double Ustar = sqrt(u*u + v*v);

        if(Ustar<=0.0 || q.law==rheo::BASAL_OFF)
        {
        QX(i,j) = hh*u;
        QY(i,j) = hh*v;
        continue;
        }

        // basal friction with the buoyant normal stress (rho_s - sub rho_w) g_n h
        s.h = hh;
        s.gn = gn_cell(p);
        s.sigma_n = (rhos - sb*rhow)*s.gn*hh;
        s.lambda = 0.0;

        if(p->W126==1)
        s.lambda = lambda(i,j);

        // advected ratio: at the middle of the water step (lambda is updated once per step)
        if(p->W126==2)
        s.lambda = rheo::pore_relax(lambda(i,j),p->W127_lambdah,p->W127_D,hh,0.5*p->dt);

        // W 164: dilatancy. tan(psi) at the speed U0 of the start of the sub-step (h0, qx0, qy0; not the trial speed of
        // the stage, which carries the driving of the whole step without the friction and held a dense layer at rest),
        // consistent with the excess pore pressure p_e = - C1 tan(psi) U0 it drives; a contracting layer (tan(psi) < 0): p_e > 0 at that speed into the pore
        // pressure ratio (sigma_n' >= 0.02 sigma_n); tan(psi) added to the friction coefficient mu_c of the law there,
        // as the ratio 1 - (1 - lambda)(1 + tan(psi)/mu_c); a dilating layer (tan(psi) > 0): the suction - p_e is linear in
        // the speed, the friction (mu_c + tan(psi)) C1 tan(psi) U a linear drag, backward Euler after the step of the law
        // (together exact backward Euler for Coulomb)
        double drag = 0.0;

        if(dil)
        {
            double C1;
            const double U0 = h0(i,j)>hmin ? sqrt(qx0(i,j)*qx0(i,j) + (p->j_dir==1 ? qy0(i,j)*qy0(i,j) : 0.0))/h0(i,j) : 0.0;
            const double tp = dil_tpsi(hh,U0,sb,s.sigma_n,s.lambda,phis(i,j),C1);
            double lt = s.lambda;

            if(tp<0.0 && s.sigma_n>0.0)
            lt = MAX(lt,MIN(lt - C1*U0*tp/s.sigma_n,0.98));

            s.lambda = lt;
            const double se = MAX(1.0 - lt,0.0)*s.sigma_n;

            if(se>0.0)
            {
            const double muc = rheo::basal_stress(q,s,Ustar)/se;
            const double fac = muc>0.0 ? MAX(1.0 + tp/muc,0.0) : 0.0;

            s.lambda = 1.0 - (1.0 - lt)*fac;

            if(tp>0.0)
            drag = MAX(muc + tp,0.0)*C1*tp/(rhos*hh);
            }
        }

        double sp = rheo::implicit_speed(q,s,Ustar,dt);

        if(drag>0.0)
        sp /= 1.0 + dt*drag;

        const double f = sp/Ustar;

        QX(i,j) = hh*u*f;
        QY(i,j) = hh*v*f;
    }
}

void slide2D::positivity(lexer *p, slice &H, slice &QX, slice &QY)
{
    // round-off below zero (the hydrostatic reconstruction keeps the depth positive under the CFL)
    SLICELOOP4
    if(H(i,j)<0.0)
    {
    H(i,j) = 0.0;
    QX(i,j) = 0.0;
    QY(i,j) = 0.0;
    }
}

// ---------------------------------------------------------------------------------------------
// pore pressure: first-order upwind advection with the slide velocity, then relaxation (as W 126 2 in SFLOW)
// ---------------------------------------------------------------------------------------------

void slide2D::pore_update(lexer *p, ghostcell *pgc)
{
    pgc->gcsl_start4(p,lambda,50);

    SLICELOOP4
    {
        lnew(i,j) = lambda(i,j);
        const double hh = h(i,j);

        if(hh<=hmin)
        continue;

        const double u = qx(i,j)/hh;
        const double v = p->j_dir==1 ? qy(i,j)/hh : 0.0;
        double adv = 0.0;

        if(u>0.0)
        adv += u*(lambda(i,j) - lambda(i-1,j))/p->DXP[IM1];
        else
        adv += u*(lambda(i+1,j) - lambda(i,j))/p->DXP[IP];

        if(p->j_dir==1)
        {
        if(v>0.0)
        adv += v*(lambda(i,j) - lambda(i,j-1))/p->DYP[JM1];
        else
        adv += v*(lambda(i,j+1) - lambda(i,j))/p->DYP[JP];
        }

        lnew(i,j) = rheo::pore_relax(lambda(i,j) - p->dt*adv,p->W127_lambdah,p->W127_D,hh,p->dt);
    }

    SLICELOOP4
    lambda(i,j) = lnew(i,j);

    pgc->gcsl_start4(p,lambda,50);
}

// ---------------------------------------------------------------------------------------------
// dilatancy (W 164)
// ---------------------------------------------------------------------------------------------

double slide2D::permeability(double ph) const
{
    if(kperm>0.0)
    return kperm;

    // Carman-Kozeny (as the CFD three-phase slide)
    ph = MIN(MAX(ph,0.3),0.7);
    return dgrain*dgrain*(1.0-ph)*(1.0-ph)*(1.0-ph)/(180.0*ph*ph);
}

// tan(psi) = K3 (phi - phi_c/(1 + J^1/2)) with J = sub eta_f gamma/p', gamma = W 164 U/h, and the effective stress that
// includes the excess pore pressure it drives, p' = max((1 - lambda) sigma_n + C1 U tan(psi), 0.02 sigma_n): f(t) = t - K3
// (phi - phi_eq(J(t))) increases with t, bisection in [-0.6, 0.6]. C1 = sub eta_f W 164 h/(2 k): the Darcy flow of the
// dilation gamma tan(psi) through the layer, drawn in (expelled) through its top, base impermeable, p_e = - C1 U tan(psi)
double slide2D::dil_tpsi(double hh, double U, double sb, double sig0, double lam, double ph, double &C1) const
{
    const double ef = sb*etaf;
    const double gam = cgam*U/MAX(hh,1.0e-12);

    C1 = ef*cgam*hh/(2.0*permeability(ph));

    if(ef*gam<=0.0 || sig0<=0.0)
    return MIN(MAX(K3*(ph - phic),-0.6),0.6);

    const double C = C1*U;
    const double s0 = MAX(1.0 - lam,0.0)*sig0;
    const double pmin = 0.02*sig0;

    auto f = [&](double t)
    {
        const double J = rheo::viscous_number(gam,MAX(s0 + C*t,pmin),ef);
        return t - rheo::tan_psi(K3,ph,phic,J);
    };

    double a = -0.6, b = 0.6;

    if(f(a)>=0.0)
    return a;

    if(f(b)<=0.0)
    return b;

    for(int n=0; n<60 && b-a>1.0e-9; ++n)
    {
    const double m = 0.5*(a + b);

        if(f(m)>0.0)
        b = m;
        else
        a = m;
    }

    return 0.5*(a + b);
}

// W 164: the mass fluxes of the faces over the water step, SSP-RK2 weights (the transport of the solid volume)
void slide2D::flux_sum(lexer *p, double w)
{
    if(!dil)
    return;

    for(i=-1; i<p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    fsx(i,j) += w*FhX(i,j);

    if(p->j_dir==1)
    for(i=0; i<p->knox; ++i)
    for(j=-1; j<p->knoy; ++j)
    fsy(i,j) += w*FhY(i,j);
}

// the solid fraction: the solid volume phi h carried with the mass fluxes of the step (upwind phi; a uniform phi stays
// uniform, the phi of a thin film enters a thicker cell with its volume), then d(phi)/dt = - phi gamma tan(psi) over the
// water step (exponential), tan(psi) of the end of the step; 0.3 <= phi <= 0.75. The advective form carried the phi of
// the thin film the layer leaves behind (dilated to 0.3) into the layer, fluidised it there and piled it up
void slide2D::dilatancy_update(lexer *p, ghostcell *pgc)
{
    pgc->gcsl_start4(p,phis,50);

    auto up = [&](double F, double pl, double pr){ return F>0.0 ? F*pl : F*pr; };

    SLICELOOP4
    {
        const double hh = h(i,j);
        lnew(i,j) = phis(i,j);

        if(hh<=hmin)
        {
        tpsi(i,j) = 0.0;
        continue;
        }

        double vs = phis(i,j)*hstep(i,j)
                  - (up(fsx(i,j),phis(i,j),phis(i+1,j)) - up(fsx(i-1,j),phis(i-1,j),phis(i,j)))/p->DXN[IP];

        if(p->j_dir==1)
        vs -= (up(fsy(i,j),phis(i,j),phis(i,j+1)) - up(fsy(i,j-1),phis(i,j-1),phis(i,j)))/p->DYN[JP];

        const double ph = MIN(MAX(vs/hh,0.3),0.75);

        const double U = sqrt(qx(i,j)*qx(i,j) + (p->j_dir==1 ? qy(i,j)*qy(i,j) : 0.0))/hh;
        const double sb = sub_cell(p,i,j);
        const double sig0 = (rhos - sb*rhow)*gn_cell(p)*hh;
        double C1;
        const double tp = dil_tpsi(hh,U,sb,sig0,p->W126>0 ? lambda(i,j) : 0.0,ph,C1);

        tpsi(i,j) = tp;
        lnew(i,j) = MIN(MAX(ph*exp(-p->dt*cgam*U/hh*tp),0.3),0.75);
    }

    SLICELOOP4
    phis(i,j) = lnew(i,j);

    pgc->gcsl_start4(p,phis,50);
    pgc->gcsl_start4(p,tpsi,50);
}

// ---------------------------------------------------------------------------------------------
// drag of the slide on the water: tau_i/rho_w = - C_d |u_w - u_s| (u_w - u_s), linearised implicit factor
// ---------------------------------------------------------------------------------------------

// NHFLOW (W 155): tau_i/rho_w = C_d |u_w - u_s| (u_w - u_s) on the bottom layer of the water (thickness dz), with the
// velocity of that layer; point-implicit as the bed friction (nhflow_bcmom)
void slide2D::drag_nhflow(lexer *p, fdm_nhf *d, int comp, double *VEL, double *F, slice &WL)
{
    if(Cd<=0.0)
    return;
    
    k = 0;
    
    SLICELOOP4
    if(p->DF[IJK]>0 && p->wet[IJ]==1 && WL(i,j)>wdry && h(i,j)>hmin)
    {
    const double us = qx(i,j)/h(i,j);
    const double vs = p->j_dir==1 ? qy(i,j)/h(i,j) : 0.0;
    const double du = d->U[IJK] - us;
    const double dv = p->j_dir==1 ? d->V[IJK] - vs : 0.0;
    const double dz = p->DZN[KP]*WL(i,j);
    const double c = sub_cell(p,i,j)*Cd*sqrt(du*du + dv*dv)/dz;
    
    F[IJK] -= c*(VEL[IJK] - (comp==0 ? us : vs))*WL(i,j)/(1.0 + p->dt*c);
    }
}

void slide2D::si_drag_reset(lexer *p)
{
    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    sgc(i,j) = suc(i,j) = svc(i,j) = 0.0;
}

double slide2D::si_drag(lexer *p, int ii, int jj, double uw, double vw, double dz)
{
    // the thin edges of the slide (below max(10 W 158, 1 mm)) carry no drag: the velocity of a column of a few W 158 is
    // not reliable, and with the drag rate of a thin water layer at the shoreline it set the water there
    const double hdrag = MAX(10.0*hmin,1.0e-3);
    
    if(h(ii,jj)<=hdrag || dz<=1.0e-12)
    return 0.0;
    
    const double us = qx(ii,jj)/h(ii,jj);
    const double vs = p->j_dir==1 ? qy(ii,jj)/h(ii,jj) : 0.0;
    const double du = uw - us, dv = (vw - vs)*p->y_dir;
    const double c = sub_cell(p,ii,jj)*Cd*sqrt(du*du + dv*dv)/dz;
    
    sgc(ii,jj) = c;
    suc(ii,jj) = us;
    svc(ii,jj) = vs;
    
    return c;
}

void slide2D::u_source(lexer *p, fdm2D *b, slice &U, slice &V, slice &WL)
{
    if(Cd<=0.0)
    return;

    SLICELOOP4
    if(p->wet[IJ]==1 && WL(i,j)>p->A244 && h(i,j)>hmin)
    {
    const double us = qx(i,j)/h(i,j);
    const double vs = p->j_dir==1 ? qy(i,j)/h(i,j) : 0.0;
    const double du = U(i,j) - us;
    const double dv = p->j_dir==1 ? V(i,j) - vs : 0.0;
    const double c = sub_cell(p,i,j)*Cd*sqrt(du*du + dv*dv);

    b->F(i,j) -= c*du/(1.0 + p->dt*c/WL(i,j));
    }
}

void slide2D::v_source(lexer *p, fdm2D *b, slice &U, slice &V, slice &WL)
{
    if(Cd<=0.0 || p->j_dir==0)
    return;

    SLICELOOP4
    if(p->wet[IJ]==1 && WL(i,j)>p->A244 && h(i,j)>hmin)
    {
    const double us = qx(i,j)/h(i,j);
    const double vs = qy(i,j)/h(i,j);
    const double du = U(i,j) - us;
    const double dv = V(i,j) - vs;
    const double c = sub_cell(p,i,j)*Cd*sqrt(du*du + dv*dv);

    b->G(i,j) -= c*dv/(1.0 + p->dt*c/WL(i,j));
    }
}

// ---------------------------------------------------------------------------------------------
// volumes and log
// ---------------------------------------------------------------------------------------------

double slide2D::volume_slide(lexer *p, ghostcell *pgc)
{
    double v = 0.0;

    SLICELOOP4
    v += h(i,j)*p->DXN[IP]*(p->j_dir==1 ? p->DYN[JP] : 1.0);

    return pgc->globalsum(v);
}

double slide2D::volume_water(lexer *p, ghostcell *pgc)
{
    double v = 0.0;

    SLICELOOP4
    v += (*wWL)(i,j)*p->DXN[IP]*(p->j_dir==1 ? p->DYN[JP] : 1.0);

    return pgc->globalsum(v);
}

void slide2D::log(lexer *p, ghostcell *pgc, int mode)
{
    const double vs = volume_slide(p,pgc);
    const double vw = volume_water(p,pgc);
    const double hf = MAX(10.0*hmin,1.0e-3);

    double hmax=0.0, umax=0.0;
    double xmin=1.0e20, xmax=-1.0e20, ymin=1.0e20, ymax=-1.0e20;
    double fmin=1.0e20, fmax=-1.0e20, tmin=1.0e20, tmax=-1.0e20;

    SLICELOOP4
    {
    hmax = MAX(hmax,h(i,j));

        if(dil && h(i,j)>hf)
        {
        fmin = MIN(fmin,phis(i,j));
        fmax = MAX(fmax,phis(i,j));
        tmin = MIN(tmin,tpsi(i,j));
        tmax = MAX(tmax,tpsi(i,j));
        }

    if(h(i,j)>hmin)
    umax = MAX(umax, sqrt(qx(i,j)*qx(i,j) + qy(i,j)*qy(i,j))/h(i,j));

        if(h(i,j)>hf)
        {
        xmin = MIN(xmin,p->XP[IP]);
        xmax = MAX(xmax,p->XP[IP]);
        ymin = MIN(ymin,p->YP[JP]);
        ymax = MAX(ymax,p->YP[JP]);
        }
    }

    hmax = pgc->globalmax(hmax);
    umax = pgc->globalmax(umax);
    xmin = pgc->globalmin(xmin);
    xmax = pgc->globalmax(xmax);
    ymin = pgc->globalmin(ymin);
    ymax = pgc->globalmax(ymax);

    if(dil)
    {
    fmin = pgc->globalmin(fmin);
    fmax = pgc->globalmax(fmax);
    tmin = pgc->globalmin(tmin);
    tmax = pgc->globalmax(tmax);
    }

    if(p->mpirank==0)
    {
    const double t = (mode==0) ? p->simtime : p->simtime + p->dt;     // start: called before simtime += dt

    logout<<setprecision(10)<<t<<" "<<p->dt<<" "<<(mode==0 ? 0 : nsub)<<" "<<vs<<" "<<(vol0>0.0 ? (vs-vol0)/vol0 : 0.0)<<" "<<vw
          <<" "<<hmax<<" "<<umax<<" "<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax;

    if(dil)
    logout<<" "<<fmin<<" "<<fmax<<" "<<tmin<<" "<<tmax;

    logout<<endl;

    if(mode==1 && (p->count%p->P12==0))
    cout<<"slide: substeps "<<nsub<<"  V "<<setprecision(8)<<vs<<"  h_max "<<setprecision(4)<<hmax<<"  |u|_max "<<umax<<"  front "<<xmax<<endl;
    }
}

// ---------------------------------------------------------------------------------------------
// output
// ---------------------------------------------------------------------------------------------

void slide2D::offset_vtp(lexer *p, int *offset, int &n)
{
    // slide_h, slide_velocity (3), slide_lambda; W 164: slide_phi, slide_tanpsi
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
    ++n;
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)*3+sizeof(int);
    ++n;
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
    ++n;

    if(dil)
    for(int kk=0; kk<2; ++kk)
    {
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
    ++n;
    }
}

void slide2D::name_vtp(std::ostream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"slide_h\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"slide_velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"slide_lambda\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;

    if(dil)
    {
    result<<"<DataArray type=\"Float32\" Name=\"slide_phi\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"slide_tanpsi\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    }
}

void slide2D::name_pvtp(std::ostream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"slide_h\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"slide_velocity\" NumberOfComponents=\"3\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"slide_lambda\"/>\n";

    if(vtp_dil)
    {
    result<<"<PDataArray type=\"Float32\" Name=\"slide_phi\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"slide_tanpsi\"/>\n";
    }
}

void slide2D::print_vtp(lexer *p, std::ostream &result)
{
    int iin;
    float ffn;

    // the cell velocities (Rx, Ry are free between the steps), ghost cells included for the interpolation
    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    {
    Rx(i,j) = h(i,j)>hmin ? qx(i,j)/h(i,j) : 0.0;
    Ry(i,j) = h(i,j)>hmin ? qy(i,j)/h(i,j) : 0.0;
    }

    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
    ffn=float(p->sl_ipol4(h));
    result.write((char*)&ffn, sizeof(float));
    }

    iin=sizeof(float)*p->pointnum2D*3;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
    ffn=float(p->sl_ipol4(Rx));
    result.write((char*)&ffn, sizeof(float));
    ffn=float(p->sl_ipol4(Ry));
    result.write((char*)&ffn, sizeof(float));
    ffn=0.0f;
    result.write((char*)&ffn, sizeof(float));
    }

    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
    ffn=float(p->sl_ipol4(lambda));
    result.write((char*)&ffn, sizeof(float));
    }

    if(dil)
    {
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
    ffn=float(p->sl_ipol4(phis));
    result.write((char*)&ffn, sizeof(float));
    }

    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
    ffn=float(p->sl_ipol4(tpsi));
    result.write((char*)&ffn, sizeof(float));
    }
    }
}

void slide2D::ckpt_fields(const std::function<void(const char*, char, double*)> &fd, const std::function<void(const char*, char, int*)> &fi)
{
    fd("slide_h",'S',h.V);
    fd("slide_qx",'S',qx.V);
    fd("slide_qy",'S',qy.V);
    fd("slide_lambda",'S',lambda.V);
    fd("slide_zb",'S',zb.V);
    
    // W 162: the filtered slide (the start of the filter solve of the next step)
    if(p_W162>0.0)
    fd("slide_hf",'S',hf.V);

    // W 164: the solid fraction and tan(psi)
    if(dil)
    {
    fd("slide_phi",'S',phis.V);
    fd("slide_tanpsi",'S',tpsi.V);
    }
}
