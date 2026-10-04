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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"CPM.h"
#include"part.h"
#include"lexer.h"
#include"fdm.h"
#include"sediment_fdm.h"
#include<algorithm>
#include<vector>
#include"ghostcell.h"

// MP-PIC parcel acceleration
//   dUp/dt = Dp (Uf - Up) - grad(p)/rho_p + g - grad(Ps)/(theta rho_p)
// returns the explicit part F,G,H = Dp Uf - grad(p)/rho_p + g - grad(Ps)/(theta rho_p) (+ solid forcing)
// and the implicit drag coefficient Dpx (=Dpy=Dpz), velocity update: Up = (Up + dt F)/(1 + dt Dp)
void CPM::advec_mppic(lexer *p, fdm *a, part &P, sediment_fdm *s, turbulence *pturb,
                        double *PX, double *PY, double *PZ, double *PU, double *PV, double *PW,
                        double &F, double &G, double &H, double dt)
{
    // fluid pressure gradient
    dPx_val = p->ccipol4a(dPx,PX[n],PY[n],PZ[n]);
    dPy_val = p->j_dir==1 ? p->ccipol4a(dPy,PX[n],PY[n],PZ[n]) : 0.0;
    dPz_val = p->ccipol4a(dPz,PX[n],PY[n],PZ[n]);

    // gravity
    Bx = p->W20;
    By = p->W21;
    Bz = p->W22;

    // fluid velocity
    uf = p->ccipol1c(a->u,PX[n],PY[n],PZ[n]);
    vf = p->j_dir==1 ? p->ccipol2c(a->v,PX[n],PY[n],PZ[n]) : 0.0;
    wf = p->ccipol3c(a->w,PX[n],PY[n],PZ[n]);
    
    // S 10 1: the bed is a solid boundary for the fluid
    liftx=lifty=liftz=0.0;
    
    if(p->S10!=2)
    {
        shelter = shelter_factor(p,s,PX[n],PY[n],PZ[n],PU[n],PV[n],PW[n]);
        expo_w = expo.size()==P.index ? expo[n] : -1.0;

        // sub-grid bedload layer (Q 58): the flow moves the bed through the layer only, the exposed
        // grains are fully sheltered (fluid at rest at the grain level, drag against it, no lift)
        if(p->Q58>0)
        shelter = 0.0;
        nearbed_velocity(p,a,PX[n],PY[n],PZ[n],P.D[n]);
        shelter = 1.0;
        expo_w = -1.0;

        // sub-grid bedload layer (Q 58): grains resting on the bed within one cell above the bed
        // level are part of the bed surface, they are moved by the layer only (fluid at rest at the
        // grain level); moving grains (suspension, settling, avalanches) see the resolved flow
        // reduced towards the bed, u delta/h
        if(p->Q58>0)
        {
            double topo = ptopo(p,a,PX[n],PY[n],PZ[n]);

            if(topo>=0.0)
            {
                i = p->posc_i(PX[n]);
                j = p->posc_j(PY[n]);
                k = p->posc_k(PZ[n]);

                double h = p->j_dir==1 ? (1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]) : 0.5*(p->DXN[IP]+p->DZN[KP]);

                if(topo<h)
                {
                    double sp = sqrt(PU[n]*PU[n] + PV[n]*PV[n] + PW[n]*PW[n]);
                    double wgt = sp < 0.1*settling_velocity(p,P.D[n]) ? 0.0 : topo/h;

                    uf *= wgt;
                    vf *= wgt;
                    wf *= wgt;
                    liftx=lifty=liftz=0.0;
                }
            }
        }
    }
    
    // two-way coupling, S 10 2 (porous bed): interstitial velocity u/eps
    if(p->Q50==1 && p->S10==2)
    {
        double eps = MAX(1.0-p->ccipol4a(Ts,PX[n],PY[n],PZ[n]), 1.0-theta_max);
        eps = MIN(eps,1.0);
        uf/=eps;
        vf/=eps;
        wf/=eps;
    }

    P.Uf[n] = uf;
    P.Vf[n] = vf;
    P.Wf[n] = wf;
    
    // relative velocity
    Urel = uf-PU[n];
    Vrel = vf-PV[n];
    Wrel = wf-PW[n];
    
    Tsval = p->ccipol4a(Ts,PX[n],PY[n],PZ[n]);
    
    // drag coefficient
    Dpx = drag_model(p,P.D[n],P.RO[n],sqrt(Urel*Urel + Vrel*Vrel + Wrel*Wrel),Tsval);
    Dpy = Dpz = Dpx;

    // explicit forces
    F = Dpx*uf - dPx_val/P.RO[n] + Bx + liftx;
    G = Dpy*vf - dPy_val/P.RO[n] + By + lifty;
    H = Dpz*wf - dPz_val/P.RO[n] + Bz + liftz;
    
    // inter-particle normal stress
    if(p->Q12>0)
    {
    dTx_val = p->ccipol4a(dTx,PX[n],PY[n],PZ[n]);
    dTy_val = p->j_dir==1 ? p->ccipol4a(dTy,PX[n],PY[n],PZ[n]) : 0.0;
    dTz_val = p->ccipol4a(dTz,PX[n],PY[n],PZ[n]);
    
    double Tdiv = P.RO[n]*MAX(Tsval, 0.5*theta_bed);
    
    F -= dTx_val/Tdiv;
    G -= dTy_val/Tdiv;
    H -= dTz_val/Tdiv;
    }
    
    if(p->j_dir==0)
    G=0.0;
    
    P.Test[n] = Tsval;

    // solid forcing
    double fx,fy,fz;
    fx = p->ccipol1c(a->fbh1,PX[n],PY[n],PZ[n])*(0.0-PU[n])/dt;
    fy = p->ccipol2c(a->fbh2,PX[n],PY[n],PZ[n])*(0.0-PV[n])/dt;
    fz = p->ccipol3c(a->fbh3,PX[n],PY[n],PZ[n])*(0.0-PW[n])/dt;

    F += fx;
    G += fy;
    H += fz;

    // relax
    if(p->Q73>0)
    {
    double r = rf(p,PX[n],PY[n]);
    F *= r;
    G *= r;
    H *= r;
    }

    // error call
    if(F!=F || G!=G || H!=H || PU[n]!=PU[n] || PV[n]!=PV[n] || PW[n]!=PW[n])
    {
        cout<<"CPM NaN detected.\nUrel: "<<Urel<<" Vrel: "<<Vrel<<" Wrel: "<<Wrel<<"\nDrag: "<<Dpx<<"\nTs: "<<Tsval<<endl;
        cout<<"F: "<<F<<" G: "<<G<<" H: "<<H<<endl;
        cout<<"dTx: "<<dTx_val<<" dTy: "<<dTy_val<<" dTz: "<<dTz_val<<endl;
        exit(1);
    }
}

/*--------------------------------------------------------------------
fluid velocity at parcels inside the bed for S 10 1, where the bed (topo<0) is a
solid boundary for the fluid and the bed surface lies on the grid scale:

  - surface layer, depth below the bed surface delta < d (the top grain layer): exposed grains,
    fluid velocity sampled 0.5 h above the bed surface along the bed normal and
    scaled to the grain with the logarithmic law of the wall,
        u_g = u_ref ln(30 d/ks)/ln(30 (0.5h)/ks),   ks = S 21 d50, d the parcel diameter
  - deeper: the pore water moves with the grains (local solid velocity)
  - linear transition between d and 2 d (the exposed volume per bed area is theta d,
    independent of the cell size)
--------------------------------------------------------------------*/
void CPM::nearbed_velocity(lexer *p, fdm *a, double xp, double yp, double zp, double dp)
{
    liftx=lifty=liftz=0.0;
    nb_w=nb_fac=0.0;
    
    double topo = ptopo(p,a,xp,yp,zp);
    
    if(topo>=0.0)
    return;
    
    i = p->posc_i(xp);
    j = p->posc_j(yp);
    k = p->posc_k(zp);
    
    double h = p->j_dir==1 ? (1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]) : 0.5*(p->DXN[IP]+p->DZN[KP]);
    double delta = -topo;
    
    // pore water with the grains
    double ug = p->ccipol4a(Us,xp,yp,zp);
    double vg = p->j_dir==1 ? p->ccipol4a(Vs,xp,yp,zp) : 0.0;
    double wg = p->ccipol4a(Ws,xp,yp,zp);
    
    nb_ug=ug; nb_vg=vg; nb_wg=wg;
    
    // exposed layer: the top grain layer, depth below the bed surface delta < d,
    // transition to the pore water velocity between d and 2 d
    double de = MIN(dp, 0.25*h);
    double w = 0.0;
    
    // exposure from the column ranking (exposure_update), else from the depth
    if(expo_w>=0.0)
    w = expo_w;
    
    else if(delta<de)
    w = 1.0;
    
    else if(delta<2.0*de)
    w = (2.0*de-delta)/de;
    
    if(w>0.0)
    {
        // bed normal from the topo level set
        double nx = (ptopo(p,a,xp+0.5*h,yp,zp) - ptopo(p,a,xp-0.5*h,yp,zp));
        double ny = p->j_dir==1 ? (ptopo(p,a,xp,yp+0.5*h,zp) - ptopo(p,a,xp,yp-0.5*h,zp)) : 0.0;
        double nz = (ptopo(p,a,xp,yp,zp+0.5*h) - ptopo(p,a,xp,yp,zp-0.5*h));
        double nm = sqrt(nx*nx + ny*ny + nz*nz);
        
        if(nm>1.0e-10)
        {
            nx/=nm;
            ny/=nm;
            nz/=nm;
        }
        else
        {
            nx=ny=0.0;
            nz=1.0;
        }
        
        double zref = 0.5*h;
        double xr = xp + nx*(delta+zref);
        double yr = yp + ny*(delta+zref);
        double zr = zp + nz*(delta+zref);
        
        double ur = p->ccipol1c(a->u,xr,yr,zr);
        double vr = p->j_dir==1 ? p->ccipol2c(a->v,xr,yr,zr) : 0.0;
        double wr = p->ccipol3c(a->w,xr,yr,zr);
        
        double ks = MAX(p->S21*p->S20, 1.0e-6);
        // velocity at the top of the grain dp in the log layer of the bed roughness ks = S 21 d50:
        // small grains of a mixture lie lower in the roughness and are hidden (WP6)
        double fac = log(MAX(30.0*dp/ks,1.0+1.0e-6))/log(MAX(30.0*zref/ks,1.0+1.0e-6));
        fac = MAX(0.0,MIN(fac,1.0));
        
        uf = w*shelter*fac*ur + (1.0-w)*ug;
        vf = w*shelter*fac*vr + (1.0-w)*vg;
        wf = w*shelter*fac*wr + (1.0-w)*wg;
        
        nb_w=w; nb_fac=fac; nb_xr=xr; nb_yr=yr; nb_zr=zr;
        
        // lift on the exposed grains (Q 54 C_L, Wiberg & Smith 1985):
        //   F_L = 0.5 rho_f C_L A (u_T^2 - u_B^2),  u_T at the grain top, u_B = 0 in the roughness,
        //   per unit mass: 0.75 C_L rho_f/rho_p u_T^2/d, along the bed normal
        if(p->Q54>0.0 && shelter>0.0)
        {
            double ut = fac*ur, vt = fac*vr, wt = fac*wr;
            double un = ut*nx + vt*ny + wt*nz;
            ut -= un*nx;
            vt -= un*ny;
            wt -= un*nz;
            
            double aL = w*0.75*p->Q54*p->W1/p->S22*(ut*ut + vt*vt + wt*wt)/MAX(dp,1.0e-9);
            
            liftx = aL*nx;
            lifty = aL*ny;
            liftz = aL*nz;
        }
    }
    else
    {
        uf = ug;
        vf = vg;
        wf = wg;
    }
}

/*--------------------------------------------------------------------
Bagnold sheltering (Q 57 1, S 10 1)

Bagnold (1956), Owen (1964): the grains in motion carry part of the bed shear stress
and pass it on to the bed in their contacts; the grains at rest only see the rest,

    tau_G = sum over the moving grains of a bed column  mu_s m g' f_m / A
    psi   = max(0, 1 - tau_G/tau_b)

tau_b: bed shear stress of the fluid from the log law at the reference point 0.5 h above
the bed (same velocity as the near-bed closure, ks = S 21 d50), f_m = min(1, |u_p|/(0.1 u*)) the
mobility of a grain. An exposed grain at rest sees the near-bed velocity reduced by psi,
a moving grain the full velocity (it travels in the flow above the bed):

    u_g = w fac u_ref (f_m + (1-f_m) psi) + (1-w) u_s

tau_G is averaged over 1 s and over the neighbouring bed columns (a parcel stands for many
grains). So entrainment stops when the moving grains carry tau_b - tau_c, while they move with the
velocity of the undisturbed near-bed flow: the transport rate follows n u_p with
n ~ (tau_b - tau_c)/(mu_s rho' g d) and u_p ~ u* - u*_c (the structure of MPM type formulas).
--------------------------------------------------------------------*/
void CPM::bagnold_update(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    if(p->Q57!=1 || p->S10==2)
    return;
    
    const int nj = p->knoy+2;
    tauG.assign((p->knox+2)*nj, 0.0);
    tauB.assign((p->knox+2)*nj, 0.0);
    
    // bed shear stress per column: log law at 0.5 h above the bed level
    const double ks = MAX(p->S21*p->S20, 1.0e-6);
    
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    {
        k = MAX(0,MIN(p->knoz-1, p->posc_k(s->bedzh(i,j))));
        double h = p->j_dir==1 ? (1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]) : 0.5*(p->DXN[IP]+p->DZN[KP]);
        double zr = s->bedzh(i,j) + 0.5*h;
        
        double ur = p->ccipol1c(a->u,p->XP[IP],p->YP[JP],zr);
        double vr = p->j_dir==1 ? p->ccipol2c(a->v,p->XP[IP],p->YP[JP],zr) : 0.0;
        
        double uplus = log(MAX(30.0*0.5*h/ks, 1.0+1.0e-6))/0.4;
        double us = sqrt(ur*ur + vr*vr)/uplus;
        
        tauB[(i+1)*nj + j+1] = p->W1*us*us;
    }
    
    const double vpar = P.ParcelFactor*Vp;
    const double gmag = sqrt(p->W20*p->W20 + p->W21*p->W21 + p->W22*p->W22);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        i = p->posc_i(P.X[n]);
        j = p->j_dir==1 ? p->posc_j(P.Y[n]) : 0;
        
        if(i<0 || i>=p->knox || j<0 || j>=p->knoy)
        continue;
        
        k = p->posc_k(P.Z[n]);
        double h = p->j_dir==1 ? (1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]) : 0.5*(p->DXN[IP]+p->DZN[KP]);
        
        // grains of the bed below the exposed layer do not count
        if(ptopo(p,a,P.X[n],P.Y[n],P.Z[n]) < -2.0*MIN(P.D[n],0.25*h))
        continue;
        
        double ustar = sqrt(tauB[(i+1)*nj + j+1]/p->W1);
        double sp = sqrt(P.U[n]*P.U[n] + P.V[n]*P.V[n] + P.W[n]*P.W[n]);
        double fm = MIN(1.0, sp/MAX(0.1*ustar,1.0e-10));
        
        tauG[(i+1)*nj + j+1] += fm*mu_s*vpar*(p->S22 - p->W1)*gmag/(p->DXN[IP]*p->DYN[JP]);
    }
    
    // a parcel stands for many grains and a bed column holds less than one moving parcel
    // at Bagnold's limit: tau_G is averaged in time (relaxation time 1 s, once per fluid step)
    // and in space (binomial filter over the neighbouring columns)
    if(bag_count!=p->count)
    {
        double r = bag_count<0 ? 1.0 : MIN(1.0, p->dt/1.0);
        bag_count = p->count;
        
        SLICELOOP4
        tauGf(i,j) += r*(tauG[(i+1)*nj + j+1] - tauGf(i,j));
        
        for(int qn=0;qn<4;++qn)
        {
            pgc->gcsl_start4(p,tauGf,1);
            
            if(perx==1)
            for(j=0;j<p->knoy;++j)
            {
                tauGf(-1,j) = tauGf(p->knox-1,j);
                tauGf(p->knox,j) = tauGf(0,j);
            }
            
            if(pery==1)
            for(i=0;i<p->knox;++i)
            {
                tauGf(i,-1) = tauGf(i,p->knoy-1);
                tauGf(i,p->knoy) = tauGf(i,0);
            }
            
            std::vector<double> tmp((p->knox+2)*nj,0.0);
            
            SLICELOOP4
            {
                double v = 0.5*tauGf(i,j) + 0.25*(tauGf(i-1,j) + tauGf(i+1,j));
                
                if(p->j_dir==1)
                v = 0.5*v + 0.25*(tauGf(i,j-1) + tauGf(i,j+1));
                
                tmp[(i+1)*nj + j+1] = v;
            }
            
            SLICELOOP4
            tauGf(i,j) = tmp[(i+1)*nj + j+1];
        }
        
        pgc->gcsl_start4(p,tauGf,1);
    }
}

double CPM::shelter_factor(lexer *p, sediment_fdm *s, double xp, double yp, double zp, double up, double vp, double wp)
{
    if(p->Q57!=1 || tauG.empty())
    return 1.0;
    
    int ii = p->posc_i(xp);
    int jj = p->j_dir==1 ? p->posc_j(yp) : 0;
    
    if(ii<0 || ii>=p->knox || jj<0 || jj>=p->knoy)
    return 1.0;
    
    double taub = MAX(tauB[(ii+1)*(p->knoy+2) + jj+1],1.0e-12);
    double psi = MAX(0.0, 1.0 - tauGf(ii,jj)/taub);
    
    double ustar = sqrt(taub/p->W1);
    double sp = sqrt(up*up + vp*vp + wp*wp);
    double fm = MIN(1.0, sp/MAX(0.1*ustar,1.0e-10));
    
    return fm + (1.0-fm)*psi;
}

/*--------------------------------------------------------------------
exposure of the bed parcels (S 10 1)

Only the top grain layer of the bed is exposed to the flow. A parcel position inside a
cell carries no meaning on the grain scale, so the exposed parcels are found by ranking:
in each bed column the topmost parcels inside the bed (topo < 0) up to the volume of one
grain layer

    V_e = theta_0 d A        (A: area of the column)

are exposed, the last one with the remaining fraction. The exposed volume per bed area is
theta_0 d for any cell size and parcel factor.
--------------------------------------------------------------------*/
void CPM::exposure_update(lexer *p, fdm *a)
{
    if(p->S10==2)
    {
        expo.clear();
        return;
    }
    
    expo.assign(P.index, 0.0);
    
    const int nj = p->knoy+2;
    const double vpar = P.ParcelFactor*Vp;
    
    std::vector<std::vector<std::pair<double,int>>> col((p->knox+2)*nj);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        i = p->posc_i(P.X[n]);
        j = p->j_dir==1 ? p->posc_j(P.Y[n]) : 0;
        
        if(i<0 || i>=p->knox || j<0 || j>=p->knoy)
        continue;
        
        if(ptopo(p,a,P.X[n],P.Y[n],P.Z[n]) >= 0.0 || P.Hop[n]>0.0)
        continue;
        
        col[(i+1)*nj + j+1].push_back(std::make_pair(P.Z[n],int(n)));
    }
    
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    {
        auto &c = col[(i+1)*nj + j+1];
        
        if(c.empty())
        continue;
        
        std::sort(c.begin(), c.end(), [](const std::pair<double,int> &a1, const std::pair<double,int> &a2){return a1.first>a2.first;});
        
        double Ve = theta_0*p->S20*p->DXN[IP]*p->DYN[JP];
        
        for(size_t q=0; q<c.size() && Ve>0.0; ++q)
        {
            int nn = c[q].second;
            double f = MIN(1.0, Ve/vpar);
            
            // exposure weight: share of the parcel in the top grain layer
            expo[nn] = f;
            Ve -= f*vpar;
        }
    }
}
