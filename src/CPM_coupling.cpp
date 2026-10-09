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

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

/*--------------------------------------------------------------------
Two-way coupling (Q 50 1)

Fluid momentum of the two-fluid equations, eps = 1 - theta the fluid fraction,
u_i the interstitial fluid velocity:

    eps rho_f Du_i/Dt = -eps grad(p) + eps rho_f g - F,    F = sum_p w_p m_p D (u_i - u_p) / V

Divided by eps it is the single-phase momentum equation of the fluid solver with
the drag reaction F/eps as source term:

    rho_f Du/Dt = -grad(p) + rho_f g - F/eps

The share theta grad(p) of the pressure gradient that is carried by the parcels
stays in the implicit pressure of the projection, no explicit pressure gradient
is fed back (that is unstable). A fixed bed gives grad(p) = -F/eps, the pressure
drop of Ergun with the Gidaspow drag (Q 51 2). The parcels feel the full pressure
gradient -V_p grad(p) (buoyancy and dynamic part) in the parcel equation.

The fluid solver velocity u is taken as u_i = u/eps in the drag, the source is
point implicit in each RK stage of the momentum equation:

    u_new = (u + alpha dt KU) / (1 + alpha dt K)

    K  = sum w m D/(rho_f V eps^2)        KU = sum w m D u_p/(rho_f V eps)

K and KU are cell centred, averaged to the velocity faces.

S 10 1 (bed as solid boundary with wall function, near-bed closure for the exposed grains):
the fluid fraction is 1, the grains above the bed act on the fluid at their position; the
exposed grains see the log-law velocity u_g = w fac u_ref + (1-w) u_s scaled from the
reference point 0.5 h above the bed, their reaction acts there:

    K = m D w fac/(rho_f V),   KU = m D (u_p - (1-w) u_s)/(rho_f V)

A grain at rest is part of the bed roughness, whose drag is the bed shear stress of the wall
function; the moving grains (mobility fm = min(1, |u_p|/(0.1 |u_g|))) carry momentum away from
the near-bed flow. This is Bagnold's mechanism: the more grains move, the slower the near-bed
flow, until the fluid stress on the bed is near critical. No mixture continuity source.

Continuity: with u = eps u_i the superficial fluid velocity, the mixture is divergence free,

    div(u) = -div(theta u_p) = d theta/dt

so the fluid gives way to the parcels (return flow of a settling suspension). The source
Dsrc = -div(theta u_p) enters the right hand side of the pressure Poisson equation,
theta u_p from the parcel deposition, zero flux through walls.
--------------------------------------------------------------------*/

void CPM::coupling_update(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    if(p->Q50!=1 || p->Q11!=2)
    return;
    
    // current solid fraction
    volfrac_update(p,pgc,s,P.X,P.Y,P.Z,P.U,P.V,P.W);
    
    // mixture continuity for the porous bed (S 10 2); with S 10 1 the bed is solid for the fluid
    if(p->S10==2)
    continuity_source(p,pgc);
    
    int qi,qj,qk;
    double w,eps,us,vs,ws,slip,D,Tsp;
    
    const double vpar = P.ParcelFactor*Vp;
    
    for(i=-1;i<p->knox+1;++i)
    for(j=-1;j<p->knoy+1;++j)
    for(k=-1;k<p->knoz+1;++k)
    {
        Kc(i,j,k) = 0.0;
        KUx(i,j,k) = 0.0;
        KUy(i,j,k) = 0.0;
        KUz(i,j,k) = 0.0;
    }
    
    const bool porous = p->S10==2;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE && P.Hop[n]<=0.0)
    {
        // S 10 1: the bed is a solid boundary for the fluid (wall function at the bed),
        // see the comment at the top: only moving grains act on the fluid
        if(!porous)
        {
            uf = p->ccipol1c(a->u,P.X[n],P.Y[n],P.Z[n]);
            vf = p->j_dir==1 ? p->ccipol2c(a->v,P.X[n],P.Y[n],P.Z[n]) : 0.0;
            wf = p->ccipol3c(a->w,P.X[n],P.Y[n],P.Z[n]);
            
            nearbed_velocity(p,a,P.X[n],P.Y[n],P.Z[n],P.D[n]);
            
            double topo = ptopo(p,a,P.X[n],P.Y[n],P.Z[n]);
            
            // grains inside the bed below the exposed layer: no coupling
            if(topo<0.0 && nb_w<=0.0)
            continue;
            
            slip = sqrt((uf-P.U[n])*(uf-P.U[n]) + (vf-P.V[n])*(vf-P.V[n]) + (wf-P.W[n])*(wf-P.W[n]));
            D = drag_model(p,P.D[n],P.RO[n],slip,0.0);
            double mD = P.RO[n]*vpar*D;
            
            double kf=1.0, kux=P.U[n], kuy=P.V[n], kuz=P.W[n];
            double xk=P.X[n], yk=P.Y[n], zk=P.Z[n];
            
            if(topo<0.0)
            {
                // exposed grain: drag with the grain-scale velocity u_g = w fac u_ref + (1-w) u_s,
                // reaction at the reference point of u_ref; a grain at rest is part of the bed
                // roughness (its drag is in the wall function), the reaction of the moving grains
                // takes the momentum they carry away from the near-bed flow (Bagnold)
                double ug2 = uf*uf + vf*vf + wf*wf;
                double up2 = P.U[n]*P.U[n] + P.V[n]*P.V[n] + P.W[n]*P.W[n];
                double fm = MIN(1.0, sqrt(up2/MAX(0.01*ug2,1.0e-20)));
                
                if(fm<=1.0e-6)
                continue;
                
                mD *= fm;
                kf = nb_w*nb_fac;
                kux = P.U[n] - (1.0-nb_w)*nb_ug;
                kuy = P.V[n] - (1.0-nb_w)*nb_vg;
                kuz = P.W[n] - (1.0-nb_w)*nb_wg;
                xk=nb_xr; yk=nb_yr; zk=nb_zr;
            }
            
            kernel(p,xk,yk,zk);
            
            for(qi=0;qi<2;++qi)
            for(qj=0;qj<2;++qj)
            for(qk=0;qk<2;++qk)
            {
                w = kw[qi][qj][qk];
                
                if(w>0.0)
                {
                Kc(ki[qi],kj[qj],kk[qk])  += w*mD*kf;
                KUx(ki[qi],kj[qj],kk[qk]) += w*mD*kux;
                KUy(ki[qi],kj[qj],kk[qk]) += w*mD*kuy;
                KUz(ki[qi],kj[qj],kk[qk]) += w*mD*kuz;
                }
            }
            
            continue;
        }
        
        Tsp = cip4a(p,Ts,P.X[n],P.Y[n],P.Z[n]);
        eps = MAX(1.0-Tsp, 1.0-theta_max);
        eps = MIN(eps,1.0);
        
        // interstitial fluid velocity
        us = p->ccipol1c(a->u,P.X[n],P.Y[n],P.Z[n])/eps;
        vs = p->j_dir==1 ? p->ccipol2c(a->v,P.X[n],P.Y[n],P.Z[n])/eps : 0.0;
        ws = p->ccipol3c(a->w,P.X[n],P.Y[n],P.Z[n])/eps;
        
        slip = sqrt((us-P.U[n])*(us-P.U[n]) + (vs-P.V[n])*(vs-P.V[n]) + (ws-P.W[n])*(ws-P.W[n]));
        
        D = drag_model(p,P.D[n],P.RO[n],slip,Tsp);
        
        // drag coefficient of the fluid source F/eps, per unit fluid velocity of the solver
        double mD = P.RO[n]*vpar*D/eps;
        
        kernel(p,P.X[n],P.Y[n],P.Z[n]);
        
        for(qi=0;qi<2;++qi)
        for(qj=0;qj<2;++qj)
        for(qk=0;qk<2;++qk)
        {
            w = kw[qi][qj][qk];
            
            if(w>0.0)
            {
            Kc(ki[qi],kj[qj],kk[qk])  += w*mD/eps;
            KUx(ki[qi],kj[qj],kk[qk]) += w*mD*P.U[n];
            KUy(ki[qi],kj[qj],kk[qk]) += w*mD*P.V[n];
            KUz(ki[qi],kj[qj],kk[qk]) += w*mD*P.W[n];
            }
        }
    }
    
    pfold(p,Kc);
    
    pgc->start4a_sum(p,Kc,1);
    pfold(p,KUx);
    pgc->start4a_sum(p,KUx,1);
    pfold(p,KUy);
    pgc->start4a_sum(p,KUy,1);
    pfold(p,KUz);
    pgc->start4a_sum(p,KUz,1);
    
    BASELOOP
    {
        double rV = p->W1*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
        
        Kc(i,j,k) /= rV;
        KUx(i,j,k) /= rV;
        KUy(i,j,k) /= rV;
        KUz(i,j,k) /= rV;
    }
    
    pgc->start4a(p,Kc,1);
    pgc->start4a(p,KUx,1);
    pgc->start4a(p,KUy,1);
    pgc->start4a(p,KUz,1);
}

// momentum source of the parcels on the fluid, called in each RK stage of the momentum equation
void CPM::fluid_forcing(lexer *p, fdm *a, ghostcell *pgc, double alpha, field &u, field &v, field &w)
{
    if(p->Q50!=1)
    return;
    
    double K,KU,adt = alpha*p->dt;
    
    ULOOP
    {
        K  = 0.5*(Kc(i,j,k) + Kc(i+1,j,k));
        KU = 0.5*(KUx(i,j,k) + KUx(i+1,j,k));
        
        u(i,j,k) = (u(i,j,k) + adt*KU)/(1.0 + adt*K);
    }
    
    if(p->j_dir==1)
    VLOOP
    {
        K  = 0.5*(Kc(i,j,k) + Kc(i,j+1,k));
        KU = 0.5*(KUy(i,j,k) + KUy(i,j+1,k));
        
        v(i,j,k) = (v(i,j,k) + adt*KU)/(1.0 + adt*K);
    }
    
    WLOOP
    {
        K  = 0.5*(Kc(i,j,k) + Kc(i,j,k+1));
        KU = 0.5*(KUz(i,j,k) + KUz(i,j,k+1));
        
        w(i,j,k) = (w(i,j,k) + adt*KU)/(1.0 + adt*K);
    }
    
    pgc->start1(p,u,10);
    pgc->start2(p,v,11);
    pgc->start3(p,w,12);
}

void CPM::continuity_source(lexer *p, ghostcell *pgc)
{
    double qe,qw,qn,qs,qt,qb;
    
    // solid volume flux theta u_p, cell centred
    // KUx, KUy, KUz serve as work fields, they are filled afterwards
    BASELOOP
    {
        KUx(i,j,k) = Ts(i,j,k)*Us(i,j,k);
        KUy(i,j,k) = Ts(i,j,k)*Vs(i,j,k);
        KUz(i,j,k) = Ts(i,j,k)*Ws(i,j,k);
    }
    
    pgc->start4a(p,KUx,1);
    pgc->start4a(p,KUy,1);
    pgc->start4a(p,KUz,1);
    
    BASELOOP
    {
        qe = wallcell(p,i+1,j,k) ? 0.0 : 0.5*(KUx(i,j,k) + KUx(i+1,j,k));
        qw = wallcell(p,i-1,j,k) ? 0.0 : 0.5*(KUx(i,j,k) + KUx(i-1,j,k));
        qt = wallcell(p,i,j,k+1) ? 0.0 : 0.5*(KUz(i,j,k) + KUz(i,j,k+1));
        qb = wallcell(p,i,j,k-1) ? 0.0 : 0.5*(KUz(i,j,k) + KUz(i,j,k-1));
        
        qn = qs = 0.0;
        
        if(p->j_dir==1)
        {
        qn = wallcell(p,i,j+1,k) ? 0.0 : 0.5*(KUy(i,j,k) + KUy(i,j+1,k));
        qs = wallcell(p,i,j-1,k) ? 0.0 : 0.5*(KUy(i,j,k) + KUy(i,j-1,k));
        }
        
        Dsrc(i,j,k) = -(qe-qw)/p->DXN[IP] - (qn-qs)/p->DYN[JP] - (qt-qb)/p->DZN[KP];
    }
    
    pgc->start4a(p,Dsrc,1);
}

// no coupling before the parcels move
void CPM::coupling_reset(lexer *p, ghostcell *pgc)
{
    if(p->Q50!=1 || p->Q11!=2)
    return;
    
    for(i=-1;i<p->knox+1;++i)
    for(j=-1;j<p->knoy+1;++j)
    for(k=-1;k<p->knoz+1;++k)
    {
        Kc(i,j,k) = 0.0;
        KUx(i,j,k) = 0.0;
        KUy(i,j,k) = 0.0;
        KUz(i,j,k) = 0.0;
        Dsrc(i,j,k) = 0.0;
    }
}
