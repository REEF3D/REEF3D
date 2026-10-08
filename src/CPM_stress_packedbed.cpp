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
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

/*--------------------------------------------------------------------
Packed bed particle stress (Q 12 2)

The Snider (2001) stress is a function of the local solid fraction only. In a packed
bed it can carry the weight of the sediment only by compaction towards close packing,
which is slow, stiff and noisy, and it gives no shear resistance, i.e. no angle of
repose and no threshold of motion. Here the packed bed is described with the
soil mechanics picture of an enduring contact network:

1) effective stress of the contact network (Terzaghi): the submerged weight of the
   bed (theta >= theta_bed, Q 26) is transmitted down through the contact network,
   integrated along gravity from the top of the bed:

       Pov(z) = int_z^top (rho_s - rho_f)|g| theta dz'

   Pov is the isotropic intergranular pressure (lateral stress coefficient K = 1).
   Its vertical gradient carries the grains of a bed at rest without compaction,
   its horizontal gradient pushes the grains of a slope downslope.

2) contact pressure against compaction beyond the packing of the bed theta_0 = 1-n,
   following Johnson & Jackson (1987):

       Pc = Fr (theta - theta_0)^eta0 / (theta_max - theta)^eta1    theta > theta_0

   continued linearly above theta_0 + 0.5 (theta_max - theta_0) to bound the stiffness

3) Coulomb friction of each parcel on the bed below it, implicit, with the normal
   load given by the support of the contact network and the friction coefficient
   mu(I) of the mu(I) rheology (Jop et al. 2006):

       mu(I) = mu_s + (mu_2 - mu_s)/(I0/I + 1),  I = d50 gamma / sqrt(P/rho_s)

   with static friction (stick) below mu_s. A slope fails at tan(beta) = mu_s and an
   exposed grain moves when the drag exceeds mu_s times its submerged weight.

4) dilatancy (Q 55): theta_0 is lowered in sheared layers, theta_0(I), see dilatancy()

5) jammed bed: the parcels below a mobile surface layer of the bed are held at rest (the forcing
   of the fluid bed, X 41 h, without the bedload layer; with it a constraint on the parcel velocity
   below the layer of thickness Q 64 h under the iso-surface of the parcels, see advec_mppic).
   The friction 3) holds the parcels of the surface layer: a 30 degree slope stays, steeper slopes
   fail by avalanches of the surface layer to 29-32 degrees (submerged wedge, mu_s = 0.63).
   A collapsing sand tower spreads further (17-20 degrees, run-out of the dynamic collapse).
   Without 5) the isotropic stress (K = 1) and the grain-scale friction do not hold the deeper bed;
   deep failures need a frictional stress on the grid (open, WP4).

Ps = Pov + Pc acts on the parcels with -grad(Ps)/(theta rho_s).
--------------------------------------------------------------------*/

// Johnson & Jackson form above theta_0, continued linearly above theta_cap = theta_0 + 0.5 (theta_max - theta_0)
// to bound the stiffness, i.e. the wave speed of the particle stress and the sub-step size
// t0: packing of the contact network, theta_0 at rest, lower in a sheared layer (dilatancy, Q 55)
double CPM::contact_pressure(double theta, double t0)
{
    if(theta<=t0)
    return 0.0;
    
    double tmax = t0 + (theta_max - theta_0);
    double theta_cap = t0 + 0.5*(tmax - t0);
    double th = MIN(theta,theta_cap);
    
    double Pc = Fr*pow(th-t0,eta0)/pow(tmax-th,eta1);
    
    if(theta>theta_cap)
    Pc += contact_pressure_deriv(theta_cap,t0)*(theta-theta_cap);
    
    return Pc;
}

double CPM::contact_pressure_deriv(double theta, double t0)
{
    if(theta<=t0)
    return 0.0;
    
    double tmax = t0 + (theta_max - theta_0);
    double theta_cap = t0 + 0.5*(tmax - t0);
    double th = MIN(theta,theta_cap);
    
    double Pc = Fr*pow(th-t0,eta0)/pow(tmax-th,eta1);
    
    return Pc*(eta0/(th-t0) + eta1/(tmax-th));
}

/*--------------------------------------------------------------------
dilatancy (Q 55): a sheared granular layer is looser than the bed at rest, the packing
follows the inertial number like the friction coefficient (phi(I) of the mu(I) rheology,
Forterre & Pouliquen 2008):

    theta_0(I) = max(theta_bed, theta_0 - Q55 I),   I = d50 gamma / sqrt(P/rho_s)

gamma: shear rate of the grid solid velocity, P: overburden effective stress, at least the
weight of one grain layer. The contact pressure and the capacity of the grid-limited step
use theta_0(I): a sheared surface layer dilates instead of jamming at the packing of the bed,
a bed at rest (gamma = 0) keeps theta_0.
--------------------------------------------------------------------*/
void CPM::dilatancy(lexer *p, ghostcell *pgc)
{
    if(p->Q55<=0.0)
    {
        BASELOOP
        T0e(i,j,k) = theta_0;
        
        pgc->start4a(p,T0e,1);
        return;
    }
    
    double rog = (p->S22 - p->W1)*sqrt(p->W20*p->W20 + p->W21*p->W21 + p->W22*p->W22);
    
    auto der = [&](field &f, int di, int dj, int dk, double h)
    {
        bool wm = wallcell(p,i-di,j-dj,k-dk);
        bool wp = wallcell(p,i+di,j+dj,k+dk);
        
        if(wm && wp)
        return 0.0;
        
        if(wm)
        return (f(i+di,j+dj,k+dk)-f(i,j,k))/h;
        
        if(wp)
        return (f(i,j,k)-f(i-di,j-dj,k-dk))/h;
        
        return (f(i+di,j+dj,k+dk)-f(i-di,j-dj,k-dk))/(2.0*h);
    };
    
    BASELOOP
    {
        double hx = p->DXN[IP], hy = p->DYN[JP], hz = p->DZN[KP];
        
        double dudx = der(Us,1,0,0,hx), dudz = der(Us,0,0,1,hz);
        double dwdx = der(Ws,1,0,0,hx), dwdz = der(Ws,0,0,1,hz);
        double dudy=0.0, dvdx=0.0, dvdy=0.0, dvdz=0.0, dwdy=0.0;
        
        if(p->j_dir==1)
        {
        dudy = der(Us,0,1,0,hy);
        dvdx = der(Vs,1,0,0,hx);
        dvdy = der(Vs,0,1,0,hy);
        dvdz = der(Vs,0,0,1,hz);
        dwdy = der(Ws,0,1,0,hy);
        }
        
        double gamma = sqrt(2.0*(dudx*dudx + dvdy*dvdy + dwdz*dwdz)
                     + (dudz+dwdx)*(dudz+dwdx) + (dudy+dvdx)*(dudy+dvdx) + (dvdz+dwdy)*(dvdz+dwdy));
        
        double P = MAX(Pov(i,j,k), rog*MAX(Ts(i,j,k),theta_bed)*p->S20);
        double I = p->S20*gamma/sqrt(P/p->S22);
        
        T0e(i,j,k) = MAX(theta_bed, theta_0 - p->Q55*I);
    }
    
    pgc->start4a(p,T0e,1);
}

void CPM::stress_packedbed(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    stress_overburden(p,pgc,s);
    
    if(p->Q67>0 && p->Q58>0 && p->S10==1)
    stress_yield(p,pgc);
    
    dilatancy(p,pgc);
    
    cmax=0.0;
    
    BASELOOP
    {
        Tau(i,j,k) = Pov(i,j,k) + contact_pressure(Ts(i,j,k),T0e(i,j,k));
        
        cmax = MAX(cmax,contact_pressure_deriv(Ts(i,j,k),T0e(i,j,k)));
    }
    
    cmax = sqrt(pgc->globalmax(cmax)/p->S22);
    
    pgc->start4a(p,Tau,1);
    
    // solid fraction gradient: bed normal for the friction
    if(p->Q13==1)
    gradient(p,pgc,Ts,dSx,dSy,dSz);
}

// effective stress of the contact network, column integration from the top
// load of the bed: rho' g theta phi(theta), phi ramps from 0 to 1 between 0.5 theta_bed and theta_bed,
// a gap (theta < 0.5 theta_bed) separates suspended sediment from the bed
// with a vertical domain decomposition the integration is repeated until the values
// from the subdomains above have arrived
void CPM::stress_overburden(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double rog = (p->S22 - p->W1)*fabs(p->W22);
    double tgap = 0.5*theta_bed;
    double change,old,Lk,Lkp;
    bool above;
    
    int itermax = zsplit==1 ? 1000 : 1;
    
    for(int qn=0; qn<itermax; ++qn)
    {
        change = 0.0;
        
        ILOOP
        JLOOP
        for(k=p->knoz-1; k>=0; --k)
        {
            PBASECHECK
            {
            old = Pov(i,j,k);
            
            if(Ts(i,j,k)<tgap)
            Pov(i,j,k) = 0.0;
            
            else
            {
                Lk = rog*Ts(i,j,k)*MIN(1.0,(Ts(i,j,k)-tgap)/tgap);
                
                above = (k+1<p->knoz || p->nb6>=0) && Ts(i,j,k+1)>=tgap;
                
                if(above)
                {
                Lkp = rog*Ts(i,j,k+1)*MIN(1.0,(Ts(i,j,k+1)-tgap)/tgap);
                Pov(i,j,k) = Pov(i,j,k+1) + 0.5*Lkp*p->DZN[KP1] + 0.5*Lk*p->DZN[KP];
                }
                
                else
                Pov(i,j,k) = 0.5*Lk*p->DZN[KP];
            }
            
            change = MAX(change, fabs(Pov(i,j,k)-old));
            }
        }
        
        pgc->start4a(p,Pov,1);
        
        if(zsplit==1)
        {
            change = pgc->globalmax(change);
            
            if(change<1.0e-8*rog*hmin)
            break;
        }
    }
}

/*--------------------------------------------------------------------
yield of the jammed bed (Q 67), Mohr-Coulomb in the column

The jammed bed (Q 64) holds every parcel below the mobile surface layer at rest, so a slope
steeper than the angle of repose fails only by avalanches of the surface layer, and an
undercut flank never fails as a whole. With Q 67 the bed below the surface layer is jammed
only where its contact network can carry the load in friction. The lateral load of the
column above a point is the integral of the horizontal gradient of the overburden stress,

    F(z) = | int_z^top grad_h(Pov) dz' |,

the frictional resistance on the plane below is mu_s Pov(z). For an infinite slope of angle
beta, F/Pov = tan(beta) at any depth, so the bed yields at tan(beta) > mu_s down to the
bottom, as in the classical infinite-slope analysis; under a flat bed F = 0. The yield ratio

    Yr = F/(mu_s Pov)

releases the jam smoothly between Yr = 0.95 and 1.05 (yield_weight); the parcels of the yielding
bed are then held only by the Coulomb/mu(I) friction on their substrate, and the bed jams
again once its slope is below the angle of repose. A free face (vertical wall of sand, the
flank of a scour hole) has a large lateral gradient and yields from its top down.

The gradient is central, one-sided next to cells without load (outside the bed, solids,
walls); the column integration follows stress_overburden (repeated with a vertical domain
decomposition until the values from above have arrived).
--------------------------------------------------------------------*/

void CPM::stress_yield(lexer *p, ghostcell *pgc)
{
    double rog = (p->S22 - p->W1)*fabs(p->W22);
    double tgap = 0.5*theta_bed;
    double change,old;
    
    auto dPdh = [&](int di, int dj, double &gx)
    {
        // horizontal derivative of Pov in direction (di,dj) at (i,j,k), one-sided at cells without load
        bool okm = p->flag4[(i-di-p->imin)*p->jmax*p->kmax + (j-dj-p->jmin)*p->kmax + k-p->kmin]>0;
        bool okp = p->flag4[(i+di-p->imin)*p->jmax*p->kmax + (j+dj-p->jmin)*p->kmax + k-p->kmin]>0;
        double h = di==1 ? p->DXN[IP] : p->DYN[JP];
        double pm = Pov(i-di,j-dj,k), pp = Pov(i+di,j+dj,k), pc = Pov(i,j,k);
        
        if(okm && okp)
        gx = (pp - pm)/(2.0*h);
        
        else if(okp)
        gx = (pp - pc)/h;
        
        else if(okm)
        gx = (pc - pm)/h;
        
        else
        gx = 0.0;
    };
    
    int itermax = zsplit==1 ? 1000 : 1;
    
    // column integral of the lateral load: Fxy holds int grad_x, Yr temporarily int grad_y
    for(int qn=0; qn<itermax; ++qn)
    {
        change = 0.0;
        
        ILOOP
        JLOOP
        {
            double fxa=0.0, fya=0.0;   // integrand at the cell above
            
            for(k=p->knoz-1; k>=0; --k)
            {
                PBASECHECK
                {
                old = Fxy(i,j,k);
                
                if(Ts(i,j,k)<tgap || Pov(i,j,k)<=0.0)
                {
                    Fxy(i,j,k) = Yr(i,j,k) = 0.0;
                    fxa = fya = 0.0;
                }
                
                else
                {
                    double gx=0.0, gy=0.0;
                    dPdh(1,0,gx);
                    
                    if(p->j_dir==1)
                    dPdh(0,1,gy);
                    
                    bool above = (k+1<p->knoz || p->nb6>=0) && Ts(i,j,k+1)>=tgap && Pov(i,j,k+1)>0.0;
                    
                    if(above)
                    {
                        // integrand of the cell above from this loop; at the top of a subdomain
                        // (vertical split) the one of this cell
                        double fx1 = k+1<p->knoz ? fxa : gx;
                        double fy1 = k+1<p->knoz ? fya : gy;
                        
                        Fxy(i,j,k) = Fxy(i,j,k+1) + 0.5*fx1*p->DZN[KP1] + 0.5*gx*p->DZN[KP];
                        Yr(i,j,k)  = Yr(i,j,k+1)  + 0.5*fy1*p->DZN[KP1] + 0.5*gy*p->DZN[KP];
                    }
                    
                    else
                    {
                        Fxy(i,j,k) = 0.5*gx*p->DZN[KP];
                        Yr(i,j,k)  = 0.5*gy*p->DZN[KP];
                    }
                    
                    fxa = gx;
                    fya = gy;
                }
                
                change = MAX(change, fabs(Fxy(i,j,k)-old));
                }
            }
        }
        
        pgc->start4a(p,Fxy,1);
        pgc->start4a(p,Yr,1);
        
        if(zsplit==1)
        {
            change = pgc->globalmax(change);
            
            if(change<1.0e-8*rog*hmin)
            break;
        }
    }
    
    // yield ratio
    BASELOOP
    {
        double F = sqrt(Fxy(i,j,k)*Fxy(i,j,k) + Yr(i,j,k)*Yr(i,j,k));
        double Pv = Pov(i,j,k);
        
        Fxy(i,j,k) = F;
        Yr(i,j,k) = Pv>1.0e-6*rog*hmin ? F/(mu_s*Pv) : 0.0;
    }
    
    pgc->start4a(p,Fxy,1);
    pgc->start4a(p,Yr,1);
}

// jam weight of the yield (Q 67): 1 within the yield (Yr <= 0.95), 0 at Yr >= 1.05, smooth in between
double CPM::yield_weight(lexer *p, double xp, double yp, double zp)
{
    double yr = p->ccipol4a(Yr,xp,yp,zp);
    
    if(yr<=0.95)
    return 1.0;
    
    if(yr>=1.05)
    return 0.0;
    
    double xi = (yr-0.95)/0.1;
    
    return 0.5*(1.0 + cos(PI*xi));
}

/*--------------------------------------------------------------------
Coulomb friction of a parcel on its substrate, implicit

  - contact normal e: at the bed surface the bed normal grad(theta)/|grad(theta)|
    (pointing into the bed), inside the bed (small grad(theta)) gravity
  - substrate: the bed one cell along e, grid solid velocity,
    zero velocity at the bottom, lateral walls and solid bodies
  - contact if the substrate belongs to the bed: theta_substrate >= theta_bed
  - normal load per unit mass: support by the contact network along e,
        a_N = max(grad(Ps).e, 0)/(theta rho_s)
    in a bed at rest a_N is the submerged weight of the grains
  - the slip tangential to e is stopped if |s_t| <= dt mu_s a_N (static friction),
    otherwise reduced by dt mu(I) a_N (kinetic friction)

dt: effective time step of the velocity update, fac: 1/(1+dt Dp) of the implicit drag
--------------------------------------------------------------------*/

void CPM::friction(lexer *p, fdm *a, double xp, double yp, double zp, double &up, double &vp, double &wp, double dt, double fac)
{
    double gmag = sqrt(p->W20*p->W20 + p->W21*p->W21 + p->W22*p->W22);
    
    if(gmag<1.0e-10)
    return;
    
    i = p->posc_i(xp);
    j = p->posc_j(yp);
    k = p->posc_k(zp);
    
    double h = p->j_dir==1 ? (1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]) : 0.5*(p->DXN[IP]+p->DZN[KP]);
    
    // contact normal: bed normal at the surface, gravity inside the bed
    double gx = p->ccipol4a(dSx,xp,yp,zp);
    double gy = p->j_dir==1 ? p->ccipol4a(dSy,xp,yp,zp) : 0.0;
    double gz = p->ccipol4a(dSz,xp,yp,zp);
    double gr = sqrt(gx*gx + gy*gy + gz*gz);
    
    double w = MIN(1.0, gr*h/(0.5*theta_0));
    
    double ex = (1.0-w)*p->W20/gmag + (gr>1.0e-10 ? w*gx/gr : 0.0);
    double ey = (1.0-w)*p->W21/gmag + (gr>1.0e-10 ? w*gy/gr : 0.0);
    double ez = (1.0-w)*p->W22/gmag + (gr>1.0e-10 ? w*gz/gr : 0.0);
    double emag = sqrt(ex*ex + ey*ey + ez*ez);
    
    if(emag<1.0e-10)
    {
        ex = p->W20/gmag;
        ey = p->W21/gmag;
        ez = p->W22/gmag;
    }
    else
    {
        ex/=emag;
        ey/=emag;
        ez/=emag;
    }
    
    // substrate
    double xs = xp + ex*h;
    double ys = yp + ey*h;
    double zs = zp + ez*h;
    
    double Usub,Vsub,Wsub,Tsub;
    bool wall=false;
    
    if(zs<p->global_zmin || zs>p->global_zmax)
    wall=true;
    
    if(p->j_dir==1 && (ys<p->global_ymin || ys>p->global_ymax))
    wall=true;
    
    if(p->ccipol4_b(a->solid,xs,ys,zs)<0.0)
    wall=true;
    
    if(wall)
    {
        Usub=Vsub=Wsub=0.0;
    }
    else
    {
        Tsub = p->ccipol4a(Ts,xs,ys,zs);
        
        if(Tsub<theta_bed)
        return;
        
        Usub = p->ccipol4a(Us,xs,ys,zs);
        Vsub = p->j_dir==1 ? p->ccipol4a(Vs,xs,ys,zs) : 0.0;
        Wsub = p->ccipol4a(Ws,xs,ys,zs);
    }
    
    // a grain resting on a packed substrate or a wall: see the static friction below
    bool packed = wall || Tsub >= p->ccipol4a(T0e,xs,ys,zs) - 0.05;
    
    // normal load per unit mass: support by the contact network
    double Tsp = MAX(p->ccipol4a(Ts,xp,yp,zp),theta_bed);
    double dTe = p->ccipol4a(dTx,xp,yp,zp)*ex + (p->j_dir==1?p->ccipol4a(dTy,xp,yp,zp)*ey:0.0) + p->ccipol4a(dTz,xp,yp,zp)*ez;
    
    // (the same solid fraction as for the stress force on the parcel in advec_mppic, MAX(Ts, 0.5 theta_bed):
    // with MAX(Ts, theta_bed) the grains in the dilute surface cells of a slope had up to half the normal load)
    double aN = MAX(dTe,0.0)/(MAX(p->ccipol4a(Ts,xp,yp,zp),0.5*theta_bed)*p->S22);
    
    if(aN<1.0e-12)
    return;
    
    // tangential slip
    double sx = up - Usub;
    double sy = vp - Vsub;
    double sz = wp - Wsub;
    
    double sn = sx*ex + sy*ey + sz*ez;
    
    sx -= sn*ex;
    sy -= sn*ey;
    sz -= sn*ez;
    
    double smag = sqrt(sx*sx + sy*sy + sz*sz);
    
    // static friction: stick
    if(smag <= fac*dt*mu_s*aN)
    {
        up -= sx;
        vp -= sy;
        wp -= sz;
        
        // a sticking grain on a packed substrate or a wall is at rest on it: also no velocity
        // into it (inelastic normal contact). The contact network supports the grain on the
        // grain scale; without this the grains sink inside a packed cell (the support from the
        // cell-centred stress gradient is exact only on average over the cell), gather at the
        // faces of the full cells below and jitter on the bottom. Sliding grains keep their
        // normal velocity, so slopes can avalanche.
        if(packed && sn>0.0)
        {
            up -= sn*ex;
            vp -= sn*ey;
            wp -= sn*ez;
        }
        
        // inside the packed bed (own cell packed as well) the sticking grain is jammed:
        // it moves with its substrate, also no drift away from it (the stress gradient
        // balances gravity only on average over a cell, near the bottom with the ghost cells)
        if(packed && p->ccipol4a(Ts,xp,yp,zp) >= p->ccipol4a(T0e,xp,yp,zp) - 0.05)
        {
            up = Usub;
            vp = Vsub;
            wp = Wsub;
        }
        return;
    }
    
    // kinetic friction, mu(I) with I = d50 gamma / sqrt(P/rho_p)
    double Peff = MAX(p->ccipol4a(Tau,xp,yp,zp), (p->S22-p->W1)*gmag*Tsp*p->S20);
    double gamma = smag/h;
    double I = p->S20*gamma/sqrt(Peff/p->S22);
    double mu = mu_s + (mu_2-mu_s)/(I0/MAX(I,1.0e-10) + 1.0);
    
    double dU = MIN(fac*dt*mu*aN, smag);
    
    up -= dU*sx/smag;
    vp -= dU*sy/smag;
    wp -= dU*sz/smag;
}
