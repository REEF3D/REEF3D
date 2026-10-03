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
    
    // normal load per unit mass: support by the contact network
    double Tsp = MAX(p->ccipol4a(Ts,xp,yp,zp),theta_bed);
    double dTe = p->ccipol4a(dTx,xp,yp,zp)*ex + (p->j_dir==1?p->ccipol4a(dTy,xp,yp,zp)*ey:0.0) + p->ccipol4a(dTz,xp,yp,zp)*ez;
    
    double aN = MAX(dTe,0.0)/(Tsp*p->S22);
    
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
    
    if(smag<1.0e-12)
    return;
    
    // static friction: stick
    if(smag <= fac*dt*mu_s*aN)
    {
        up -= sx;
        vp -= sy;
        wp -= sz;
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
