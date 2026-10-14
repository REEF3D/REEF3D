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

#include"bedshear.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"turbulence.h"
#include"sliceint.h"


bedshear::bedshear(lexer *p, turbulence *ppturb) : norm_vec(p), kappa(0.4), uw_mean(p), vw_mean(p), uw_old(p), vw_old(p), uw_hat(p)
{
    tau=0.0;
    tauc=0.0;
    pturb=ppturb;
    tw_old=0.0;
    tw_dt=0.0;
    tw_count=0;
}

// S 16 9 ------------------------------------------------------------------------------------------
// The log law at the sample height assumes that the boundary layer reaches the sample. Under waves the
// wave boundary layer is only a few mm thick, the sample (1.6 DZN, or the first cell) lies in the
// free stream, and the log law then gives a far too small shear stress. The wave friction stress
//     tau_w = 0.5 rho fw(A/ks) |u_w|^2,   fw = 1.39 (A/z0)^-0.52 <= 0.3, z0 = ks/30   (Soulsby 1997)
// uses the oscillatory part u_w = u - u_mean of the same velocity sample (u_mean: running mean over
// one wave period, the current); A = Uhat T/2pi with the orbital amplitude of the sample
//     Uhat^2 = |u_w|^2 + |du_w/dt|^2/omega^2   (exact for a harmonic motion)
// and T = S 58 (if > 0) or the period of the wave generation. S 16 9 takes max(tau_log, tau_w):
// in a current (no wave period, or u_w -> 0) it is S 16 1. Wave-current interaction is not included.
void bedshear::tau_wave_step(lexer *p)
{
    // time since the previous call: the derivative of the sample needs a new time level
    tw_dt = p->simtime - tw_old;
    
    if(tw_dt>1.0e-12)
    {
    tw_old = p->simtime;
    ++tw_count;
    }
}

double bedshear::tau_wave(lexer *p, double uvel, double vvel, double ks, double density)
{
    double T = p->S58>0.0 ? p->S58 : p->wT;
    
    if(T<=1.0e-10)
    return 0.0;
    
    double omega = 2.0*PI/T;
    
    // current: running mean over one wave period; oscillatory part and its amplitude
    // (kept when called twice at the same time level)
    if(tw_dt>1.0e-12)
    {
        if(tw_count>1)
        {
        double r = MIN(tw_dt/T, 1.0);
        uw_mean(i,j) += r*(uvel - uw_mean(i,j));
        vw_mean(i,j) += r*(vvel - vw_mean(i,j));
        }
        
        double uw = uvel - uw_mean(i,j);
        double vw = vvel - vw_mean(i,j);
        
        double uhat2 = uw*uw + vw*vw;
        
        if(tw_count>1)
        {
        double dudt = (uw - uw_old(i,j))/tw_dt;
        double dvdt = (vw - vw_old(i,j))/tw_dt;
        
        uhat2 += (dudt*dudt + dvdt*dvdt)/(omega*omega);
        }
        
    uw_hat(i,j) = sqrt(uhat2);
    uw_old(i,j) = uw;
    vw_old(i,j) = vw;
    }
    
    double A = uw_hat(i,j)/omega;
    double z0 = MAX(ks,1.0e-12)/30.0;
    
    double fw = A>1.0e-12 ? MIN(1.39*pow(A/z0,-0.52), 0.3) : 0.3;
    
    return 0.5*density*fw*(uw_old(i,j)*uw_old(i,j) + vw_old(i,j)*vw_old(i,j));
}

bedshear::~bedshear()
{
}

// CFD
void bedshear::taubed(lexer *p, fdm * a, ghostcell *pgc, sediment_fdm *s)
{
	int count;
	double zval,fac,topoval,taukin,tauvel,density;
    double tau_eff,tau_i;
    
    if(p->S16==9)
    tau_wave_step(p);
    
    SEDSLICELOOP
    {
    k = s->bedk(i,j);
    
    xip= p->XP[IP];
	yip= p->YP[JP];
    dist = p->DZN[KP];
    tau_eff = 0.0;
		
    density = s->ro(i,j);
    
	
    if(p->S16==1 || p->S16==9)
    {
    dist = p->DZN[KP];
    
    zval = s->bedzh(i,j) + 1.6*dist;
    

    uvel=p->ccipol1(a->u,xip,yip,zval);
    vvel=p->ccipol2(a->v,xip,yip,zval);
    wvel=p->ccipol3(a->w,xip,yip,zval);

    u_abs = sqrt(uvel*uvel + vvel*vvel  + wvel*wvel);

    // log law at the height of the velocity sample
    u_plus = uplus(1.6*dist,s->ks_eff(i,j));

    tau_eff = density*(u_abs*u_abs)/(u_plus*u_plus);
    
    // tau_wave updates the running mean: called once (MAX is a macro)
    if(p->S16==9)
    {
    double tau_w = tau_wave(p,uvel,vvel,s->ks_eff(i,j),density);
    tau_eff = MAX(tau_eff, tau_w);
    }
    }
    
    
    if(p->S16==2)
    {
    dist = 1.1*p->DZN[KP];
    
    zval = s->bedzh(i,j) + dist;
    
        uvel=p->ccipol1(a->u,xip,yip,zval);
        vvel=p->ccipol2(a->v,xip,yip,zval);
        wvel=p->ccipol3(a->w,xip,yip,zval);
        
    u_abs = sqrt(uvel*uvel + vvel*vvel  + wvel*wvel);

    u_plus = uplus(dist,s->ks_eff(i,j));
    
    
    zval = s->bedzh(i,j) + 0.5*p->DZN[KP];

    tau_eff = MAX(density*(u_abs*u_abs)/(u_plus*u_plus), density*pturb->ccipol_a_kinval(p,pgc,xip,yip,zval)*0.3);
    }
    
    if(p->S16==3)
    {
	double v_t,v_d;
		
	xip = p->XP[IP];
	yip = p->YP[JP];

    zval = s->bedzh(i,j) + 0.6*p->DZN[KP];
	

    uvel=p->ccipol1(a->u,xip,yip,zval);
    vvel=p->ccipol2(a->v,xip,yip,zval);
    wvel=p->ccipol3(a->w,xip,yip,zval);
        
    v_d=p->ccipol4(a->visc,xip,yip,zval);
    v_t=p->ccipol4(a->eddyv,xip,yip,zval);

    u_abs = sqrt(uvel*uvel + vvel*vvel);
    
    // velocity gradient between the bed and the sample at 0.6 DZN
    tau_eff = density*(v_d + v_t)*(u_abs/(0.6*p->DZN[KP]));
    }
    
	if(p->S16==4)
    {
    zval = s->bedzh(i,j) + 0.5*p->DZN[KP];
    
    tau_eff = density*pturb->ccipol_a_kinval(p,pgc,xip,yip,zval)*0.3;
    }
	
	if(p->S16==5)
    {
    zval = s->bedzh(i,j) + 0.5*p->DZN[KP];
    
        uvel=p->ccipol1(a->u,xip,yip,zval);
        vvel=p->ccipol2(a->v,xip,yip,zval);
        wvel=p->ccipol3(a->w,xip,yip,zval);

        u_abs = sqrt(uvel*uvel + vvel*vvel);
			   
    u_plus = uplus(0.5*p->DZN[KP],s->ks(i,j));


    tauvel=density*(u_abs*u_abs)/(u_plus*u_plus);
	
	taukin=density*pturb->ccipol_a_kinval(p,pgc,xip,yip,zval)*0.3;
	
	tau_eff = sqrt(fabs(tauvel))*sqrt(fabs(taukin));
    }
    
    if(p->S16==6)
    {
        double Cval,wh;
        double bedlevel,waterlevel;
        count=0;
        uvel=vvel=wh=0.0;
        KLOOP
        PCHECK
        {
            if(a->phi(i,j,k)>=0.0)
            {
            uvel+=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
            vvel+=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
            ++count;
            wh+=p->DZN[KP];
            }
        }

        // dry column: no shear
        tau_eff = 0.0;
        
        if(count>0)
        {
        uvel=uvel/double(count);
        vvel=vvel/double(count);

        Cval=18.0*log10(MAX((12.0*wh)/s->ks(i,j),10.0));

        u_abs = sqrt(uvel*uvel + vvel*vvel);
	
        tau_eff = density*pow(sqrt(9.81)*(u_abs/Cval),2.0);
        }
        
        k = s->bedk(i,j);
    }
    
    if(p->S16==7)
    {
    // horizontal velocity at half the first cell height (zval was not set in this branch)
    zval = s->bedzh(i,j) + 0.5*p->DZN[KP];

    uvel=p->ccipol1(a->u,xip,yip,zval);
    vvel=p->ccipol2(a->v,xip,yip,zval);


    u_abs = sqrt(uvel*uvel + vvel*vvel);

    u_plus = uplus(0.5*p->DZN[KP],s->ks(i,j));

    tau_eff = density*(u_abs*u_abs)/(u_plus*u_plus);
    }
    
    
    if(p->S16==8)
    {
    double E,B,dB,ustar,y_plus,visc,ks_plus,tau0;
    double v_d,v_t;
    
    visc=p->W2;
    
    zval = s->bedzh(i,j) + 1.6*p->DZN[KP];
    
        uvel=p->ccipol1(a->u,xip,yip,zval);
        vvel=p->ccipol2(a->v,xip,yip,zval);
        wvel=p->ccipol3(a->w,xip,yip,zval);
        
        v_d=p->ccipol4(a->visc,xip,yip,zval);
        v_t=p->ccipol4(a->eddyv,xip,yip,zval);

        
    // predictor    
    dist = 1.6*p->DZN[KP];
    u_abs = sqrt(uvel*uvel + vvel*vvel);
    u_plus = uplus(dist,s->ks(i,j));
    tau0=tau=density*(u_abs*u_abs)/(u_plus*u_plus);
    ustar=sqrt(tau/density);
    
    visc = v_d;
    
    
        // corrector
        for(int qn=0; qn<5;++qn)
        {
        y_plus = ustar*dist;
        
        ks_plus = visc + 0.246*ustar*s->ks(i,j);
        
        u_plus = (1.0/kappa)*log(MAX(y_plus/ks_plus,1.0)) + 5.0;

        
        tau_eff = MIN(density*(u_abs*u_abs)/pow((u_plus>1.0e-4?u_plus:1.0e20),2.0), tau0*3.5);

        ustar=sqrt(tau_eff/density);
        }
    }
    
    // S 16 10: the bed shear stress of the CPM bedload layer (bedshear::column_tau): log law at the first cell
    // centre at least half a cell above the bed level, or S 63 cells above it, and max(log law, 0.3 rho k) with
    // k of the first fluid cell (validation 16: in separated flow the log law gives a fraction of the skin friction)
    if(p->S16==10)
    {
        double tx=0.0, ty=0.0, hc=0.0;
        bool rest=false;
        
        bedshear::column_tau(p,a,pturb,i,j,s->bedzh(i,j),s->ks_eff(i,j),p->S63,true,density,tx,ty,hc,rest);
        
        tau_eff = sqrt(tx*tx + ty*ty);
    }
    
    // tau_i for Ti: grain roughness, log law at the height of the velocity sample
        dist = 1.6*p->DZN[KP];
        
        zval = s->bedzh(i,j) + dist;
        
        uvel=p->ccipol1(a->u,xip,yip,zval);
        vvel=p->ccipol2(a->v,xip,yip,zval);
        wvel=p->ccipol3(a->w,xip,yip,zval);
        
        u_abs = sqrt(uvel*uvel + vvel*vvel  + wvel*wvel);
        
        u_plus = uplus(dist,s->ks(i,j));

        tau_i = density*(u_abs*u_abs)/(u_plus*u_plus);
    
    s->tau_eff(i,j) = tau_eff;
    s->tau_i(i,j) = tau_i;
    s->shearvel_eff(i,j) = sqrt(tau_eff/density);
    s->shields_eff(i,j) = tau_eff/((p->S22-density)*fabs(p->W22)*p->S20);
    }
}

void bedshear::taucritbed(lexer *p, fdm * a, ghostcell *pgc, sediment_fdm *s)
{
	double density = p->W1;
    
    SEDSLICELOOP
    {
    // water density (a->ro(i,j,0) was the density of the lowest grid layer, not of the bed cell)
    density = p->W1;
    
    tauc = s->reduce(i,j) * (p->S30*fabs(p->W22)*(p->S22-density))*p->S20;
  
    s->tau_crit(i,j) = tauc;
    s->shearvel_crit(i,j) = sqrt(tauc/density);
    s->shields_crit(i,j) = p->S30*s->reduce(i,j);
    }
    
    SEDSLICELOOP
    s->MOB(i,j) = s->shields_eff(i,j)/(fabs(s->shields_crit(i,j))>1.0e-10?s->shields_crit(i,j):1.0e10);
}


/*--------------------------------------------------------------------
bed shear stress of a bed column from the flow above the bed level set (S 16 10 and the CPM bedload layer)

  reference cell: the first cell centre at least half a cell above the bed level (distance from the topo
  level set, with a rounding tolerance: a flat bed on a cell face puts the centre exactly half a cell above
  it); a solid body above the bed (a pipeline with a small gap) takes the fluid cell below it if there is one,
  else the column is at rest
  href > 0: the velocity of the log law is taken href cells above the bed, linear between the cell centres of
  the column, at most half way to a solid body or the free surface above (the first cell next to the bed is
  under-resolved in an accelerated flow, e.g. the gap below a pipeline)
  log law:  tau = rho (u / u+)^2,  u+ = ln(30 z_r/k_s)/kappa
  kterm:    tau = max(tau, 0.3 rho k) with k of the reference cell: in a separated or reattaching flow the
            near-bed turbulence moves the grains while the mean velocity next to the bed is small
            (backward-facing step: the log law gives 0-13 % of the measured skin friction in the
            recirculation, 0.3 rho k its size, validation 16)
  k is the value of the first fluid cell in both sediment models: the Eulerian S 16 2 interpolates k half a
  cell above the bed level, which takes in the cell below the bed where the bed cuts through the cells
--------------------------------------------------------------------*/
bool bedshear::column_tau(lexer *p, fdm *a, turbulence *pt, int i, int j, double zbed, double ks, double href, bool kterm,
                          double rho, double &tx, double &ty, double &hcell, bool &rest)
{
    int kb = MAX(0,MIN(p->knoz-1, p->posc_k(zbed)));
    double ur=0.0, vr=0.0, zr=0.0;
    bool found=false;
    int k;
    
    rest=false;
    ks = MAX(ks, 1.0e-6);

    for(k=kb;k<p->knoz;++k)
    {
        double h = p->DZN[k+marge];

        if(a->topo(i,j,k) >= 0.5*h*(1.0-1.0e-6))
        {
            bool body = p->solidread>0 && a->solid(i,j,k)<0.0;
            
            if(body && k>kb && a->topo(i,j,k-1)>0.0 && a->solid(i,j,k-1)>0.0)
            {
                --k;
                body = false;
            }
            
            rest = body;
            ur = body ? 0.0 : 0.5*(a->u(i-1,j,k) + a->u(i,j,k));
            vr = (body || p->j_dir==0) ? 0.0 : 0.5*(a->v(i,j-1,k) + a->v(i,j,k));
            zr = a->topo(i,j,k);
            hcell = h;
            found = true;
            break;
        }
    }

    // no fluid cell above the bed in this subdomain (vertical decomposition)
    if(!found)
    {
        tx = ty = 0.0;
        hcell = p->DZN[kb+marge];
        return false;
    }

    if(href>0.0 && !rest)
    {
        int k0 = k;
        double h = hcell;
        double zb = p->ZP[k0+marge] - a->topo(i,j,k0);
        double zref = zb + href*h;
        
        // top of the fluid column above the bed: solid body or free surface
        for(int kk=k0+1; kk<p->knoz; ++kk)
        {
            double dtop = 1.0e20;
            
            if(p->solidread>0 && a->solid(i,j,kk)<0.0)
            dtop = p->ZP[kk+marge] + a->solid(i,j,kk) - zb;
            
            if(a->phi(i,j,kk)<0.0)
            dtop = MIN(dtop, p->ZP[kk+marge] + a->phi(i,j,kk) - zb);
            
            if(dtop<1.0e19)
            {
                zref = MIN(zref, zb + 0.5*dtop);
                break;
            }
        }
        
        zref = MAX(zref, zb + zr);
        
        // linear between the cell centres around zref
        int kr = k0;
        while(kr+1<p->knoz && p->ZP[kr+1+marge]<=zref)
        ++kr;
        
        if(kr>k0 && kr+1<p->knoz && !(p->solidread>0 && a->solid(i,j,kr+1)<0.0))
        {
            double f = (zref - p->ZP[kr+marge])/(p->ZP[kr+1+marge] - p->ZP[kr+marge]);
            double u0 = 0.5*(a->u(i-1,j,kr) + a->u(i,j,kr)), u1 = 0.5*(a->u(i-1,j,kr+1) + a->u(i,j,kr+1));
            ur = (1.0-f)*u0 + f*u1;
            
            if(p->j_dir==1)
            {
                double v0 = 0.5*(a->v(i,j-1,kr) + a->v(i,j,kr)), v1 = 0.5*(a->v(i,j-1,kr+1) + a->v(i,j,kr+1));
                vr = (1.0-f)*v0 + f*v1;
            }
            
            zr = zref - zb;
        }
        else if(kr>k0)
        {
            ur = 0.5*(a->u(i-1,j,kr) + a->u(i,j,kr));
            vr = p->j_dir==1 ? 0.5*(a->v(i,j-1,kr) + a->v(i,j,kr)) : 0.0;
            zr = p->ZP[kr+marge] - zb;
        }
    }

    double um = sqrt(ur*ur + vr*vr);
    bool kt = kterm && pt!=nullptr && !rest;
    double tau = tau_logk(um, zr, ks, kt ? pt->kinval(i,j,k) : 0.0, kt, rho);

    tx = um>1.0e-12 ? tau*ur/um : 0.0;
    ty = um>1.0e-12 ? tau*vr/um : 0.0;
    
    return true;
}

// the closure of column_tau, shared by the CFD and the NHFLOW column (S 16 10) and the CPM bedload layer
double bedshear::tau_logk(double um, double zr, double ks, double kin, bool kterm, double rho)
{
    double up = log(MAX(30.0*zr/ks, 1.0+1.0e-6))/0.4;
    double us = um/up;
    double tau = rho*us*us;
    
    if(kterm)
    tau = MAX(tau, 0.3*rho*MAX(kin,0.0));
    
    return tau;
}
