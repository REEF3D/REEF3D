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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_waves.h"
#include"sediment_fdm.h"
#include"seastate_sflow.h"
#include"seastate_dispersion.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<cmath>
#include<vector>

sediment_waves::sediment_waves(lexer *p) : Hs(p),Tm10(p),Tp(p),Urms(p),Tbot(p),ex(p),ey(p),Dr(p),
                                           Tw(p),Dbr(p),kb(p),Ure(p),Uw(p),Sk(p),As(p),ua(p),
                                           tauc(p),tauw(p),taum(p),taumax(p),fw(p),stamp(-1)
{
    SLICEBASELOOP
    {
    Hs(i,j)=Tm10(i,j)=Tp(i,j)=Urms(i,j)=Tbot(i,j)=ex(i,j)=ey(i,j)=0.0;
    Dr(i,j)=-1.0;
    Tw(i,j)=Dbr(i,j)=kb(i,j)=Ure(i,j)=Uw(i,j)=Sk(i,j)=As(i,j)=ua(i,j)=0.0;
    tauc(i,j)=tauw(i,j)=taum(i,j)=taumax(i,j)=fw(i,j)=0.0;
    }
}

void sediment_waves::start(lexer *p, fdm2D *b, ghostcell *pgc, sediment_fdm *s)
{
    if(p->S110==1)
    seastate(p,b);
    
    if(p->S110==2)
    prescribed(p,s);
    
    derived(p,pgc,s);
}

void sediment_waves::seastate(lexer *p, fdm2D *b)
{
    // the spectra change only with a wave step
    if(b->wave==nullptr || b->wave->wave_steps()==stamp)
    return;
    
    b->wave->bed_waves(p,Hs,Tm10,Tp,Urms,Tbot,ex,ey,Dr);
    stamp = b->wave->wave_steps();
}

void sediment_waves::prescribed(lexer *p, sediment_fdm *s)
{
    const double pi = 3.14159265358979323846;
    const double sig = (p->S112>0.0) ? 2.0*pi/p->S112 : 0.0;
    const double cx = cos(p->S113*pi/180.0), cy = sin(p->S113*pi/180.0);
    
    SLICELOOP4
    {
    const double h = s->waterlevel(i,j);
    
    Hs(i,j)=Tm10(i,j)=Tp(i,j)=Urms(i,j)=Tbot(i,j)=0.0;
    ex(i,j) = cx;
    ey(i,j) = cy;
    Dr(i,j) = -1.0;
    
        if(h<p->S128 || !(sig>0.0) || !(p->S111>0.0))
        continue;
    
    const double H  = MIN(p->S111, p->S114*h);
    const double kh = seastate_wavenumber(sig,h)*h;
    
    Hs(i,j)   = H;
    Tp(i,j)   = p->S112;
    Tm10(i,j) = 0.9*p->S112;
    Tbot(i,j) = p->S112;
    Urms(i,j) = (kh<25.0) ? sig*H/(4.0*sinh(kh)) : 0.0;
    }
}

void sediment_waves::derived(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    const double pi = 3.14159265358979323846;
    const double g = fabs(p->W22);
    double B,psi;
    
    SLICELOOP4
    {
    const double h = s->waterlevel(i,j);
    const double Hrms = Hs(i,j)/sqrt(2.0);
    
    Tw(i,j)=Dbr(i,j)=kb(i,j)=Ure(i,j)=Uw(i,j)=Sk(i,j)=As(i,j)=ua(i,j)=0.0;
    
        if(h<p->S128 || !(Hs(i,j)>0.0))
        continue;
    
    Tw(i,j) = (p->S110==1 && Tbot(i,j)>0.0) ? Tbot(i,j) : Tp(i,j);
    
        // breaking dissipation per area and density: roller (SEASTATE surfbeat) or Baldock et al. (1998)
        if(Dr(i,j)>=0.0)
        Dbr(i,j) = Dr(i,j);
        
        if(Dr(i,j)<0.0 && Tp(i,j)>0.0 && Hrms>0.0)
        {
        const double Hb = 0.73*h;
        Dbr(i,j) = 0.25*g/Tp(i,j)*exp(-(Hb/Hrms)*(Hb/Hrms))*(Hb*Hb + Hrms*Hrms);
        }
    
        // stirring by breaking (S 118 1): turbulence at the bed
        if(p->S118==1 && Hrms>0.0 && Dbr(i,j)>0.0)
        {
        const double x = MIN(h/Hrms,50.0);
        kb(i,j) = pow(Dbr(i,j),2.0/3.0)/MAX(exp(x) - 1.0,1.0e-10);
        }
    
    Ure(i,j) = sqrt(Urms(i,j)*Urms(i,j) + 1.45*kb(i,j));
    Uw(i,j)  = sqrt(2.0)*Ure(i,j);
    
        // skewness and asymmetry (Ruessink et al. 2012)
        if(p->S116==1)
        {
        const double T = (Tm10(i,j)>0.0) ? Tm10(i,j) : Tp(i,j);
        
            if(T>0.0)
            {
            const double kh = seastate_wavenumber(2.0*pi/T,h)*h;
            const double Ur = 0.375*Hs(i,j)*(kh/h)/MAX(kh*kh*kh,1.0e-12);
            
            ruessink(Ur,B,psi);
            
            Sk(i,j) = B*cos(psi);
            As(i,j) = B*sin(psi);
            }
        }
    
    ua(i,j) = p->S117*(Sk(i,j) - As(i,j))*Ure(i,j);
    }
    
    // 1D (one cell in y): the component of the wave direction along x
    if(p->j_dir==0)
    SLICELOOP4
    ey(i,j) = 0.0;
    
    pgc->gcsl_start4(p,Ure,1);
    pgc->gcsl_start4(p,Uw,1);
    pgc->gcsl_start4(p,ua,1);
    pgc->gcsl_start4(p,ex,1);
    pgc->gcsl_start4(p,ey,1);
}

void sediment_waves::ruessink(double Ur, double &B, double &psi)
{
    const double pi = 3.14159265358979323846;
    
    Ur = MAX(Ur,1.0e-12);
    
    B   = 0.857/(1.0 + exp((-0.471 - log10(Ur))/0.297));
    psi = (-90.0 + 90.0*tanh(0.815/pow(Ur,0.672)))*pi/180.0;
}

void sediment_waves::abreu(double B, double psi, int nph, double *u)
{
    const double pi = 3.14159265358979323846;
    const double b = sqrt(2.0)*B/sqrt(9.0 + 2.0*B*B);
    const double r = 2.0*b/(1.0 + b*b);
    const double phi = -psi - 0.5*pi;
    const double f = sqrt(1.0 - r*r);
    double mean=0.0, var=0.0;
    
    for(int n=0; n<nph; ++n)
    {
    const double wt = 2.0*pi*double(n)/double(nph);
    u[n] = f*(sin(wt) + r*sin(phi)/(1.0 + f))/(1.0 - r*cos(wt + phi));
    mean += u[n];
    }
    
    mean /= double(nph);
    
    for(int n=0; n<nph; ++n)
    {
    u[n] -= mean;
    var += u[n]*u[n];
    }
    
    var /= double(nph);
    
    const double sc = (var>0.0) ? 1.0/sqrt(var) : 0.0;
    
    for(int n=0; n<nph; ++n)
    u[n] *= sc;
}

void sediment_waves::cycle(lexer *p, int ii, int jj, int nph, double *uw)
{
    const double pi = 3.14159265358979323846;
    const double U = Ure(ii,jj);
    const double B = sqrt(Sk(ii,jj)*Sk(ii,jj) + As(ii,jj)*As(ii,jj));
    
    if(p->S116==1 && B>1.0e-6)
    {
    abreu(B,atan2(As(ii,jj),Sk(ii,jj)),nph,uw);
    
        for(int n=0; n<nph; ++n)
        uw[n] *= U;
    
    return;
    }
    
    for(int n=0; n<nph; ++n)
    uw[n] = sqrt(2.0)*U*sin(2.0*pi*double(n)/double(nph));
}

void sediment_waves::bedshear(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    const double pi = 3.14159265358979323846;
    const double rho = p->W1;
    double u,v,uabs,cphi,sphi,A,z0;
    
    SEDSLICELOOP
    {
    u = 0.5*(s->P(i,j)+s->P(i-1,j));
    v = 0.5*(s->Q(i,j)+s->Q(i,j-1));
    uabs = sqrt(u*u + v*v);
    
    tauc(i,j) = s->tau_eff(i,j);
    tauw(i,j) = fw(i,j) = 0.0;
    taum(i,j) = taumax(i,j) = tauc(i,j);
    
        if(Uw(i,j)>0.0 && Tw(i,j)>0.0)
        {
        A  = Uw(i,j)*Tw(i,j)/(2.0*pi);
        z0 = MAX(s->ks_eff(i,j),1.0e-6)/30.0;
        
        fw(i,j) = MIN(1.39*pow(MAX(A/z0,1.0e-6),-0.52),0.3);
        tauw(i,j) = 0.5*rho*fw(i,j)*Uw(i,j)*Uw(i,j);
        
            if(tauc(i,j)+tauw(i,j)>0.0)
            taum(i,j) = tauc(i,j)*(1.0 + 1.2*pow(tauw(i,j)/(tauc(i,j)+tauw(i,j)),3.2));
        
        cphi = (uabs>1.0e-10) ? (u*ex(i,j) + v*ey(i,j))/uabs : 0.0;
        sphi = sqrt(MAX(1.0 - cphi*cphi,0.0));
        
        taumax(i,j) = sqrt(pow(taum(i,j) + tauw(i,j)*fabs(cphi),2.0) + pow(tauw(i,j)*sphi,2.0));
        }
    
    const double tau = (p->S115==2 || p->S11>=6) ? taumax(i,j) : taum(i,j);
    
    s->tau_eff(i,j) = tau;
    s->tau_i(i,j) = tau;
    s->shearvel_eff(i,j) = sqrt(tau/rho);
    s->shields_eff(i,j) = tau/((p->S22-rho)*fabs(p->W22)*p->S20);
    }
    
    pgc->gcsl_start4(p,s->tau_eff,1);
}

void sediment_waves::eulerian(lexer *p, fdm2D *b, ghostcell *pgc, sediment_fdm *s)
{
    if(p->S110!=1 || p->A751!=1 || p->S119!=0 || b->wave==nullptr)
    return;
    
    const double r = b->wave->stokes_ramp(p);
    
    SLICELOOP1
    if(b->wave->wave_active(i,j)==1 && b->wave->wave_active(i+1,j)==1 && b->hx(i,j)>=p->S128)
    s->P(i,j) -= r*0.5*(b->wave->stokes_x()(i,j) + b->wave->stokes_x()(i+1,j))/b->hx(i,j);
    
    if(p->j_dir==1)
    SLICELOOP2
    if(b->wave->wave_active(i,j)==1 && b->wave->wave_active(i,j+1)==1 && b->hy(i,j)>=p->S128)
    s->Q(i,j) -= r*0.5*(b->wave->stokes_y()(i,j) + b->wave->stokes_y()(i,j+1))/b->hy(i,j);
    
    pgc->gcsl_start1(p,s->P,10);
	pgc->gcsl_start2(p,s->Q,11);
}
