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

#include"sediment_susp2D.h"
#include"sediment_inflow.h"
#include"sediment_waves.h"
#include"bedload_waves.h"
#include"sediment_fdm.h"
#include"seastate_sflow.h"
#include"ckpt_registry.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<cmath>
#include<iomanip>

sediment_susp2D::sediment_susp2D(lexer *p) : M(p),C(p),Ceq(p),ED(p),hn(p),rate(p),beta(p),Fx(p),Gx(p),Fy(p),Gy(p),clip(0.0)
{
    SLICEBASELOOP
    M(i,j)=C(i,j)=Ceq(i,j)=ED(i,j)=hn(i,j)=rate(i,j)=beta(i,j)=0.0;
    
    tacc[0] = 0.0;
    tacc[1] = 0.0;
    
    ckid = ckpt_registry::add_raw(p,tacc,2,false);
}

sediment_susp2D::~sediment_susp2D()
{
    ckpt_registry::remove(ckid);
}

void sediment_susp2D::ini(lexer *p, fdm2D *b, ghostcell *pgc, sediment_fdm *s)
{
    if(b->FEacc_x==nullptr)
    {
    b->FEacc_x = new slice1(p);
    b->FEacc_y = new slice2(p);
    }
    
    s->edrate = &rate;
}

double sediment_susp2D::minmod(double a, double b)
{
    if(a*b<=0.0)
    return 0.0;
    
    return (fabs(a)<fabs(b)) ? a : b;
}

double sediment_susp2D::cin(lexer *p, int ii, int jj)
{
    const double c = (pin!=nullptr) ? pin->conc(p) : p->S124;
    
    return (c>=0.0) ? c : Ceq(ii,jj);
}

double sediment_susp2D::rouse(double h, double a, double Z) const
{
    // F = c_a/C: h over the integral of the Rouse profile from a to h (log-spaced Simpson rule)
    if(Z<1.0e-6 || !(a<h) || a<=0.0)
    return 1.0;
    
    const int n = 40;
    const double L = log(h/a);
    double sum = 0.0;
    
    for(int q=0; q<=n; ++q)
    {
    const double sv = double(q)/double(n);
    const double z  = a*exp(sv*L);
    const double f  = (z<h) ? pow((a/(h-a))*((h-z)/z),Z)*z*L : 0.0;
    const double wq = (q==0 || q==n) ? 1.0 : ((q%2==1) ? 4.0 : 2.0);
    sum += wq*f;
    }
    
    const double I = sum/(3.0*double(n));
    
    if(!(I>0.0))
    return 200.0;
    
    return MIN(MAX(h/I,1.0),200.0);
}

void sediment_susp2D::start(lexer *p, fdm2D *b, ghostcell *pgc, sediment_fdm *s, bedload_waves *pbw)
{
    pin = s->pinflow;
    
    // first step: the water depth of the previous step
    if(tacc[1]==0.0)
    {
    SLICELOOP4
    hn(i,j) = b->WL(i,j);
    
    pgc->gcsl_start4(p,hn,1);
    tacc[1] = 1.0;
    }
    
    equilibrium(p,b,pgc,s,pbw);
    
    advection(p,b,pgc,s);
    
    // output |q C| of the advective fluxes (before the diffusion uses Fx, Fy)
    SLICELOOP4
    {
    const double fx = 0.5*(Fx(i,j) + Fx(i-1,j));
    const double fy = (p->j_dir==1) ? 0.5*(Fy(i,j) + Fy(i,j-1)) : 0.0;
    s->qbs(i,j) = sqrt(fx*fx + fy*fy);
    }
    
    diffusion(p,b,pgc);
    
    exchange(p,b,pgc,s);
    
    tacc[0] += p->dt;
    
    // output: C, C_eq
    SLICELOOP4
    {
    s->cb(i,j)  = C(i,j);
    s->cbe(i,j) = Ceq(i,j);
    }
    
    pgc->gcsl_start4(p,s->cb,1);
    pgc->gcsl_start4(p,s->cbe,1);
    pgc->gcsl_start4(p,s->qbs,1);
    
    if(p->count%p->P12==0)
    {
    const double vol = total(p,pgc);
    const double cl  = pgc->globalsum(clip);
    
    // the sediment of the bed (solids, since the start) and of the exchange not yet given to the bed
    double bed=0.0, pend=0.0;
    
    SEDSLICELOOP
    {
    bed  += (1.0-p->S24)*(s->bedzh(i,j)-s->bedzh0(i,j))*p->DXN[IP]*p->DYN[JP];
    pend -= ED(i,j)*p->DXN[IP]*p->DYN[JP];
    }
    
    bed  = pgc->globalsum(bed);
    pend = pgc->globalsum(pend);
    
    if(p->mpirank==0)
    cout<<"suspended load (S 120): volume "<<setprecision(8)<<vol<<" m^3, bed change "<<bed<<" m^3, pending exchange "<<pend
        <<" m^3, balance "<<vol+bed+pend<<" m^3, clipped "<<cl<<" m^3"<<endl;
    }
}

void sediment_susp2D::equilibrium(lexer *p, fdm2D *b, ghostcell *pgc, sediment_fdm *s, bedload_waves *pbw)
{
    const double rho = p->W1;
    const double g = fabs(p->W22);
    const double d = p->S20;
    const double Delta = p->S22/p->W1 - 1.0;
    const double Dstar = d*pow(Delta*g/(p->W2*p->W2),1.0/3.0);
    const double kappa = 0.4;
    sediment_waves *w = s->pwave;
    
    SLICEBASELOOP
    {
    Ceq(i,j) = 0.0;
    beta(i,j) = 0.0;
    }
    
    SLICELOOP4
    {
    const double h = b->WL(i,j);
    
        if(h<p->S128 || p->wet[IJ]==0)
        continue;
    
    const double u = 0.5*(s->P(i,j)+s->P(i-1,j));
    const double v = 0.5*(s->Q(i,j)+s->Q(i,j-1));
    const double U = sqrt(u*u + v*v);
    const double ks = MAX(s->ks(i,j),1.0e-6);
    const double a  = MIN(MAX(0.5*s->ks_eff(i,j),0.01*h),0.5*h);
    
    // u* of the profile: bed shear stress with waves
    double tb = s->tau_eff(i,j);
    
    if(w!=nullptr)
    tb = MAX(tb,w->taumax(i,j));
    
    const double ustar = sqrt(MAX(tb,0.0)/rho);
    const double Z = (ustar>1.0e-8) ? s->ws/(kappa*ustar) : 1.0e3;
    const double F = (p->S121>0.0) ? p->S121 : rouse(h,a,Z);
    
        if(p->S120==1)
        {
        const double fc  = 0.24/pow(log10(MAX(12.0*h/ks,2.0)),2.0);
        double tcw = 0.125*rho*fc*U*U;
        
            if(w!=nullptr && w->Uw(i,j)>0.0 && w->Tw(i,j)>0.0)
            {
            const double A  = w->Uw(i,j)*w->Tw(i,j)/(2.0*3.14159265358979323846);
            const double fw = MIN(exp(-6.0 + 5.2*pow(MAX(A/ks,1.0e-6),-0.19)),0.3);
            const double muw = MIN(MAX(0.7/Dstar,0.14),0.35);     // wave efficiency factor (van Rijn 2007b)
            tcw += muw*0.25*rho*fw*w->Uw(i,j)*w->Uw(i,j);
            }
        
        const double tcr = MAX(s->tau_crit(i,j),1.0e-10);
        const double T = MAX((tcw - tcr)/tcr,0.0);
        const double ca = MIN(0.015*d*pow(Dstar,-0.3)*pow(T,1.5)/a,0.05);
        
        Ceq(i,j) = MIN(ca/F,p->S126);
        }
        
        if(p->S120==2 && pbw!=nullptr)
        Ceq(i,j) = MIN(pbw->Ls(i,j)/h,p->S126);
    
    // exchange with the bed only on erodible bed cells inside the erosion window (S 71, S 72) and outside held cells (S 78)
    const bool bed = p->DFBED[IJ]>0 && s->guard(i,j)>0.0 && p->XP[IP]>=p->S71 && p->XP[IP]<=p->S72;
    
        // no exchange (C_eq is kept for the inflow concentration, S 124 < 0)
        if(!bed)
        continue;
    
    double ws = s->ws;
    
    if(p->S127==1)
    ws *= pow(MAX(1.0-C(i,j),0.0),4.7);
    
    beta(i,j) = ws*F;
    }
    
    pgc->gcsl_start4(p,Ceq,1);
}

void sediment_susp2D::advection(lexer *p, fdm2D *b, ghostcell *pgc, sediment_fdm *s)
{
    const bool acc = (b->FEacc_x!=nullptr && b->FEacc_step==p->count);
    // Stokes transport of the waves (SEASTATE): part of the water flux with the vortex force (A 751 2), already in
    // the host flux with the radiation stress (A 751 1); the sediment moves with the Eulerian flux (S 119 0) or the
    // Lagrangian one (S 119 1)
    const bool stokes = (b->wave!=nullptr);
    const double r = stokes ? b->wave->stokes_ramp(p) : 0.0;
    double ms;
    const int muscl = (p->S125==2) ? 1 : 0;
    sediment_waves *w = s->pwave;
    double q,Cf;
    
    // x faces i+1/2
    for(i=-1; i<p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    {
    q = acc ? (*b->FEacc_x)(i,j) : b->P(i,j)*b->hx(i,j);
    ms = 0.0;
    
        if(stokes && b->wave->wave_active(i,j)==1 && b->wave->wave_active(i+1,j)==1)
        ms = r*0.5*(b->wave->stokes_x()(i,j) + b->wave->stokes_x()(i+1,j));
    
        if(p->A751==2)
        q += ms;
    
    const bool edgeW = (i+p->origin_i==-1);
    const bool edgeE = (i+p->origin_i==p->gknox-1);
    
        // walls and solids (inside the domain)
        if((!edgeW && p->flagslice4[IJ]<0) || (!edgeE && p->flagslice4[Ip1J]<0))
        q = 0.0;
    
    Gx(i,j) = q;
    
        // Eulerian sediment flux: without the Stokes transport
        if(p->S119==0)
        q -= ms;
    
        // asymmetry transport velocity of the waves (S 117, XBeach): the sediment moves with u + ua, the water with u
        if(w!=nullptr && ((edgeW||edgeE) ? q!=0.0 : (p->flagslice4[IJ]>0 && p->flagslice4[Ip1J]>0)) && p->wet[IJ]==1 && p->wet[Ip1J]==1)
        q += 0.5*(w->ua(i,j)*w->ex(i,j)*b->WL(i,j) + w->ua(i+1,j)*w->ex(i+1,j)*b->WL(i+1,j));
    
    Cf = 0.0;
    
        if(q>0.0)
        {
        if(edgeW)
        Cf = cin(p,i+1,j);
        else
        Cf = C(i,j) + muscl*0.5*minmod(C(i,j)-C(i-1,j), C(i+1,j)-C(i,j));
        }
        
        if(q<0.0)
        {
        if(edgeE)
        Cf = cin(p,i,j);
        else
        Cf = C(i+1,j) - muscl*0.5*minmod(C(i+1,j)-C(i,j), C(i+2,j)-C(i+1,j));
        }
    
    Fx(i,j) = q*Cf;
    }
    
    // y faces j+1/2
    for(i=0; i<p->knox; ++i)
    for(j=-1; j<p->knoy; ++j)
    {
    Gy(i,j) = Fy(i,j) = 0.0;
    
        if(p->j_dir==0)
        continue;
    
    q = acc ? (*b->FEacc_y)(i,j) : b->Q(i,j)*b->hy(i,j);
    ms = 0.0;
    
        if(stokes && b->wave->wave_active(i,j)==1 && b->wave->wave_active(i,j+1)==1)
        ms = r*0.5*(b->wave->stokes_y()(i,j) + b->wave->stokes_y()(i,j+1));
    
        if(p->A751==2)
        q += ms;
    
    const bool edgeS = (j+p->origin_j==-1);
    const bool edgeN = (j+p->origin_j==p->gknoy-1);
    
        if((!edgeS && p->flagslice4[IJ]<0) || (!edgeN && p->flagslice4[IJp1]<0))
        q = 0.0;
    
    Gy(i,j) = q;
    
        if(p->S119==0)
        q -= ms;
    
        if(w!=nullptr && ((edgeS||edgeN) ? q!=0.0 : (p->flagslice4[IJ]>0 && p->flagslice4[IJp1]>0)) && p->wet[IJ]==1 && p->wet[IJp1]==1)
        q += 0.5*(w->ua(i,j)*w->ey(i,j)*b->WL(i,j) + w->ua(i,j+1)*w->ey(i,j+1)*b->WL(i,j+1));
    
    Cf = 0.0;
    
        if(q>0.0)
        {
        if(edgeS)
        Cf = cin(p,i,j+1);
        else
        Cf = C(i,j) + muscl*0.5*minmod(C(i,j)-C(i,j-1), C(i,j+1)-C(i,j));
        }
        
        if(q<0.0)
        {
        if(edgeN)
        Cf = cin(p,i,j);
        else
        Cf = C(i,j+1) - muscl*0.5*minmod(C(i,j+1)-C(i,j), C(i,j+2)-C(i,j+1));
        }
    
    Fy(i,j) = q*Cf;
    }
    
    // conservative update, with the residual of the continuity carried by the cell concentration
    SLICELOOP4
    {
    const double divF = (Fx(i,j)-Fx(i-1,j))/p->DXN[IP] + (Fy(i,j)-Fy(i,j-1))/p->DYN[JP];
    const double divG = (Gx(i,j)-Gx(i-1,j))/p->DXN[IP] + (Gy(i,j)-Gy(i,j-1))/p->DYN[JP];
    const double h1 = b->WL(i,j);
    
    M(i,j) += -p->dt*divF + C(i,j)*(h1 - hn(i,j) + p->dt*divG);
    
        if(M(i,j)<0.0)
        {
        clip += -M(i,j)*p->DXN[IP]*p->DYN[JP];
        M(i,j) = 0.0;
        }
    
    C(i,j) = (h1>=p->S128 && p->wet[IJ]==1) ? M(i,j)/h1 : 0.0;
    }
    
    pgc->gcsl_start4(p,C,1);
}

void sediment_susp2D::diffusion(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // K = S 122 + nu_t/S 123, explicit, in sub-steps
    double kmax = 0.0, rmax = 0.0;
    
    SLICELOOP4
    {
    const double K = p->S122 + MAX(b->eddyv(i,j),0.0)/MAX(p->S123,1.0e-10);
    const double inv = 1.0/(p->DXN[IP]*p->DXN[IP]) + (p->j_dir==1 ? 1.0/(p->DYN[JP]*p->DYN[JP]) : 0.0);
    kmax = MAX(kmax,K);
    rmax = MAX(rmax,K*inv);
    }
    
    kmax = pgc->globalmax(kmax);
    rmax = pgc->globalmax(rmax);
    
    if(!(kmax>0.0))
    return;
    
    const int nsub = MAX(1,int(ceil(p->dt*rmax/0.45)));
    const double dts = p->dt/double(nsub);
    
    for(int n=0; n<nsub; ++n)
    {
        // x faces
        for(i=-1; i<p->knox; ++i)
        for(j=0; j<p->knoy; ++j)
        {
        Fx(i,j) = 0.0;
        
        const bool edge = (i+p->origin_i==-1 || i+p->origin_i==p->gknox-1);
        
            if(edge || p->flagslice4[IJ]<0 || p->flagslice4[Ip1J]<0 || p->wet[IJ]==0 || p->wet[Ip1J]==0)
            continue;
        
        const double hf = MIN(b->WL(i,j),b->WL(i+1,j));
        const double K  = p->S122 + 0.5*(MAX(b->eddyv(i,j),0.0)+MAX(b->eddyv(i+1,j),0.0))/MAX(p->S123,1.0e-10);
        Fx(i,j) = hf*K*(C(i+1,j)-C(i,j))/p->DXP[IP];
        }
        
        // y faces
        for(i=0; i<p->knox; ++i)
        for(j=-1; j<p->knoy; ++j)
        {
        Fy(i,j) = 0.0;
        
        const bool edge = (j+p->origin_j==-1 || j+p->origin_j==p->gknoy-1);
        
            if(p->j_dir==0 || edge || p->flagslice4[IJ]<0 || p->flagslice4[IJp1]<0 || p->wet[IJ]==0 || p->wet[IJp1]==0)
            continue;
        
        const double hf = MIN(b->WL(i,j),b->WL(i,j+1));
        const double K  = p->S122 + 0.5*(MAX(b->eddyv(i,j),0.0)+MAX(b->eddyv(i,j+1),0.0))/MAX(p->S123,1.0e-10);
        Fy(i,j) = hf*K*(C(i,j+1)-C(i,j))/p->DYP[JP];
        }
        
        SLICELOOP4
        {
        const double h1 = b->WL(i,j);
        
        M(i,j) += dts*((Fx(i,j)-Fx(i-1,j))/p->DXN[IP] + (Fy(i,j)-Fy(i,j-1))/p->DYN[JP]);
        
            if(M(i,j)<0.0)
            {
            clip += -M(i,j)*p->DXN[IP]*p->DYN[JP];
            M(i,j) = 0.0;
            }
        
        C(i,j) = (h1>=p->S128 && p->wet[IJ]==1) ? M(i,j)/h1 : 0.0;
        }
        
        pgc->gcsl_start4(p,C,1);
    }
}

void sediment_susp2D::exchange(lexer *p, fdm2D *b, ghostcell *pgc, sediment_fdm *s)
{
    SLICELOOP4
    {
    const double h1 = b->WL(i,j);
    
        // dry: the sediment of the column is deposited on an erodible bed; cells without one (structures, held
        // cells S 78, outside S 71 - S 72) keep it until they are wet again
        if(h1<p->S128 || p->wet[IJ]==0)
        {
        const bool bed = p->DFBED[IJ]>0 && s->guard(i,j)>0.0 && p->XP[IP]>=p->S71 && p->XP[IP]<=p->S72;
        
            if(bed)
            {
            ED(i,j) -= M(i,j);
            M(i,j) = 0.0;
            }
        
        C(i,j) = 0.0;
        hn(i,j) = h1;
        continue;
        }
    
    const double bt = beta(i,j)*p->dt;
    const double Mn = (M(i,j) + bt*Ceq(i,j))/(1.0 + bt/h1);
    
    ED(i,j) += Mn - M(i,j);
    M(i,j) = Mn;
    C(i,j) = M(i,j)/h1;
    hn(i,j) = h1;
    }
    
    pgc->gcsl_start4(p,M,1);
    pgc->gcsl_start4(p,C,1);
    pgc->gcsl_start4(p,hn,1);
}

void sediment_susp2D::bed_rate(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    SLICELOOP4
    rate(i,j) = (tacc[0]>0.0) ? ED(i,j)/tacc[0] : 0.0;
    
    pgc->gcsl_start4(p,rate,1);
    
    s->edrate = &rate;
}

void sediment_susp2D::bed_done(lexer *p)
{
    SLICEBASELOOP
    ED(i,j) = 0.0;
    
    tacc[0] = 0.0;
}

double sediment_susp2D::total(lexer *p, ghostcell *pgc)
{
    double vol = 0.0;
    
    SLICELOOP4
    vol += M(i,j)*p->DXN[IP]*p->DYN[JP];
    
    return pgc->globalsum(vol);
}
