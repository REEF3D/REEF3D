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

#include"bedload_waves.h"
#include"sediment_waves.h"
#include"sediment_fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include<cmath>

bedload_waves::bedload_waves(lexer *p) : tx(p),ty(p),Ls(p)
{
    SLICEBASELOOP
    tx(i,j)=ty(i,j)=Ls(i,j)=0.0;
}

void bedload_waves::start(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double qx,qy,qa;
    
    g = fabs(p->W22);
    d = s->dk;
    Delta = p->S22/p->W1 - 1.0;
    Dstar = d*pow(Delta*g/(p->W2*p->W2),1.0/3.0);
    
    SLICEBASELOOP
    {
    s->qbe(i,j) = 0.0;
    Ls(i,j) = 0.0;
    }
    
    SEDSLICELOOP
    {
    u = 0.5*(s->P(i,j)+s->P(i-1,j));
    v = 0.5*(s->Q(i,j)+s->Q(i,j-1));
    h = s->waterlevel(i,j);
    
    qx = qy = 0.0;
    
        if(s->active(i,j)==1 && h>=p->S128)
        {
        if(p->S11==6)
        soulsby_vanrijn(p,s,qx,qy);
        
        if(p->S11==7)
        vanrijn2007(p,s,qx,qy);
        
        if(p->S11==8)
        bailard(p,s,qx,qy);
        
        if(p->S11==9)
        bijker(p,s,qx,qy);
        }
    
    qa = sqrt(qx*qx + qy*qy);
    
    s->qbe(i,j) = qa;
    
        // direction of the transport; the current where there is none
        if(qa>0.0)
        {
        tx(i,j) = qx/qa;
        ty(i,j) = qy/qa;
        }
        else
        {
        const double ua = sqrt(u*u + v*v);
        tx(i,j) = (ua>1.0e-10) ? u/ua : 0.0;
        ty(i,j) = (ua>1.0e-10) ? v/ua : 0.0;
        }
    }
    
    pgc->gcsl_start4(p,s->qbe,1);
    pgc->gcsl_start4(p,tx,1);
    pgc->gcsl_start4(p,ty,1);
    pgc->gcsl_start4(p,Ls,1);
}

void bedload_waves::soulsby_vanrijn(lexer *p, sediment_fdm *s, double &qx, double &qy)
{
    sediment_waves *w = s->pwave;
    const double ua  = (w!=nullptr) ? w->ua(i,j) : 0.0;
    const double ex  = (w!=nullptr) ? w->ex(i,j) : 0.0;
    const double ey  = (w!=nullptr) ? w->ey(i,j) : 0.0;
    const double Ure = (w!=nullptr) ? w->Ure(i,j) : 0.0;
    
    const double utx = u + ua*ex;
    const double uty = v + ua*ey;
    const double Ut  = sqrt(utx*utx + uty*uty);
    
    const double z0 = 0.006;
    const double CD = pow(0.40/(log(MAX(h/z0,7.38905609893065)) - 1.0),2.0);
    const double d90 = 2.0*d;
    const double lg = log10(MAX(4.0*h/d90,1.0));
    const double Ucr = (d<=0.0005) ? 0.19*pow(d,0.1)*lg : 8.5*pow(d,0.6)*lg;
    
    const double stir = sqrt(Ut*Ut + 0.018/CD*Ure*Ure) - Ucr;
    
    if(stir<=0.0)
    return;
    
    const double den = pow(Delta*g*d,1.2);
    const double Asb = 0.005*h*pow(d/h,1.2)/den;
    const double Ass = 0.012*d*pow(Dstar,-0.6)/den;
    const double st  = pow(stir,2.4);
    
    Ls(i,j) = Ass*st;
    
    // bed part, plus the suspended part without a suspended load model (SFLOW S 120, NHFLOW S 12)
    const double A = (p->S120>0 || (p->A10==5 && p->S12>0)) ? Asb : Asb + Ass;
    
    qx = A*st*utx;
    qy = A*st*uty;
}

void bedload_waves::vanrijn2007(lexer *p, sediment_fdm *s, double &qx, double &qy)
{
    sediment_waves *w = s->pwave;
    const int nph = sediment_waves::nphase;
    double uw[sediment_waves::nphase];
    
    const double ks = MAX(s->ks(i,j),1.0e-6);
    const double fc = 0.24/pow(log10(MAX(12.0*h/ks,2.0)),2.0);
    const double taucr = MAX(s->tau_crit(i,j),1.0e-10);
    const double rho = p->W1;
    double ex=0.0, ey=0.0, fw=0.0;
    
    if(w!=nullptr && w->Uw(i,j)>0.0 && w->Tw(i,j)>0.0)
    {
    const double A = w->Uw(i,j)*w->Tw(i,j)/(2.0*3.14159265358979323846);
    fw = MIN(exp(-6.0 + 5.2*pow(MAX(A/ks,1.0e-6),-0.19)),0.3);
    ex = w->ex(i,j);
    ey = w->ey(i,j);
    w->cycle(p,i,j,nph,uw);
    }
    else
    for(int n=0; n<nph; ++n)
    uw[n] = 0.0;
    
    const double sc = sqrt(fc/8.0), sw = sqrt(fw/2.0);
    const double fac = 0.5*d*pow(Dstar,-0.3);
    
    for(int n=0; n<nph; ++n)
    {
    const double usx = sc*u + sw*uw[n]*ex;
    const double usy = sc*v + sw*uw[n]*ey;
    const double us  = sqrt(usx*usx + usy*usy);
    const double T   = (rho*us*us - taucr)/taucr;
    
        if(T>0.0 && us>0.0)
        {
        const double q = fac*us*pow(T,1.5);
        qx += q*usx/us;
        qy += q*usy/us;
        }
    }
    
    qx /= double(nph);
    qy /= double(nph);
}

void bedload_waves::bailard(lexer *p, sediment_fdm *s, double &qx, double &qy)
{
    sediment_waves *w = s->pwave;
    const int nph = sediment_waves::nphase;
    double uw[sediment_waves::nphase];
    
    const double ks = MAX(s->ks(i,j),1.0e-6);
    const double fc = 0.24/pow(log10(MAX(12.0*h/ks,2.0)),2.0);
    const double U  = sqrt(u*u + v*v);
    double ex=0.0, ey=0.0, fw=0.0, Uwa=0.0;
    
    if(w!=nullptr && w->Uw(i,j)>0.0 && w->Tw(i,j)>0.0)
    {
    Uwa = w->Uw(i,j);
    const double A = Uwa*w->Tw(i,j)/(2.0*3.14159265358979323846);
    fw = MIN(exp(-6.0 + 5.2*pow(MAX(A/ks,1.0e-6),-0.19)),0.3);
    ex = w->ex(i,j);
    ey = w->ey(i,j);
    w->cycle(p,i,j,nph,uw);
    }
    else
    for(int n=0; n<nph; ++n)
    uw[n] = 0.0;
    
    const double a  = (U+Uwa>0.0) ? U/(U+Uwa) : 1.0;
    const double cf = a*fc/8.0 + (1.0-a)*fw/2.0;
    const double eb = 0.1, es = 0.02, tanphi = 0.63;
    const double ws = MAX(s->ws,1.0e-10);
    
    double m3x=0.0, m3y=0.0, m4x=0.0, m4y=0.0, a3=0.0;
    
    for(int n=0; n<nph; ++n)
    {
    const double ux = u + uw[n]*ex;
    const double uy = v + uw[n]*ey;
    const double ua = sqrt(ux*ux + uy*uy);
    
    m3x += ua*ua*ux;
    m3y += ua*ua*uy;
    m4x += ua*ua*ua*ux;
    m4y += ua*ua*ua*uy;
    a3  += ua*ua*ua;
    }
    
    m3x/=double(nph); m3y/=double(nph); m4x/=double(nph); m4y/=double(nph); a3/=double(nph);
    
    const double cb = cf*eb/(Delta*g*tanphi);
    const double cs = cf*es/(Delta*g*ws);
    
    Ls(i,j) = cs*a3;
    
    qx = cb*m3x;
    qy = cb*m3y;
    
    // suspended part without a suspended load model (SFLOW S 120, NHFLOW S 12)
    if(p->S120==0 && !(p->A10==5 && p->S12>0))
    {
    qx += cs*m4x;
    qy += cs*m4y;
    }
}

void bedload_waves::bijker(lexer *p, sediment_fdm *s, double &qx, double &qy)
{
    sediment_waves *w = s->pwave;
    const double U = sqrt(u*u + v*v);
    
    if(U<=1.0e-10)
    return;
    
    const double ks  = MAX(s->ks_eff(i,j),1.0e-6);
    const double d90 = 2.0*d;
    const double C   = 18.0*log10(MAX(12.0*h/ks,2.0));
    const double C90 = 18.0*log10(MAX(12.0*h/d90,2.0));
    const double mu  = MIN(pow(C/C90,1.5),1.0);
    double xi=0.0, Uw=0.0;
    
    if(w!=nullptr && w->Uw(i,j)>0.0 && w->Tw(i,j)>0.0)
    {
    Uw = w->Uw(i,j);
    const double A  = Uw*w->Tw(i,j)/(2.0*3.14159265358979323846);
    const double fw = MIN(exp(-5.977 + 5.213*pow(MAX(A/ks,1.0e-6),-0.194)),0.3);
    xi = C*sqrt(fw/(2.0*g));
    }
    
    const double b = 2.0;
    const double q = b*d*sqrt(g)*U/C*exp(-0.27*Delta*d*C*C/(mu*(U*U + 0.5*xi*xi*Uw*Uw)));
    
    qx = q*u/U;
    qy = q*v/U;
}
