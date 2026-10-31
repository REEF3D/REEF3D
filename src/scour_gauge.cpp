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

#include"scour_gauge.h"
#include"lexer.h"
#include"reefckpt.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<cmath>
#include<cstdio>
#include<iomanip>
#include<iostream>
#include<algorithm>

namespace
{
    const double scour_pi = 3.14159265358979323846;
    const double scour_kappa = 0.4;
}

scour_gauge::scour_gauge(lexer *p, int type_, const char *module) : type(type_)
{
    num = (type==0) ? p->P83 : p->P84;
    x   = (type==0) ? p->P83_x : p->P84_x;
    y   = (type==0) ? p->P83_y : p->P84_y;
    D   = (type==0) ? p->P83_D : p->P84_D;

    // sediment: d50 S 20, density S 22, porosity S 24, critical Shields S 30, angle of repose S 81
    g      = 9.81;
    s      = p->S22/p->W1;
    d50    = p->S20;
    poro   = p->S24;
    thcr   = p->S30;
    cotphi = 1.0/tan(p->S81*scour_pi/180.0);
    ks     = 2.5*d50;
    q0     = sqrt(g*(s-1.0)*d50*d50*d50);

    // P 86: statistics window (0: 50 Tp, set at the first step), Tp, backfilling factor, output interval
    Tw    = p->P86_Tw;
    Tp_in = p->P86_Tp;
    cb    = p->P86_cb>0.0 ? p->P86_cb : 1.0;
    dtout = p->P86_dt;

    // P 87: sediment balance; beta such that the live-bed equilibrium is 1.3 D for KC -> infinity
    amax = std::max(p->P87_amax, 1.0+1.0e-6);
    beta = p->P87_beta>0.0 ? p->P87_beta : log(amax)/1.3;
    kB   = p->P87_k;

    const std::vector<double> zero(num,0.0);
    wsum=mub=mvb=muu=mvv=muc=mvc=maa=mh=ubp=vbp=hasp=tact=S=V=zero;
    nextprint = -1.0e20;
    warned = false;

    out = nullptr;

    if(p->mpirank==0 && num>0)
    {
        char folder[200], name[300];
        snprintf(folder,sizeof(folder),"./REEF3D_%s_Scour",module);
        mkdir(folder,0777);

        out = new std::ofstream[num];

        for(int n=0;n<num;++n)
        {
            snprintf(name,sizeof(name),"%s/REEF3D-%s-Scour-%s-%i.dat",folder,module,(type==0?"SF":"SB"),n+1);

            const bool head = reefckpt_open(out[n],name,p->I45);   // I 45: continued

            if(head)
            {
            out[n]<<(type==0 ? "Scour gauge P 83 (Sumer-Fredsoe incremental) ID: " : "Scour gauge P 84 (sediment balance) ID: ")<<n+1<<std::endl<<std::endl;
            out[n]<<"x_coord \t y_coord \t D"<<std::endl;
            out[n]<<p->Xout(x[n],y[n])<<" \t "<<p->Yout(x[n],y[n])<<" \t "<<D[n]<<std::endl<<std::endl;
            out[n]<<"d50 = "<<d50<<"  s = "<<s<<"  n = "<<poro<<"  theta_cr = "<<thcr<<"  phi = "<<p->S81<<std::endl;
            if(type==1)
            out[n]<<"amax = "<<amax<<"  beta = "<<beta<<"  k = "<<kB<<std::endl;
            out[n]<<std::endl;

            if(type==0)
            out[n]<<"t \t S \t S/D \t S_eq/D \t T \t KC \t theta_w \t theta_c \t Ucw \t Um \t Uc \t Tp"<<std::endl;
            else
            out[n]<<"t \t S \t S/D \t alpha \t theta(t) \t KC \t theta_w \t theta_c \t Ucw \t Um \t Uc \t Tp"<<std::endl;
            }
        }

        std::cout<<"scour gauges "<<(type==0?"P 83":"P 84")<<": "<<num<<std::endl;
    }
}

scour_gauge::~scour_gauge()
{
    if(out!=nullptr)
    {
        for(int n=0;n<num;++n)
        out[n].close();

        delete [] out;
    }
}

// the shape of the equilibrium depth in KC and Ucw (Sumer & Fredsoe 2001), the current limit blended in
// for Ucw > 0.75 (their formula goes to zero for KC -> 0 also in a pure current)
double scour_gauge::seq_shape(double KC, double Ucw) const
{
    const double A = 0.03 + 0.75*pow(Ucw,2.6);
    const double B = 6.0*exp(-4.7*Ucw);

    const double Fw = KC>B ? 1.0-exp(-A*(KC-B)) : 0.0;
    const double wc = std::min(std::max((Ucw-0.75)/0.2,0.0),1.0);

    return std::max(Fw,wc);
}

void scour_gauge::statistics(lexer *p, int n, stats &st) const
{
    const double W  = wsum[n]>0.0 ? wsum[n] : 1.0;
    const double ub = mub[n]/W, vb = mvb[n]/W;
    const double var = std::max(muu[n]/W - ub*ub + mvv[n]/W - vb*vb, 0.0);
    const double cu = muc[n]/W, cv = mvc[n]/W;

    st.h  = mh[n]/W;
    st.Um = sqrt(2.0*var);
    st.Uc = sqrt(cu*cu + cv*cv);
    st.ecx = st.Uc>1.0e-12 ? cu/st.Uc : 0.0;
    st.ecy = st.Uc>1.0e-12 ? cv/st.Uc : 0.0;

    // Tp: P 86, the wave input (peak period, regular waves their period), else from the bed velocity
    // (Tm02 = 2 pi sigma_u/sigma_a, Tp = 1.28 Tm02, JONSWAP)
    st.Tp = 0.0;
    if(Tp_in>0.0)
    st.Tp = Tp_in;
    else if(p->wTp>0.0)
    st.Tp = p->wTp;
    else if(p->wT>0.0)
    st.Tp = p->wT;
    else if(maa[n]>0.0 && var>0.0)
    st.Tp = 1.28*2.0*scour_pi*sqrt(var/(maa[n]/W));

    st.KC = st.Tp>0.0 ? st.Um*st.Tp/D[n] : 0.0;

    // wave friction factor (Sumer et al. 1992), capped for a/ks -> 0
    const double a = st.Um*st.Tp/(2.0*scour_pi);
    st.fw = a>0.0 ? std::min(0.04*pow(a/ks,-0.25),0.3) : 0.0;
    st.theta_w = 0.5*st.fw*st.Um*st.Um/(g*(s-1.0)*d50);

    // current: log law, the velocity at D/2 above the bed (at most the depth)
    const double z0 = ks/30.0;
    const double zc = std::max(std::min(0.5*D[n],st.h),2.0*z0);
    const double ufc = scour_kappa*st.Uc/log(zc/z0);
    st.ufc2 = ufc*ufc;
    st.theta_c = st.ufc2/(g*(s-1.0)*d50);

    st.Ucw = (st.Uc+st.Um)>1.0e-12 ? st.Uc/(st.Uc+st.Um) : 0.0;
    st.F = seq_shape(st.KC,st.Ucw);
    st.r = sqrt((st.theta_w+st.theta_c)/thcr);
}

// the scour depth from the volume of the inverted cone around the pile, V = pi cot (R S^2 + cot S^3/3):
// Newton from an upper bound (V convex in S) converges from above
double scour_gauge::depth_from_volume(double Vol, double R, double S0) const
{
    if(Vol<=0.0)
    return 0.0;

    const double c = cotphi;
    double Sv = std::min(sqrt(Vol/(scour_pi*c*R)), cbrt(3.0*Vol/(scour_pi*c*c)));

    if(S0>0.0 && S0<Sv && scour_pi*c*S0*S0*(R+c*S0/3.0)>=Vol)
    Sv = S0;

    for(int it=0; it<50; ++it)
    {
        const double f  = scour_pi*c*Sv*Sv*(R + c*Sv/3.0) - Vol;
        const double df = scour_pi*c*(2.0*R*Sv + c*Sv*Sv);

        if(df<=0.0)
        break;

        const double dS = f/df;
        Sv -= dS;

        if(fabs(dS)<1.0e-12*(Sv+R))
        break;
    }

    return std::max(Sv,0.0);
}

void scour_gauge::step(lexer *p, const double *smp)
{
    if(num==0)
    return;

    const double dt = p->dt;

    if(dt<=0.0)
    return;

    // the window: 50 Tp (P 86, else the wave input, else 10 s)
    if(Tw<=0.0)
    {
        const double Tpr = Tp_in>0.0 ? Tp_in : (p->wTp>0.0 ? p->wTp : (p->wT>0.0 ? p->wT : 10.0));
        Tw = 50.0*Tpr;

        if(p->mpirank==0)
        std::cout<<"scour gauges "<<(type==0?"P 83":"P 84")<<": statistics window Tw = "<<Tw<<" s"<<std::endl;
    }

    const double wa = std::min(dt/Tw,1.0);

    // output: every step, or every P86_dt seconds
    bool print = dtout<=0.0;

    if(dtout>0.0)
    {
        if(nextprint<-1.0e19)
        nextprint = p->simtime;

        if(p->simtime>=nextprint)
        {
            print = true;
            while(nextprint<=p->simtime)
            nextprint += dtout;
        }
    }

    for(int n=0;n<num;++n)
    {
        const double ub = smp[n*nsmp+0];
        const double vb = smp[n*nsmp+1];
        const double uc = smp[n*nsmp+2];
        const double vc = smp[n*nsmp+3];
        const double h  = smp[n*nsmp+4];

        if(ub<-1.0e19)
        {
            if(!warned && p->mpirank==0)
            std::cout<<"scour gauge "<<(type==0?"P 83 ":"P 84 ")<<n+1<<" lies outside the domain"<<std::endl;
            continue;
        }

        // running (exponentially weighted) means; wsum is the weight of the incomplete window
        auto ew = [wa](double &m, double v) { m = (1.0-wa)*m + wa*v; };

        ew(wsum[n],1.0);
        ew(mub[n],ub);
        ew(mvb[n],vb);
        ew(muu[n],ub*ub);
        ew(mvv[n],vb*vb);
        ew(muc[n],uc);
        ew(mvc[n],vc);
        ew(mh[n],h);

        if(hasp[n]>0.5)
        {
            const double ax = (ub-ubp[n])/dt, ay = (vb-vbp[n])/dt;
            ew(maa[n],ax*ax+ay*ay);
        }
        else
        maa[n] = 0.0;

        ubp[n] = ub;
        vbp[n] = vb;
        hasp[n] = 1.0;

        tact[n] += dt;

        stats st;
        statistics(p,n,st);

        double outa = 0.0, outb = 0.0;

        const bool active = tact[n]>=Tw && h>1.0e-6 && D[n]>0.0;

        // P 83: Sumer-Fredsoe incremental
        if(type==0)
        {
            const double KI = std::min(st.r,1.0);
            const double Seq = 1.3*st.F*KI*D[n];

            const double Tsw = (st.KC>0.0 && st.theta_w>0.0) ? 1.0e-6*pow(st.KC/st.theta_w,3.0) : -1.0;
            const double Tsc = (st.Uc>0.0 && st.theta_c>0.0) ? (1.0/2000.0)*(st.h/D[n])*pow(st.theta_c,-2.2) : -1.0;

            double Ts = -1.0;
            if(Tsw>0.0 && Tsc>0.0)
            Ts = exp((1.0-st.Ucw)*log(Tsw) + st.Ucw*log(Tsc));
            else if(Tsw>0.0)
            Ts = Tsw;
            else if(Tsc>0.0)
            Ts = Tsc;

            const double T = Ts>0.0 ? Ts*D[n]*D[n]/q0 : -1.0;

            // no sediment motion at the pile (U/Ucr < 0.5) or no time scale: the scour hole stays
            if(active && st.r>=0.5 && T>0.0)
            {
                const double Tu = S[n]>Seq ? cb*T : T;
                S[n] = Seq + (S[n]-Seq)*exp(-dt/Tu);
            }

            outa = Seq/D[n];
            outb = T;
        }

        // P 84: sediment balance
        if(type==1)
        {
            const double alpha = pow(amax,st.F)*exp(-beta*S[n]/D[n]);

            // instantaneous bed shear stress: the oscillating part of the bed velocity (wave friction)
            // plus the current (log law), as vectors
            const double ut = ub - mub[n]/wsum[n];
            const double vt = vb - mvb[n]/wsum[n];
            const double um = sqrt(ut*ut + vt*vt);

            const double tx = 0.5*st.fw*um*ut + st.ufc2*st.ecx;
            const double ty = 0.5*st.fw*um*vt + st.ufc2*st.ecy;

            const double theta = sqrt(tx*tx + ty*ty)/(g*(s-1.0)*d50);

            auto Phi = [this](double th) { return th>thcr ? 8.0*pow(th-thcr,1.5) : 0.0; };

            if(active)
            {
                const double dV = dt*kB*D[n]*q0*(Phi(alpha*theta) - Phi(theta))/(1.0-poro);
                V[n] = std::max(V[n]+dV,0.0);
                S[n] = depth_from_volume(V[n],0.5*D[n],S[n]);
            }

            outa = alpha;
            outb = theta;
        }

        if(print && p->mpirank==0)
        {
            out[n]<<std::setprecision(9)<<p->simtime<<" \t "<<std::setprecision(7)<<S[n]<<" \t "<<S[n]/D[n]
                  <<" \t "<<outa<<" \t "<<outb<<" \t "<<st.KC<<" \t "<<st.theta_w<<" \t "<<st.theta_c
                  <<" \t "<<st.Ucw<<" \t "<<st.Um<<" \t "<<st.Uc<<" \t "<<st.Tp<<std::endl;
        }
    }

    warned = true;
}

// checkpoints: num, Tw, nextprint, then 15 values per gauge
void scour_gauge::ckpt_state(std::vector<double> &v) const
{
    v.push_back(num);
    v.push_back(Tw);
    v.push_back(nextprint);

    for(int n=0;n<num;++n)
    {
        const double a[15] = {wsum[n],mub[n],mvb[n],muu[n],mvv[n],muc[n],mvc[n],maa[n],mh[n],ubp[n],vbp[n],hasp[n],tact[n],S[n],V[n]};
        v.insert(v.end(),a,a+15);
    }
}

void scour_gauge::ckpt_restore(const std::vector<double> &v, size_t &pos)
{
    auto nx = [&]() { return pos<v.size() ? v[pos++] : 0.0; };

    const int nc = int(nx());
    const double Twc = nx();
    const double npc = nx();

    if(nc>0)
    {
        if(Twc>0.0)
        Tw = Twc;
        nextprint = npc;
    }

    for(int n=0;n<nc;++n)
    {
        double a[15];
        for(int q=0;q<15;++q)
        a[q] = nx();

        if(n<num)
        {
            wsum[n]=a[0]; mub[n]=a[1]; mvb[n]=a[2]; muu[n]=a[3]; mvv[n]=a[4]; muc[n]=a[5]; mvc[n]=a[6];
            maa[n]=a[7]; mh[n]=a[8]; ubp[n]=a[9]; vbp[n]=a[10]; hasp[n]=a[11]; tact[n]=a[12]; S[n]=a[13]; V[n]=a[14];
        }
    }
}
