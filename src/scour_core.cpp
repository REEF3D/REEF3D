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

#include"scour_core.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<sys/stat.h>
#include<sys/types.h>
#include<cstdio>

scour_core::scour_core(lexer *p, ghostcell *pgc, const char *solver_, slice *depth_, bool computed_current) : ub(p),vb(p),
                       Um(p),ubmax(p),KC(p),Uc(p),Ucw(p),SD(p),S(p),
                       thmax(p),thcr(p),mob(p),armSR(p),armdreq(p),armdmax(p),
                       ucx(p),ucy(p),
                       Ucx(p),Ucy(p),su(p),sv(p),suu(p),svv(p),suv(p),
                       suc(p),svc(p),uhi(p),ulo(p),vhi(p),vlo(p),
                       solver(solver_),depth(depth_),has_current(computed_current)
{
    g   = fabs(p->W22)>1.0e-10 ? fabs(p->W22) : 9.81;
    rho = p->W1;
    nu  = p->W2;

    p_aamp = p->P197_aamp;
    p_kt   = p->P197_kt;
    p_gE   = p->P197_gE;
    p_nr   = p->P197_nr;
    p_rhos = p->P197_rho;

    wsum = 0.0;
    warned_T = false;

    SLICELOOP4
    {
    ub(i,j)=vb(i,j)=0.0;
    Um(i,j)=ubmax(i,j)=KC(i,j)=Uc(i,j)=Ucw(i,j)=SD(i,j)=S(i,j)=0.0;
    thmax(i,j)=thcr(i,j)=mob(i,j)=armSR(i,j)=armdreq(i,j)=armdmax(i,j)=0.0;
    su(i,j)=sv(i,j)=suu(i,j)=svv(i,j)=suv(i,j)=0.0;
    ucx(i,j)=ucy(i,j)=suc(i,j)=svc(i,j)=0.0;
    uhi(i,j)=vhi(i,j)=-1.0e20;
    ulo(i,j)=vlo(i,j)=1.0e20;
    }

    // ambient current: constant (P 189 U_c direction), the direction it flows to, counterclockwise
    // from the x-axis in degrees
    const double dir = p->P189_dir*(PI/180.0);

    SLICELOOP4
    {
    Ucx(i,j) = p->P189_U*cos(dir);
    Ucy(i,j) = p->P189_U*sin(dir);
    }

    // or from the file scour_current.dat (P 196 1)
    if(p->P196==1)
    current_file(p,pgc);

    // P 198: scour gauges
    gauge_ini(p);

    // P 193: scour protection gauges
    pgauge_ini(p);

    if(p->mpirank==0)
    {
    cout<<solver<<" scour (P 187): D = "<<p->P188_D<<" m, statistics from t = "<<p->P188_ts<<" s, ";
    if(p->P188_T>0.0)
    cout<<"T = "<<p->P188_T<<" s"<<endl;
    else
    cout<<"T from the wave input"<<endl;

    if(p->P196==1)
    cout<<solver<<" scour (P 196): ambient current from scour_current.dat"<<endl;
    else
    cout<<solver<<" scour (P 189): ambient current "<<p->P189_U<<" m/s towards "<<p->P189_dir<<" deg"<<endl;

    if(has_current)
    cout<<solver<<" scour: plus the computed current (time mean of the depth-averaged velocity)"<<endl;

    cout<<solver<<" scour: sediment d50 = "<<p->S20<<" m (S 20), density "<<p->S22<<" kg/m3 (S 22); "
        <<"armour d = "<<p->P197_d<<" m, density "<<p->P197_rho<<" kg/m3 (P 197)"<<endl;
    cout<<solver<<" scour protection (P 197, Soulsby as REEF_Tools): alpha "<<p->P197_aamp<<", kt "<<p->P197_kt
        <<", gamma_E "<<p->P197_gE<<", ks = "<<p->P197_nr<<" D50, Shields 0.055"<<endl;

    if(p->P180!=1 && p->P187==1)
    cout<<solver<<" scour: P 187 writes the bed with the free-surface output, which needs P 180 1"<<endl;
    }
}

void scour_core::current_file(lexer *p, ghostcell *pgc)
{
    // x y Ux Uy per line, in the coordinates of the output (as the bed file); # comments.
    // Inverse distance weighting (power 2) of all points onto the cells.
    std::vector<double> X,Y,U,V;

    std::ifstream in("scour_current.dat");

    if(!in.is_open())
    {
        if(p->mpirank==0)
        cout<<solver<<" scour: scour_current.dat not found, P 189 current used"<<endl;
        return;
    }

    std::string line;
    while(std::getline(in,line))
    {
        const size_t q = line.find_first_not_of(" \t\r");
        if(q==std::string::npos || line[q]=='#')
        continue;

        std::istringstream ls(line);
        double x,y,u,v;
        if(ls>>x>>y>>u>>v)
        {
        X.push_back(x);
        Y.push_back(y);
        U.push_back(u);
        V.push_back(v);
        }
    }
    in.close();

    const int np = int(X.size());

    if(np==0)
    {
        if(p->mpirank==0)
        cout<<solver<<" scour: scour_current.dat holds no points, P 189 current used"<<endl;
        return;
    }

    SLICELOOP4
    {
        const double xc = p->Xout(p->XP[IP],p->YP[JP]);
        const double yc = p->Yout(p->XP[IP],p->YP[JP]);

        double ws=0.0, us=0.0, vs=0.0;
        int hit=-1;

        for(int n=0; n<np; ++n)
        {
            const double r2 = (X[n]-xc)*(X[n]-xc) + (Y[n]-yc)*(Y[n]-yc);

            if(r2<1.0e-12)
            {
            hit=n;
            break;
            }

            const double w = 1.0/r2;
            ws += w;
            us += w*U[n];
            vs += w*V[n];
        }

        if(hit>=0)
        {
        Ucx(i,j) = U[hit];
        Ucy(i,j) = V[hit];
        }
        else
        {
        Ucx(i,j) = us/ws;
        Ucy(i,j) = vs/ws;
        }
    }

    if(p->mpirank==0)
    cout<<solver<<" scour: "<<np<<" current points read from scour_current.dat"<<endl;
}

void scour_core::sample(lexer *p, ghostcell *pgc)
{
    if(p->simtime<p->P188_ts || p->dt<=0.0)
    return;

    const double w = p->dt;

    // P 199 Tw > 0: exponentially weighted statistics with the time scale Tw, else all from t_s on
    if(p->P199_Tw>0.0)
    {
    const double f = exp(-p->dt/p->P199_Tw);

    wsum *= f;

    SLICELOOP4
    {
    su(i,j)  *= f;
    sv(i,j)  *= f;
    suu(i,j) *= f;
    svv(i,j) *= f;
    suv(i,j) *= f;
    suc(i,j) *= f;
    svc(i,j) *= f;
    }
    }

    wsum += w;

    SLICELOOP4
    WETDRY
    {
    const double u = ub(i,j);
    const double v = vb(i,j);

    su(i,j)  += w*u;
    sv(i,j)  += w*v;
    suu(i,j) += w*u*u;
    svv(i,j) += w*v*v;
    suv(i,j) += w*u*v;

    if(has_current)
    {
    suc(i,j) += w*ucx(i,j);
    svc(i,j) += w*ucy(i,j);
    }

    ubmax(i,j) = MAX(ubmax(i,j),sqrt(u*u + v*v));

    uhi(i,j) = MAX(uhi(i,j),u);
    ulo(i,j) = MIN(ulo(i,j),u);
    vhi(i,j) = MAX(vhi(i,j),v);
    vlo(i,j) = MIN(vlo(i,j),v);
    }
}

double scour_core::period(lexer *p)
{
    double T = p->P188_T;

    if(T<=0.0)
    T = p->wTp>0.0 ? p->wTp : p->wT;

    if(T<=0.0 && !warned_T && p->mpirank==0)
    {
    cout<<solver<<" scour: no wave period (P 188 T or the wave input), KC is zero"<<endl;
    warned_T = true;
    }

    return MAX(T,0.0);
}

double scour_core::thetacr(double d, double s) const
{
    const double Dstar = pow(g*(s-1.0)/(nu*nu),1.0/3.0)*d;

    return 0.30/(1.0 + 1.2*Dstar) + 0.055*(1.0 - exp(-0.020*Dstar));
}

double scour_core::taumax(double Uc_, double Um_, double T, double h, double cosphi, double d) const
{
    const double z0 = d/12.0;   // k_s = 2.5 d50, z0 = k_s/30

    // current: depth-averaged velocity, logarithmic profile
    double lg = log(MAX(h,1.0e-3)/z0) - 1.0;
    lg = MAX(lg,1.0);
    const double CD = pow(0.40/lg,2.0);
    const double tc = rho*CD*Uc_*Uc_;

    // waves: rough and smooth (laminar or turbulent) friction factor, the larger one
    double tw = 0.0;

    if(Um_>1.0e-10 && T>0.0)
    {
    const double A   = Um_*T/(2.0*PI);
    const double fwr = 1.39*pow(A/z0,-0.52);
    const double Rew = Um_*A/nu;
    const double fws = Rew<=5.0e5 ? 2.0*pow(Rew,-0.5) : 0.0521*pow(Rew,-0.187);

    tw = 0.5*rho*MAX(fwr,fws)*Um_*Um_;
    }

    if(tc+tw<=0.0)
    return 0.0;

    const double tm = tc*(1.0 + 1.2*pow(tw/(tc+tw),3.2));
    const double cp = MIN(fabs(cosphi),1.0);
    const double sp = sqrt(MAX(1.0-cp*cp,0.0));

    return sqrt(pow(tm + tw*cp,2.0) + pow(tw*sp,2.0));
}

// the scour parameters of the cell (i,j) in the order of names()
void scour_core::cell(lexer *p, int ii, int jj, double T, double *q)
{
    i=ii;
    j=jj;

    const double D  = MAX(p->P188_D,1.0e-10);

    const double d50 = MAX(p->S20,1.0e-6);
    const double s   = p->S22/rho;
    const double tcs = thetacr(d50,s);

    const double da  = MAX(p->P197_d,1.0e-6);

    // ambient current, plus the time mean of the computed depth-averaged current (NHFLOW)
    double ux = Ucx(i,j);
    double uy = Ucy(i,j);

    if(has_current && wsum>0.0)
    {
    ux += suc(i,j)/wsum;
    uy += svc(i,j)/wsum;
    }

    const double uc = sqrt(ux*ux + uy*uy);

    for(int n=0; n<NQ; ++n)
    q[n]=0.0;

    q[0]  = ub(i,j);
    q[1]  = vb(i,j);
    q[3]  = ubmax(i,j);
    q[5]  = uc;
    q[10] = tcs;
    q[14] = armdmax(i,j);

    // dry, or no statistics yet: before P 188 t_s and during the first wave period after it
    if(p->wet[IJ]==0 || wsum<=0.0 || wsum<T)
    return;

    // orbital velocity amplitude and principal direction of the bed velocity
    double um=0.0, cphw=1.0, sphw=0.0;

    if(wsum>0.0)
    {
    const double mu = su(i,j)/wsum;
    const double mv = sv(i,j)/wsum;
    const double vu = MAX(suu(i,j)/wsum - mu*mu,0.0);
    const double vv = MAX(svv(i,j)/wsum - mv*mv,0.0);
    const double cv = suv(i,j)/wsum - mu*mv;

    um = sqrt(2.0*(vu + vv));

    // P 199 m 1: the maximum orbital velocity (regular waves), the larger excursion about the mean
    if(p->P199_m==1 && uhi(i,j)>-1.0e19)
    {
    const double eu = MAX(uhi(i,j)-mu, mu-ulo(i,j));
    const double ev = MAX(vhi(i,j)-mv, mv-vlo(i,j));
    um = sqrt(eu*eu + ev*ev);
    }

    const double phw = 0.5*atan2(2.0*cv,vu-vv);
    cphw = cos(phw);
    sphw = sin(phw);
    }

    // angle between the current and the wave axis
    const double cosphi = uc>1.0e-10 ? (ux*cphw + uy*sphw)/uc : 1.0;

    const double kc  = um*T/D;
    const double ucw = (uc+um)>1.0e-10 ? uc/(uc+um) : 0.0;

    // Sumer & Fredsoe (2002), live bed
    double sd = 0.0;

    if(ucw>=0.99 && uc>1.0e-10)
    sd = 1.3;
    else if(um>1.0e-10)
    {
    const double A = 0.03 + 0.75*pow(ucw,2.6);
    const double B = 6.0*exp(-4.7*ucw);

    sd = kc>B ? 1.3*(1.0 - exp(-A*(kc-B))) : 0.0;
    }

    // bed sediment: maximum Shields parameter of waves and current, mobility
    const double h = MAX((*depth)(i,j),1.0e-3);   // still water depth

    const double th = taumax(uc,um,T,h,cosphi,d50)/(rho*(s-1.0)*g*d50);

    // scour protection (Soulsby, as REEF_Tools): stability ratio of the stone d (> 1: not stable),
    // the required D50 and its maximum over time
    const protection pa = protect(um,uc,T,h,da);
    const double dq = protect_dreq(um,uc,T,h,da);

    q[2]  = um;
    q[4]  = kc;
    q[6]  = ucw;
    q[7]  = sd;
    q[8]  = sd*D;
    q[9]  = th;
    q[11] = th/tcs;
    q[12] = pa.sr;
    q[13] = dq;
    q[14] = wsum>=NACC*T ? MAX(armdmax(i,j),dq) : armdmax(i,j);
}

// scour protection after Soulsby (1997), the procedure of REEF_Tools: the stresses for the stone size d
scour_core::protection scour_core::protect(double Um_, double Uc_, double T, double h, double d) const
{
    protection r = {0.0,0.0,0.0,0.0,0.0,0.0};

    const double amp = p_aamp*p_kt;
    const double ub  = Um_*amp;
    const double uc  = Uc_*amp;
    const double ks  = p_nr*d;

    // waves: rough friction factor, A_m/z0 with z0 = ks/30
    if(ub>1.0e-10 && T>0.0)
    {
    const double omega = 2.0*PI/T;
    const double fw = 1.39*pow(30.0*ub/(omega*ks),-0.52);
    r.tw = 0.5*rho*fw*ub*ub;
    }

    // current: Manning n from D50, Chezy, depth-averaged velocity
    if(uc>1.0e-10)
    {
    const double n = pow(2.5*d,1.0/6.0)/20.0;
    const double C = pow(h,1.0/6.0)/n;
    r.tc = rho*(g/(C*C))*uc*uc;
    }

    if(r.tw + r.tc > 0.0)
    r.tm = r.tc*(1.0 + 1.2*pow(r.tw/(r.tc+r.tw),3.2));

    // collinear waves and current (phi = 0), environmental factor
    r.tmax = p_gE*(r.tm + r.tw);
    r.tcr  = 0.055*(p_rhos-rho)*g*d;
    r.sr   = r.tcr>0.0 ? r.tmax/r.tcr : 0.0;

    return r;
}

// the D50 for which gamma_E tau_max = tau_cr (fixed point: tau_max depends weakly on D50 through ks and n)
double scour_core::protect_dreq(double Um_, double Uc_, double T, double h, double d0) const
{
    double d = MAX(d0,1.0e-4);

    for(int it=0; it<60; ++it)
    {
    const protection r = protect(Um_,Uc_,T,h,d);

    if(r.tmax<=0.0)
    return 0.0;

    const double dn = MIN(MAX(r.tmax/(0.055*(p_rhos-rho)*g),1.0e-6),50.0);

    if(fabs(dn-d)<1.0e-7*d)
    return dn;

    d = dn;
    }

    return d;
}

const char *scour_core::name(int n)
{
    static const char *const names[NQ] = {"bedvel_u","bedvel_v","Um","ubmax","KC","Uc","Ucw","S_D","S",
                                          "theta_max","theta_cr","mobility","armour_SR","armour_d_req",
                                          "armour_d_req_max"};
    return names[n];
}

slice *scour_core::field(int n)
{
    slice *const f[NQ] = {&ub,&vb,&Um,&ubmax,&KC,&Uc,&Ucw,&SD,&S,&thmax,&thcr,&mob,&armSR,&armdreq,&armdmax};
    return f[n];
}

void scour_core::evaluate(lexer *p, ghostcell *pgc)
{
    const double T = period(p);
    double q[NQ];

    SLICELOOP4
    {
        cell(p,i,j,T,q);

        for(int n=2; n<NQ; ++n)   // the bed velocity and its maximum are the fields themselves
        if(n!=3)
        (*field(n))(i,j) = q[n];
    }

    for(int n=0; n<NQ; ++n)
    pgc->gcsl_start4(p,*field(n),1);
}

// P 198: the scour gauges, the cell that contains the gauge (as the wave gauges P 51)
void scour_core::gauge_ini(lexer *p)
{
    ngauge = p->P198;

    if(ngauge<=0)
    return;

    gi.resize(ngauge);
    gj.resize(ngauge);
    gflag.resize(ngauge);
    gval.resize(size_t(ngauge)*NQ);
    gmax.assign(ngauge,0.0);

    for(int n=0; n<ngauge; ++n)
    {
    gi[n] = p->posc_i(p->P198_x[n]);
    gj[n] = p->j_dir==1 ? p->posc_j(p->P198_y[n]) : 0;

    gflag[n] = (gi[n]>=0 && gi[n]<p->knox && gj[n]>=0 && gj[n]<p->knoy) ? 1 : 0;
    }

    nextgauge = p->simtime;

    if(p->mpirank!=0)
    return;

    char gdir[128];
    snprintf(gdir,sizeof(gdir),"./REEF3D_%s_Scour",solver);
    mkdir(gdir,0777);

    gout.resize(ngauge);

    for(int n=0; n<ngauge; ++n)
    {
    char gname[256];
    snprintf(gname,sizeof(gname),"./REEF3D_%s_Scour/REEF3D-%s-Scour-Gauge-%03i.dat",solver,solver,n+1);

    gout[n] = new ofstream(gname);
    ofstream &o = *gout[n];

    o<<solver<<" scour gauge "<<n+1<<" (P 198)"<<endl;
    o<<"x_coord y_coord: "<<p->Xout(p->P198_x[n],p->P198_y[n])<<" "<<p->Yout(p->P198_x[n],p->P198_y[n])<<endl;
    o<<"D = "<<p->P188_D<<" m, statistics from t = "<<p->P188_ts<<" s";
    if(p->P199_Tw>0.0)
    o<<", time window "<<p->P199_Tw<<" s (exponential)";
    o<<endl;
    o<<"d50 = "<<p->S20<<" m, sediment density "<<p->S22<<", armour d = "<<p->P197_d<<" m, density "<<p->P197_rho
     <<", alpha "<<p->P197_aamp<<", kt "<<p->P197_kt<<", gamma_E "<<p->P197_gE<<", ks/D50 "<<p->P197_nr<<endl<<endl;

    o<<"time";
    for(int q=0; q<NQ; ++q)
    o<<"\t "<<name(q);
    o<<endl;
    }
}

void scour_core::gauges(lexer *p, ghostcell *pgc)
{
    // the bed output accumulates armour_d_req_max at these times too (P 187 with P 199 dt > 0)
    const bool acc = p->P187==1 && p->P199_dt>0.0;

    if(ngauge<=0 && npg<=0 && !acc)
    return;

    if(p->P199_dt>0.0)
    {
        if(p->simtime<nextgauge)
        return;

        nextgauge += p->P199_dt*MAX(1.0,floor((p->simtime-nextgauge)/p->P199_dt)+1.0);
    }

    if(acc)
    accumulate(p);

    pgauges(p,pgc);

    if(ngauge<=0)
    return;

    const double T = period(p);

    for(size_t q=0; q<gval.size(); ++q)
    gval[q] = -1.0e20;

    for(int n=0; n<ngauge; ++n)
    if(gflag[n]==1)
    cell(p,gi[n],gj[n],T,&gval[size_t(n)*NQ]);

    pgc->globalmax(gval.data(),int(gval.size()));

    // the maximum required D50 over the gauge outputs
    for(int n=0; n<ngauge; ++n)
    {
    if(wsum>=NACC*T)
    gmax[n] = MAX(gmax[n],gval[size_t(n)*NQ+13]);
    gval[size_t(n)*NQ+14] = gmax[n];
    }

    if(p->mpirank==0)
    for(int n=0; n<ngauge; ++n)
    {
        ofstream &o = *gout[n];

        o<<p->simtime;
        for(int q=0; q<NQ; ++q)
        o<<"\t "<<gval[size_t(n)*NQ+q];
        o<<"\n";

        if(++gflush>=100*ngauge)
        {
        for(ofstream *f : gout)
        f->flush();
        gflush=0;
        }
    }
}

// the maximum required D50 of the bed output, between its outputs at the times of P 199 dt
void scour_core::accumulate(lexer *p)
{
    const double T = period(p);

    if(wsum<NACC*T)
    return;

    double q[NQ];

    SLICELOOP4
    {
    cell(p,i,j,T,q);
    armdmax(i,j) = q[14];
    }
}

// P 193: the scour protection gauges, the cell that contains the gauge, at the times of P 199 dt
void scour_core::pgauge_ini(lexer *p)
{
    npg = p->P193;

    if(npg<=0)
    return;

    pi_.resize(npg);
    pj_.resize(npg);
    pflag.resize(npg);
    pval.resize(size_t(npg)*NP);
    pmax.assign(npg,0.0);

    for(int n=0; n<npg; ++n)
    {
    pi_[n] = p->posc_i(p->P193_x[n]);
    pj_[n] = p->j_dir==1 ? p->posc_j(p->P193_y[n]) : 0;

    pflag[n] = (pi_[n]>=0 && pi_[n]<p->knox && pj_[n]>=0 && pj_[n]<p->knoy) ? 1 : 0;
    }

    if(p->mpirank!=0)
    return;

    char gdir[128];
    snprintf(gdir,sizeof(gdir),"./REEF3D_%s_ScourProtection",solver);
    mkdir(gdir,0777);

    pout.resize(npg);

    for(int n=0; n<npg; ++n)
    {
    char gname[256];
    snprintf(gname,sizeof(gname),"./REEF3D_%s_ScourProtection/REEF3D-%s-ScourProtection-Gauge-%03i.dat",solver,solver,n+1);

    pout[n] = new ofstream(gname);
    ofstream &o = *pout[n];

    o<<solver<<" scour protection gauge "<<n+1<<" (P 193), Soulsby (1997) as REEF_Tools"<<endl;
    o<<"x_coord y_coord: "<<p->Xout(p->P193_x[n],p->P193_y[n])<<" "<<p->Yout(p->P193_x[n],p->P193_y[n])<<endl;
    o<<"stone d = "<<p->P197_d<<" m, density "<<p->P197_rho<<", water density "<<rho
     <<", alpha "<<p->P197_aamp<<", kt "<<p->P197_kt<<", gamma_E "<<p->P197_gE<<", ks/D50 "<<p->P197_nr
     <<", Shields 0.055, statistics from t = "<<p->P188_ts<<" s";
    if(p->P199_Tw>0.0)
    o<<", time window "<<p->P199_Tw<<" s (exponential)";
    o<<endl;
    o<<"Um, Uc: without the amplification alpha kt; tau_w, tau_c, tau_m, tau_max (with gamma_E), tau_cr and SR "
     <<"for the stone d; D50_req: gamma_E tau_max = tau_cr; D50_req_max: its maximum over time"<<endl<<endl;

    o<<"time	 Um	 Uc	 tau_w	 tau_c	 tau_m	 tau_max	 tau_cr	 SR	 D50_req	 D50_req_max"<<endl;
    }
}

void scour_core::pgauges(lexer *p, ghostcell *pgc)
{
    if(npg<=0)
    return;

    const double T = period(p);
    const double da = MAX(p->P197_d,1.0e-6);
    double q[NQ];

    for(size_t m=0; m<pval.size(); ++m)
    pval[m] = -1.0e20;

    for(int n=0; n<npg; ++n)
    if(pflag[n]==1)
    {
    cell(p,pi_[n],pj_[n],T,q);

    double *v = &pval[size_t(n)*NP];

    const double h = MAX((*depth)(pi_[n],pj_[n]),1.0e-3);
    const protection r = protect(q[2],q[5],T,h,da);

    v[0] = q[2];
    v[1] = q[5];
    v[2] = r.tw;
    v[3] = r.tc;
    v[4] = r.tm;
    v[5] = r.tmax;
    v[6] = r.tcr;
    v[7] = r.sr;
    v[8] = q[13];
    v[9] = 0.0;
    }

    pgc->globalmax(pval.data(),int(pval.size()));

    for(int n=0; n<npg; ++n)
    {
    if(wsum>=NACC*T)
    pmax[n] = MAX(pmax[n],pval[size_t(n)*NP+8]);
    pval[size_t(n)*NP+9] = pmax[n];
    }

    if(p->mpirank==0)
    for(int n=0; n<npg; ++n)
    {
        ofstream &o = *pout[n];

        o<<p->simtime;
        for(int m=0; m<NP; ++m)
        o<<"\t "<<pval[size_t(n)*NP+m];
        o<<"\n";
    }
}

scour_core::~scour_core()
{
    for(ofstream *f : gout)
    {
    f->flush();
    delete f;
    }

    for(ofstream *f : pout)
    {
    f->flush();
    delete f;
    }
}

void scour_core::gauge_flush()
{
    for(ofstream *f : gout)
    f->flush();

    for(ofstream *f : pout)
    f->flush();
}
