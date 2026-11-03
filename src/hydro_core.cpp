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

#include"hydro_core.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<iostream>
#include<sstream>

using namespace std;

const char *hydro_core::parname[hydro_core::NPAR] =
{"TT","CFMAX","SFCF","PCORR","CFR","CWH","CFGLAC","FC","LP","BETA","PERC","UZL","K0","K1","K2","MAXBAS","PCALT","TCALT",
 "CH_N","CH_WA","CH_WB","INI_SWE","INI_SM","INI_SUZ","INI_SLZ"};

// defaults: typical HBV-light values
static const double pardefault[hydro_core::NPAR] =
{0.0,       // TT     threshold temperature [C]
 3.5,       // CFMAX  degree-day factor [mm/C/day]
 1.0,       // SFCF   snowfall correction [-]
 1.0,       // PCORR  precipitation correction [-]
 0.05,      // CFR    refreezing coefficient [-]
 0.1,       // CWH    water holding capacity of the snow [-]
 1.5,       // CFGLAC glacier melt factor relative to CFMAX [-]
 250.0,     // FC     field capacity [mm]
 0.7,       // LP     soil moisture above which evaporation is potential, fraction of FC [-]
 2.0,       // BETA   shape of the recharge [-]
 1.5,       // PERC   percolation to the lower zone [mm/day]
 20.0,      // UZL    threshold of the quick flow [mm]
 0.2,       // K0     quick flow recession [1/day]
 0.1,       // K1     upper zone recession [1/day]
 0.01,      // K2     lower zone recession [1/day]
 1.5,       // MAXBAS base of the triangular weights [day]
 10.0,      // PCALT  precipitation increase with elevation [%/100 m]
 0.6,       // TCALT  temperature decrease with elevation [C/100 m]
 0.035,     // CH_N   Manning n of the channel [s/m^(1/3)]
 2.5,       // CH_WA  channel width W = CH_WA (upstream area/km2)^CH_WB [m]
 0.5,       // CH_WB
 0.0,       // INI_SWE initial snow water equivalent [mm]
 0.5,       // INI_SM  initial soil moisture, fraction of FC [-]
 0.0,       // INI_SUZ initial upper zone [mm]
 0.0};      // INI_SLZ initial lower zone [mm]

bool hydro_core::active(lexer *p)
{
    return p->A10==2 && p->R100>0;
}

hydro_core::hydro_core(lexer *p, ghostcell *pgc) : t(0.0), dt(p->R101), N(0), nb(1), ret(-1), Qret(0.0), Vto(0.0), Vfrom(0.0),
    ni(0), nj(0), x0(0.0), y0(0.0), D(1.0), Vp(0.0), Vea(0.0), Vgm(0.0), Vout(0.0), S0(0.0), n0(0), n1(0), fhasE(false), nmb(1), tout(0.0), tvtk(0.0), tlast(-1.0), nvtk(0), kw(p->R106==1), first(true)
{
}

hydro_core::~hydro_core()
{
}

void hydro_core::stop(lexer *p, ghostcell *pgc, const string &msg) const
{
    if(p->mpirank==0)
    cout<<endl<<"!!! HYDRO: "<<msg<<" !!!"<<endl<<endl;

    pgc->final(true);
}

// the next line that is not empty and not a // comment
static bool next_line(ifstream &f, string &line)
{
    while(getline(f,line))
    {
        const size_t a=line.find_first_not_of(" \t\r");
        if(a==string::npos || line.compare(a,2,"//")==0)
        continue;
        return true;
    }
    return false;
}

void hydro_core::read_grid(lexer *p, ghostcell *pgc)
{
    ifstream f("hydro_grid.dat");
    if(!f)
    stop(p,pgc,"R 100: no hydro_grid.dat (DIVEMesh W 10)");

    string line, key;
    int version=0;
    next_line(f,line);
    {
        istringstream s(line);
        s>>key>>version;
        if(key!="REEF3D_HYDRO_GRID" || version<1 || version>2)
        stop(p,pgc,"hydro_grid.dat: unknown format");
    }

    auto header = [&](const char *name, istringstream &s)
    {
        if(!next_line(f,line))
        stop(p,pgc,string("hydro_grid.dat: no '")+name+"'");
        s.clear();
        s.str(line);
        s>>key;
        if(key!=name)
        stop(p,pgc,string("hydro_grid.dat: '")+name+"' expected");
    };

    istringstream s;
    header("origin",s);   s>>x0>>y0;
    header("cellsize",s); s>>D;
    header("raster",s);   s>>ni>>nj;
    header("bands",s);    s>>nb;
    header("cells",s);    s>>N;

    if(N<=0 || nb<=0)
    stop(p,pgc,"hydro_grid.dat: no cells");

    I.resize(N); J.resize(N); down.resize(N); sub.resize(N); zone.resize(N);
    xc.resize(N); yc.resize(N); area.resize(N); zmean.resize(N); len.resize(N); slope.resize(N); aup.resize(N); lakef.resize(N);
    bz.resize(N*nb); bf.resize(N*nb); bg.resize(N*nb);

    for(int n=0; n<N; ++n)
    {
        if(!next_line(f,line))
        stop(p,pgc,"hydro_grid.dat: too few cells");

        istringstream c(line);
        int id;
        c>>id>>I[n]>>J[n]>>xc[n]>>yc[n]>>area[n]>>zmean[n]>>down[n]>>len[n]>>slope[n]>>aup[n]>>sub[n]>>zone[n]>>lakef[n];
        for(int b=0; b<nb; ++b)
        c>>bz[n*nb+b]>>bf[n*nb+b]>>bg[n*nb+b];

        if(c.fail() || id!=n)
        stop(p,pgc,"hydro_grid.dat: could not read cell "+to_string(n));

        if(down[n]>=0 && down[n]<=n)
        stop(p,pgc,"hydro_grid.dat: the cells are not in topological order");
    }

    int ng=0;
    header("gauges",s); s>>ng;
    gid.resize(ng);
    gcell.resize(ng);
    for(int g=0; g<ng; ++g)
    {
        if(!next_line(f,line))
        stop(p,pgc,"hydro_grid.dat: too few gauges");
        istringstream c(line);
        double gx,gy;
        c>>gid[g]>>gx>>gy>>gcell[g];
        if(c.fail() || gcell[g]<0 || gcell[g]>=N)
        stop(p,pgc,"hydro_grid.dat: could not read gauge "+to_string(g));
    }

    // version 2: the fine flow paths (x y of the pixel centres)
    poff.assign(N+1,0);
    if(version>=2)
    {
        header("paths",s);
        for(int n=0; n<N; ++n)
        {
            if(!next_line(f,line))
            stop(p,pgc,"hydro_grid.dat: too few paths");
            istringstream c(line);
            int id, np;
            c>>id>>np;
            for(int k=0; k<np; ++k)
            {
                double a,b;
                c>>a>>b;
                px.push_back(a);
                py.push_back(b);
            }
            if(c.fail() || id!=n || np<1)
            stop(p,pgc,"hydro_grid.dat: could not read the path of cell "+to_string(n));
            poff[n+1]=int(px.size());
        }
    }

    // partition: contiguous ranges of cells
    const int P=p->mpi_size;
    cnt.assign(P,0);
    displ.assign(P,0);
    for(int r=0; r<P; ++r)
    {
        cnt[r] = N/P + (r<N%P ? 1 : 0);
        displ[r] = (r==0) ? 0 : displ[r-1]+cnt[r-1];
    }
    n0=displ[p->mpirank];
    n1=n0+cnt[p->mpirank];
}

// hydro_param.dat: lines 'NAME value' apply to the cells of the last selector line
// ('all', 'zone <id>', 'subcatchment <id>', 'cell <id>'), later lines override earlier ones
void hydro_core::read_param(lexer *p, ghostcell *pgc)
{
    par.resize(N*NPAR);
    for(int n=0; n<N; ++n)
    for(int k=0; k<NPAR; ++k)
    par[n*NPAR+k]=pardefault[k];

    ifstream f("hydro_param.dat");
    if(!f)
    {
        if(p->mpirank==0)
        cout<<"HYDRO: no hydro_param.dat, default parameters"<<endl;
    }

    string line;
    int sel=0, selid=0;     // 0 all, 1 zone, 2 subcatchment, 3 cell
    while(f && next_line(f,line))
    {
        istringstream s(line);
        string key;
        s>>key;

        if(key=="all")
        {sel=0; continue;}
        if(key=="zone")
        {sel=1; s>>selid; continue;}
        if(key=="subcatchment")
        {sel=2; s>>selid; continue;}
        if(key=="cell")
        {sel=3; s>>selid; continue;}

        int k=0;
        while(k<NPAR && key!=parname[k])
        ++k;

        double v;
        s>>v;
        if(k==NPAR || s.fail())
        stop(p,pgc,"hydro_param.dat: could not read '"+line+"'");

        for(int n=0; n<N; ++n)
        if(sel==0 || (sel==1 && zone[n]==selid) || (sel==2 && sub[n]==selid) || (sel==3 && n==selid))
        par[n*NPAR+k]=v;
    }

    for(int n=0; n<N; ++n)
    {
        if(P(n,FC)<=0.0 || P(n,LP)<=0.0 || P(n,LP)>1.0)
        stop(p,pgc,"hydro_param.dat: FC > 0 and 0 < LP <= 1");
        for(int k : {K0,K1,K2})
        if(P(n,k)<0.0 || P(n,k)>1.0)
        stop(p,pgc,"hydro_param.dat: 0 <= K0, K1, K2 <= 1 [1/day]");
        if(P(n,MAXBAS)<0.0 || P(n,CH_N)<=0.0 || P(n,CH_WA)<=0.0)
        stop(p,pgc,"hydro_param.dat: MAXBAS >= 0, CH_N > 0, CH_WA > 0");
    }
}

/*
hydro_forcing.dat
    stations <ns>
    <id> <x> <y> <z>                 ns lines
    pet <0|1>                        1: PET columns given
    <time[s]> then per station P[mm/h] T[C] (PET[mm/h])
a row holds until the time of the next row, the last row holds to the end
*/
void hydro_core::read_forcing(lexer *p, ghostcell *pgc)
{
    ifstream f("hydro_forcing.dat");
    if(!f)
    stop(p,pgc,"R 100: no hydro_forcing.dat");

    string line, key;
    int ns=0, pet=0;

    next_line(f,line);
    {
        istringstream s(line);
        s>>key>>ns;
        if(key!="stations" || ns<=0)
        stop(p,pgc,"hydro_forcing.dat: 'stations <n>' expected");
    }

    sx.resize(ns); sy.resize(ns); sz.resize(ns);
    for(int k=0; k<ns; ++k)
    {
        next_line(f,line);
        istringstream s(line);
        int id;
        s>>id>>sx[k]>>sy[k]>>sz[k];
        if(s.fail())
        stop(p,pgc,"hydro_forcing.dat: could not read station "+to_string(k+1));
    }

    next_line(f,line);
    {
        istringstream s(line);
        s>>key>>pet;
        if(key!="pet")
        stop(p,pgc,"hydro_forcing.dat: 'pet <0|1>' expected");
    }
    fhasE=(pet==1);

    while(next_line(f,line))
    {
        istringstream s(line);
        double tt;
        s>>tt;
        vector<double> P(ns), T(ns), E(ns,0.0);
        for(int k=0; k<ns; ++k)
        {
            s>>P[k]>>T[k];
            if(fhasE)
            s>>E[k];
        }
        if(s.fail())
        stop(p,pgc,"hydro_forcing.dat: could not read the row at t = "+to_string(tt));
        if(!ft.empty() && tt<=ft.back())
        stop(p,pgc,"hydro_forcing.dat: the times must increase");

        ft.push_back(tt);
        fP.push_back(P);
        fT.push_back(T);
        fE.push_back(E);
    }

    if(ft.empty())
    stop(p,pgc,"hydro_forcing.dat: no rows");

    // inverse distance weights (power 2) of the owned cells
    wgt.assign(size_t(n1-n0)*ns,0.0);
    for(int n=n0; n<n1; ++n)
    {
        double *w=&wgt[size_t(n-n0)*ns];
        int hit=-1;
        double sum=0.0;
        for(int k=0; k<ns; ++k)
        {
            const double d2=pow(xc[n]-sx[k],2)+pow(yc[n]-sy[k],2);
            if(d2<1.0e-6)
            hit=k;
            w[k] = 1.0/max(d2,1.0e-6);
            sum+=w[k];
        }
        for(int k=0; k<ns; ++k)
        w[k] = (hit>=0) ? (k==hit ? 1.0 : 0.0) : w[k]/sum;
    }

    mP.assign(ns,0.0);
    mT.assign(ns,0.0);
    mE.assign(ns,0.0);
}

// mean of the piecewise constant forcing over [ta,tb]
void hydro_core::forcing_mean(lexer *p, double ta, double tb)
{
    const int ns=int(sx.size());
    fill(mP.begin(),mP.end(),0.0);
    fill(mT.begin(),mT.end(),0.0);
    fill(mE.begin(),mE.end(),0.0);

    const int nr=int(ft.size());
    int r = int(upper_bound(ft.begin(),ft.end(),ta)-ft.begin())-1;
    r=max(r,0);

    double tc=ta;
    while(tc<tb)
    {
        const double te = (r+1<nr) ? min(tb,ft[r+1]) : tb;
        const double w = (te-tc)/(tb-ta);
        if(te>tc)
        for(int k=0; k<ns; ++k)
        {
            mP[k]+=w*fP[r][k];
            mT[k]+=w*fT[r][k];
            mE[k]+=w*fE[r][k];
        }
        tc=te;
        ++r;
        if(r>=nr)
        r=nr-1;
    }
}

// Oudin et al. (2005): PET = Re/(lambda rho) (T+5)/100 [mm/day] for T > -5 C, Re: extraterrestrial radiation (FAO-56)
double hydro_core::pet_oudin(lexer *p, double T, double tday) const
{
    if(T<=-5.0)
    return 0.0;

    const double PI=3.14159265358979323846;
    const double phi=p->R104*PI/180.0;
    const double J=fmod(p->R105-1.0+tday,365.0)+1.0;
    const double dr=1.0+0.033*cos(2.0*PI*J/365.0);
    const double decl=0.409*sin(2.0*PI*J/365.0-1.39);
    const double ws=acos(max(-1.0,min(1.0,-tan(phi)*tan(decl))));
    const double Re=24.0*60.0/PI*0.0820*dr*(ws*sin(phi)*sin(decl)+cos(phi)*cos(decl)*sin(ws));   // MJ/m2/day

    return max(0.0,Re)/2.45*(T+5.0)/100.0;   // mm/day
}

void hydro_core::ini(lexer *p, ghostcell *pgc)
{
    if(dt<=0.0)
    stop(p,pgc,"R 101: the time step must be > 0");

    read_grid(p,pgc);
    read_param(p,pgc);
    read_forcing(p,pgc);

    const int no=n1-n0;
    SP.assign(size_t(no)*nb,0.0);
    WC.assign(size_t(no)*nb,0.0);
    SM.assign(size_t(no)*nb,0.0);
    SUZ.assign(no,0.0);
    SLZ.assign(no,0.0);

    for(int n=n0; n<n1; ++n)
    {
        for(int b=0; b<nb; ++b)
        {
            SP[(n-n0)*nb+b]=P(n,INI_SWE);
            SM[(n-n0)*nb+b]=P(n,INI_SM)*P(n,FC);
        }
        SUZ[n-n0]=P(n,INI_SUZ);
        SLZ[n-n0]=P(n,INI_SLZ);
    }

    // MAXBAS: triangle of base B = MAXBAS/dt (steps) and area 1, weight of step k = integral over [k,k+1]
    const double dtd=dt/86400.0;
    nmb=1;
    for(int n=0; n<N; ++n)
    nmb=max(nmb,int(ceil(P(n,MAXBAS)/dtd-1.0e-9)));

    mbw.assign(size_t(N)*nmb,0.0);
    for(int n=0; n<N; ++n)
    {
        const double B=P(n,MAXBAS)/dtd;
        double *w=&mbw[size_t(n)*nmb];
        if(B<=1.0)
        {
            w[0]=1.0;
            continue;
        }

        auto F = [&](double x)      // cumulative area of the triangle up to x
        {
            x=max(0.0,min(B,x));
            return (x<=0.5*B) ? 2.0*x*x/(B*B) : 1.0-2.0*(B-x)*(B-x)/(B*B);
        };
        double s=0.0;
        for(int k=0; k<nmb; ++k)
        {
            w[k]=F(k+1.0)-F(double(k));
            s+=w[k];
        }
        for(int k=0; k<nmb; ++k)
        w[k]/=s;
    }
    buf.assign(size_t(no)*nmb,0.0);

    qgen.assign(N,0.0);
    share.assign(N,0.0);
    enter.assign(N,0);
    addv.assign(N,0.0);
    qsf.assign(N,0.0);
    Qent.assign(N,0.0);
    cP.assign(no,0.0); cT.assign(no,0.0); cE.assign(no,0.0); cEa.assign(no,0.0); cRun.assign(no,0.0);

    // channels: A = alpha Q^0.6 for a wide channel, alpha = (n W^(2/3) / sqrt(S))^0.6
    Q.assign(N,0.0);
    Qold.assign(N,0.0);
    Qin.assign(N,0.0);
    alpha.assign(N,0.0);
    for(int n=0; n<N; ++n)
    {
        const double W=P(n,CH_WA)*pow(max(aup[n],area[n])*1.0e-6,P(n,CH_WB));
        alpha[n]=pow(P(n,CH_N)*pow(W,2.0/3.0)/sqrt(max(slope[n],1.0e-6)),0.6);
        len[n]=max(len[n],1.0);
    }

    if(p->R108==1)
    read_state(p,pgc);

    if(p->mpirank==0)
    {
        double A=0.0;
        for(int n=0; n<N; ++n)
        A+=area[n];
        cout<<"HYDRO: "<<N<<" cells ("<<A*1.0e-6<<" km2), "<<nb<<" bands, "<<sx.size()<<" stations, "
            <<gid.size()<<" gauges, dt "<<dt<<" s, "<<p->mpi_size<<" ranks"<<endl;
    }
}
