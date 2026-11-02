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

#include"hydro_sflow.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"reefckpt.h"
#include<algorithm>
#include<cmath>
#include<iomanip>
#include<iostream>
#include<sstream>
#include<sys/stat.h>
#include<sys/types.h>

using namespace std;

hydro_sflow::hydro_sflow(lexer *p, fdm2D *b, ghostcell *pgc) : S(p), R(p), F(p), Ra(p), Rb(p)
{
    SLICEBASELOOP
    {
    S(i,j) = 0.0;
    R(i,j) = 0.0;
    F(i,j) = 0.0;
    Ra(i,j) = 0.0;
    Rb(i,j) = 0.0;
    }
}

void hydro_sflow::stop(lexer *p, ghostcell *pgc, const string &msg) const
{
    if(p->mpirank==0)
    cout<<endl<<"!!! HYDRO: "<<msg<<" !!!"<<endl<<endl;

    pgc->final(true);
}

double hydro_sflow::xc(lexer *p, int ii, int jj) const
{
    return p->geo_curv==1 ? p->XC2D[(ii-p->imin)*p->jmax + (jj-p->jmin)] : p->XP[ii+p->marge];
}

double hydro_sflow::yc(lexer *p, int ii, int jj) const
{
    return p->geo_curv==1 ? p->YC2D[(ii-p->imin)*p->jmax + (jj-p->jmin)] : p->YP[jj+p->marge];
}

bool hydro_sflow::read_series(const string &file, vector<double> &t, vector<double> &v)
{
    ifstream in(file);

    if(!in)
    return false;

    t.clear();
    v.clear();
    string line;

    while(getline(in,line))
    {
        const size_t c = line.find_first_not_of(" \t\r");

        if(c==string::npos || line[c]=='#' || line[c]=='/')
        continue;

        istringstream ls(line);
        double a,q;

        if(ls>>a>>q)
        {
        t.push_back(a);
        v.push_back(q);
        }
    }

    for(size_t n=1; n<t.size(); ++n)
    if(!(t[n]>t[n-1]))
    return false;

    return !t.empty();
}

double hydro_sflow::ipol(const vector<double> &t, const vector<double> &v, double tt, int mode)
{
    if(tt<=t.front())
    return v.front();

    if(tt>=t.back())
    return v.back();

    size_t n = size_t(upper_bound(t.begin(),t.end(),tt) - t.begin()) - 1;

    if(mode==0)
    return v[n];

    const double w = (tt-t[n])/(t[n+1]-t[n]);
    return (1.0-w)*v[n] + w*v[n+1];
}

// Green-Ampt: the potential infiltration of a step from the cumulative infiltration F0, the root F1 of
// g(F1) = F1 - F0 - P ln((F1 + P)/(F0 + P)) - Ks dt (ponded cell), g increasing and convex: Newton from an
// upper bound, safeguarded by bisection
double hydro_sflow::ga_step(double F0, double P, double Ksdt)
{
    if(Ksdt<=0.0)
    return 0.0;

    if(P<=0.0)
    return Ksdt;

    auto g = [&](double f1) { return f1 - F0 - P*log((f1+P)/(F0+P)) - Ksdt; };

    double lo = F0, hi = F0 + Ksdt + 2.0*sqrt(2.0*P*Ksdt) + Ksdt;

    while(g(hi)<0.0)
    hi = F0 + 2.0*(hi-F0);

    double f1 = hi;

    for(int it=0; it<60; ++it)
    {
        const double gv = g(f1);

        if(fabs(gv)<=1.0e-14*(Ksdt + f1))
        break;

        if(gv>0.0)
        hi = f1;
        else
        lo = f1;

        double fn = f1 - gv*(f1+P)/MAX(f1,1.0e-300);

        if(!(fn>lo && fn<hi))
        fn = 0.5*(lo+hi);

        f1 = fn;
    }

    return MAX(f1-F0,0.0);
}

void hydro_sflow::ini(lexer *p, fdm2D *b, ghostcell *pgc)
{
    int ok;

    if(p->mpirank==0)
    {
    cout<<"HYDRO (SFLOW): rain R 10 "<<p->R10<<", evaporation R 20 "<<p->R20<<", infiltration R 30 "<<p->R30
        <<", point sources R 40 "<<p->R40<<", roughness raster R 50 "<<p->R50<<endl;

        if(wetting(p) && p->A243==1)
        cout<<"!!! HYDRO: A 243 1 resets the rain on dry cells with the film, use A 243 2"<<endl;

        if(p->R50>0 && p->A218==0)
        cout<<"!!! HYDRO: R 50 roughness raster without bed friction (A 218 0)"<<endl;
    }

    if(p->R10<0 || p->R10>3)
    stop(p,pgc,"R 10 must be 0, 1, 2 or 3");

    if(p->R30<0 || p->R30>2)
    stop(p,pgc,"R 30 must be 0, 1 or 2");

    if(p->R30==2 && (p->R32_Ks<=0.0 || p->R32_psi<0.0 || p->R32_dth<0.0))
    stop(p,pgc,"R 30 2 needs R 32 Ks psi dtheta with Ks > 0, psi >= 0, dtheta >= 0");

    // rain time series
    if(p->R10==2)
    {
        ok = read_series("rain.dat",rt,rv) ? 1 : 0;
        if(pgc->globalimin(ok)==0)
        stop(p,pgc,"R 10 2: no readable rain.dat (lines t [s] rate [mm/h], times increasing)");
    }

    // rain rasters
    if(p->R10==3)
    {
        ifstream in("rain_series.dat");
        string line;
        ok = in ? 1 : 0;

        while(ok && getline(in,line))
        {
            const size_t c = line.find_first_not_of(" \t\r");

            if(c==string::npos || line[c]=='#' || line[c]=='/')
            continue;

            istringstream ls(line);
            double t;
            string f;

            if(ls>>t>>f)
            {
                if(!st.empty() && !(t>st.back()))
                ok = 0;
            st.push_back(t);
            sf.push_back(f);
            }
        }

        if(st.empty())
        ok = 0;

        if(pgc->globalimin(ok)==0)
        stop(p,pgc,"R 10 3: no readable rain_series.dat (lines t [s] file.asc, times increasing)");
    }

    // point sources: the cell that holds (x,y), one rank
    src.resize(p->R40);

    for(int n=0; n<p->R40; ++n)
    {
        source &s = src[n];
        s.ID = p->R40_ID[n];
        s.x = p->R40_x[n];
        s.y = p->R40_y[n];
        s.Q = p->R40_Q[n];

        int ii=0, jj=0;
        int found = p->geo_cell(s.x,s.y,ii,jj) ? 1 : 0;

        if(found==1 && p->flagslice4[(ii-p->imin)*p->jmax + (jj-p->jmin)]<0)
        found = 0;

        const int first = pgc->globalimin(found==1 ? p->mpirank : 1<<30);

        if(first==(1<<30))
        {
            ostringstream os;
            os<<"point source R 40 "<<s.ID<<" at ("<<s.x<<", "<<s.y<<") is not in a fluid cell of the grid";
            stop(p,pgc,os.str());
        }

        s.own = (p->mpirank==first);

        if(s.own)
        {
        s.ci = ii;
        s.cj = jj;
        s.area = p->cellarea(ii,jj);
        }
    }

    for(int n=0; n<p->R41; ++n)
    {
        int k = -1;

        for(int q=0; q<int(src.size()); ++q)
        if(src[q].ID==p->R41_ID[n])
        k = q;

        if(k<0)
        {
            ostringstream os;
            os<<"R 41 "<<p->R41_ID[n]<<": no point source R 40 with this ID";
            stop(p,pgc,os.str());
        }

        ostringstream name;
        name<<"hydrograph_src_"<<p->R41_ID[n]<<".dat";
        ok = read_series(name.str(),src[k].ht,src[k].hq) ? 1 : 0;

        if(pgc->globalimin(ok)==0)
        stop(p,pgc,"R 41: no readable "+name.str()+" (lines t [s] Q [m3/s], times increasing)");
    }

    // Manning n raster: ks = (20 n)^6, n = ks^(1/6)/20 as sflow_rough_manning
    if(p->R50==1)
    {
        hydro_raster mr;
        string err;
        ok = mr.read("manning.asc",err) ? 1 : 0;

        if(pgc->globalimin(ok)==0)
        stop(p,pgc,"R 50 1: "+(err.empty() ? string("manning.asc") : err));

        int inside=0, bad=0;

        SLICELOOP4
        {
            const double x = xc(p,i,j), y = yc(p,i,j);
            double n;

            if(mr.at(p->Xout(x,y),p->Yout(x,y),n))
            {
                if(n>0.0)
                {
                b->ks(i,j) = pow(20.0*n,6.0);
                ++inside;
                }
                else
                ++bad;
            }
        }

        pgc->gcsl_start4(p,b->ks,50);

        inside = pgc->globalisum(inside);
        bad = pgc->globalisum(bad);

        if(p->mpirank==0)
        cout<<"HYDRO: Manning raster manning.asc on "<<inside<<" of "<<p->cellnumtot2D<<" cells (the others keep B 50)"
            <<(bad>0 ? ", n <= 0 ignored on " : "")<<(bad>0 ? to_string(bad)+" cells" : "")<<endl;

        if(p->mpirank==0 && p->S10>0)
        cout<<"!!! HYDRO: with sediment (S 10) the roughness follows the sediment, not R 50"<<endl;
    }

    Vstore0 = Vstore = storage(p,b,pgc);
    tnext = p->simtime;

    // budget file
    if(p->mpirank==0)
    {
        mkdir("./REEF3D_SFLOW_Log",0777);

        if(reefckpt_open(bout,"./REEF3D_SFLOW_Log/REEF3D-SFLOW-Hydro-Budget.dat",p->I45))
        {
        bout<<"# REEF3D::HYDRO volume budget, cumulative volumes [m3] since the start; residual = dS - (rain + sources - abstraction - evaporation - infiltration)"<<endl;
        bout<<"# time [s] \t rain \t sources \t abstraction \t evaporation \t infiltration \t storage \t dS \t residual \t Qi_int \t Qo_int \t rain_rate [m3/s] \t wet_cells"<<endl;
        }
    }
}

void hydro_sflow::load_record(lexer *p, ghostcell *pgc, int n, slice4 &dst)
{
    hydro_raster r;
    string err;
    const int ok = r.read(sf[n],err) ? 1 : 0;

    if(pgc->globalimin(ok)==0)
    stop(p,pgc,"R 10 3: "+(err.empty() ? sf[n] : err));

    SLICEBASELOOP
    dst(i,j) = 0.0;

    SLICELOOP4
    {
        const double x = xc(p,i,j), y = yc(p,i,j);
        double v;

        if(r.at(p->Xout(x,y),p->Yout(x,y),v) && v>0.0)
        dst(i,j) = v/3.6e6;
    }
}

void hydro_sflow::rain_rasters(lexer *p, ghostcell *pgc, double t)
{
    const int N = int(st.size());
    int n = 0;
    double w = 0.0;

    if(t>=st.back())
    n = N-1;

    else if(t>st.front())
    {
    n = int(upper_bound(st.begin(),st.end(),t) - st.begin()) - 1;
    w = (t-st[n])/(st[n+1]-st[n]);
    }

    const bool lin = (p->R12==1 && n<N-1);

    if(rk!=n)
    {
        if(lin && rk==n-1 && n>0)
        {
            SLICEBASELOOP
            Ra(i,j) = Rb(i,j);
        }
        else
        load_record(p,pgc,n,Ra);

        if(n<N-1 && p->R12==1)
        load_record(p,pgc,n+1,Rb);

        rk = n;

        if(p->mpirank==0)
        cout<<"HYDRO: rain raster "<<sf[n]<<" (t = "<<st[n]<<" s)"<<endl;
    }

    SLICELOOP4
    R(i,j) = lin ? (1.0-w)*Ra(i,j) + w*Rb(i,j) : Ra(i,j);
}

void hydro_sflow::prepare(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const double t = p->simtime;

    // rain [m/s]
    if(p->R10==1)
    SLICELOOP4
    R(i,j) = MAX(p->R11,0.0)/3.6e6;

    if(p->R10==2)
    {
    const double r = MAX(ipol(rt,rv,t,p->R12==1 ? 1 : 0),0.0)/3.6e6;

        SLICELOOP4
        R(i,j) = r;
    }

    if(p->R10==3)
    rain_rasters(p,pgc,t);

    double q = 0.0;

    SLICELOOP4
    {
    S(i,j) = R(i,j);
    q += R(i,j)*p->cellarea(i,j);
    }

    qrain = pgc->globalsum(q);

    // point sources, Q > 0 (the same discharge on every rank)
    qsrc = 0.0;

    for(source &s : src)
    {
        s.q = s.ht.empty() ? s.Q : ipol(s.ht,s.hq,t,1);

        if(s.q>0.0)
        {
        qsrc += s.q;

            if(s.own)
            S(s.ci,s.cj) += s.q/s.area;
        }
    }
}

void hydro_sflow::mass_source(lexer *p, fdm2D *b, slice &K)
{
    SLICELOOP4
    K(i,j) += S(i,j);
}

void hydro_sflow::sinks(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const double dt = p->dt;
    const double wd = p->A244;
    const double pe = MAX(p->R20,0.0)/8.64e7*dt;
    const double fc = p->R30==1 ? MAX(p->R31,0.0)/3.6e6*dt : 0.0;
    const double Ksdt = p->R32_Ks/3.6e6*dt;
    const double P = p->R32_psi*p->R32_dth;

    double ve=0.0, vi=0.0, va=0.0;

    const bool losses = (pe>0.0 || p->R30>0);
    bool abstraction = false;

    for(const source &s : src)
    if(s.q<0.0)
    abstraction = true;

    if(losses || abstraction)
    SLICELOOP4
    {
        const double avail = MAX(b->WL(i,j) - wd, 0.0);

        if(avail<=0.0)
        continue;

        const double A = p->cellarea(i,j);

        double e = pe;
        double f = p->R30==1 ? fc : (p->R30==2 ? ga_step(F(i,j),P,Ksdt) : 0.0);
        double a = 0.0;

        if(abstraction)
        for(const source &s : src)
        if(s.own && s.q<0.0 && s.ci==i && s.cj==j)
        a += -s.q*dt/A;

        const double tot = e + f + a;

        if(tot<=0.0)
        continue;

        const double take = MIN(tot,avail);
        const double sc = take/tot;

        b->WL(i,j) -= take;
        b->eta(i,j) = b->WL(i,j) - b->depth(i,j);

        F(i,j) += f*sc;

        ve += e*sc*A;
        vi += f*sc*A;
        va += a*sc*A;
    }

    Vevap += pgc->globalsum(ve);
    Vinf += pgc->globalsum(vi);
    Vabs += pgc->globalsum(va);

    Vrain += dt*qrain;
    Vsrc += dt*qsrc;
    Vqi += dt*p->Qi;
    Vqo += dt*p->Qo;

    // the storage of the water column the next step starts from (before sflow_eta::depth_update: the
    // semi-implicit scheme takes back the film it gives dry cells there, sflow_momentum_SI::si_restore)
    Vstore = storage(p,b,pgc);
}

double hydro_sflow::storage(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double v = 0.0;

    SLICELOOP4
    v += b->WL(i,j)*p->cellarea(i,j);

    return pgc->globalsum(v);
}

void hydro_sflow::budget(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const bool due = (p->R90<=0.0 || p->simtime>=tnext-1.0e-9*MAX(1.0,fabs(tnext)));

    if(!due && !(p->count%p->P12==0))
    return;

    int wet = 0;

    SLICELOOP4
    if(p->wet[IJ]==1)
    ++wet;

    wet = pgc->globalisum(wet);

    const double dS = Vstore - Vstore0;
    const double res = dS - (Vrain + Vsrc - Vabs - Vevap - Vinf);

    if(p->mpirank==0 && due)
    {
    bout<<setprecision(10)<<p->simtime<<" \t "<<Vrain<<" \t "<<Vsrc<<" \t "<<Vabs<<" \t "<<Vevap<<" \t "<<Vinf
        <<" \t "<<Vstore<<" \t "<<dS<<" \t "<<res<<" \t "<<Vqi<<" \t "<<Vqo<<" \t "<<qrain<<" \t "<<wet<<endl;
    }

    if(due && p->R90>0.0)
    while(tnext<=p->simtime+1.0e-9*MAX(1.0,fabs(tnext)))
    tnext += p->R90;

    if(p->mpirank==0 && p->count%p->P12==0)
    cout<<"hydro: rain "<<setprecision(5)<<qrain<<" m3/s, storage "<<Vstore<<" m3, residual "<<res<<" m3, wet cells "<<wet<<endl;
}

void hydro_sflow::ckpt_fields(const function<void(const char*, char, double*)> &fd, const function<void(const char*, char, int*)> &fi)
{
    fd("hydro_F",'S',F.V);
}

void hydro_sflow::ckpt_state(vector<double> &v) const
{
    v = {Vrain,Vsrc,Vabs,Vevap,Vinf,Vqi,Vqo,Vstore0,tnext};
}

bool hydro_sflow::ckpt_restore(const vector<double> &v)
{
    if(v.size()!=9)
    return false;

    Vrain = v[0];
    Vsrc = v[1];
    Vabs = v[2];
    Vevap = v[3];
    Vinf = v[4];
    Vqi = v[5];
    Vqo = v[6];
    Vstore0 = v[7];
    tnext = v[8];

    rk = -2;
    return true;
}
