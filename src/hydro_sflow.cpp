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

hydro_sflow::hydro_sflow(lexer *p, fdm2D *b, ghostcell *pgc, const hydro_sflow *lv0) : S(p), R(p), F(p), Rc(p), fc(p), Pc(p), sc(p), L(p), C(p), Bq(p), Fp(p)
{
    in = (lv0==nullptr) ? this : lv0;

    SLICEBASELOOP
    {
    S(i,j) = 0.0;
    R(i,j) = 0.0;
    F(i,j) = 0.0;
    Rc(i,j) = 0.0;
    fc(i,j) = 0.0;
    Pc(i,j) = 0.0;
    sc(i,j) = 0.0;
    L(i,j) = 0.0;
    C(i,j) = 0.0;
    Bq(i,j) = 0.0;
    Fp(i,j) = 0.0;
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
    ifstream f(file);

    if(!f)
    return false;

    t.clear();
    v.clear();
    string line;

    while(getline(f,line))
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

// the rain rate [m/s] of a raster [mm/h] at a point, 0 outside and at NODATA
double hydro_sflow::rate(const hydro_raster &r, double x, double y)
{
    double v;

    if(r.at(x,y,v) && v>0.0)
    return v/3.6e6;

    return 0.0;
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

    if(p->R30==2 && p->R33==0 && (p->R32_Ks<=0.0 || p->R32_psi<0.0 || p->R32_dth<0.0))
    stop(p,pgc,"R 30 2 needs R 32 Ks psi dtheta with Ks > 0, psi >= 0, dtheta >= 0");

    if(p->R34<0.0 || p->R34>1.0)
    stop(p,pgc,"R 34: the sealed fraction is between 0 and 1");

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
        ifstream f("rain_series.dat");
        string line;
        ok = f ? 1 : 0;

        while(ok && getline(f,line))
        {
            const size_t c = line.find_first_not_of(" \t\r");

            if(c==string::npos || line[c]=='#' || line[c]=='/')
            continue;

            istringstream ls(line);
            double t;
            string fn;

            if(ls>>t>>fn)
            {
                if(!st.empty() && !(t>st.back()))
                ok = 0;
            st.push_back(t);
            sf.push_back(fn);
            }
        }

        if(st.empty())
        ok = 0;

        if(pgc->globalimin(ok)==0)
        stop(p,pgc,"R 10 3: no readable rain_series.dat (lines t [s] file.asc, times increasing)");
    }

    // point sources
    src.resize(p->R40);

    for(int n=0; n<p->R40; ++n)
    {
    src[n].ID = p->R40_ID[n];
    src[n].x = p->R40_x[n];
    src[n].y = p->R40_y[n];
    src[n].Q = p->R40_Q[n];
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

    // rasters of the roughness and the soil: read by every rank, kept for the patches and blocks
    auto raster = [&](const char *file, hydro_raster &r, bool needed, const string &key)
    {
        string err;
        ifstream test(file);
        int there = test ? 1 : 0;
        test.close();

        if(pgc->globalimin(there)==0)
        {
            if(needed)
            stop(p,pgc,key+": no "+file);
            return false;
        }

        const int rok = r.read(file,err) ? 1 : 0;

        if(pgc->globalimin(rok)==0)
        stop(p,pgc,key+": "+(err.empty() ? string(file) : err));

        if(p->mpirank==0)
        cout<<"HYDRO: "<<file<<" ("<<r.nc<<" x "<<r.nr<<" cells)"<<endl;

        return true;
    };

    if(p->R50==1)
    hman = raster("manning.asc",rman,true,"R 50 1");

    if(p->R33==1 && p->R30==1)
    hinf = raster("inf_rate.asc",rinf,true,"R 33 1 with R 30 1");

    if(p->R33==1 && p->R30==2)
    {
    hks = raster("ga_ks.asc",rks,false,"R 33 1");
    hpsi = raster("ga_psi.asc",rpsi,false,"R 33 1");
    hdth = raster("ga_dtheta.asc",rdth,false,"R 33 1");

        if(!hks && !hpsi && !hdth)
        stop(p,pgc,"R 33 1 with R 30 2: none of ga_ks.asc, ga_psi.asc, ga_dtheta.asc");
    }

    if(p->R35==1)
    hseal = raster("sealed.asc",rseal,true,"R 35 1");

    // the grid of level 0: rasters at its cells and the cells of the point sources, one rank per source
    grid_ini(p,b);

    for(int n=0; n<p->R40; ++n)
    {
        int here = 0;

        for(const cell &c : cells)
        if(c.n==n)
        here = 1;

        const int first = pgc->globalimin(here==1 ? p->mpirank : 1<<30);

        if(here==1 && p->mpirank!=first)
        cells.erase(remove_if(cells.begin(),cells.end(),[n](const cell &c){ return c.n==n; }),cells.end());

        if(first==(1<<30) && p->mpirank==0)
        cout<<"!!! HYDRO: point source R 40 "<<src[n].ID<<" at ("<<src[n].x<<", "<<src[n].y<<") is not in a fluid cell of level 0"
            <<" (a block of a multi-block grid may hold it)"<<endl;
    }

    if(hman)
    {
        pgc->gcsl_start4(p,b->ks,50);

        if(p->mpirank==0 && p->S10>0)
        cout<<"!!! HYDRO: with sediment (S 10) the roughness follows the sediment, not R 50"<<endl;
    }

    Vstore0 = Vstore = grid_storage(p,b,[](int,int){return true;});
    Vstore0 = Vstore = pgc->globalsum(Vstore0);
    tnext = p->simtime;

    // budget file
    if(p->mpirank==0)
    {
        mkdir("./REEF3D_SFLOW_Log",0777);

        if(reefckpt_open(bout,"./REEF3D_SFLOW_Log/REEF3D-SFLOW-Hydro-Budget.dat",p->I45))
        {
        bout<<"# REEF3D::HYDRO volume budget, cumulative volumes [m3] since the start; residual = dS - (rain + sources - abstraction - evaporation - infiltration); error = residual - (inflow - outflow)"<<endl;
        bout<<"# time [s] \t rain \t sources \t abstraction \t evaporation \t infiltration \t storage \t dS \t residual \t inflow \t outflow \t error \t rain_rate [m3/s] \t wet_cells"<<endl;
        }
    }
}

// a point in a cell of this grid (Cartesian: the nodes; curvilinear: lexer::geo_cell)
bool hydro_sflow::locate(lexer *p, double x, double y, int &ci, int &cj)
{
    if(in==this || p->geo_curv==1)
    {
        if(!p->geo_cell(x,y,ci,cj))
        return false;
    }

    else
    {
        bool found = false;

        for(int ii=0; ii<p->knox && !found; ++ii)
        if(x>=p->XN[ii+p->marge] && x<p->XN[ii+1+p->marge])
        for(int jj=0; jj<p->knoy && !found; ++jj)
        if(p->j_dir==0 || (y>=p->YN[jj+p->marge] && y<p->YN[jj+1+p->marge]))
        {
        ci = ii;
        cj = jj;
        found = true;
        }

        if(!found)
        return false;
    }

    return p->flagslice4[(ci-p->imin)*p->jmax + (cj-p->jmin)]>0;
}

void hydro_sflow::grid_ini(lexer *p, fdm2D *b)
{
    const hydro_sflow &I = *in;

    // the infiltration parameters and the sealed fraction
    SLICEBASELOOP
    {
    fc(i,j) = p->R30==1 ? MAX(p->R31,0.0)/3.6e6 : p->R32_Ks/3.6e6;
    Pc(i,j) = p->R32_psi*p->R32_dth;
    sc(i,j) = p->R34;
    }

    if(I.hman || I.hinf || I.hks || I.hpsi || I.hdth || I.hseal)
    SLICELOOP4
    {
        const double xw = p->Xout(xc(p,i,j),yc(p,i,j)), yw = p->Yout(xc(p,i,j),yc(p,i,j));
        double v;

        // Manning n: ks = (20 n)^6
        if(I.hman && I.rman.at(xw,yw,v) && v>0.0)
        b->ks(i,j) = pow(20.0*v,6.0);

        if(I.hinf && I.rinf.at(xw,yw,v) && v>=0.0)
        fc(i,j) = v/3.6e6;

        if(I.hks && I.rks.at(xw,yw,v) && v>=0.0)
        fc(i,j) = v/3.6e6;

        double psi = p->R32_psi, dth = p->R32_dth;

        if(I.hpsi && I.rpsi.at(xw,yw,v) && v>=0.0)
        psi = v;

        if(I.hdth && I.rdth.at(xw,yw,v) && v>=0.0)
        dth = v;

        if(I.hpsi || I.hdth)
        Pc(i,j) = psi*dth;

        if(I.hseal && I.rseal.at(xw,yw,v))
        sc(i,j) = MIN(MAX(v,0.0),1.0);
    }

    // the point sources in cells of this grid
    cells.clear();

    for(int n=0; n<int(I.src.size()); ++n)
    {
        cell c;
        c.n = n;

        if(locate(p,I.src[n].x,I.src[n].y,c.ci,c.cj))
        {
        c.area = p->cellarea(c.ci,c.cj);
        cells.push_back(c);
        }
    }
}

void hydro_sflow::load_record(lexer *p, ghostcell *pgc, int n, hydro_raster &dst)
{
    string err;
    const int ok = dst.read(sf[n],err) ? 1 : 0;

    if(pgc->globalimin(ok)==0)
    stop(p,pgc,"R 10 3: "+(err.empty() ? sf[n] : err));
}

void hydro_sflow::inputs_step(lexer *p, ghostcell *pgc)
{
    const double t = p->simtime;

    if(p->R10==1)
    rw = MAX(p->R11,0.0)/3.6e6;

    if(p->R10==2)
    rw = MAX(ipol(rt,rv,t,p->R12==1 ? 1 : 0),0.0)/3.6e6;

    if(p->R10==3)
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

        rlin = (p->R12==1 && n<N-1);

        if(rk!=n)
        {
            if(rlin && rk==n-1 && n>0)
            ra = rb;
            else
            load_record(p,pgc,n,ra);

            if(n<N-1 && p->R12==1)
            load_record(p,pgc,n+1,rb);

            rk = n;

            if(p->mpirank==0)
            cout<<"HYDRO: rain raster "<<sf[n]<<" (t = "<<st[n]<<" s)"<<endl;
        }

        rw = w;
    }

    // point sources (the same discharge on every rank), Q > 0
    qsrc = 0.0;

    for(source &s : src)
    {
        s.q = s.ht.empty() ? s.Q : ipol(s.ht,s.hq,t,1);

        if(s.q>0.0)
        qsrc += s.q;
    }
}

void hydro_sflow::grid_prepare(lexer *p)
{
    const hydro_sflow &I = *in;

    // rain [m/s]
    if(p->R10==1 || p->R10==2)
    SLICELOOP4
    R(i,j) = I.rw;

    if(p->R10==3)
    SLICELOOP4
    {
        const double xw = p->Xout(xc(p,i,j),yc(p,i,j)), yw = p->Yout(xc(p,i,j),yc(p,i,j));
        R(i,j) = I.rlin ? (1.0-I.rw)*rate(I.ra,xw,yw) + I.rw*rate(I.rb,xw,yw) : rate(I.ra,xw,yw);
    }

    SLICELOOP4
    S(i,j) = R(i,j);

    // the catchment model (R 100 2)
    if(cpl)
    SLICELOOP4
    S(i,j) += C(i,j);

    // point sources, Q > 0
    for(const cell &c : cells)
    if(I.src[c.n].q>0.0)
    S(c.ci,c.cj) += I.src[c.n].q/c.area;
}

double hydro_sflow::grid_rain(lexer *p, const function<bool(int,int)> &leaf)
{
    double q = 0.0;

    SLICELOOP4
    if(leaf(i,j))
    q += R(i,j)*p->cellarea(i,j);

    return q;
}

void hydro_sflow::mass_source(lexer *p, fdm2D *b, slice &K)
{
    SLICELOOP4
    K(i,j) += S(i,j);
}

void hydro_sflow::grid_sinks(lexer *p, fdm2D *b, double dt, const function<bool(int,int)> &leaf, double *v)
{
    const hydro_sflow &I = *in;
    const double wd = p->A244;
    const double pe = MAX(p->R20,0.0)/8.64e7*dt;

    // the cumulative rain of the cells
    SLICELOOP4
    {
    Rc(i,j) += R(i,j)*dt;
    L(i,j) = 0.0;
    }

    bool abstraction = false;

    for(const cell &c : cells)
    if(I.src[c.n].q<0.0)
    abstraction = true;

    if(!(pe>0.0 || p->R30>0 || abstraction))
    return;

    SLICELOOP4
    {
        const double avail = MAX(b->WL(i,j) - wd, 0.0);

        if(avail<=0.0)
        continue;

        const double A = p->cellarea(i,j);

        // the infiltration of the soil (the pervious part of the cell, F its cumulative infiltration) and of the cell
        double e = pe;
        const double fs = p->R30==1 ? fc(i,j)*dt : (p->R30==2 ? ga_step(F(i,j),Pc(i,j),fc(i,j)*dt) : 0.0);
        const double f = (1.0 - sc(i,j))*fs;
        double a = 0.0;

        if(abstraction)
        for(const cell &c : cells)
        if(I.src[c.n].q<0.0 && c.ci==i && c.cj==j)
        a += -I.src[c.n].q*dt/A;

        const double tot = e + f + a;

        if(tot<=0.0)
        continue;

        const double take = MIN(tot,avail);
        const double s = take/tot;

        b->WL(i,j) -= take;
        b->eta(i,j) = b->WL(i,j) - b->depth(i,j);
        L(i,j) = take;

        F(i,j) += fs*s;

        if(leaf(i,j))
        {
        v[0] += e*s*A;
        v[1] += f*s*A;
        v[2] += a*s*A;
        }
    }
}

double hydro_sflow::grid_storage(lexer *p, fdm2D *b, const function<bool(int,int)> &leaf)
{
    double v = 0.0;

    SLICELOOP4
    if(leaf(i,j))
    v += b->WL(i,j)*p->cellarea(i,j);

    return v;
}

void hydro_sflow::step_volumes(lexer *p, ghostcell *pgc, double ve, double vi, double va, double store)
{
    const double dt = p->dt;

    Vevap += pgc->globalsum(ve);
    Vinf += pgc->globalsum(vi);
    Vabs += pgc->globalsum(va);

    Vrain += dt*qrain;
    Vsrc += dt*qsrc;
    Vcpl += dt*qcpl;
    Vqi += dt*(bq_set ? bq_in : p->Qi);
    Vqo += dt*(bq_set ? bq_out : p->Qo);
    bq_set = false;

    // the storage of the water column the next step starts from (before sflow_eta::depth_update: the
    // semi-implicit scheme takes back the film it gives dry cells there, sflow_momentum_SI::si_restore)
    Vstore = pgc->globalsum(store);
}

void hydro_sflow::prepare(lexer *p, fdm2D *b, ghostcell *pgc)
{
    inputs_step(p,pgc);
    qsrc += qcpl;
    grid_prepare(p);

    qrain = pgc->globalsum(grid_rain(p,[](int,int){return true;}));
}

void hydro_sflow::sinks(lexer *p, fdm2D *b, ghostcell *pgc)
{
    auto all = [](int,int){return true;};
    double v[3] = {0.0,0.0,0.0};

    grid_sinks(p,b,p->dt,all,v);
    step_volumes(p,pgc,v[0],v[1],v[2],grid_storage(p,b,all));
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
    const double err = res - (Vqi - Vqo);

    if(p->mpirank==0 && due)
    {
    bout<<setprecision(10)<<p->simtime<<" \t "<<Vrain<<" \t "<<Vsrc<<" \t "<<Vabs<<" \t "<<Vevap<<" \t "<<Vinf
        <<" \t "<<Vstore<<" \t "<<dS<<" \t "<<res<<" \t "<<Vqi<<" \t "<<Vqo<<" \t "<<err<<" \t "<<qrain<<" \t "<<wet<<endl;
    }

    if(due && p->R90>0.0)
    while(tnext<=p->simtime+1.0e-9*MAX(1.0,fabs(tnext)))
    tnext += p->R90;

    if(p->mpirank==0 && p->count%p->P12==0)
    cout<<"hydro: rain "<<setprecision(5)<<qrain<<" m3/s, storage "<<Vstore<<" m3, residual "<<res<<" m3, error "<<err<<" m3, wet cells "<<wet<<endl;
}

void hydro_sflow::offset_vtp(lexer *p, int *offset, int &n)
{
    for(int q=0; q<3; ++q)
    {
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
    ++n;
    }
}

void hydro_sflow::name_vtp(std::ostream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"rain\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"rain_depth\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"infiltration\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
}

void hydro_sflow::name_pvtp(std::ostream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"rain\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"rain_depth\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"infiltration\"/>\n";
}

void hydro_sflow::print_vtp(lexer *p, std::ostream &result)
{
    int iin;
    float ffn;

    // rain [mm/h], cumulative rain and infiltration [mm]
    for(int q=0; q<3; ++q)
    {
        slice4 &f = (q==0) ? R : (q==1) ? Rc : F;
        const double fac = (q==0) ? 3.6e6 : 1.0e3;

        iin=sizeof(float)*p->pointnum2D;
        result.write((char*)&iin, sizeof(int));
        TPSLICELOOP
        {
        ffn=float(fac*p->sl_ipol4(f));
        result.write((char*)&ffn, sizeof(float));
        }
    }
}

void hydro_sflow::ckpt_fields(const function<void(const char*, char, double*)> &fd, const function<void(const char*, char, int*)> &fi)
{
    fd("hydro_F",'S',F.V);
    fd("hydro_Rc",'S',Rc.V);

    // coupled with the catchment model (R 100 2)
    if(cpl)
    {
    fd("hydro_C",'S',C.V);
    fd("hydro_Bq",'S',Bq.V);
    fd("hydro_Fp",'S',Fp.V);
    }
}

void hydro_sflow::ckpt_state(vector<double> &v) const
{
    v = {Vrain,Vsrc,Vabs,Vevap,Vinf,Vqi,Vqo,Vstore0,tnext,Vstore,Vcpl,qcpl};
}

bool hydro_sflow::ckpt_restore(const vector<double> &v)
{
    // 9 values: a checkpoint before R 100 2 (no storage of the step, no coupling)
    if(v.size()!=9 && v.size()!=12)
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
    store0_leaf = true;

    if(v.size()==12)
    {
    Vstore = v[9];
    Vcpl = v[10];
    qcpl = v[11];
    }

    rk = -2;
    return true;
}
