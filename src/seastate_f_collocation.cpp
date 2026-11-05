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


#include"seastate_f.h"
#include"seastate_collocation.h"
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_param.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<iostream>
#include<sstream>
#include<string>
#include<vector>
#include<mpi.h>
#include<sys/stat.h>
#include<sys/types.h>

/*--------------------------------------------------------------------
Online collocation (A 787 1): the observations of seastate-collocation.csv
(seastate_collocation.h) sampled during the run on the rank that owns
the cell of each observation.

Every call at the end of a time step (model time tn):
  - an observation at t <= tn is finished: linear in time between the
    sample taken at the end of the previous step (t_s < t) and the
    model at tn; without a sample (the first step) the model at tn;
  - an observation at tn < t <= tn + 1.5 dt is sampled at tn for the
    next call.
The finished matchups are gathered on rank 0, written to
REEF3D_SEASTATE_Collocation/seastate-collocation.csv and counted in the
statistics (from model time A 788 on, in a stationary run all; wet stencils only); the summary
in the log every A 789 s of model time and at the end, the table in
REEF3D_SEASTATE_Collocation/seastate-collocation-stats.dat (one numeric row per source and
one (-1) for all).

The model at a point: Hs, Tm01, Tp and the mean direction of the four
cell centres around it (level 0; the spectra of the ghost cells of the
neighbour ranks included), bilinear (the direction as a unit vector);
wet only if all four cells with weight are wet.
A stationary run (A 700 2) takes the final field for all observations.
A restart (I 45 1) appends to the file of the run it continues and
counts its matchups again; observations before the time of the restart
are not repeated.
--------------------------------------------------------------------*/

void seastate_f::colloc_ini(lexer *p, ghostcell *pgc)
{
    pcol = new seastate_collocation;
    std::string err;

    if(!pcol->read("seastate-collocation.csv",tref,err))
    {
        if(p->mpirank==0)
        cout<<endl<<"SEASTATE input error  --  A 787 1: "<<err<<endl<<endl;
    pgc->final(true);
    }

    const size_t n = pcol->obs.size();
    col_ci.assign(n,-1);
    col_cj.assign(n,-1);
    col_ts.assign(n,-1.0e300);
    col_s.assign(4*n,0.0);
    col_sw.assign(n,0);
    col_done.assign(n,0);
    col_next = 0;

    // the cell of every observation on this rank (the faces XN, YN of the local cells, ascending)
    long owned = 0;
    for(size_t k=0; k<n; ++k)
    {
    const double xp = pcol->obs[k].x, yp = pcol->obs[k].y;
    const double *xn = &p->XN[marge], *yn = &p->YN[marge];
        if(!(xp>=xn[0] && xp<xn[p->knox]))
        continue;
        if(p->j_dir==1 && !(yp>=yn[0] && yp<yn[p->knoy]))
        continue;
        if(p->j_dir==0 && p->knoy!=1)
        continue;
    const int ci = int(std::upper_bound(xn,xn+p->knox+1,xp) - xn) - 1;
    const int cj = (p->j_dir==0) ? 0 : int(std::upper_bound(yn,yn+p->knoy+1,yp) - yn) - 1;
        if(ci<0 || ci>=p->knox || cj<0 || cj>=p->knoy)
        continue;
    col_ci[k] = ci;
    col_cj[k] = cj;
    ++owned;
    }
    const long owned_all = long(pgc->globalisum(int(owned)));

    // observations before the start of this run (the start, or the time of the restart) are not taken (nonstationary)
    const double t0 = p->simtime;
    long before = 0;
    for(size_t k=0; k<n; ++k)
        if(p->A700!=2 && pcol->obs[k].t < t0 - 1.0e-6)
        {
        col_done[k] = 1;
        ++before;
        }

    // a restart appends to the file of the run it continues and counts its matchups again
    const bool append = restart && std::ifstream("./REEF3D_SEASTATE_Collocation/seastate-collocation.csv").good();

    if(p->mpirank==0)
    {
    mkdir("./REEF3D_SEASTATE_Collocation",0777);

        if(append)
        {
        std::ifstream in("./REEF3D_SEASTATE_Collocation/seastate-collocation.csv");
        std::string line;
        std::getline(in,line);
            while(std::getline(in,line))
            {
            std::vector<std::string> c;
            std::string v;
            std::istringstream s(line);
            while(std::getline(s,v,','))
            c.push_back(v);
                if(c.size()<13 || std::strtod(c[1].c_str(),nullptr)>=t0-1.0e-6)
                continue;
                if(c[12]=="1")
                {
                const int src = int(std::find(pcol->sources.begin(),pcol->sources.end(),c[4]) - pcol->sources.begin());
                pcol->count(src<int(pcol->sources.size()) ? src : -1,std::strtod(c[6].c_str(),nullptr),std::strtod(c[7].c_str(),nullptr));
                }
            }
        }
        else
        {
        ofstream out("./REEF3D_SEASTATE_Collocation/seastate-collocation.csv");
        out<<"time,t,x,y,source,id,hs_obs,hs,tm01,tp,dir,wet,stat"<<endl;
        }

    cout<<"SEASTATE collocation (A 787 1): "<<n<<" observations of seastate-collocation.csv ("<<pcol->sources.size()<<" sources), "
        <<owned_all<<" inside the grid";
        if(before>0)
        cout<<", "<<before<<" before the start";
    cout<<"; statistics from t = "<<p->A788<<" s (A 788), summary every "<<p->A789<<" s (A 789)"
        <<" -> REEF3D_SEASTATE_Collocation/seastate-collocation.csv"<<endl;
    }

    col_tlog = (p->A789>0.0) ? (std::floor(t0/p->A789 + 1.0e-9) + 1.0)*p->A789 : 1.0e300;

    // the samples at the start for the observations of the first step
    if(p->A700!=2)
    colloc(p,pgc,false);
}

// the model at observation k on this rank: Hs, Tm01, Tp, direction [deg]; false if a cell of the stencil is dry
bool seastate_f::colloc_sample(lexer *p, int k, double *v)
{
    const seastate_grid &g = *e->grid;
    const double xp = pcol->obs[k].x, yp = pcol->obs[k].y;
    const int ci = col_ci[k], cj = col_cj[k];

    // the four cell centres around the point (ghost cells of the neighbour ranks included)
    const int i0 = (xp>=p->XP[ci+marge]) ? ci : ci-1;
    const int j0 = (p->j_dir==0 || yp>=p->YP[cj+marge]) ? cj : cj-1;
    const double wx = std::min(std::max((xp - p->XP[i0+marge])/(p->XP[i0+1+marge] - p->XP[i0+marge]),0.0),1.0);
    const double wy = (p->j_dir==0) ? 0.0 : std::min(std::max((yp - p->YP[j0+marge])/(p->YP[j0+1+marge] - p->YP[j0+marge]),0.0),1.0);

    double hs = 0.0, tm = 0.0, tp = 0.0, dx = 0.0, dy = 0.0;
    bool wet = true;
    seastate_param sp;

    for(int a=0; a<2; ++a)
    for(int c=0; c<2; ++c)
    {
    const double w = (a==0 ? 1.0-wx : wx)*(c==0 ? 1.0-wy : wy);
        if(w<=0.0)
        continue;
    const int ii = i0+a, jj = j0+c;
        if(e->wet(ii,jj)!=1 || e->N->spec(ii,jj)==nullptr)
        {
        wet = false;
        continue;
        }
    sp.compute(g,e->N->spec(ii,jj));
    hs += w*sp.Hs;
    tm += w*sp.Tm01;
    tp += w*sp.Tp;
    dx += w*sp.Hs*std::cos(sp.dir*3.14159265358979323846/180.0);
    dy += w*sp.Hs*std::sin(sp.dir*3.14159265358979323846/180.0);
    }

    v[0] = hs;
    v[1] = tm;
    v[2] = tp;
    double d = std::atan2(dy,dx)*180.0/3.14159265358979323846;
    v[3] = (d<0.0) ? d+360.0 : d;
    return wet;
}

// the observations due at model time tn (fin: all that are left, the stationary end of the run)
void seastate_f::colloc(lexer *p, ghostcell *pgc, bool fin)
{
    const double tn = p->simtime;
    const double look = 1.5*std::max(dtw,p->dt);
    const size_t n = pcol->obs.size();
    std::vector<double> send;

    for(size_t k=col_next; k<n; ++k)
    {
    const seastate_obs &o = pcol->obs[k];
        if(!fin && o.t > tn + look)
        break;
        if(col_done[k])
        continue;

    const bool own = (col_ci[k]>=0);

        if(fin || o.t <= tn + 1.0e-6)
        {
            if(own)
            {
            double v[4];
            const bool wet = colloc_sample(p,int(k),v);
            const double ts = col_ts[k];

                // linear in time from the sample at the end of the previous step (both wet)
                if(!fin && ts>-1.0e299 && ts<o.t && tn>ts && wet && col_sw[k])
                {
                const double w = (o.t - ts)/(tn - ts);
                const double *s = &col_s[4*k];
                for(int q=0; q<3; ++q)
                v[q] = s[q] + w*(v[q] - s[q]);
                double dd = v[3] - s[3];
                dd -= 360.0*std::floor((dd + 180.0)/360.0);
                v[3] = std::fmod(s[3] + w*dd + 360.0,360.0);
                }

            send.push_back(double(k));
            for(int q=0; q<4; ++q)
            send.push_back(v[q]);
            send.push_back(wet ? 1.0 : 0.0);
            }
        col_done[k] = 1;
        }
        else if(own)
        {
        double v[4];
        col_sw[k] = colloc_sample(p,int(k),v) ? 1 : 0;
        col_ts[k] = tn;
        for(int q=0; q<4; ++q)
        col_s[4*k+q] = v[q];
        }
    }

    while(col_next<n && col_done[col_next])
    ++col_next;

    // the finished matchups on rank 0
    int cnt = int(send.size());
    std::vector<int> cnts(p->mpi_size,0), disp(p->mpi_size,0);
    MPI_Gather(&cnt,1,MPI_INT,cnts.data(),1,MPI_INT,0,pgc->mpi_comm);

    std::vector<double> all;
    if(p->mpirank==0)
    {
    for(int r=1; r<p->mpi_size; ++r)
    disp[r] = disp[r-1] + cnts[r-1];
    all.resize(size_t(disp[p->mpi_size-1] + cnts[p->mpi_size-1]));
    }
    MPI_Gatherv(send.data(),cnt,MPI_DOUBLE,all.data(),cnts.data(),disp.data(),MPI_DOUBLE,0,pgc->mpi_comm);

    if(p->mpirank!=0)
    return;

    const size_t m = all.size()/6;
    std::vector<size_t> ord(m);
    for(size_t q=0; q<m; ++q)
    ord[q] = q;
    std::sort(ord.begin(),ord.end(),[&](size_t a, size_t b){return all[6*a]<all[6*b];});

    if(m>0)
    {
    ofstream out("./REEF3D_SEASTATE_Collocation/seastate-collocation.csv",ios_base::app);
    out<<std::fixed;
        for(size_t q : ord)
        {
        const double *r = &all[6*q];
        const seastate_obs &o = pcol->obs[size_t(r[0])];
        const bool wet = (r[5]>0.5);
        const bool stat = wet && (p->A700==2 || o.t>=p->A788 - 1.0e-6);
            if(stat)
            pcol->count(o.src,o.hs,r[1]);
        out<<seastate_collocation::stamp(tref + o.t)<<","<<setprecision(1)<<o.t<<","<<o.x<<","<<o.y<<","<<pcol->sources[o.src]<<","
           <<pcol->ids[o.id]<<","<<setprecision(3)<<o.hs<<","<<r[1]<<","<<r[2]<<","<<r[3]<<","<<setprecision(1)<<r[4]<<","
           <<(wet ? 1 : 0)<<","<<(stat ? 1 : 0)<<"\n";
        }
    }

    if(fin || tn>=col_tlog - 1.0e-6)
    {
        while(col_tlog<=tn + 1.0e-6)
        col_tlog += (p->A789>0.0 ? p->A789 : 1.0e300);
    cout<<"SEASTATE collocation at t = "<<tn<<" s ("<<seastate_collocation::stamp(tref + tn)<<"): Hs, model - observation [m]"<<endl
        <<pcol->summary();
    }
}

// the end of the run: a stationary run takes the final field for all observations; the table of the statistics
void seastate_f::colloc_end(lexer *p, ghostcell *pgc)
{
    if(pcol==nullptr)
    return;

    colloc(p,pgc,p->A700==2);

    long left = 0;
    for(size_t k=0; k<col_done.size(); ++k)
    if(!col_done[k] && col_ci[k]>=0)
    ++left;
    left = long(pgc->globalisum(int(left)));

    if(p->mpirank!=0)
    return;

    ofstream out("./REEF3D_SEASTATE_Collocation/seastate-collocation-stats.dat");
    out<<"# REEF3D::SEASTATE online collocation (A 787 1): Hs, model - observation [m], matchups from t = "<<p->A788<<" s, wet"<<endl;
    out<<"# model time at the end "<<p->simtime<<" s ("<<seastate_collocation::stamp(tref + p->simtime)<<"); "<<left
       <<" observations inside the grid after the end of the run"<<endl;
    out<<"# source: N, mean observation, bias, RMSE, scatter index, correlation; -1 all sources"<<endl;
    for(size_t k=0; k<pcol->sources.size(); ++k)
    out<<"#   "<<k<<" "<<pcol->sources[k]<<endl;
    out<<setprecision(8);
    for(size_t k=0; k<=pcol->sources.size(); ++k)
    {
    const seastate_colloc_stats &a = (k<pcol->sources.size()) ? pcol->per_source[k] : pcol->total;
    out<<(k<pcol->sources.size() ? long(k) : -1L)<<" \t "<<a.n<<" \t "<<a.mean_obs()<<" \t "<<a.bias()<<" \t "<<a.rmse()<<" \t "<<a.si()<<" \t "<<a.r()<<endl;
    }

    cout<<"SEASTATE collocation, end of the run: Hs, model - observation [m]"<<endl<<pcol->summary();
}
