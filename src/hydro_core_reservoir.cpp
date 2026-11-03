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
#include"reefckpt.h"
#include<algorithm>
#include<cmath>
#include<iomanip>
#include<iostream>
#include<sstream>
#include<sys/stat.h>

/*
Reservoirs (R 110 1): hydro_reservoir.dat, one line per reservoir
    id x y area[m2] bottom[m] level[m] crest[m] weir_width[m] weir_coefficient[m^0.5/s] release
x, y: a point in the hydrology cell of the dam (DEM coordinates); release 1: the release Q(t) from
hydro_release_<id>.dat (lines t [s] Q [m3/s], linear, clamped), 0: none.
The reservoir is a level pool with vertical walls in place of the channel of its cell: storage area (h - bottom),
inflow = the discharge of the upstream cells and the runoff of the cell, outflow = the release (at most the water
above the bottom) + the spill over the weir cw width (h - crest)^1.5, implicit in h (Newton with bisection):
    area (h - h_old) = dt (inflow - release - cw width max(h - crest, 0)^1.5)
*/
void hydro_core::read_reservoirs(lexer *p, ghostcell *pgc)
{
    rescell.assign(N,-1);
    resv.clear();

    if(p->R110!=1)
    return;

    ifstream f("hydro_reservoir.dat");
    if(!f)
    stop(p,pgc,"R 110 1: no hydro_reservoir.dat");

    vector<int> map(size_t(ni)*nj,-1);
    for(int n=0; n<N; ++n)
    map[I[n]+size_t(ni)*J[n]]=n;

    string line;
    while(getline(f,line))
    {
        const size_t a=line.find_first_not_of(" \t\r");
        if(a==string::npos || line.compare(a,2,"//")==0)
        continue;

        istringstream s(line);
        reservoir r;
        double x, y;
        int rel=0;
        s>>r.id>>x>>y>>r.area>>r.zb>>r.h>>r.crest>>r.width>>r.cw>>rel;
        if(s.fail() || r.area<=0.0 || r.h<r.zb)
        stop(p,pgc,"hydro_reservoir.dat: could not read '"+line+"'");

        const int II=int(floor((x-x0)/D)), JJ=int(floor((y-y0)/D));
        r.cell = (II>=0 && JJ>=0 && II<ni && JJ<nj) ? map[II+size_t(ni)*JJ] : -1;
        if(r.cell<0)
        stop(p,pgc,"hydro_reservoir.dat: reservoir "+to_string(r.id)+" is not in a hydrology cell");
        if(rescell[r.cell]>=0)
        stop(p,pgc,"hydro_reservoir.dat: two reservoirs in cell "+to_string(r.cell));

        if(rel==1)
        {
            ifstream g("hydro_release_"+to_string(r.id)+".dat");
            double t1, q1;
            while(g>>t1>>q1)
            {
                r.rt.push_back(t1);
                r.rq.push_back(q1);
            }
            if(r.rt.empty())
            stop(p,pgc,"R 110 1: no readable hydro_release_"+to_string(r.id)+".dat");
        }

        r.qin=r.qrel=r.qspill=0.0;
        rescell[r.cell]=int(resv.size());
        resv.push_back(r);
    }

    if(p->mpirank==0)
    cout<<"HYDRO: "<<resv.size()<<" reservoirs"<<endl;
}

double hydro_core::reservoir_step(lexer *p, reservoir &r, double qin)
{
    // the release of the step (time mean of the linear series, by the midpoint), at most the water above the bottom
    double qr=0.0;
    if(!r.rt.empty())
    {
        const double tm=t+0.5*dt;
        if(tm<=r.rt.front())
        qr=r.rq.front();
        else if(tm>=r.rt.back())
        qr=r.rq.back();
        else
        {
            const size_t k=size_t(upper_bound(r.rt.begin(),r.rt.end(),tm)-r.rt.begin());
            const double w=(tm-r.rt[k-1])/(r.rt[k]-r.rt[k-1]);
            qr=(1.0-w)*r.rq[k-1]+w*r.rq[k];
        }
        qr=max(0.0,min(qr,(r.area*(r.h-r.zb)+dt*qin)/dt));
    }

    // the level: area (h - h0) - dt (qin - qr - spill(h)) = 0, increasing in h
    const double h0=r.h, cww=r.cw*r.width;
    auto spill=[&](double h){ return h>r.crest ? cww*pow(h-r.crest,1.5) : 0.0; };
    auto F=[&](double h){ return r.area*(h-h0)-dt*(qin-qr-spill(h)); };

    // F(bottom) <= 0 (the release is limited to the water above the bottom), F(level without spill) >= 0
    double lo=r.zb, hi=max(r.zb,h0+dt*(qin-qr)/r.area);
    double h=hi;
    if(F(hi)>0.0)
    for(int it=0; it<200; ++it)
    {
        const double f=F(h);
        if(f>0.0) hi=h; else lo=h;
        const double df=r.area+dt*(h>r.crest ? 1.5*cww*sqrt(h-r.crest) : 0.0);
        double hn=h-f/df;
        if(!(hn>lo && hn<hi))
        hn=0.5*(lo+hi);
        if(fabs(hn-h)<=1.0e-13*max(1.0,fabs(hn)) || hi-lo<=1.0e-13*max(1.0,fabs(hi)))
        {
            h=hn;
            break;
        }
        h=hn;
    }

    r.h=h;
    r.qin=qin;
    r.qrel=qr;
    r.qspill=spill(h);
    // the outflow from the volume balance (exact for the budget)
    return qin-r.area*(h-h0)/dt;
}

double hydro_core::reservoir_storage() const
{
    double s=0.0;
    for(const reservoir &r : resv)
    s+=r.area*(r.h-r.zb);
    return s;
}

void hydro_core::reservoir_output(lexer *p)
{
    if(resv.empty() || p->mpirank!=0)
    return;

    if(!fr.is_open())
    {
        mkdir("./REEF3D_HYDRO",0777);
        if(reefckpt_open(fr,"./REEF3D_HYDRO/REEF3D-HYDRO-Reservoirs.dat",p->R100==2 ? p->I45 : 0))
        {
        fr<<"// REEF3D::HYDRO reservoirs: per reservoir level [m], inflow, release, spill [m3/s] of the last step"<<endl<<"// time[s]";
        for(const reservoir &r : resv)
        fr<<" h_"<<r.id<<" qin_"<<r.id<<" qrel_"<<r.id<<" qspill_"<<r.id;
        fr<<endl;
        }
        fr<<setprecision(12);
    }

    fr<<t;
    for(const reservoir &r : resv)
    fr<<" "<<r.h<<" "<<r.qin<<" "<<r.qrel<<" "<<r.qspill;
    fr<<endl;
}
