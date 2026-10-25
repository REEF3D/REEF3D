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

#include"met_provider.h"
#include"met_source.h"
#include"lexer.h"
#include"ghostcell.h"
#include<cmath>
#include<iostream>

met_provider *met_provider::instance = nullptr;

met_provider* met_provider::get(lexer *p, ghostcell *pgc)
{
    if(instance==nullptr)
    instance = new met_provider(p,pgc);

    return instance;
}

void met_provider::stop(lexer *p, ghostcell *pgc, const std::string &msg) const
{
    std::cout<<std::endl<<"REEF3D met forcing, input error (rank "<<p->mpirank<<")  --  "<<msg<<std::endl<<std::endl;

    pgc->final(true);
}

met_provider::met_provider(lexer *p, ghostcell *pgc)
{
    std::string err;

    if(p->W205==1)
    frame = met_frame(p->W205_lon,p->W205_lat,p->W205_rot,p->W206_x,p->W206_y);

    if(p->W200==1)
    {
    src = new met_source_uniform(p->W201_u,p->W201_dir,p->W202);

    tref = (p->W203>0.0) ? met_datetime(p->W203) : 0.0;

        if(p->mpirank==0)
        std::cout<<"REEF3D met forcing (W 200 1): uniform wind U10 "<<p->W201_u<<" m/s to "<<p->W201_dir
                 <<" deg, air pressure "<<p->W202<<" Pa"<<std::endl;
    }

    if(p->W200==2)
    {
    met_source_file *fs = new met_source_file;

        if(!fs->open("met-forcing.r3m",err))
        stop(p,pgc,"W 200 2: "+err);

    src = fs;

    const met_file &f = fs->f;

    tref = (p->W203>0.0) ? met_datetime(p->W203) : f.time(0);

        if(p->mpirank==0)
        {
        std::cout<<"REEF3D met forcing (W 200 2): met-forcing.r3m, "<<f.nx<<" x "<<f.ny<<" nodes from ("<<f.x0<<", "<<f.y0
                 <<") m, spacing "<<f.dx<<" x "<<f.dy<<" m, "<<f.nrec<<" records "<<met_date(f.time(0))<<" - "<<met_date(f.time(f.nrec-1))
                 <<(f.nvar==3 ? ", u10 v10 pa" : ", u10 v10 (pa 101325 Pa)")<<(f.moving==1 ? ", moving grid" : "")<<std::endl;

            if(f.avg>0.0)
            std::cout<<"  averaging period of u10: "<<f.avg<<" s"<<std::endl;

            if(tref + p->N41 < f.time(0) || tref > f.time(f.nrec-1))
            std::cout<<"  warning: the run lies outside the records, the forcing is held at the first or last record"<<std::endl;

            // frame of the file against W 205 / W 206
            if(f.frame==1 && p->W205==1)
            {
            const double d = std::fabs(f.lon0-p->W205_lon) + std::fabs(f.lat0-p->W205_lat) + std::fabs(f.rot-p->W205_rot);
            const double o = std::fabs(f.X0-p->W206_x) + std::fabs(f.Y0-p->W206_y);

                if(d>1.0e-9 || o>1.0e-6)
                std::cout<<"  warning: the frame of the file (lon0 "<<f.lon0<<", lat0 "<<f.lat0<<", rot "<<f.rot<<", X0 "<<f.X0<<", Y0 "<<f.Y0
                         <<") differs from W 205 / W 206"<<std::endl;
            }
        }
    }

    if(p->mpirank==0)
    {
    std::cout<<"  model time 0 = "<<met_date(tref)<<(p->W203>0.0 ? " (W 203)" : (p->W200==2 ? " (first record)" : ""))
             <<", drag law W 207 "<<p->W207;

        if(p->W204>0.0)
        std::cout<<", ramp "<<p->W204<<" s";

    std::cout<<std::endl;

    std::cout<<"  air pressure gradient (W 208) "<<(p->W208==1 ? "on" : "off")<<", inverse barometer on the open edges (W 210) "
             <<(p->W210==1 ? "on" : "off")<<", reference pressure "<<p->W209<<" Pa (W 209)"<<std::endl;

        if(frame.set)
        std::cout<<"  frame (W 205, W 206): transverse Mercator about lon "<<frame.lon0<<", lat "<<frame.lat0<<" deg, x axis "<<frame.rot
                 <<" deg from east, origin at ("<<frame.X0<<", "<<frame.Y0<<") m"<<std::endl;
    }
}

void met_provider::window(double xs, double xe, double ys, double ye)
{
    src->window(xs,xe,ys,ye);

    // a grown window reads the records again at the next update
    has_t = false;
}

void met_provider::update(lexer *p, ghostcell *pgc, double simtime)
{
    if(has_t && simtime==tnow)
    return;

    std::string err;

    if(!src->advance(tref + simtime,err))
    stop(p,pgc,"W 200: "+err);

    tnow = simtime;
    has_t = true;
}

void met_provider::at(double x, double y, double &u, double &v, double &pa) const
{
    src->at(tref + tnow,x,y,u,v,pa);
}

double met_provider::ramp(lexer *p, double simtime) const
{
    if(p->W204>0.0 && simtime<p->W204)
    return 0.5*(1.0 - std::cos(3.14159265358979323846*simtime/p->W204));

    return 1.0;
}

double met_provider::ib(lexer *p, ghostcell *pgc, double x, double y)
{
    update(p,pgc,p->simtime);

    double u, v, pa;
    at(x,y,u,v,pa);

    return -ramp(p,p->simtime)*(pa - p->W209)/(p->W1*std::fabs(p->W22));
}

void met_provider::check(lexer *p, ghostcell *pgc)
{
    std::string msg;

    if(p->W200<0 || p->W200>2)
    msg = "W 200: unknown met forcing; use 0 (off), 1 (uniform, W 201, W 202) or 2 (met file met-forcing.r3m).";

    else if(p->W200>0 && p->A10!=2 && p->A10!=5)
    msg = "W 200: the met forcing is for REEF3D::SFLOW (A 10 2) and REEF3D::NHFLOW (A 10 5).";

    else if(p->W200>0 && p->A10==5 && p->A570>0)
    msg = "W 200 and A 570: use either the met forcing (W 200) or the NHFLOW wind forcing (A 570-578), not both.";

    else if(p->W200>0 && p->G50>0)
    msg = "W 200: the met forcing does not support curvilinear grids (G 50) yet.";

    else if(p->W200>0 && (p->W207<1 || p->W207>4))
    msg = "W 207: unknown drag law; use 1 (Garratt), 2 (Wu), 3 (Smith and Banke) or 4 (Zijlema et al.).";

    else if(p->W200==1 && !(p->W201_u>=0.0))
    msg = "W 201: the wind speed must be zero or positive.";

    else if(p->W200>0 && (p->W208<0 || p->W208>1 || p->W210<0 || p->W210>1))
    msg = "W 208, W 210: 0 (off) or 1 (on).";

    else if(p->W211<0 || p->W211>2)
    msg = "W 211: unknown Coriolis option; use 0 (off), 1 (f-plane, W 212 or the latitude of W 205) or 2 (f per cell from the frame W 205).";

    else if(p->W211>0 && p->A10!=2 && p->A10!=5)
    msg = "W 211: the Coriolis force is for REEF3D::SFLOW (A 10 2) and REEF3D::NHFLOW (A 10 5).";

    else if(p->W211>0 && p->G50>0)
    msg = "W 211: the Coriolis force does not support curvilinear grids (G 50) yet.";

    else if(p->W211==1 && p->W212==0.0 && p->W205==0)
    msg = "W 211 1: give the Coriolis parameter (W 212) or the frame with the latitude (W 205).";

    else if(p->W211==2 && p->W205==0)
    msg = "W 211 2: the latitude of every cell needs the frame W 205.";

    if(msg.empty() && p->W211>0 && p->mpirank==0)
    {
    const double f0 = p->W212!=0.0 ? p->W212 : 2.0*7.2921159e-5*std::sin(p->W205_lat*3.14159265358979323846/180.0);

        if(p->W211==1)
        std::cout<<"REEF3D Coriolis force (W 211 1): f-plane, f = "<<f0<<" 1/s"<<(p->W212!=0.0 ? " (W 212)" : " (latitude of W 205)")<<std::endl;

        if(p->W211==2)
        std::cout<<"REEF3D Coriolis force (W 211 2): f = 2 Omega sin(latitude) of every cell from the frame W 205"<<std::endl;
    }

    if(!msg.empty())
    {
        if(p->mpirank==0)
        std::cout<<std::endl<<msg<<std::endl<<std::endl;

    pgc->final(true);
    }
}
