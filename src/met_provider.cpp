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
#include"met_blend.h"
#include"met_drag.h"
#include"lexer.h"
#include"ghostcell.h"
#include<cmath>
#include<iostream>

met_provider *met_provider::instance = nullptr;

namespace
{
    // frame of a met file against W 205 / W 206 (rank 0 log)
    void frame_check(lexer *p, const met_file &f, const char *name)
    {
        if(f.frame==1 && p->W205==1)
        {
        const double d = std::fabs(f.lon0-p->W205_lon) + std::fabs(f.lat0-p->W205_lat) + std::fabs(f.rot-p->W205_rot);
        const double o = std::fabs(f.X0-p->W206_x) + std::fabs(f.Y0-p->W206_y);

            if(d>1.0e-9 || o>1.0e-6)
            std::cout<<"  warning: the frame of "<<name<<" (lon0 "<<f.lon0<<", lat0 "<<f.lat0<<", rot "<<f.rot<<", X0 "<<f.X0<<", Y0 "<<f.Y0
                     <<") differs from W 205 / W 206"<<std::endl;
        }
    }
}

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

            frame_check(p,f,"met-forcing.r3m");
        }
    }

    if(p->W200==3)
    {
    met_source_vortex *vs = new met_source_vortex(frame);
    vs->v.par.model = p->W215;
    vs->v.par.wrf = p->W216;
    vs->v.par.one2ten = p->W217;
    vs->v.par.inflow = p->W218;
    vs->v.par.trans = p->W219;

        if(!vs->open("met-vortex.trk",err))
        stop(p,pgc,"W 200 3: "+err);

    src = vs;

    const met_vortex &v = vs->v;

    tref = (p->W203>0.0) ? met_datetime(p->W203) : v.rec.front().t;

        if(p->mpirank==0)
        {
        const char *mname[4] = {"", "Holland (1980)", "GAHM", "Holland et al. (2010)"};
        double lon, lat;
        v.centre(tref,lon,lat);

        std::cout<<"REEF3D met forcing (W 200 3): vortex "<<mname[p->W215]<<" (W 215) from met-vortex.trk, "<<v.rec.size()<<" records "
                 <<met_date(v.rec.front().t)<<" - "<<met_date(v.rec.back().t)<<", rho_a "<<v.rho<<std::endl;
        std::cout<<"  Wrf "<<p->W216<<" (W 216), 1-min to 10-min "<<p->W217<<" (W 217), inflow angle "<<(p->W218==1 ? "ADCIRC" : "none")
                 <<" (W 218), translation "<<(p->W219==1 ? "scaled by Vg/Vmax" : (p->W219==2 ? "full" : "none"))<<" (W 219)"<<std::endl;
        std::cout<<"  centre at model time 0: lon "<<lon<<", lat "<<lat<<std::endl;

            if(tref + p->N41 < v.rec.front().t || tref > v.rec.back().t)
            std::cout<<"  warning: the run lies outside the records, the vortex is held at the first or last record"<<std::endl;
        }
    }

    if(p->W200==4)
    {
    met_blend_par bp;
    bp.nest = p->W220_n;
    bp.vortex = p->W220_v;
    bp.nest_width = p->W221_w;
    bp.nest_time = p->W221_t;
    bp.basis = p->W222;
    bp.r1 = p->W222_r1;
    bp.r2 = p->W222_r2;
    bp.pressure = p->W223;
    bp.vortex_time = p->W224;
    bp.shape = p->W225;
    bp.env = p->W226;
    bp.nest_relative = p->W227;

    met_source_blend *bs = new met_source_blend(frame,bp);
    bs->vx.par.model = p->W215;
    bs->vx.par.wrf = p->W216;
    bs->vx.par.one2ten = p->W217;
    bs->vx.par.inflow = p->W218;
    bs->vx.par.trans = p->W219;

        if(!bs->open("met-forcing.r3m","met-nest.r3m","met-vortex.trk",err))
        stop(p,pgc,"W 200 4: "+err);

    src = bs;

    const met_file &f = bs->bg;

    tref = (p->W203>0.0) ? met_datetime(p->W203) : f.time(0);

        if(p->mpirank==0)
        {
        std::cout<<"REEF3D met forcing (W 200 4): blend of"<<std::endl;
        std::cout<<"  background met-forcing.r3m, "<<f.nx<<" x "<<f.ny<<" nodes from ("<<f.x0<<", "<<f.y0<<") m, spacing "<<f.dx<<" x "<<f.dy
                 <<" m, "<<f.nrec<<" records "<<met_date(f.time(0))<<" - "<<met_date(f.time(f.nrec-1))
                 <<(f.nvar==3 ? ", u10 v10 pa" : ", u10 v10 (pa 101325 Pa)")<<std::endl;
        frame_check(p,f,"met-forcing.r3m");

            if(bp.nest==1)
            {
            const met_file &n = bs->nest;
            std::cout<<"  nest met-nest.r3m, "<<n.nx<<" x "<<n.ny<<" nodes, spacing "<<n.dx<<" x "<<n.dy<<" m"<<(n.moving==1 ? ", moving" : "")
                     <<", "<<n.nrec<<" records "<<met_date(n.time(0))<<" - "<<met_date(n.time(n.nrec-1))
                     <<(n.nvar==3 ? ", u10 v10 pa" : ", u10 v10 (pressure of the background)")
                     <<"; taper "<<bp.nest_width<<" m at the edges, "<<bp.nest_time<<" s at the ends (W 221)"<<std::endl;
            std::cout<<"  nest in time (W 227): "<<(bp.nest_relative==1 ? "storm-relative, the records shifted with the track centre (window margin "
                                                                          +std::to_string(int(bs->nest_margin))+" m)" : "linear at fixed points")<<std::endl;
            frame_check(p,n,"met-nest.r3m");
            }

            if(bp.vortex==1)
            {
            const char *mname[4] = {"", "Holland (1980)", "GAHM", "Holland et al. (2010)"};
            const char *bname[4] = {"", "x rmw", "x mean 34 kt radius", "m"};
            const met_vortex &v = bs->vx;
            std::cout<<"  vortex "<<mname[p->W215]<<" (W 215) from met-vortex.trk, "<<v.rec.size()<<" records "
                     <<met_date(v.t_first())<<" - "<<met_date(v.t_last())<<", rho_a "<<v.rho
                     <<", Wrf "<<p->W216<<", 1-min to 10-min "<<p->W217<<", inflow "<<p->W218<<", translation "<<p->W219<<std::endl;
            std::cout<<"  vortex inside "<<bp.r1<<", background beyond "<<bp.r2<<" "<<bname[bp.basis]<<" (W 222), "
                     <<(bp.shape==1 ? "cosine" : "linear")<<" taper (W 225), "<<bp.vortex_time<<" s at the ends of the track (W 224), pressure "
                     <<(bp.pressure==1 ? "deficit rescaled to the background on the outer ring (W 223 1)" : "as fitted (W 223 0)")<<std::endl;
            std::cout<<"  environment (W 226): "<<(bp.env==1 ? "the background minus the vortex on the outer ring (mean wind, pressure gradient) added to the vortex" : "off, the vortex alone inside r1")<<std::endl;

                if(tref + p->N41 < v.t_first() || tref > v.t_last())
                std::cout<<"  warning: the run lies outside the track, the vortex is off"<<std::endl;
            }

            if(tref + p->N41 < f.time(0) || tref > f.time(f.nrec-1))
            std::cout<<"  warning: the run lies outside the background records, the background is held at the first or last record"<<std::endl;
        }
    }

    if(p->mpirank==0)
    {
    const int law = met_drag_law(p->W207,p->W200);

    std::cout<<"  model time 0 = "<<met_date(tref)<<(p->W203>0.0 ? " (W 203)" : (p->W200>=2 ? (p->W200==4 ? " (first background record)" : " (first record)") : ""))
             <<", drag law W 207 "<<law<<(p->W207==0 ? " (auto)" : "");

        if(law==5)
        std::cout<<" (cap "<<p->W213<<", W 213)";

        if(law==7)
        std::cout<<" (wave-dependent, Taylor and Yelland 2001 with the SEASTATE Hs, Tp"<<(p->W229==1 ? ", Charnock floor" : "")<<", cap "<<p->W228
                 <<" (W 228); law 5 with the cap "<<p->W213<<" where there are no waves)";

        if(law==6)
        std::cout<<" (floor "<<p->W214<<", W 214)";

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

    if(p->W200<0 || p->W200>4)
    msg = "W 200: unknown met forcing; use 0 (off), 1 (uniform, W 201, W 202), 2 (met file met-forcing.r3m), 3 (vortex from met-vortex.trk) or 4 (blend of met-forcing.r3m, met-nest.r3m and met-vortex.trk, W 220-225).";

    else if(p->W200>0 && p->A10!=2 && p->A10!=5 && !(p->A10==7 && p->A730==2))
    msg = "W 200: the met forcing is for REEF3D::SFLOW (A 10 2), REEF3D::NHFLOW (A 10 5) and the wind of REEF3D::SEASTATE (A 10 7 with A 730 2).";

    else if(p->W200>0 && p->A10==5 && p->A570>0)
    msg = "W 200 and A 570: use either the met forcing (W 200) or the NHFLOW wind forcing (A 570-578), not both.";

    else if(p->W200>0 && p->G50==2)
    msg = "W 200: the met forcing supports the curvilinear grids of DIVEMesh (G 50 1), not curvgrid.dat (G 50 2) yet.";

    else if(p->W200>0 && (p->W207<0 || p->W207>7))
    msg = "W 207: unknown drag law; use 0 (auto), 1 (Garratt), 2 (Wu), 3 (Smith and Banke), 4 (Zijlema et al.), 5 (Garratt capped at W 213), 6 (Zijlema et al. with the floor W 214) or 7 (wave-dependent, Taylor and Yelland, with the SEASTATE coupling A 750).";

    else if(p->W200>0 && p->W207==7 && p->A10!=7 && p->A750!=1 && p->A750!=2)
    msg = "W 207 7: the wave-dependent drag needs the sea state of the SEASTATE coupling (A 750 1 or 2) of SFLOW or NHFLOW.";

    else if(p->W200>0 && p->W207==7 && (!(p->W228>0.0) || (p->W229!=0 && p->W229!=1)))
    msg = "W 228, W 229: the cap of W 207 7 must be positive, the floor 0 (none) or 1 (Charnock).";

    else if(p->W200>0 && (!(p->W213>0.0) || !(p->W214>=0.0)))
    msg = "W 213, W 214: the cap must be positive, the floor zero or positive.";

    else if((p->W200==3 || (p->W200==4 && (p->W220_v==1 || p->W227==1))) && p->W205==0)
    msg = "W 200 3/4: the vortex is given in longitude and latitude; the frame W 205 (and W 206) places it on the grid.";

    else if(p->W200>=3 && (p->W215<1 || p->W215>3))
    msg = "W 215: unknown vortex model; use 1 (Holland 1980), 2 (GAHM) or 3 (Holland et al. 2010).";

    else if(p->W200==4 && (p->W220_n<0 || p->W220_n>1 || p->W220_v<0 || p->W220_v>1))
    msg = "W 220 n v: the nest layer (met-nest.r3m) and the vortex layer (met-vortex.trk), each 0 (off) or 1 (on).";

    else if(p->W200==4 && p->W220_n==1 && (!(p->W221_w>0.0) || !(p->W221_t>=0.0)))
    msg = "W 221 width T: the taper width at the nest edges must be positive, the taper time zero or positive.";

    else if(p->W200==4 && p->W220_v==1 && (p->W222<1 || p->W222>3 || !(p->W222_r1>=0.0) || !(p->W222_r2>p->W222_r1)))
    msg = "W 222 basis r1 r2: basis 1 (rmw), 2 (mean 34 kt radius) or 3 (metres), 0 <= r1 < r2.";

    else if(p->W200==4 && (p->W223<0 || p->W223>1 || !(p->W224>=0.0) || p->W225<1 || p->W225>2 || p->W226<0 || p->W226>1 || p->W227<0 || p->W227>1))
    msg = "W 223 0/1, W 224 >= 0, W 225 1 (cosine) or 2 (linear), W 226 0/1, W 227 0/1.";

    else if(p->W200>=3 && (p->W218<0 || p->W218>1 || p->W219<0 || p->W219>2 || !(p->W216>0.0) || !(p->W217>0.0)))
    msg = "W 216 - W 219: W 216, W 217 positive, W 218 0 or 1, W 219 0, 1 or 2.";

    else if(p->W200==1 && !(p->W201_u>=0.0))
    msg = "W 201: the wind speed must be zero or positive.";

    else if(p->W200>0 && (p->W208<0 || p->W208>1 || p->W210<0 || p->W210>1))
    msg = "W 208, W 210: 0 (off) or 1 (on).";

    else if(p->W211<0 || p->W211>2)
    msg = "W 211: unknown Coriolis option; use 0 (off), 1 (f-plane, W 212 or the latitude of W 205) or 2 (f per cell from the frame W 205).";

    else if(p->W211>0 && p->A10!=2 && p->A10!=5)
    msg = "W 211: the Coriolis force is for REEF3D::SFLOW (A 10 2) and REEF3D::NHFLOW (A 10 5).";

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
