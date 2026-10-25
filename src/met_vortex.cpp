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

#include"met_vortex.h"
#include"met_file.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<sstream>

/*--------------------------------------------------------------------
Parametric tropical cyclone vortex, see met_vortex.h. Every formula
follows tools/met_track.py (evaluate, at_time) operation by operation.
--------------------------------------------------------------------*/

namespace
{
    const double pi = 3.14159265358979323846;
    const double deg = pi/180.0;
    const double Omega = 7.2921159e-5;
    const double Rearth = 6371000.0;
    const double quad[4] = {45.0, 135.0, 225.0, 315.0};      // NE SE SW NW, bearing clockwise from north
}

bool met_vortex::open(const std::string &file, std::string &error)
{
    std::ifstream in(file.c_str());

    if(!in.is_open())
    {
    error = "cannot open "+file;
    return false;
    }

    std::string line, magic;
    int version = 0, nrec = 0;

    while(std::getline(in,line))
    if(!line.empty() && line[0]!='#')
    break;

    std::istringstream hs(line);
    hs>>magic>>version>>nrec>>rho;

    if(magic!="REEF3D-MET-TRACK" || version!=1 || nrec<1 || !(rho>0.0))
    {
    error = file+": not a REEF3D vortex track (REEF3D-MET-TRACK 1 nrec rho_a)";
    return false;
    }

    while((int)rec.size()<nrec && std::getline(in,line))
    {
        if(line.empty() || line[0]=='#')
        continue;

    std::istringstream ls(line);
    double stamp;
    met_vortex_rec r;

    ls>>stamp>>r.lon>>r.lat>>r.ut>>r.vt>>r.pc>>r.pn>>r.vmaxg>>r.rmw>>r.B80>>r.vms>>r.bs>>r.xn>>r.rn>>r.Bg0>>r.phi0;

        for(int q=0; q<4; ++q)
        for(int k=0; k<3; ++k)
        ls>>r.R[q][k]>>r.Rm[q][k]>>r.Bg[q][k]>>r.phi[q][k];

        if(!ls)
        {
        error = file+": record "+std::to_string(rec.size()+1)+" incomplete (64 numbers)";
        return false;
        }

    r.t = met_datetime(stamp);

        if(!rec.empty() && !(r.t>rec.back().t))
        {
        error = file+": the record times must increase (record "+std::to_string(rec.size()+1)+")";
        return false;
        }

        if(!(r.pn>r.pc) || !(r.vmaxg>0.0) || !(r.rmw>0.0))
        {
        error = file+": record "+std::to_string(rec.size()+1)+" needs pn > pc, vmaxg > 0, rmw > 0";
        return false;
        }

    rec.push_back(r);
    }

    if((int)rec.size()!=nrec)
    {
    error = file+": "+std::to_string(rec.size())+" records, the header gives "+std::to_string(nrec);
    return false;
    }

    return true;
}

void met_vortex::dist_bearing(double lon1, double lat1, double lon2, double lat2, double &d, double &brg)
{
    const double p1 = lat1*deg, p2 = lat2*deg;
    const double dl = (lon2 - lon1)*deg;
    const double a = std::pow(std::sin((p2 - p1)/2.0),2.0) + std::cos(p1)*std::cos(p2)*std::pow(std::sin(dl/2.0),2.0);

    d = 2.0*Rearth*std::asin(std::sqrt(std::min(1.0,a)));

    double b = std::atan2(std::sin(dl)*std::cos(p2), std::cos(p1)*std::sin(p2) - std::sin(p1)*std::cos(p2)*std::cos(dl))*(180.0/pi);
    b = std::fmod(b,360.0);
    brg = b<0.0 ? b + 360.0 : b;
}

double met_vortex::fang(double r, double rm)
{
    if(r<rm)
    return 10.0*r/rm;

    if(r<1.2*rm)
    return 10.0 + 75.0*(r/rm - 1.0);

    return 25.0;
}

double met_vortex::holland80(double r, double rm, double B, double dp, double rho, double f)
{
    const double x = std::pow(rm/r,B);
    const double af = std::fabs(f);

    return std::sqrt(B*dp/rho*x*std::exp(-x) + std::pow(r*af/2.0,2.0)) - r*af/2.0;
}

double met_vortex::gahm(double r, double vmaxg, double rm, double Bg, double phi, double f)
{
    const double af = std::fabs(f);
    const double ro = af>0.0 ? vmaxg/(af*rm) : 1.0e30;
    const double x = std::pow(rm/r,Bg);

    return std::sqrt(std::pow(vmaxg,2.0)*(1.0 + 1.0/ro)*std::exp(phi*(1.0 - x))*x + std::pow(r*af/2.0,2.0)) - r*af/2.0;
}

double met_vortex::holland10(double r, double vms, double rm, double bs, double xn, double rn)
{
    const double g = std::pow(rm/r,bs)*std::exp(1.0 - std::pow(rm/r,bs));
    const double x = (r<=rm || rn<=rm) ? 0.5 : 0.5 + (r - rm)*(xn - 0.5)/(rn - rm);

    return vms*std::pow(g,x);
}

void met_vortex::quad_params(const met_vortex_rec &r, int q, double rr, double &rm, double &bg, double &ph) const
{
    // the isotachs of the quadrant, outermost first
    int idx[3], n = 0;

    for(int k=0; k<3; ++k)
    if(r.R[q][k]>0.0)
    idx[n++] = k;

    if(n==0)
    {
    rm = r.rmw;
    bg = r.Bg0;
    ph = r.phi0;
    return;
    }

    std::stable_sort(idx,idx+n,[&](int a, int b) {return r.R[q][a]>r.R[q][b];});

    if(rr>=r.R[q][idx[0]])
    {
    rm = r.Rm[q][idx[0]];  bg = r.Bg[q][idx[0]];  ph = r.phi[q][idx[0]];
    return;
    }

    if(rr<=r.R[q][idx[n-1]])
    {
    rm = r.Rm[q][idx[n-1]];  bg = r.Bg[q][idx[n-1]];  ph = r.phi[q][idx[n-1]];
    return;
    }

    for(int m=0; m<n-1; ++m)
    {
    const int a = idx[m], b = idx[m+1];

        if(r.R[q][b]<=rr && rr<=r.R[q][a])
        {
        const double w = (rr - r.R[q][b])/(r.R[q][a] - r.R[q][b]);
        rm = w*r.Rm[q][a]  + (1.0 - w)*r.Rm[q][b];
        bg = w*r.Bg[q][a]  + (1.0 - w)*r.Bg[q][b];
        ph = w*r.phi[q][a] + (1.0 - w)*r.phi[q][b];
        return;
        }
    }
}

void met_vortex::azimuth_params(const met_vortex_rec &r, double theta, double rr, double &rm, double &bg, double &ph) const
{
    double a = std::fmod(theta - 45.0,360.0);

    if(a<0.0)
    a += 360.0;

    if(a>=360.0)
    a = 0.0;

    const int i = std::min(int(std::floor(a/90.0)),3);
    const int j = (i + 1)%4;
    const double d = a - 90.0*i;

    double rmi, bgi, phi_i, rmj, bgj, phj;
    quad_params(r,i,rr,rmi,bgi,phi_i);
    quad_params(r,j,rr,rmj,bgj,phj);

    if(d<1.0)
    {
    rm = rmi;  bg = bgi;  ph = phi_i;
    return;
    }

    if(d>89.0)
    {
    rm = rmj;  bg = bgj;  ph = phj;
    return;
    }

    const double wi = 1.0/std::pow(d,2.0), wj = 1.0/std::pow(90.0 - d,2.0);

    rm = (wi*rmi + wj*rmj)/(wi + wj);
    bg = (wi*bgi + wj*bgj)/(wi + wj);
    ph = (wi*phi_i + wj*phj)/(wi + wj);
}

void met_vortex::eval(const met_vortex_rec &r, double lonc, double latc, double ut, double vt,
                      double lon, double lat, double &ue, double &vn, double &pa) const
{
    double rr, theta;
    dist_bearing(lonc,latc,lon,lat,rr,theta);

    if(rr<1000.0)
    {
    ue = 0.0;
    vn = 0.0;
    pa = r.pc;
    return;
    }

    const double f = 2.0*Omega*std::sin(latc*deg);
    const double s = f>=0.0 ? 1.0 : -1.0;
    const double dp = r.pn - r.pc;
    double rm, V, gam, sp;

    if(par.model==1)
    {
    rm = r.rmw;
    V = holland80(rr,rm,r.B80,dp,rho,f);
    pa = r.pc + dp*std::exp(-std::pow(rm/rr,r.B80));
    gam = V/r.vmaxg;
    sp = par.wrf*V;
    }
    else if(par.model==2)
    {
    double bg, ph;
    azimuth_params(r,theta,rr,rm,bg,ph);
    V = gahm(rr,r.vmaxg,rm,bg,ph,f);
    pa = r.pc + dp*std::exp(-ph*std::pow(rm/rr,bg));
    gam = V/r.vmaxg;
    sp = par.wrf*V;
    }
    else
    {
    rm = r.rmw;
    V = holland10(rr,r.vms,rm,r.bs,r.xn,r.rn);
    pa = r.pc + dp*std::exp(-std::pow(rm/rr,r.bs));
    gam = V/r.vms;
    sp = V;
    }

    const double beta = par.inflow==1 ? fang(rr,rm) : 0.0;
    const double ang = ((90.0 - theta) + s*(90.0 + beta))*deg;

    ue = sp*std::cos(ang);
    vn = sp*std::sin(ang);

    if(par.trans==1)
    {
    ue = ue + gam*ut;
    vn = vn + gam*vt;
    }
    else if(par.trans==2)
    {
    ue = ue + ut;
    vn = vn + vt;
    }

    ue = par.one2ten*ue;
    vn = par.one2ten*vn;
}

int met_vortex::bracket(double t, double &w) const
{
    w = 0.0;
    const int n = (int)rec.size();

    if(n==1 || t<=rec[0].t)
    return 0;

    if(t>=rec[n-1].t)
    return n-1;

    // the last record at or before t (numpy searchsorted, side right, minus one)
    int k = int(std::upper_bound(rec.begin(),rec.end(),t,[](double v, const met_vortex_rec &r) {return v<r.t;}) - rec.begin()) - 1;
    w = (t - rec[k].t)/(rec[k+1].t - rec[k].t);
    return k;
}

void met_vortex::centre(double t, double &lon, double &lat) const
{
    double w;
    const int k = bracket(t,w);
    const met_vortex_rec &a = rec[k];
    const met_vortex_rec &b = w>0.0 ? rec[k+1] : rec[k];

    lon = a.lon + w*(b.lon - a.lon);
    lat = a.lat + w*(b.lat - a.lat);
}

void met_vortex::at(double t, double lon, double lat, double &ue, double &vn, double &pa) const
{
    double w;
    const int k = bracket(t,w);
    const met_vortex_rec &a = rec[k];
    const met_vortex_rec &b = w>0.0 ? rec[k+1] : rec[k];

    const double lonc = a.lon + w*(b.lon - a.lon), latc = a.lat + w*(b.lat - a.lat);
    const double ut = a.ut + w*(b.ut - a.ut), vt = a.vt + w*(b.vt - a.vt);

    eval(a,lonc,latc,ut,vt,lon,lat,ue,vn,pa);

    if(w==0.0)
    return;

    double ub, vb, pb;
    eval(b,lonc,latc,ut,vt,lon,lat,ub,vb,pb);

    ue = ue + w*(ub - ue);
    vn = vn + w*(vb - vn);
    pa = pa + w*(pb - pa);
}
