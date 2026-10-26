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

#include"met_blend.h"
#include<algorithm>
#include<cmath>

/*--------------------------------------------------------------------
Composite met source W 200 4, see met_blend.h. The formulas follow
tools/met_forcing.py (class Composite) operation by operation.
--------------------------------------------------------------------*/

namespace
{
    const double pi = 3.14159265358979323846;
    const int nring = 16;

    // f0 + w (f1 - f0); w = 0 and w = 1 give f0 and f1 bitwise
    inline double mix(double f0, double f1, double w)
    {
        if(w<=0.0)
        return f0;

        if(w>=1.0)
        return f1;

        return f0 + w*(f1 - f0);
    }
}

double met_source_blend::taper(double s, int shape)
{
    if(s<=0.0)
    return 0.0;

    if(s>=1.0)
    return 1.0;

    return shape==2 ? s : 0.5*(1.0 - std::cos(pi*s));
}

double met_source_blend::taper_time(double t, double t0, double t1, double T, int shape)
{
    if(t<t0 || t>t1)
    return 0.0;

    if(!(T>0.0))
    return 1.0;

    return taper(std::min(t - t0, t1 - t)/T, shape);
}

bool met_source_blend::open(const std::string &bgfile, const std::string &nestfile, const std::string &trkfile, std::string &error)
{
    if(!bg.open(bgfile,error))
    return false;

    if(par.nest==1)
    {
        if(!nest.open(nestfile,error))
        return false;

        if(nest.nx<2 || nest.ny<2)
        {
        error = nestfile+": a nest needs a grid of at least 2 x 2 nodes";
        return false;
        }
    }

    if(par.vortex==1 || (par.nest==1 && par.nest_relative==1))
    if(!vx.open(trkfile,error))
    return false;

    // W 227 1: the largest move of the storm between two nest records (window margin of the nest)
    if(par.nest==1 && par.nest_relative==1)
    {
    double xp = 0.0, yp = 0.0;

        for(int n=0; n<nest.nrec; ++n)
        {
        double lon, lat, x, y;
        vx.centre(nest.time(n),lon,lat);
        frame.to_model(lon,lat,x,y);

            if(n>0)
            nest_margin = std::max(nest_margin, std::sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp)));

        xp = x;
        yp = y;
        }

    nest_margin += 2.0*std::max(nest.dx,nest.dy);
    }

    return true;
}

void met_source_blend::window(double xs, double xe, double ys, double ye)
{
    bg.window(xs,xe,ys,ye);

    if(par.nest==1)
    nest.window(xs-nest_margin,xe+nest_margin,ys-nest_margin,ye+nest_margin);
}

void met_source_blend::radii(double t, double &a, double &b) const
{
    double lon, lat, p0, p1, rm, r3;
    vx.state(t,lon,lat,p0,p1,rm,r3);

    const double scale = par.basis==1 ? rm : (par.basis==2 ? r3 : 1.0);

    a = par.r1*scale;
    b = par.r2*scale;
}

bool met_source_blend::ring(int n, met_blend_ring &rg, std::string &error)
{
    const double t = bg.time(n);
    double lon, lat, p0, p1, rm, r3, a, b, xc, yc;
    vx.state(t,lon,lat,p0,p1,rm,r3);
    radii(t,a,b);
    frame.to_model(lon,lat,xc,yc);

    double X[nring], Y[nring], R[nring];
    double su = 0.0, sv = 0.0, sp = 0.0, spv = 0.0, sx = 0.0, sy = 0.0, sr = 0.0;

    for(int m=0; m<nring; ++m)
    {
    double lon2, lat2, x, y, ub = 0.0, vb = 0.0, pb;
    met_vortex::destination(lon,lat,360.0*double(m)/double(nring),b,lon2,lat2);
    frame.to_model(lon2,lat2,x,y);

        if(!bg.sample(n,2,x,y,pb,error))
        return false;

        if(par.env==1)
        if(!bg.sample(n,0,x,y,ub,error) || !bg.sample(n,1,x,y,vb,error))
        return false;

    // the vortex at the same point and time, wind along the model axes
    double ue, vn, pv;
    vx.at(t,lon2,lat2,ue,vn,pv);
    const double chi = frame.north_geo(lon2,lat2);
    const double uv = ue*std::sin(chi) + vn*std::cos(chi);
    const double vv = -ue*std::cos(chi) + vn*std::sin(chi);

    X[m] = x - xc;
    Y[m] = y - yc;
    R[m] = pb - pv;
    sp += pb;  spv += pv;  sr += R[m];  sx += X[m];  sy += Y[m];
    su += ub - uv;  sv += vb - vv;
    }

    rg.p  = sp/double(nring);
    rg.pv = spv/double(nring);
    rg.u  = par.env==1 ? su/double(nring) : 0.0;
    rg.v  = par.env==1 ? sv/double(nring) : 0.0;

    // gradient of the background minus the vortex: least-squares plane through the ring
    const double mx = sx/double(nring), my = sy/double(nring), mr = sr/double(nring);
    double sxx = 0.0, syy = 0.0, sxy = 0.0, sxp = 0.0, syp = 0.0;

    for(int m=0; m<nring; ++m)
    {
    const double dx = X[m] - mx, dy = Y[m] - my, dr = R[m] - mr;
    sxx += dx*dx;  syy += dy*dy;  sxy += dx*dy;  sxp += dx*dr;  syp += dy*dr;
    }

    const double det = sxx*syy - sxy*sxy;
    rg.gx = det>0.0 ? (sxp*syy - syp*sxy)/det : 0.0;
    rg.gy = det>0.0 ? (syp*sxx - sxp*sxy)/det : 0.0;

    ++ring_samples;

    return true;
}

bool met_source_blend::advance(double t, std::string &error)
{
    tnow = t;

    if(!bg.advance(t,error))
    return false;

    // nest: time weight and the grid at time t
    wn_t = 0.0;

    if(par.nest==1)
    {
        if(!nest.advance(t,error))
        return false;

    wn_t = taper_time(t,nest.time(0),nest.time(nest.nrec-1),par.nest_time,par.shape);
    nest.origin(t,nest_x0,nest_y0);
    nest_lx = double(nest.nx-1)*nest.dx;
    nest_ly = double(nest.ny-1)*nest.dy;

        if(par.nest_relative==1)
        {
        // shifts of the two records with the storm centre, the grid at time t from the shifted records
        double lon, lat, xc, yc;
        vx.centre(t,lon,lat);
        frame.to_model(lon,lat,xc,yc);

        nest_w = nest.time_weight(t);

        double o[2][2];

            for(int q=0; q<2; ++q)
            {
            double xq, yq;
            vx.centre(nest.rec_time(q),lon,lat);
            frame.to_model(lon,lat,xq,yq);
            nest_d[q][0] = xc - xq;
            nest_d[q][1] = yc - yq;
            nest.rec_origin(q,o[q][0],o[q][1]);
            o[q][0] += nest_d[q][0];
            o[q][1] += nest_d[q][1];
            }

        nest_x0 = nest_w>0.0 ? o[0][0] + nest_w*(o[1][0] - o[0][0]) : o[0][0];
        nest_y0 = nest_w>0.0 ? o[0][1] + nest_w*(o[1][1] - o[0][1]) : o[0][1];
        }
    }

    // vortex: time weight, centre, radii, the background on the ring
    wv_t = 0.0;

    if(par.vortex==1)
    {
    wv_t = taper_time(t,vx.t_first(),vx.t_last(),par.vortex_time,par.shape);
    vx.state(t,lonc,latc,pc,pn,rmw,r34);
    radii(t,r1m,r2m);
    frame.to_model(lonc,latc,xcm,ycm);

    env = met_blend_ring();
    env.p = pn;
    env.pv = pn;

        if(par.pressure==1 || par.env==1)
        {
        // background records around t, as met_file::advance and weight
        int k = 0;
        double w = 0.0;

            if(bg.nrec>1)
            {
                while(k+1<bg.nrec && bg.time(k+1)<=t)
                ++k;

            k = std::min(std::max(k,0),bg.nrec-2);
            w = (t - bg.time(k))/(bg.time(k+1) - bg.time(k));
            w = std::min(std::max(w,0.0),1.0);
            }

        const int k1 = bg.nrec>1 ? k+1 : k;
        const int need[2] = {k, k1};

            for(int q=0; q<2; ++q)
            {
                if(ring_n[q]==need[q])
                continue;

                // reuse the other slot
                if(ring_n[1-q]==need[q])
                {
                std::swap(ring_n[0],ring_n[1]);
                std::swap(ring_s[0],ring_s[1]);

                    if(ring_n[q]==need[q])
                    continue;
                }

                if(!ring(need[q],ring_s[q],error))
                return false;

            ring_n[q] = need[q];
            }

        const met_blend_ring &a = ring_s[0], &b = ring_s[1];
        auto lin = [w](double f0, double f1) {return w>0.0 ? f0 + w*(f1 - f0) : f0;};

        env.p  = lin(a.p,b.p);
        env.pv = lin(a.pv,b.pv);
        env.u  = lin(a.u,b.u);
        env.v  = lin(a.v,b.v);
        env.gx = lin(a.gx,b.gx);
        env.gy = lin(a.gy,b.gy);
        }
    }

    return true;
}

void met_source_blend::at(double t, double x, double y, double &u, double &v, double &pa) const
{
    bg.at(t,x,y,u,v,pa);

    // nest
    if(wn_t>0.0)
    {
    const double d = std::min(std::min(x - nest_x0, nest_x0 + nest_lx - x), std::min(y - nest_y0, nest_y0 + nest_ly - y));
    const double w = wn_t*taper(d/par.nest_width,par.shape);

        if(w>0.0)
        {
        double un, vn, pn_;

            if(par.nest_relative==1)
            {
            nest.at_record(0,x-nest_d[0][0],y-nest_d[0][1],un,vn,pn_);

                if(nest_w>0.0)
                {
                double u1, v1, p1;
                nest.at_record(1,x-nest_d[1][0],y-nest_d[1][1],u1,v1,p1);
                un  = un  + nest_w*(u1 - un);
                vn  = vn  + nest_w*(v1 - vn);
                pn_ = pn_ + nest_w*(p1 - pn_);
                }
            }
            else
            nest.at(t,x,y,un,vn,pn_);

        u = mix(u,un,w);
        v = mix(v,vn,w);

            if(nest.nvar==3)
            pa = mix(pa,pn_,w);
        }
    }

    // vortex
    if(wv_t>0.0)
    {
    double lon, lat, r, brg;
    frame.to_geo(x,y,lon,lat);
    met_vortex::dist_bearing(lonc,latc,lon,lat,r,brg);

    const double w = wv_t*taper((r2m - r)/(r2m - r1m),par.shape);

        if(w>0.0)
        {
        double ue, vn, pv;
        vx.at(t,lon,lat,ue,vn,pv);

        // true north at angle chi from the model x axis (as met_source_vortex)
        const double chi = frame.north_geo(lon,lat);
        const double uv = ue*std::sin(chi) + vn*std::cos(chi);
        const double vv = -ue*std::cos(chi) + vn*std::sin(chi);

            if(par.pressure==1)
            pv = pc + (pv - pc)*std::max((env.p - pc)/(env.pv - pc),0.0);

        double uw = uv, vw = vv;

            // environment: the mean background wind and the pressure gradient on the ring
            if(par.env==1)
            {
            uw = uv + env.u;
            vw = vv + env.v;
            pv = pv + env.gx*(x - xcm) + env.gy*(y - ycm);
            }

        u = mix(u,uw,w);
        v = mix(v,vw,w);
        pa = mix(pa,pv,w);
        }
    }
}
