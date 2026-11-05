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

#include"6DOF_deform.h"
#include"6DOF_geometry.h"
#include<cmath>
#include<algorithm>
#include<utility>

namespace
{
    const double fish_pi = 3.14159265358979323846;
    const int fish_bins = 256;
}

sixdof_deform::sixdof_deform(int type, double a, double lambda, double f, double tramp, int lateral, bool recoil)
              : type_(type), lat_(lateral), a_(a), lambda_(lambda>0.0 ? lambda : 1.0), f_(f), tramp_(tramp), recoil_(recoil)
{
}

void sixdof_deform::setup(const sixdof_geometry &g)
{
    // body length and head (largest x of the body frame)
    xhead_ = -1.0e20;
    xtail_ =  1.0e20;

    for(int n=0; n<g.tricount; ++n)
    for(int q=0; q<3; ++q)
    {
        xhead_ = std::max(xhead_, g.tri_x0[n][q]);
        xtail_ = std::min(xtail_, g.tri_x0[n][q]);
    }

    L_ = (xhead_>xtail_) ? xhead_ - xtail_ : 1.0;

    // volume left of the plane x = xi:  V(xi) = int_V H(xi - x) dV = oint min(x,xi) n_x dS
    // (divergence theorem with the field (min(x,xi),0,0)), with a midpoint rule on a uniform
    // subdivision of every triangle into M^2 triangles of an x-extent of at most half a bin
    const double dxb = L_/double(fish_bins);
    std::vector<std::pair<double,double> > qp;     // x of the quadrature point, n_x dS

    for(int n=0; n<g.tricount; ++n)
    {
        const double ax = g.tri_x0[n][0], ay = g.tri_y0[n][0], az = g.tri_z0[n][0];
        const double bx = g.tri_x0[n][1], by = g.tri_y0[n][1], bz = g.tri_z0[n][1];
        const double cx = g.tri_x0[n][2], cy = g.tri_y0[n][2], cz = g.tri_z0[n][2];

        // x component of the area vector (outward orientation of the hull)
        const double Sx = 0.5*((by-ay)*(cz-az) - (bz-az)*(cy-ay));

        if(Sx==0.0)
        continue;

        const double ext = std::max(ax,std::max(bx,cx)) - std::min(ax,std::min(bx,cx));
        const int M = std::min(64, std::max(1, int(std::ceil(ext/(0.5*dxb)))));
        const double w = Sx/double(M*M);

        for(int i=0; i<M; ++i)
        for(int j=0; j<M-i; ++j)
        {
            // upward triangle, centroid at barycentric (i+1/3, j+1/3)/M
            double u = (i + 1.0/3.0)/double(M);
            double v = (j + 1.0/3.0)/double(M);
            qp.push_back(std::make_pair(ax + u*(bx-ax) + v*(cx-ax), w));

            // downward triangle
            if(i+j<M-1)
            {
            u = (i + 2.0/3.0)/double(M);
            v = (j + 2.0/3.0)/double(M);
            qp.push_back(std::make_pair(ax + u*(bx-ax) + v*(cx-ax), w));
            }
        }
    }

    std::sort(qp.begin(), qp.end());

    const int nq = int(qp.size());
    std::vector<double> sw(nq+1,0.0), swx(nq+1,0.0);    // prefix sums of w and w x

    for(int q=0; q<nq; ++q)
    {
        sw[q+1]  = sw[q]  + qp[q].second;
        swx[q+1] = swx[q] + qp[q].second*qp[q].first;
    }

    // V(xi) = sum_{x<xi} w x + xi sum_{x>=xi} w
    auto Vleft = [&](double xi)
    {
        const int m = int(std::lower_bound(qp.begin(), qp.end(), std::make_pair(xi,-1.0e300)) - qp.begin());
        return swx[m] + xi*(sw[nq] - sw[m]);
    };

    // closed hull: V(x_head) = volume; an inward oriented hull gives a negative volume
    vol_ = Vleft(xhead_) - Vleft(xtail_);
    const double sgn = (vol_<0.0) ? -1.0 : 1.0;
    vol_ = std::fabs(vol_);

    xb_.assign(fish_bins,0.0);
    mb_.assign(fish_bins,0.0);
    S0_ = S1_ = S2_ = 0.0;

    double Vl = Vleft(xtail_);

    for(int b=0; b<fish_bins; ++b)
    {
        const double xr = xtail_ + (b+1)*dxb;
        const double Vr = Vleft(xr);

        xb_[b] = xtail_ + (b+0.5)*dxb;
        mb_[b] = std::max(0.0, sgn*(Vr - Vl));
        Vl = Vr;

        S0_ += mb_[b];
        S1_ += mb_[b]*xb_[b];
        S2_ += mb_[b]*xb_[b]*xb_[b];
    }
}

double sixdof_deform::envelope(double s) const
{
    if(type_==2)
    return (a_/0.1)*(0.02 - 0.0825*s + 0.1625*s*s);

    if(type_==3)
    return a_*std::exp(s - 1.0);

    return a_*(s + 0.03125)/1.03125;
}

void sixdof_deform::raw(double x, double &h, double &ht) const
{
    double s = (xhead_ - x)/L_;
    s = std::min(1.0, std::max(0.0, s));

    const double A = L_*envelope(s);
    const double om = 2.0*fish_pi*f_;
    const double ph = 2.0*fish_pi*s/lambda_ - om*t_;

    h  = r_*A*std::sin(ph);
    ht = -r_*A*om*std::cos(ph) + rdot_*A*std::sin(ph);
}

void sixdof_deform::fit(const std::vector<double> &g, double &k0, double &k1) const
{
    // k0 + k1 x: least squares fit of g with the weights of the mass distribution
    double G0 = 0.0, G1 = 0.0;

    for(int b=0; b<int(xb_.size()); ++b)
    {
        G0 += mb_[b]*g[b];
        G1 += mb_[b]*xb_[b]*g[b];
    }

    const double det = S0_*S2_ - S1_*S1_;

    if(S0_<=0.0 || std::fabs(det)<=1.0e-30*S0_*S2_)
    {
    k0 = k1 = 0.0;
    return;
    }

    k1 = (S0_*G1 - S1_*G0)/det;
    k0 = (G0 - k1*S1_)/S0_;
}

void sixdof_deform::set_time(double t)
{
    t_ = t;

    // ramp r(t) = 0.5 (1 - cos(pi t/t_ramp)) and its rate
    r_ = 1.0;
    rdot_ = 0.0;

    if(tramp_>0.0 && t_<tramp_)
    {
        r_    = (t_>0.0) ? 0.5*(1.0 - std::cos(fish_pi*t_/tramp_)) : 0.0;
        rdot_ = (t_>0.0) ? 0.5*(fish_pi/tramp_)*std::sin(fish_pi*t_/tramp_) : 0.0;
    }

    c0_ = c1_ = d0_ = d1_ = 0.0;

    if(!recoil_ || xb_.empty())
    return;

    // recoil: no net lateral momentum and no angular momentum of the deformation about the CoG
    std::vector<double> hb(xb_.size()), htb(xb_.size());

    for(int b=0; b<int(xb_.size()); ++b)
    raw(xb_[b],hb[b],htb[b]);

    fit(hb,c0_,c1_);
    fit(htb,d0_,d1_);
}

double sixdof_deform::disp(double x) const
{
    double h, ht;
    raw(x,h,ht);

    return h - c0_ - c1_*x;
}

double sixdof_deform::vel(double x) const
{
    double h, ht;
    raw(x,h,ht);

    return ht - d0_ - d1_*x;
}
