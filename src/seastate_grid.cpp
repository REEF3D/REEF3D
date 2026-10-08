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

#include"seastate_grid.h"
#include<algorithm>
#include<cmath>

seastate_grid::seastate_grid(int nsig_, double fmin_, double fmax_, int ndir_)
    : nsig(nsig_), ndir(ndir_), nbin(0), fmin(fmin_), fmax(fmax_), ratio(1.0), dtheta(0.0)
{
    const double pi = 3.14159265358979323846;

    if(nsig<2)
    error = "A 701: at least 2 frequencies are needed";
    else if(!(fmin>0.0) || !(fmax>fmin))
    error = "A 702: frequency range needs 0 < fmin < fmax";
    else if(ndir<4)
    error = "A 703: at least 4 directions are needed";

    if(!error.empty())
    {
    nsig = ndir = 0;
    return;
    }

    nbin  = nsig*ndir;
    ratio = std::pow(fmax/fmin, 1.0/double(nsig-1));

    f.resize(nsig);
    sig.resize(nsig);
    dsig.resize(nsig);

    for(int l=0; l<nsig; ++l)
    {
    f[l]   = fmin*std::pow(ratio,double(l));
    sig[l] = 2.0*pi*f[l];
    }
    f[nsig-1]   = fmax;
    sig[nsig-1] = 2.0*pi*fmax;

    // bin widths: faces at the geometric midpoints, end faces mirrored,
    // so that dsig_l = sig_l * (sqrt(r) - 1/sqrt(r)) for every bin
    const double sr = std::sqrt(ratio);

    for(int l=0; l<nsig; ++l)
    {
    const double lo = (l==0)      ? sig[l]/sr : std::sqrt(sig[l-1]*sig[l]);
    const double hi = (l==nsig-1) ? sig[l]*sr : std::sqrt(sig[l]*sig[l+1]);
    dsig[l] = hi - lo;
    }

    directions();
}

seastate_grid::seastate_grid(double frep, int ndir_)
    : nsig(1), ndir(ndir_), nbin(0), fmin(frep), fmax(frep), ratio(1.0), dtheta(0.0)
{
    const double pi = 3.14159265358979323846;

    if(!(frep>0.0))
    error = "A 771: the representative frequency must be positive";
    else if(ndir<4)
    error = "A 703: at least 4 directions are needed";

    if(!error.empty())
    {
    nsig = ndir = 0;
    return;
    }

    nbin = ndir;
    f.assign(1,frep);
    sig.assign(1,2.0*pi*frep);
    dsig.assign(1,1.0);

    directions();
}

void seastate_grid::directions()
{
    const double pi = 3.14159265358979323846;

    dtheta = 2.0*pi/double(ndir);

    theta.resize(ndir);
    costh.resize(ndir);
    sinth.resize(ndir);
    costhf.resize(ndir);
    sinthf.resize(ndir);
    quad.resize(ndir);

    for(int m=0; m<ndir; ++m)
    {
    theta[m]  = double(m)*dtheta;
    costh[m]  = std::cos(theta[m]);
    sinth[m]  = std::sin(theta[m]);
    costhf[m] = std::cos(theta[m]+0.5*dtheta);
    sinthf[m] = std::sin(theta[m]+0.5*dtheta);
    quad[m]   = std::min(int(4.0*double(m)/double(ndir)),3);
    }

    dth.assign(ndir,dtheta);
    wth.assign(ndir,1.0);
}

void seastate_grid::sector(double th1, double th2, int k)
{
    const double pi = 3.14159265358979323846;

    if(k<=1 || ndir==0)
    return;

    auto wrap = [&](double a)
    {
        a = std::fmod(a,2.0*pi);
        return (a<0.0) ? a+2.0*pi : a;
    };

    const double width = wrap(th2-th1);

    if(!(width>0.0))
    {
    error = "A 715: the direction sector is empty";
    return;
    }

    // bins (centre, lower face, upper face) of the uniform grid, the ones in the sector divided
    std::vector<double> c, lo, hi;

    for(int m=0; m<ndir; ++m)
    {
    const double t = theta[m];

        if(wrap(t-th1)<=width)
        for(int n=0; n<k; ++n)
        {
        const double a = t - 0.5*dtheta + dtheta*double(n)/double(k);
        const double b = t - 0.5*dtheta + dtheta*double(n+1)/double(k);
        c.push_back(wrap(0.5*(a+b)));
        lo.push_back(a);
        hi.push_back(b);
        }
        else
        {
        c.push_back(t);
        lo.push_back(t-0.5*dtheta);
        hi.push_back(t+0.5*dtheta);
        }
    }

    std::vector<int> o(c.size());
    for(size_t n=0; n<o.size(); ++n)
    o[n] = int(n);
    std::sort(o.begin(),o.end(),[&](int a, int b){return c[a]<c[b];});

    ndir = int(c.size());
    nbin = nsig*ndir;
    uniform = false;

    theta.resize(ndir);
    costh.resize(ndir);
    sinth.resize(ndir);
    costhf.resize(ndir);
    sinthf.resize(ndir);
    quad.resize(ndir);
    dth.resize(ndir);
    wth.resize(ndir);

    for(int m=0; m<ndir; ++m)
    {
    const int n = o[m];
    theta[m]  = c[n];
    costh[m]  = std::cos(c[n]);
    sinth[m]  = std::sin(c[n]);
    costhf[m] = std::cos(hi[n]);
    sinthf[m] = std::sin(hi[n]);
    dth[m]    = hi[n]-lo[n];
    wth[m]    = dth[m]/dtheta;
    quad[m]   = std::min(int(c[n]/(0.5*pi)),3);
    }
}

int seastate_grid::direction_bin(double th) const
{
    const double pi = 3.14159265358979323846;

    double t = std::fmod(th,2.0*pi);
    if(t<0.0)
    t += 2.0*pi;

    if(uniform)
    return int(std::lround(t/dtheta))%ndir;

    // the bin whose faces enclose t (upper face of bin m: theta_m + dth_m/2)
    for(int m=0; m<ndir; ++m)
    {
    double d = t - theta[m];
    if(d>pi)  d -= 2.0*pi;
    if(d<-pi) d += 2.0*pi;
        if(std::fabs(d)<=0.5*dth[m])
        return m;
    }
    return 0;
}
