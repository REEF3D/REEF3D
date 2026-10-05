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

#include"spectral_grid.h"
#include<cmath>

spectral_grid::spectral_grid(int nsig_, double fmin_, double fmax_, int ndir_)
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

    dtheta = 2.0*pi/double(ndir);

    theta.resize(ndir);
    costh.resize(ndir);
    sinth.resize(ndir);

    for(int m=0; m<ndir; ++m)
    {
    theta[m] = double(m)*dtheta;
    costh[m] = std::cos(theta[m]);
    sinth[m] = std::sin(theta[m]);
    }
}
