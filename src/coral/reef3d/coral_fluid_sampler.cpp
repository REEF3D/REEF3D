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

#include"coral/reef3d/coral_fluid_sampler.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"fdm_nhf.h"
#include<algorithm>

namespace
{
    // the column of (x, y) lies in this rank's cells
    bool own(lexer *p, double x, double y)
    {
        const int i = p->posc_i(x);
        if(i < 0 || i >= p->knox)
            return false;
        if(p->j_dir == 0)
            return true;
        const int j = p->posc_j(y);
        return j >= 0 && j < p->knoy;
    }

    template<class F>
    void sample_points(lexer *p, slice &WL, slice &bed, const std::vector<double> &xyz, std::vector<double> &out, F vel)
    {
        const int n = (int)xyz.size()/3;
        out.assign(5*n, 0.0);
        for(int q=0; q<n; ++q)
        {
            const double x = xyz[3*q], y = xyz[3*q+1], z = xyz[3*q+2];
            if(!own(p, x, y))
                continue;
            const double wl = p->ccslipol4(WL, x, y), zb = p->ccslipol4(bed, x, y);
            const double zs = zb + wl;
            double u = 0.0, v = 0.0, w = 0.0;
            if(z >= zb && wl > 0.0)
            {
                // above the surface: the value just below it; the interpolation of the sigma grid may give
                // nothing in its top cell (0): then 1 % and 3 % of the depth lower
                const double zc = std::min(z, zb + wl*(1.0 - 1.0e-6));
                vel(x, y, zc, u, v, w);
                for(double f : {0.01, 0.03})
                    if(u == 0.0 && v == 0.0 && w == 0.0 && zc > zb + wl*(1.0 - f))
                        vel(x, y, zb + wl*(1.0 - f), u, v, w);
            }
            out[5*q] = u;
            out[5*q+1] = v;
            out[5*q+2] = w;
            out[5*q+3] = zs;
            out[5*q+4] = 1.0;
        }
    }
}

void coral_sampler_fnpf::sample(const std::vector<double> &xyz, std::vector<double> &out)
{
    lexer *p = p_;
    fdm_fnpf *c = c_;
    sample_points(p, c->WL, c->bed, xyz, out, [&](double x, double y, double z, double &u, double &v, double &w)
    {
        u = p->ccipol7V(c->U, c->WL, c->bed, x, y, z);
        v = p->ccipol7V(c->V, c->WL, c->bed, x, y, z);
        w = p->ccipol7V(c->W, c->WL, c->bed, x, y, z);
    });
}

void coral_sampler_nhflow::sample(const std::vector<double> &xyz, std::vector<double> &out)
{
    lexer *p = p_;
    fdm_nhf *d = d_;
    sample_points(p, d->WL, d->bed, xyz, out, [&](double x, double y, double z, double &u, double &v, double &w)
    {
        u = p->ccipol4V(d->U, d->WL, d->bed, x, y, z);
        v = p->ccipol4V(d->V, d->WL, d->bed, x, y, z);
        w = p->ccipol4V(d->W, d->WL, d->bed, x, y, z);
    });
}
