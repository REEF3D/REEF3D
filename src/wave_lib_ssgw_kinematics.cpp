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
Interior kinematics for SSGW (needed by NHFLOW and CFD wave generation).

The SSGW solution is a permanent-form wave travelling at c_e, so the
fixed-frame velocity field only depends on xi = x - c_e t and z.
It is evaluated once with the Cauchy integral (computeVelocityField)
on a periodic xi-grid times a sigma-grid s = (z+d)/(eta(xi)+d) in [0,1]
and bilinearly interpolated afterwards. The top level s=1 uses the
exact surface velocities from the conformal solution, which avoids the
near-singular Cauchy kernel on the surface.

Frame: c_e (zero Eulerian mean current below the trough), consistent
with wave_eta/wave_fi, which propagate with c_e as well.
--------------------------------------------------------------------*/

#include"wave_lib_ssgw.h"
#include<algorithm>

double wave_lib_ssgw::waveFrameX(double x, double t)
{
    // periodic position in the wave frame, xi in [-L/2, L/2)
    double Lw = ParameterValue.waveLength;
    double xi = modulo(x - ParameterValue.phaseVelocity*t, Lw);

    if(xi>=0.5*Lw)
    xi -= Lw;

    return xi;
}

double wave_lib_ssgw::surfaceInterp(const std::vector<double> &f, double xi)
{
    // periodic linear interpolation on the (non-uniform) surface nodes xs
    const int n = xs.size();
    const double Lw = ParameterValue.waveLength;

    auto it = std::upper_bound(xs.begin(),xs.end(),xi);
    int ir = std::distance(xs.begin(),it);
    int il = ir-1;

    double x0,x1,f0,f1;

    if(il<0)
    {
        x0 = xs[n-1]-Lw;  f0 = f[n-1];
        x1 = xs[0];       f1 = f[0];
    }
    else if(ir>=n)
    {
        x0 = xs[n-1];     f0 = f[n-1];
        x1 = xs[0]+Lw;    f1 = f[0];
    }
    else
    {
        x0 = xs[il];      f0 = f[il];
        x1 = xs[ir];      f1 = f[ir];
    }

    return f0 + (xi-x0)/(x1-x0)*(f1-f0);
}

void wave_lib_ssgw::buildVelocityTable()
{
    const double Lw = ParameterValue.waveLength;
    const double dw = ParameterValue.waterDepth;

    Nxt = 512;
    Nst = 41;
    dxt = Lw/double(Nxt);

    etat.resize(Nxt);
    Ut.resize(Nxt*Nst);
    Wt.resize(Nxt*Nst);

    std::vector<double> xq, zq, uq, wq;
    xq.reserve(Nxt*(Nst-1));
    zq.reserve(Nxt*(Nst-1));

    for(int i=0; i<Nxt; ++i)
    {
        double xi = -0.5*Lw + i*dxt;
        etat[i] = surfaceInterp(ys,xi);

        for(int k=0; k<Nst-1; ++k)
        {
            double s = double(k)/double(Nst-1);
            xq.push_back(xi);
            zq.push_back(-dw + s*(etat[i]+dw));
        }
    }

    uq.resize(xq.size());
    wq.resize(xq.size());
    computeVelocityField(xq,zq,uq,wq);

    int q=0;
    for(int i=0; i<Nxt; ++i)
    {
        double xi = -0.5*Lw + i*dxt;

        for(int k=0; k<Nst-1; ++k)
        {
            Ut[i*Nst+k] = uq[q];
            Wt[i*Nst+k] = wq[q];
            ++q;
        }

        // free surface: exact conformal-map velocities (fixed frame)
        Ut[i*Nst+Nst-1] = surfaceInterp(us,xi);
        Wt[i*Nst+Nst-1] = surfaceInterp(vs,xi);
    }

    kinematicsReady = true;
}

double wave_lib_ssgw::tableInterp(const std::vector<double> &tab, double xi, double z)
{
    const double Lw = ParameterValue.waveLength;
    const double dw = ParameterValue.waterDepth;

    double fx = (xi + 0.5*Lw)/dxt;
    int i0 = int(std::floor(fx));
    double wx = fx - i0;
    i0 = ((i0%Nxt)+Nxt)%Nxt;
    int i1 = (i0+1)%Nxt;

    auto column = [&](int i)
    {
        // sigma coordinate of this column; points above the SSGW surface
        // (possible in the relaxation zone) take the surface value
        double s = (z+dw)/(etat[i]+dw);
        s = std::max(0.0,std::min(1.0,s));

        double fk = s*(Nst-1);
        int k0 = std::min(int(fk),Nst-2);
        double wk_ = fk - k0;

        return (1.0-wk_)*tab[i*Nst+k0] + wk_*tab[i*Nst+k0+1];
    };

    return (1.0-wx)*column(i0) + wx*column(i1);
}
