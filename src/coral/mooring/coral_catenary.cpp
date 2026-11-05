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

#include"coral/mooring/coral_catenary.h"
#include<algorithm>
#include<cmath>
#include<utility>
#include<vector>

namespace coral
{

double catenary_seabed_length(const catenary_line &c, double, double VF)
{
    return std::max(0.0, c.L - VF/c.w);
}

void catenary_ends(const catenary_line &c, double HF, double VF, double &xF, double &zF)
{
    const double L = c.L, w = c.w, EA = c.EA, CB = c.CB;
    const double LB = L - VF/w;
    if(HF <= 0.0)
    {
        // no horizontal force: the line hangs vertically from the fairlead, the rest lies on the seabed
        xF = std::max(LB, 0.0);
        zF = VF/w + VF*VF/(2.0*EA*w);
        return;
    }
    if(LB <= 0.0 || CB < 0.0)
    {
        // hanging free (CB < 0: no seabed)
        const double a = VF/HF, b = (VF - w*L)/HF;
        xF = HF/w*(std::asinh(a) - std::asinh(b)) + HF*L/EA;
        zF = HF/w*(std::sqrt(1.0 + a*a) - std::sqrt(1.0 + b*b)) + (VF*L - 0.5*w*L*L)/EA;
    }
    else
    {
        const double a = VF/HF;
        const double lf = LB - HF/(CB*w);    // length of the seabed portion with friction
        xF = LB + HF/w*std::asinh(a) + HF*L/EA;
        if(CB > 0.0)
            xF += CB*w/(2.0*EA)*(-LB*LB + (lf > 0.0 ? lf*lf : 0.0));
        zF = HF/w*(std::sqrt(1.0 + a*a) - 1.0) + VF*VF/(2.0*EA*w);
    }
}

// V/EA of a vertical hanging part of height z: z = V/w + V^2/(2 EA w)
static double w_root(const catenary_line &c, double z)
{
    // V^2/(2 EA w) + V/w - z = 0  ->  V = EA (sqrt(1 + 2 w z/EA) - 1)
    return std::sqrt(1.0 + 2.0*c.w*z/c.EA) - 1.0;
}

// Newton in (ln HF, VF) with a backtracking line search on the residual; false if not converged
static bool newton(const catenary_line &c, double xF, double zF, double &HF, double &VF)
{
    const double scale = c.L + std::fabs(xF) + std::fabs(zF);
    auto res = [&](double lh, double V, double &rx, double &rz)
    {
        double x, z;
        catenary_ends(c, std::exp(lh), V, x, z);
        rx = x - xF;
        rz = z - zF;
        return std::sqrt(rx*rx + rz*rz);
    };
    double lh = std::log(HF), V = VF;
    double rx, rz;
    double r = res(lh, V, rx, rz);
    for(int it=0; it<200; ++it)
    {
        if(!std::isfinite(r))
            return false;
        if(r < 1.0e-11*scale)
        {
            HF = std::exp(lh);
            VF = V;
            return true;
        }
        // Jacobian by central differences
        const double dl = 1.0e-7, dV = 1.0e-7*std::max(std::fabs(V), 1.0e-3*c.w*c.L);
        double a1, b1, a2, b2, a3, b3, a4, b4;
        res(lh + dl, V, a1, b1);
        res(lh - dl, V, a2, b2);
        res(lh, V + dV, a3, b3);
        res(lh, V - dV, a4, b4);
        const double ja = (a1 - a2)/(2.0*dl), jb = (a3 - a4)/(2.0*dV);
        const double jc = (b1 - b2)/(2.0*dl), jd = (b3 - b4)/(2.0*dV);
        const double det = ja*jd - jb*jc;
        if(det==0.0 || !std::isfinite(det))
            return false;
        const double sl = -( jd*rx - jb*rz)/det;
        const double sV = -(-jc*rx + ja*rz)/det;
        // backtracking: at most a factor e^3 in HF per iteration, residual must fall
        double f = std::min(1.0, 3.0/std::max(std::fabs(sl), 1.0e-300));
        bool moved = false;
        for(int k=0; k<40; ++k)
        {
            double tx, tz;
            const double rn = res(lh + f*sl, V + f*sV, tx, tz);
            if(std::isfinite(rn) && rn < (1.0 - 1.0e-4*f)*r)
            {
                lh += f*sl;
                V += f*sV;
                r = rn;
                rx = tx;
                rz = tz;
                moved = true;
                break;
            }
            f *= 0.5;
        }
        if(!moved)
            return false;
    }
    return false;
}

bool catenary_solve(const catenary_line &c, double xF, double zF, double &HF, double &VF)
{
    xF = std::max(xF, 1.0e-6*c.L);
    // slack line on the seabed: it reaches the fairlead hanging vertically with length to spare; no
    // horizontal force, the vertical part carries its own weight (z = V/w + V^2/(2 EA w))
    if(c.CB >= 0.0 && zF >= 0.0 && xF + zF <= c.L)
    {
        HF = 0.0;
        VF = c.EA*w_root(c, zF);
        return true;
    }
    // start values: the given ones, Peyrot & Goulois, and a range for slack lines on the seabed
    std::vector<std::pair<double,double>> starts;
    if(HF > 0.0 && std::isfinite(HF) && std::isfinite(VF))
        starts.push_back({HF, VF});
    {
        double lam;
        if(c.L <= std::sqrt(xF*xF + zF*zF))
            lam = 0.2;
        else
            lam = std::sqrt(3.0*((c.L*c.L - zF*zF)/(xF*xF) - 1.0));
        starts.push_back({std::max(std::fabs(0.5*c.w*xF/lam), 1.0e-6*c.w*c.L), 0.5*c.w*(zF/std::tanh(lam) + c.L)});
    }
    for(double f : {1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0})
        starts.push_back({f*c.w*c.L, c.w*std::max(zF, 1.0e-3*c.L)});
    for(const auto &s0 : starts)
    {
        double H = s0.first, V = s0.second;
        if(newton(c, xF, zF, H, V))
        {
            HF = H;
            VF = V;
            return true;
        }
    }
    return false;
}

void catenary_profile(const catenary_line &c, double HF, double VF, double s, double &x, double &z)
{
    const double L = c.L, w = c.w, EA = c.EA, CB = c.CB;
    if(HF <= 0.0)
    {
        // slack: on the seabed up to LB, then vertical (stretch of the vertical part neglected)
        const double LB = std::max(0.0, L - VF/w);
        x = std::min(s, LB);
        z = std::max(0.0, s - LB);
        return;
    }
    const double LB = (CB < 0.0) ? 0.0 : std::max(0.0, L - VF/w);
    if(LB <= 0.0)
    {
        const double VA = VF - w*L;
        x = HF/w*(std::asinh((VA + w*s)/HF) - std::asinh(VA/HF)) + HF*s/EA;
        z = HF/w*(std::sqrt(1.0 + std::pow((VA + w*s)/HF, 2)) - std::sqrt(1.0 + std::pow(VA/HF, 2))) + (VA*s + 0.5*w*s*s)/EA;
        return;
    }
    if(s <= LB)
    {
        // on the seabed: tension HF - CB w (LB - s), not below 0
        const double lf = LB - HF/(CB > 0.0 ? CB*w : 1.0e-300);
        double stretch;
        if(CB <= 0.0)
            stretch = HF*s/EA;
        else if(s <= lf)
            stretch = 0.0;
        else
        {
            const double s0 = std::max(lf, 0.0);
            // int_s0^s (HF - CB w (LB - u)) du / EA
            stretch = (HF*(s - s0) - CB*w*(LB*(s - s0) - 0.5*(s*s - s0*s0)))/EA;
        }
        x = s + stretch;
        z = 0.0;
        return;
    }
    double xb, zb;
    catenary_profile(c, HF, VF, LB, xb, zb);
    const double u = s - LB;
    x = xb + HF/w*std::asinh(w*u/HF) + HF*u/EA;
    z = HF/w*(std::sqrt(1.0 + std::pow(w*u/HF, 2)) - 1.0) + 0.5*w*u*u/EA;
}

void catenary_anchor(const catenary_line &c, double HF, double VF, double &HA, double &VA)
{
    const double LB = (c.CB < 0.0) ? 0.0 : std::max(0.0, c.L - VF/c.w);
    if(LB <= 0.0)
    {
        HA = HF;
        VA = VF - c.w*c.L;
    }
    else
    {
        HA = std::max(HF - c.CB*c.w*LB, 0.0);
        VA = 0.0;
    }
}

}
