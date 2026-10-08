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
Author: Ahmet Soydan
--------------------------------------------------------------------*/

#include"acoustics_fwh_kernel.h"
#include<cmath>
#include<utility>

fwh_permeable::fwh_permeable(double c0_, double rho0_, double dto_) : c0(c0_), rho0(rho0_), dto(dto_),
                                                                      nlev(0), k0(0), tau_first_mid(0.0), tau_last_mid(0.0), spread_any(false)
{
    tau[0]=tau[1]=tau[2]=0.0;
    Um[0]=Um[1]=Um[2]=0.0;
}

void fwh_permeable::set_medium_velocity(const double *U)
{
    Um[0]=U[0];
    Um[1]=U[1];
    Um[2]=U[2];
}

void fwh_permeable::add_panel(const fwh_panel &pan)
{
    panel.push_back(pan);
}

int fwh_permeable::add_observer(const double *x)
{
    observer ob;
    ob.pt.push_back({{x[0],x[1],x[2]},1.0});
    obs.push_back(ob);

    return int(obs.size())-1;
}

void fwh_permeable::add_image(int o, const double *x, double weight)
{
    obs[size_t(o)].pt.push_back({{x[0],x[1],x[2]},weight});
}

void fwh_permeable::step(double t, const double *p, const double *u)
{
    const size_t np = panel.size();

    if(nlev==0)
    {
        k0 = long(floor(t/dto));

        for(int l=0; l<3; ++l)
        lev[l].assign(4*np,0.0);
    }

    // shift the time levels, the newest goes to lev[2]
    std::swap(lev[0],lev[1]);
    std::swap(lev[1],lev[2]);
    tau[0]=tau[1];
    tau[1]=tau[2];
    tau[2]=t;

    for(size_t i=0; i<np; ++i)
    {
        lev[2][4*i]   = p[i];
        lev[2][4*i+1] = u[3*i];
        lev[2][4*i+2] = u[3*i+1];
        lev[2][4*i+3] = u[3*i+2];
    }

    if(nlev<3)
    ++nlev;

    if(nlev==3)
    {
        integrand();

        if(!spread_any)
        tau_first_mid = tau[1];

        tau_last_mid = tau[1];
        spread_any = true;
    }
}

void fwh_permeable::integrand()
{
    const size_t np = panel.size();
    const double c2 = c0*c0;
    const double fac = 1.0/(4.0*M_PI);

    // central difference at tau[1] on non-uniform steps
    const double ha = tau[1]-tau[0];
    const double hb = tau[2]-tau[1];
    const double ca = -hb/(ha*(ha+hb));
    const double cm = (hb-ha)/(ha*hb);
    const double cb = ha/(hb*(ha+hb));

    for(size_t i=0; i<np; ++i)
    {
        const fwh_panel &pan = panel[i];
        const double *A = &lev[0][4*i];
        const double *M = &lev[1][4*i];
        const double *B = &lev[2][4*i];

        double d[4];
        for(int q=0; q<4; ++q)
        d[q] = ca*A[q] + cm*M[q] + cb*B[q];

        const double pm = M[0];
        const double pd = d[0];
        const double *um = &M[1];
        const double *ud = &d[1];

        const double rho  = rho0 + pm/c2;
        const double rhod = pd/c2;

        const double un  = um[0]*pan.n[0] + um[1]*pan.n[1] + um[2]*pan.n[2];
        const double und = ud[0]*pan.n[0] + ud[1]*pan.n[1] + ud[2]*pan.n[2];

        const double Qd = rhod*un + rho*und;
        const double Q  = rho*un - rho0*(Um[0]*pan.n[0] + Um[1]*pan.n[1] + Um[2]*pan.n[2]);

        double L[3], Ld[3];
        for(int q=0; q<3; ++q)
        {
            L[q]  = pm*pan.n[q] + rho*(um[q]-Um[q])*un;
            Ld[q] = pd*pan.n[q] + rhod*(um[q]-Um[q])*un + rho*(ud[q]*un + (um[q]-Um[q])*und);
        }

        for(observer &ob : obs)
        for(const point &pt : ob.pt)
        {
            const double rx = pt.x[0]-pan.x[0];
            const double ry = pt.x[1]-pan.x[1];
            const double rz = pt.x[2]-pan.x[2];
            const double r  = sqrt(rx*rx + ry*ry + rz*rz);

            const double Lr  = (L[0]*rx  + L[1]*ry  + L[2]*rz)/r;
            const double Ldr = (Ld[0]*rx + Ld[1]*ry + Ld[2]*rz)/r;
            const double vr  = -(Um[0]*rx + Um[1]*ry + Um[2]*rz)/r;

            const double f = pt.w*fac*pan.dS*(Qd/r + Ldr/(c0*r) + Lr/(r*r) + Q*vr/(r*r));

            spread(ob,f,r/c0,tau[0],tau[1],tau[2]);
        }
    }
}

void fwh_permeable::spread(observer &ob, double f, double d, double ta, double tm, double tb)
{
    // hat function of the sample at tm, support (ta,tb) in source time, shifted by the delay d;
    // half-open parts (ta,tm] and (tm,tb) so that each observer time gets the weights of exactly
    // the two samples around it
    const long kmin = long(floor((ta+d)/dto));
    const long kmax = long(ceil((tb+d)/dto));

    for(long k=kmin; k<=kmax; ++k)
    {
        const double s = double(k)*dto - d;

        if(s<=ta || s>=tb)
        continue;

        const double w = s<=tm ? (s-ta)/(tm-ta) : (tb-s)/(tb-tm);

        const size_t m = size_t(k-k0);

        if(m>=ob.sig.size())
        ob.sig.resize(m+1,0.0);

        ob.sig[m] += f*w;
    }
}

void fwh_permeable::delays(int o, double &dmin, double &dmax) const
{
    dmin = 1.0e300;
    dmax = -1.0e300;

    for(const fwh_panel &pan : panel)
    for(const point &pt : obs[size_t(o)].pt)
    {
        const double rx = pt.x[0]-pan.x[0];
        const double ry = pt.x[1]-pan.x[1];
        const double rz = pt.x[2]-pan.x[2];
        const double d  = sqrt(rx*rx + ry*ry + rz*rz)/c0;

        dmin = d<dmin ? d : dmin;
        dmax = d>dmax ? d : dmax;
    }
}

bool fwh_permeable::complete_range(double dmin, double dmax, double &t0, double &t1) const
{
    // every panel has spread the samples tau_first_mid..tau_last_mid, so the observer time t is
    // complete if t-d lies in that interval for all delays d
    t0 = tau_first_mid + dmax;
    t1 = tau_last_mid + dmin;

    return spread_any && t1>=t0;
}
