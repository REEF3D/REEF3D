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

#ifndef WAVE_LIB_ZGROUP_H_
#define WAVE_LIB_ZGROUP_H_

#include<vector>
#include<map>
#include<cmath>
#include<algorithm>

/*--------------------------------------------------------------------
Cached-point evaluation of the irregular-wave libraries for the points of
one column and for points with the same phases (wave_lib_irregular_1st,
wave_lib_irregular_2nd_a / _b).

Points with the same projection on every wave direction (all points of a
row across the crests of a unidirectional sea, e.g. the cells j of a 3D
generation zone for waves along x) have the same phases: they form a
group. Per group and time:
- eta and the potential at one z are kept (the second point of the group
  takes them, bitwise equal);
- u, v, w as functions of z are interpolated in Chebyshev nodes over
  [-h, b] once the group asked for more values in each of the two
  previous steps than the interpolation needs nodes (the many levels of
  the cells of a group, or of one column of a fine CFD grid; two steps,
  as the initialisation asks the levels several times at t = 0).
  The nodes J follow from the largest vertical wave number kz of the
  profiles exp(+-K z): the interpolation error of exp(K z) over a half length r is at most
  2 (K r / 2)^J / J! of its largest value, J the smallest with 1e-14.
  b = max(z, zmax, 0) + h/10, zmax the highest level the group asked for
  in the step before (the levels of a column follow the surface); a value
  above b extends b by half the range (new nodes).
--------------------------------------------------------------------*/

struct wave_lib_zgroup
{
    struct grp
    {
        double tq=-1.0e300;             // time of the counters
        int nreq=0, nprev=0, npp=0;     // u, v, w values asked at that time and at the two times before
        double tt=-1.0e300;             // time of the table
        double tn=-1.0e300;             // time at which the group was found not worth a table
        double zhi=-1.0e300, phi=-1.0e300;  // highest level asked at tq and at the time before
        int mode=0, J=0;
        double a=0.0, b=0.0;
        std::vector<double> cu, cv, cw; // Chebyshev coefficients
        double te=-1.0e300, eta=0.0;    // eta at time te
        double tf=-1.0e300, zf=0.0, fi=0.0;  // potential at (tf, zf)
    };

    std::vector<int> gq;                // point -> group
    std::vector<grp> G;
    double kz=0.0, h=0.0;
    std::vector<double> zn, fu, fv, fw; // nodes, values at the nodes

    static constexpr double tol = 1.0e-14;
    static constexpr int Jmax = 64;

    // groups of the points x, y for the wave directions (cosb, sinb) of the M components;
    // kz_ the largest vertical wave number of the profiles, h_ the depth
    void build(const std::vector<double> &x, const std::vector<double> &y, int M, const double *cosb, const double *sinb,
               double kz_, double h_)
    {
        kz = kz_;
        h = h_;

        std::vector<std::pair<double,double>> dir;
        for(int n=0; n<M; ++n)
        {
            bool found=false;
            for(auto &d : dir)
            if(d.first==cosb[n] && d.second==sinb[n])
            {
            found=true;
            break;
            }

            if(!found)
            dir.push_back({cosb[n],sinb[n]});
        }

        const int N = int(x.size());
        std::map<std::vector<double>,int> id;
        std::vector<double> key(dir.size());
        gq.assign(N,0);

        for(int q=0; q<N; ++q)
        {
            for(size_t d=0; d<dir.size(); ++d)
            key[d] = dir[d].first*x[q] + dir[d].second*y[q];

            auto it = id.find(key);
            if(it==id.end())
            it = id.emplace(key,int(id.size())).first;
            gq[q] = it->second;
        }

        G.assign(id.size(),grp());
    }

    // the wave state changed (B 530): no value of the current time is kept
    void invalidate()
    {
        for(grp &g : G)
        {
        g.tt = g.te = g.tf = -1.0e300;
        }
    }
    
    // nodes for the range [a, b]; 0: too many
    int nodes(double a, double b) const
    {
        const double x = 0.25*kz*(b - a);
        double e = 2.0;

        for(int J=1; J<=Jmax; ++J)
        {
            e *= x/double(J);
            if(e<=tol && J>=4)
            return J;
        }
        return 0;
    }

    // eta of point q at time t: ev() once per group and time
    template<class EV>
    double eta(int q, double t, EV &&ev)
    {
        grp &g = G[gq[q]];

        if(g.te!=t)
        {
            g.eta = ev();
            g.te = t;
        }
        return g.eta;
    }

    // the potential of point q at (t, z): ev(z) once per group, time and z
    template<class EV>
    double fi(int q, double t, double z, EV &&ev)
    {
        grp &g = G[gq[q]];

        if(g.tf!=t || g.zf!=z)
        {
            g.fi = ev(z);
            g.tf = t;
            g.zf = z;
        }
        return g.fi;
    }

    // u, v, w of point q at (t, z) from the table of its group (ev(z, u, v, w) at the nodes);
    // false: not worth it (the caller evaluates directly)
    template<class EV>
    bool uvw(int q, double t, double z, int mode, EV &&ev, double &u, double &v, double &w)
    {
        grp &g = G[gq[q]];

        if(g.tq!=t)
        {
            g.npp = g.nprev;
            g.nprev = g.nreq;
            g.nreq = 0;
            g.phi = g.zhi;
            g.zhi = -1.0e300;
            g.tq = t;
        }
        ++g.nreq;
        g.zhi = fmax(g.zhi,z);

        if(g.tn==t)
        return false;

        if(g.tt!=t || g.mode!=mode || z<g.a || z>g.b)
        {
            const double a = fmin(-h,z);
            double b = fmax(fmax(z,g.phi),0.0) + 0.1*h;

            // a value above the table of this time: half its range more, so that the levels of a
            // column above the still water level (CFD) extend the table a few times only
            if(g.tt==t && g.mode==mode)
            b = fmax(b, z + 0.5*(g.b - g.a));

            const int J = nodes(a,b);

            if(J==0 || 3*std::min(g.nprev,g.npp)<4*J)
            {
            g.tn = t;
            return false;
            }

            table(g,a,b,J,mode,ev);
            g.tt = t;
        }

        // Clenshaw
        const double s = (2.0*z - g.a - g.b)/(g.b - g.a), s2 = 2.0*s;
        double u1=0.0, u2=0.0, v1=0.0, v2=0.0, w1=0.0, w2=0.0;

        for(int j=g.J-1; j>=1; --j)
        {
            const double un = s2*u1 - u2 + g.cu[j];
            const double vn = s2*v1 - v2 + g.cv[j];
            const double wn = s2*w1 - w2 + g.cw[j];
            u2 = u1; u1 = un;
            v2 = v1; v1 = vn;
            w2 = w1; w1 = wn;
        }

        u = s*u1 - u2 + 0.5*g.cu[0];
        v = s*v1 - v2 + 0.5*g.cv[0];
        w = s*w1 - w2 + 0.5*g.cw[0];
        return true;
    }

    template<class EV>
    void table(grp &g, double a, double b, int J, int mode, EV &&ev)
    {
        zn.resize(J); fu.resize(J); fv.resize(J); fw.resize(J);

        for(int j=0; j<J; ++j)
        {
            zn[j] = 0.5*(a + b) + 0.5*(b - a)*cos(M_PI*(j + 0.5)/J);
            ev(zn[j],fu[j],fv[j],fw[j]);
        }

        g.cu.assign(J,0.0); g.cv.assign(J,0.0); g.cw.assign(J,0.0);

        for(int k=0; k<J; ++k)
        {
            double su=0.0, sv=0.0, sw=0.0;

            for(int j=0; j<J; ++j)
            {
                const double c = cos(M_PI*k*(j + 0.5)/J);
                su += fu[j]*c;
                sv += fv[j]*c;
                sw += fw[j]*c;
            }

            g.cu[k] = 2.0*su/J;
            g.cv[k] = 2.0*sv/J;
            g.cw[k] = 2.0*sw/J;
        }

        g.a = a;
        g.b = b;
        g.J = J;
        g.mode = mode;
    }
};

#endif
