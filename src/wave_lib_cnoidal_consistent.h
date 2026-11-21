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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef WAVE_LIB_CNOIDAL_CONSISTENT_H_
#define WAVE_LIB_CNOIDAL_CONSISTENT_H_

// Self-consistent cnoidal wave (Fenton 1979, 1999), used by wave_lib_cnoidal_5th.
//
// - h (trough depth) is solved so that the phase-averaged surface equals the still water depth d
// - m is solved from the given period or length with the full relation L = 2 K h / alpha
// - c is set from the same velocity polynomial so that the Eulerian mean current is zero
// - w is the y-integral of -du/dx of the u polynomial, so the field is exactly divergence free
// - K, E and sn, cn, dn by AGM to machine precision
//
// Order 5 uses the Iwagaki (m -> 1) coefficients, valid for m >= 0.96; below that the third-order
// expansion with the full m dependence is used.

#include<array>
#include<cmath>
#include<stdexcept>
#include<algorithm>

namespace cnoidal_consistent
{
using Row = std::array<double,6>;     // coefficients of cn^(2k), k = 0..5
using Matrix = std::array<Row,5>;     // velocity: row j multiplies (Y/h)^(2j)

struct Model
{
    double d=0.0, H=0.0, g=9.81;
    double m=0.0, K=0.0, E=0.0, h=0.0, eps=0.0, alpha=0.0, delta=0.0, c=0.0, L=0.0, T=0.0;
    int order=3;
    Row moments{}, surface{};
    Matrix velocity{};

    // complete elliptic integrals K(m), E(m) by AGM; m is the parameter (k^2)
    static std::array<double,2> integrals(double m)
    {
        double a=1.0, b=std::sqrt(1.0-m), sum=0.5*m, weight=1.0;

        for(int i=0; i<50; ++i)
        {
            double cc = 0.5*(a-b);
            sum += weight*cc*cc;
            double next = 0.5*(a+b);
            b = std::sqrt(a*b);
            a = next;
            weight *= 2.0;

            if(std::fabs(cc)<2.0e-16*a)
            break;
        }

        double Kv = std::acos(-1.0)/(2.0*a);
        return {Kv, Kv*(1.0-sum)};
    }

    // Jacobi elliptic functions sn, cn, dn by AGM and back substitution
    static std::array<double,3> jacobi(double u, double m)
    {
        double mu[24]{1.0}, nu[24]{std::sqrt(1.0-m)}, a[24]{}, b[24]{};
        int n=0;

        while(std::fabs((mu[n]-nu[n])/(mu[n]+nu[n]))>1.0e-14 && n<22)
        {
            mu[n+1] = 0.5*(mu[n]+nu[n]);
            nu[n+1] = std::sqrt(mu[n]*nu[n]);
            ++n;
        }

        double su = std::sin(u*mu[n]);
        double cu = std::cos(u*mu[n]);

        if(std::fabs(su)<std::fabs(cu))
        {
            a[n] = mu[n]*su/cu;
            b[n] = 1.0;

            while(n>0)
            {
                --n;
                a[n] = b[n+1]*a[n+1];
                double r = a[n+1]*a[n+1]/mu[n+1];
                b[n] = (r+nu[n])/(r+mu[n]);
            }

            double dn = std::sqrt(1.0-m)/b[0];
            double cn = dn*std::copysign(1.0,cu)/std::sqrt(1.0+a[0]*a[0]);
            return {cn*a[0]/std::sqrt(1.0-m), cn, dn};
        }

        a[n] = mu[n]*cu/su;
        b[n] = 1.0;

        while(n>0)
        {
            --n;
            a[n] = b[n+1]*a[n+1];
            double r = a[n+1]*a[n+1]/mu[n+1];
            b[n] = (r+nu[n])/(r+mu[n]);
        }

        double sn = std::copysign(1.0,su)/std::sqrt(1.0+a[0]*a[0]);
        return {sn, a[0]*sn, b[0]};
    }

    static double poly(const Row &a, double q)
    {
        double s=0.0;
        for(int k=5; k>=0; --k)
        s = s*q + a[k];
        return s;
    }

    // phase average of sum a_k cn^(2k)
    double average(const Row &a) const
    {
        double s=0.0;
        for(int k=0; k<6; ++k)
        s += a[k]*moments[k];
        return s;
    }

    void set_surface(double hh)
    {
        h = hh;
        eps = H/h;
        surface = {};
        const double e2=eps*eps, e3=e2*eps;

        surface[1] = eps - 0.75*e2;
        surface[2] = 0.75*e2;

        if(order==3)
        {
            surface[1] += e3*(-61.0/(80.0*m) + 111.0/80.0);
            surface[2] += e3*( 61.0/(80.0*m) - 53.0/20.0);
            surface[3]  = e3*101.0/80.0;
        }
        else
        {
            const double e4=e3*eps, e5=e4*eps;
            surface[1] += e3*5.0/8.0 - e4*8209.0/6000.0 + e5*364671.0/196000.0;
            surface[2] += -e3*151.0/80.0 + e4*11641.0/3000.0 - e5*2920931.0/392000.0;
            surface[3]  = e3*101.0/80.0 - e4*112393.0/24000.0 + e5*2001361.0/156800.0;
            surface[4]  = e4*17367.0/8000.0 - e5*17906339.0/1568000.0;
            surface[5]  = e5*1331817.0/313600.0;
        }
    }

    void set_velocity()
    {
        velocity = {};
        velocity[0][0] = -1.0;

        auto add = [&](int power, int iy, Row a)
        {
            double f = std::pow(delta,power);
            for(int k=0; k<6; ++k)
            velocity[iy][k] += f*a[k];
        };

        if(order==3)
        {
            const double m2=m*m, m3=m2*m;
            add(1,0,{0.5-m, m});
            add(2,0,{-19.0/40.0 + 79.0/40.0*m - 79.0/40.0*m2, -1.5*m + 3.0*m2, -m2});
            add(2,1,{-0.75*m + 0.75*m2, 1.5*m - 3.0*m2, 2.25*m2});
            add(3,0,{55.0/112.0 - 3471.0/1120.0*m + 7113.0/1120.0*m2 - 2371.0/560.0*m3,
                     71.0/40.0*m - 339.0/40.0*m2 + 339.0/40.0*m3, 27.0/10.0*m2 - 27.0/5.0*m3, 6.0/5.0*m3});
            add(3,1,{9.0/8.0*m - 27.0/8.0*m2 + 9.0/4.0*m3, -9.0/4.0*m + 27.0/2.0*m2 - 27.0/2.0*m3,
                     -75.0/8.0*m2 + 75.0/4.0*m3, -15.0/2.0*m3});
            add(3,2,{-3.0/16.0*m + 9.0/16.0*m2 - 3.0/8.0*m3, 3.0/8.0*m - 51.0/16.0*m2 + 51.0/16.0*m3,
                     45.0/16.0*m2 - 45.0/8.0*m3, 45.0/16.0*m3});
        }
        else
        {
            add(1,0,{-0.5, 1.0});
            add(2,0,{-19.0/40.0, 1.5, -1.0});
            add(2,1,{0.0, -1.5, 2.25});
            add(3,0,{-55.0/112.0, 71.0/40.0, -27.0/10.0, 6.0/5.0});
            add(3,1,{0.0, -9.0/4.0, 75.0/8.0, -15.0/2.0});
            add(3,2,{0.0, 3.0/8.0, -45.0/16.0, 45.0/16.0});
            add(4,0,{-11813.0/22400.0, 53327.0/42000.0, -13109.0/3000.0, 1763.0/375.0, -197.0/125.0});
            add(4,1,{0.0, -213.0/80.0, 3231.0/160.0, -729.0/20.0, 189.0/10.0});
            add(4,2,{0.0, 9.0/16.0, -327.0/32.0, 915.0/32.0, -315.0/16.0});
            add(4,3,{0.0, -3.0/80.0, 189.0/160.0, -63.0/16.0, 189.0/64.0});
            add(5,0,{-57159.0/98560.0, -144821.0/156800.0, -1131733.0/294000.0, 757991.0/73500.0, -298481.0/36750.0, 13438.0/6125.0});
            add(5,1,{0.0, -53327.0/28000.0, 1628189.0/56000.0, -192481.0/2000.0, 11187.0/100.0, -5319.0/125.0});
            add(5,2,{0.0, 213.0/320.0, -13563.0/640.0, 68643.0/640.0, -5481.0/32.0, 1701.0/20.0});
            add(5,3,{0.0, -9.0/160.0, 267.0/64.0, -987.0/32.0, 7875.0/128.0, -567.0/16.0});
            add(5,4,{0.0, 9.0/4480.0, -459.0/1792.0, 567.0/256.0, -1215.0/256.0, 729.0/256.0});
        }

        // zero Eulerian mean current, from the same velocity polynomial (bed level)
        c = -std::sqrt(g*h)*average(velocity[0]);
    }

    void at_m(double mm, int requested_order)
    {
        m = mm;
        order = requested_order;
        auto ke = integrals(m);
        K = ke[0];
        E = ke[1];

        // <cn^(2k)> over one period
        moments[0] = 1.0;
        moments[1] = (E/K + m - 1.0)/m;
        for(int j=0; j<4; ++j)
        moments[j+2] = (2.0*j+2.0)/(2.0*j+3.0)*(2.0-1.0/m)*moments[j+1] + (2.0*j+1.0)/(2.0*j+3.0)*(1.0/m-1.0)*moments[j];

        // trough depth h from mean surface = d
        double lo = std::max(d-H,0.1*d), hi = d;
        set_surface(lo);
        double flo = h*(1.0+average(surface)) - d;
        set_surface(hi);
        double fhi = h*(1.0+average(surface)) - d;

        if(!(flo<=0.0 && fhi>=0.0))
        throw std::runtime_error("cnoidal: mean-depth root not bracketed, wave too high for this expansion");

        for(int i=0; i<70; ++i)
        {
            double mid = 0.5*(lo+hi);
            set_surface(mid);
            if(h*(1.0+average(surface))>d)
            hi = mid;
            else
            lo = mid;
        }
        set_surface(0.5*(lo+hi));

        if(order==3)
        {
            double r = eps/m;
            alpha = std::sqrt(0.75*r)*(1.0 + r*(0.25 - 7.0*m/8.0) + r*r*(1.0/32.0 - 11.0*m/32.0 + 111.0*m*m/128.0));
        }
        else
        alpha = std::sqrt(0.75*eps)*(1.0 - 5.0*eps/8.0 + 71.0*eps*eps/128.0 - 100627.0*std::pow(eps,3.0)/179200.0
                                     + 16259737.0*std::pow(eps,4.0)/28672000.0);

        delta = 4.0*alpha*alpha/3.0;
        set_velocity();
        L = 2.0*K*h/alpha;
        T = L/c;

        if(!(alpha>0.0 && c>0.0 && std::isfinite(T)))
        throw std::runtime_error("cnoidal: invalid expansion (alpha, c or T)");
    }

    // target is T (period=true) or L; m is bisected in [0.5, 1-1e-12]
    static Model solve(double d, double H, double g, double target, bool period)
    {
        if(!(d>0.0 && H>0.0 && H<0.9*d && g>0.0 && target>0.0))
        throw std::runtime_error("cnoidal: invalid input, need 0 < H < 0.9 d");

        auto branch = [&](int n)
        {
            Model a;
            a.d=d; a.H=H; a.g=g;
            double lo=0.5, hi=1.0-1.0e-12;
            a.at_m(lo,n);
            double flo = (period?a.T:a.L) - target;
            a.at_m(hi,n);
            double fhi = (period?a.T:a.L) - target;

            if(flo*fhi>0.0)
            throw std::runtime_error("cnoidal: no modulus in m = [0.5,1) for this H, T/L (Ursell number too low: use Stokes 5th, B 92 5 or 11)");

            for(int j=0; j<65; ++j)
            {
                double mid = 0.5*(lo+hi);
                a.at_m(mid,n);
                if(((period?a.T:a.L) - target)*flo>0.0)
                lo = mid;
                else
                hi = mid;
            }
            a.at_m(0.5*(lo+hi),n);
            return a;
        };

        Model a = branch(3);

        if(a.m>=0.96)
        {
            try
            {
                Model b = branch(5);
                if(b.m>=0.96)
                return b;
            }
            catch(const std::exception&)
            {
            }
        }

        return a;
    }

    std::array<double,3> phase(double x, double t, double shift) const
    {
        // reduce the argument for long runs
        double u = std::remainder(2.0*K*(x/L - t/T) + shift, 4.0*K);
        return jacobi(u,m);
    }

    double eta(double x, double t, double shift) const
    {
        auto j = phase(x,t,shift);
        return h*(1.0 + poly(surface,j[1]*j[1])) - d;
    }

    double u(double x, double z, double t, double shift) const
    {
        auto j = phase(x,t,shift);
        double q=j[1]*j[1], y=(z+d)/h, v=0.0;

        for(int k=order-1; k>=0; --k)
        v = v*y*y + poly(velocity[k],q);

        return c + std::sqrt(g*h)*v;
    }

    double w(double x, double z, double t, double shift) const
    {
        auto j = phase(x,t,shift);
        double q=j[1]*j[1], y=(z+d)/h, v=0.0;

        for(int iy=0; iy<order; ++iy)
        {
            double dq=0.0;
            for(int k=5; k>=1; --k)
            dq = dq*q + k*velocity[iy][k];

            v += std::pow(y,2*iy+1)/(2*iy+1)*dq;
        }

        return 2.0*alpha*std::sqrt(g*h)*j[0]*j[1]*j[2]*v;
    }
};
}

#endif
