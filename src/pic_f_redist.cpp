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

#include"pic_f.h"
#include"lexer.h"
#include"fdm.h"
#include<algorithm>
#include<cmath>
#include<vector>

// The signed distance in the band of the surface (F 81 1, 2), from the zero of phi = h d(c)
// (pic_f::phi_build). h d(c) is the distance of a flat surface only: where the particles near the
// surface are not on their lattice, its slope across the surface is not 1 and varies along it,
// and the density of the smoothed interface (F 45) with it. In still water the varying density on
// the air side of the interface drove air currents, these the water at the surface and the
// particles, and the particles the density again: the velocity grew exponentially from the side
// walls (0.09 m/s after 5 s); 20 iterations of the reinitialisation instead of the 3 of F 44
// removed that, at four times the cost of phi.
// Here, directly:
//  - the cells at the surface (a neighbour with the other sign) get the distance from the
//    crossings of the zero along the axes, 1/d^2 = sum over the axes of 1/s^2, with s the distance
//    of the nearest crossing along an axis, from the linear interpolation of phi to the
//    neighbour (exact for a plane), and their closest point on the surface, x - phi n, with n the
//    normal from the central differences of phi;
//  - the cells of the band (0.01 < c < 0.99) the distance of the nearest of these closest points
//    within 2 cells, with their sign.
// The zero of phi stays where the particles put it. The surface cells are found on the cells of
// this subdomain and of the margin (the values of the neighbours), so the result does not depend
// on the decomposition. The reinitialisation of the level set afterwards (F 44) continues the
// distance beyond the band.

namespace
{
inline int fid(lexer *p, int i, int j, int k)
{
    return (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin;
}
}

void pic_f::redistance(lexer *p, fdm *a, const std::vector<char> &inband)
{
    const int S = 2;
    const int fm = fmarg;
    const bool d3 = p->j_dir==1;
    const int size = p->imax*p->jmax*p->kmax;
    const double *V = a->phi.data();

    // surface cells: closest point and distance, over the cells of the subdomain and S around
    std::vector<char> iface(size,0);
    std::vector<double> cp(3*size,0.0);
    std::vector<double> dist(size,0.0);

    const int ilo = -S, ihi = p->knox+S-1;
    const int jlo = d3 ? -S : 0, jhi = d3 ? p->knoy+S-1 : 0;
    const int klo = -S, khi = p->knoz+S-1;

    auto fluid = [&](int i, int j, int k)
    {
        return i>=-fm && i<p->knox+fm && k>=-fm && k<p->knoz+fm && (!d3 || (j>=-fm && j<p->knoy+fm))
            && p->flag4[fid(p,i,j,k)]>0;
    };

    for(int i=ilo; i<=ihi; ++i)
    for(int j=jlo; j<=jhi; ++j)
    for(int k=klo; k<=khi; ++k)
    {
        if(!fluid(i,j,k))
        continue;

        const int f = fid(p,i,j,k);
        const double v0 = V[f];
        const bool w0 = v0>=0.0;
        const double X[3] = {p->XP[i+marge], p->YP[j+marge], p->ZP[k+marge]};
        double inv2=0.0;
        double g[3] = {0.0,0.0,0.0};

        for(int ax=0; ax<3; ++ax)
        {
            if(ax==1 && !d3)
            continue;

            double smin = 1.0e30;
            double vm=v0, vp=v0, xm=X[ax], xp=X[ax];

            for(int sd=-1; sd<=1; sd+=2)
            {
                const int ni = i + (ax==0 ? sd : 0);
                const int nj = j + (ax==1 ? sd : 0);
                const int nk = k + (ax==2 ? sd : 0);

                if(!fluid(ni,nj,nk))
                continue;

                const double vn = V[fid(p,ni,nj,nk)];
                const double xn = ax==0 ? p->XP[ni+marge] : (ax==1 ? p->YP[nj+marge] : p->ZP[nk+marge]);

                if(sd<0)
                {
                vm = vn;
                xm = xn;
                }
                else
                {
                vp = vn;
                xp = xn;
                }

                if((vn>=0.0)!=w0)
                smin = std::min(smin, std::fabs(xn-X[ax])*v0/(v0-vn));
            }

            if(smin<1.0e29)
            inv2 += 1.0/std::max(smin*smin,1.0e-30);

            if(xp!=xm)
            g[ax] = (vp-vm)/(xp-xm);
        }

        if(inv2<=0.0)
        continue;

        const double d = 1.0/std::sqrt(inv2);
        const double gn = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
        const double s = w0 ? d : -d;

        iface[f] = 1;
        dist[f] = d;

        for(int ax=0; ax<3; ++ax)
        cp[3*f+ax] = X[ax] - (gn>0.0 ? s*g[ax]/gn : 0.0);
    }

    // the band: the nearest closest point within S cells
    std::vector<double> out(p->knox*p->knoy*p->knoz,0.0);
    std::vector<char> set(out.size(),0);

    for(int i=0; i<p->knox; ++i)
    for(int j=0; j<p->knoy; ++j)
    for(int k=0; k<p->knoz; ++k)
    {
        const int l = lid(p,i,j,k);

        if(!inband[l])
        continue;

        const int f = fid(p,i,j,k);
        const double X[3] = {p->XP[i+marge], p->YP[j+marge], p->ZP[k+marge]};
        double dmin2 = 1.0e30;

        for(int ni=i-S; ni<=i+S; ++ni)
        for(int nj=(d3 ? j-S : 0); nj<=(d3 ? j+S : 0); ++nj)
        for(int nk=k-S; nk<=k+S; ++nk)
        {
            const int g = fid(p,ni,nj,nk);

            if(!iface[g])
            continue;

            double d2;

            if(g==f)
            d2 = dist[f]*dist[f];
            else
            {
            const double dx = X[0]-cp[3*g], dy = d3 ? X[1]-cp[3*g+1] : 0.0, dz = X[2]-cp[3*g+2];
            d2 = dx*dx + dy*dy + dz*dz;
            }

            dmin2 = std::min(dmin2,d2);
        }

        if(dmin2<1.0e29)
        {
        out[l] = (V[f]>=0.0 ? 1.0 : -1.0)*std::sqrt(dmin2);
        set[l] = 1;
        }
    }

    for(int i=0; i<p->knox; ++i)
    for(int j=0; j<p->knoy; ++j)
    for(int k=0; k<p->knoz; ++k)
    {
        const int l = lid(p,i,j,k);

        if(set[l])
        a->phi(i,j,k) = out[l];
    }
}
