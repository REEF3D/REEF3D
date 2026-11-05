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

#include"cfd_amr.h"
#include"lexer.h"
#include"6DOF_cfd.h"
#include"6DOF_obj_cfd.h"
#include<cmath>

//  The refinement zone of the floating bodies (X 10 1) of the CFD mesh refinement, G 12 r: as the zone of
//  the column models (reefamr_zone.cpp), with the vertical extent of the body.
//
//  In plan the zone of a body is a rectangle in the frame of its heading (s along, n across): the
//  extent of its triangles grown by the margin r, ahead of the bow by the distance the body travels
//  until the next regrid (1.5 u dt G 2).  Without G 13 the heading is x and the rectangle is grown by
//  that distance on every side (moored and oscillating bodies).  With G 13 L a the heading is the
//  direction of motion (the yaw angle at rest) and a wake wedge of half angle a reaches back from the
//  bow by L at most, only where the bow has been (from the start of the body).  In the vertical the zone
//  is the extent of the triangles grown by r and the vertical travel until the next regrid.  The tiles
//  of the finest level with a cell centre in the zone are flagged (seed); the nesting of reefamr3d
//  adds the coarser levels.  All ranks hold the bodies, so every rank flags the same tiles.

void cfd_amr::zone_seed(vector<vector<unsigned char>> &seed)
{
    if(!bzone || b6==nullptr)
    return;

    lexer *p = p0;
    const int l = maxlev;
    const double zr = p->G12_r;
    const bool wake = (p->G13_L>0.0);
    const double ta = tan(p->G13_a*PI/180.0);
    const double trav_dt = 1.5*p->dt*double(MAX(p->G2,1));
    const int nob = b6->objects();

    if((int)bz_x0.size()!=nob)
    {
        bz_x0.resize(nob);
        bz_y0.resize(nob);
        for(int nb=0; nb<nob; ++nb)
        {
            bz_x0[nb] = b6->object(nb)->amr_c(0);
            bz_y0[nb] = b6->object(nb)->amr_c(1);
        }
    }

    // the cells of level l with the centre in [xs,xe] along d
    auto range = [&](int d, double xs, double xe, int &lo, int &hi)
    {
        lo = gn(l,d);
        hi = -1;
        if(rr[d]==1)
        {
            lo = 0;
            hi = gn(l,d)-1;
            return;
        }
        for(int I=0; I<gn(l,d); ++I)
        {
            const double xc = 0.5*(node(l,d,I)+node(l,d,I+1));
            if(xc>=xs && xc<=xe)
            {
                lo = MIN(lo,I);
                hi = MAX(hi,I);
            }
        }
    };

    for(int nb=0; nb<nob; ++nb)
    {
        sixdof_obj_cfd *o = b6->object(nb);
        double **tx = o->amr_tri(0), **ty = o->amr_tri(1), **tz = o->amr_tri(2);
        const int nt = o->amr_tricount();
        if(nt<=0)
        continue;

        const double cx = o->amr_c(0), cy = o->amr_c(1);
        const double ux = o->amr_u(0), uy = o->amr_u(1), uz = o->amr_u(2);
        const double sp = sqrt(ux*ux + uy*uy);

        // the heading: x (no G 13), the direction of motion, the yaw angle at rest
        double ex = 1.0, ey = 0.0;
        if(wake)
        {
            if(sp>1.0e-6)
            {
                ex = ux/sp;
                ey = uy/sp;
            }
            else
            {
                ex = cos(p->psi_fb);
                ey = sin(p->psi_fb);
            }
        }

        double smin = 1.0e20, smax = -1.0e20, nmin = 1.0e20, nmax = -1.0e20, zmin = 1.0e20, zmax = -1.0e20;
        for(int n=0; n<nt; ++n)
        for(int q=0; q<3; ++q)
        {
            const double dx = tx[n][q]-cx, dy = ty[n][q]-cy;
            const double s = dx*ex + dy*ey;
            const double r = -dx*ey + dy*ex;
            smin = MIN(smin,s); smax = MAX(smax,s);
            nmin = MIN(nmin,r); nmax = MAX(nmax,r);
            zmin = MIN(zmin,tz[n][q]); zmax = MAX(zmax,tz[n][q]);
        }

        const double ahead = sp*trav_dt;
        double sfront = smax + zr + ahead, sback = smin - zr, nlo = nmin - zr, nhi = nmax + zr;
        double wl = 0.0;
        if(wake)
        {
            const double trav = sqrt(pow(cx-bz_x0[nb],2.0) + pow(cy-bz_y0[nb],2.0));
            wl = MIN(p->G13_L,(smax-smin) + trav + zr);
        }
        else
        {
            sback -= ahead;
            nlo -= ahead;
            nhi += ahead;
        }
        const double zlo = zmin - zr - fabs(uz)*trav_dt, zhi = zmax + zr + fabs(uz)*trav_dt;

        // the zone at (x,y,z)
        auto inside = [&](double x, double y, double z)
        {
            if(z<zlo || z>zhi)
            return false;
            const double dx = x-cx, dy = y-cy;
            const double s = dx*ex + dy*ey;
            const double n = -dx*ey + dy*ex;
            if(s>=sback && s<=sfront && n>=nlo && n<=nhi)
            return true;
            if(wl>0.0)
            {
                const double d = smax - s;
                if(d>=0.0 && d<=wl)
                {
                    const double half = 0.5*(nmax-nmin) + zr + d*ta;
                    if(fabs(n - 0.5*(nmin+nmax))<=half)
                    return true;
                }
            }
            return false;
        };

        // its box in plan: the corners of the rectangle and of the wake
        double bx0 = 1.0e20, bx1 = -1.0e20, by0 = 1.0e20, by1 = -1.0e20;
        auto corner = [&](double s, double n)
        {
            const double x = cx + s*ex - n*ey, y = cy + s*ey + n*ex;
            bx0 = MIN(bx0,x); bx1 = MAX(bx1,x);
            by0 = MIN(by0,y); by1 = MAX(by1,y);
        };
        corner(sback,nlo); corner(sback,nhi); corner(sfront,nlo); corner(sfront,nhi);
        if(wl>0.0)
        {
            const double half = 0.5*(nmax-nmin) + zr + wl*ta, nc = 0.5*(nmin+nmax);
            corner(smax-wl,nc-half); corner(smax-wl,nc+half);
        }

        int lo[3], hi[3];
        range(0,bx0,bx1,lo[0],hi[0]);
        range(1,by0,by1,lo[1],hi[1]);
        range(2,zlo,zhi,lo[2],hi[2]);
        if(hi[0]<lo[0] || hi[1]<lo[1] || hi[2]<lo[2])
        continue;

        int I[3];
        for(I[0]=lo[0]; I[0]<=hi[0]; ++I[0])
        for(I[1]=lo[1]; I[1]<=hi[1]; ++I[1])
        for(I[2]=lo[2]; I[2]<=hi[2]; ++I[2])
        {
            const size_t t = tid(l,I[0]/tile[0],I[1]/tile[1],I[2]/tile[2]);
            if(seed[l][t])
            continue;
            const double x = 0.5*(node(l,0,I[0])+node(l,0,I[0]+1));
            const double y = 0.5*(node(l,1,I[1])+node(l,1,I[1]+1));
            const double z = 0.5*(node(l,2,I[2])+node(l,2,I[2]+1));
            if(inside(x,y,z))
            seed[l][t] = 1;
        }
    }
}
