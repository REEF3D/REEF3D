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

#include"reefamr.h"
#include"lexer.h"
#include"6DOF_obj.h"
#include<cmath>

//  Refinement zone of the moving bodies (zone_bodies of the module): a margin zr around the
//  wetted hull, reaching ahead of the bow by the distance travelled until the next regrid,
//  and a wake wedge behind the bow with half angle za and length zL at most (only where the
//  bow has been).  With zalign the rectangle is aligned with x and y and grown on every side
//  by the distance travelled until the next regrid, without wake.  zone_setup also lists the hull triangles below the still water level
//  (ztri), which the module can use to evaluate the body on its patches.

// triangles below the still water level and the refinement zone of every body
void reefamr::zone_setup(lexer *p)
{
    vector<sixdof_obj*> obj;
    zone_bodies(obj);

    const double psi = 1.0e-8*p->DXM;
    const int nbody = (int)obj.size();

    ztri.assign(nbody,vector<int>());
    zones.assign(nbody,reefamr_zone());

    if((int)zx0.size()!=nbody)
    {
        zx0.resize(nbody);
        zy0.resize(nbody);
        for(int nb=0; nb<nbody; ++nb)
        {
            zx0[nb] = obj[nb]->amr_c(0);
            zy0[nb] = obj[nb]->amr_c(1);
        }
    }

    for(int nb=0; nb<nbody; ++nb)
    {
        sixdof_obj *o = obj[nb];
        double **tx = o->amr_tri(0), **ty = o->amr_tri(1), **tz = o->amr_tri(2);

        reefamr_zone &z = zones[nb];
        z.cx = o->amr_c(0);
        z.cy = o->amr_c(1);

        // heading: direction of motion, the yaw angle for a body at rest; zalign: the x axis
        // (a moored body has no meaningful heading, and every swing of it changed the zone)
        double ux = o->amr_u(0), uy = o->amr_u(1);
        double sp = sqrt(ux*ux + uy*uy);
        if(par.zalign)
        {
            z.ex = 1.0;
            z.ey = 0.0;
        }
        else
        if(sp>1.0e-6)
        {
            z.ex = ux/sp;
            z.ey = uy/sp;
        }
        else
        {
            z.ex = cos(p->psi_fb);
            z.ey = sin(p->psi_fb);
        }

        z.smin = z.nmin = 1.0e20;
        z.smax = z.nmax = -1.0e20;
        z.bx0 = z.by0 = 1.0e20;
        z.bx1 = z.by1 = -1.0e20;

        for(int n=0; n<o->amr_tricount(); ++n)
        {
            if(tz[n][0]>p->wd+psi && tz[n][1]>p->wd+psi && tz[n][2]>p->wd+psi)
            continue;

            ztri[nb].push_back(n);

            for(int q=0; q<3; ++q)
            {
                z.bx0 = MIN(z.bx0,tx[n][q]); z.bx1 = MAX(z.bx1,tx[n][q]);
                z.by0 = MIN(z.by0,ty[n][q]); z.by1 = MAX(z.by1,ty[n][q]);
            }

            for(int q=0; q<3; ++q)
            if(tz[n][q]<=p->wd+psi)
            {
                double dx = tx[n][q]-z.cx, dy = ty[n][q]-z.cy;
                double s = dx*z.ex + dy*z.ey;
                double r = -dx*z.ey + dy*z.ex;
                z.smin = MIN(z.smin,s); z.smax = MAX(z.smax,s);
                z.nmin = MIN(z.nmin,r); z.nmax = MAX(z.nmax,r);
            }
        }

        // the zone reaches ahead of the bow by the distance travelled until the next regrid
        const double ahead = 1.5*sp*p->dt*MAX(regrid_int,1);
        z.sfront = z.smax + par.zr + ahead;
        z.sback = z.smin - par.zr;
        z.nlo = z.nmin - par.zr;
        z.nhi = z.nmax + par.zr;

        // the wake wedge reaches back to where the bow has been
        double trav = sqrt(pow(z.cx-zx0[nb],2.0) + pow(z.cy-zy0[nb],2.0));
        z.wake = MIN(par.zL, (z.smax-z.smin) + trav + par.zr);

        // aligned zone: the travel distance on every side, no wake
        if(par.zalign)
        {
            z.sback -= ahead;
            z.nlo -= ahead;
            z.nhi += ahead;
            z.wake = 0.0;
        }

        if(z.smin>z.smax)       // no part of the body in the water
        z.smin = z.smax = z.nmin = z.nmax = z.sfront = z.sback = z.nlo = z.nhi = z.wake = 0.0;
    }
}

bool reefamr::zone_test(double x, double y)
{
    const double r = par.zr;
    const double ta = tan(par.za*PI/180.0);

    for(auto &z : zones)
    {
        if(z.smin>=z.smax)
        continue;

        double dx = x-z.cx, dy = y-z.cy;
        double s = dx*z.ex + dy*z.ey;
        double n = -dx*z.ey + dy*z.ex;

        // hull
        if(s>=z.sback && s<=z.sfront && n>=z.nlo && n<=z.nhi)
        return true;

        // wake: wedge from the bow with half angle za, length zL at most (only where the bow has been)
        if(z.wake>0.0)
        {
            double d = z.smax - s;
            if(d>=0.0 && d<=z.wake)
            {
                double half = 0.5*(z.nmax-z.nmin) + r + d*ta;
                if(fabs(n - 0.5*(z.nmin+z.nmax))<=half)
                return true;
            }
        }
    }
    return false;
}
