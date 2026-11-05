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

#include"nhflow_particle_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<cmath>

// L 21: point  x y z n ts te ws cw mode
// L 22: line   xs ys zs xe ye ze n ts te ws cw mode
// L 23: box    xs xe ys ye zs ze n ts te ws cw mode
// Horizontal coordinates are given in world coordinates and converted here,
// after the grid (and a possible grid rotation) is known.

void nhflow_particle_f::setup_releases(lexer *p)
{
    long long idoffset=0;
    release r;

    auto finish = [&](const double *v)
    {
        r.num  = (long long)(v[0]+0.5);
        r.ts   = v[1];
        r.te   = v[2];
        r.ws   = v[3];
        r.cw   = v[4];
        r.mode = int(v[5]+0.5);
        r.idoffset = idoffset;
        r.released = 0;
        idoffset += r.num;

        double x,y;
        x = p->Xin(r.xs,r.ys);
        y = p->Yin(r.xs,r.ys);
        r.xs=x; r.ys=y;
        x = p->Xin(r.xe,r.ye);
        y = p->Yin(r.xe,r.ye);
        r.xe=x; r.ye=y;

        if(r.num>0)
        R.push_back(r);
    };

    for(int q=0; q<p->L21; ++q)
    {
        const double *v = &p->L21_val[9*q];
        r.type=1;
        r.xs=r.xe=v[0];
        r.ys=r.ye=v[1];
        r.zs=r.ze=v[2];
        finish(&v[3]);
    }

    for(int q=0; q<p->L22; ++q)
    {
        const double *v = &p->L22_val[12*q];
        r.type=2;
        r.xs=v[0]; r.ys=v[1]; r.zs=v[2];
        r.xe=v[3]; r.ye=v[4]; r.ze=v[5];
        finish(&v[6]);
    }

    for(int q=0; q<p->L23; ++q)
    {
        const double *v = &p->L23_val[12*q];
        r.type=3;
        r.xs=v[0]; r.xe=v[1];
        r.ys=v[2]; r.ye=v[3];
        r.zs=v[4]; r.ze=v[5];
        finish(&v[6]);
    }
}

double nhflow_particle_f::release_time(const release &r, long long q) const
{
    if(r.te<=r.ts)
    return r.ts;

    return r.ts + (r.te-r.ts)*double(q)/double(r.num);
}

// splitmix64, reproducible uniform numbers in [0,1) per particle id
double nhflow_particle_f::hash_uniform(uint64_t id, uint64_t c)
{
    uint64_t z = id*0x9E3779B97F4A7C15ULL + (c+1)*0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z>>30))*0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z>>27))*0x94D049BB133111EBULL;
    z = z ^ (z>>31);

    return double(z>>11)*(1.0/9007199254740992.0);
}

void nhflow_particle_f::release_position(const release &r, long long q, double &xp, double &yp, double &zp) const
{
    if(r.type==1)
    {
    xp=r.xs;
    yp=r.ys;
    zp=r.zs;
    }

    if(r.type==2)
    {
    double s = (double(q)+0.5)/double(r.num);
    xp = r.xs + s*(r.xe-r.xs);
    yp = r.ys + s*(r.ye-r.ys);
    zp = r.zs + s*(r.ze-r.zs);
    }

    if(r.type==3)
    {
    uint64_t id = uint64_t(r.idoffset + q);
    xp = r.xs + hash_uniform(id,0)*(r.xe-r.xs);
    yp = r.ys + hash_uniform(id,1)*(r.ye-r.ys);
    zp = r.zs + hash_uniform(id,2)*(r.ze-r.zs);
    }
}

// every rank runs through the same release schedule and keeps the particles inside its subdomain
void nhflow_particle_f::seed(lexer *p, fdm_nhf *d, ghostcell *pgc, double tlimit, bool inclusive)
{
    const double eps = 1.0e-12;

    for(size_t rn=0; rn<R.size(); ++rn)
    {
        release &r = R[rn];

        while(r.released<r.num)
        {
            double tr = release_time(r,r.released);

            if(inclusive ? (tr>tlimit+eps) : (tr>=tlimit-eps))
            break;

            double xp,yp,zp;
            release_position(r,r.released,xp,yp,zp);

            if(p->j_dir==0)
            yp = p->YP[0+marge];

            if(owned(p,xp,yp))
            {
                nhflow_particle_data a{};
                a.x=xp;
                a.y=yp;
                a.z=zp;
                a.ws=r.ws;
                a.cw=r.cw;
                a.mode=r.mode;
                a.state=NHFP_WATER;
                a.id=double(r.idoffset+r.released);
                a.src=double(rn+1);
                a.t0=tr;

                vertical_bounds(p,d,a);
                P.push_back(a);
            }

            ++r.released;
        }
    }
}
