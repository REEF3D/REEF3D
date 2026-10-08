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

#include"seastate_obstacle.h"
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_param.h"
#include"lexer.h"
#include<algorithm>
#include<cmath>

seastate_obstacle::seastate_obstacle(lexer *p) : nobs(p->A722), imin(p->imin), jmin(p->jmin), ni(p->imax), nj(p->jmax)
{
    for(int n=0; n<nobs; ++n)
    {
    xs.push_back(p->A722_xs[n]);
    ys.push_back(p->A722_ys[n]);
    xe.push_back(p->A722_xe[n]);
    ye.push_back(p->A722_ye[n]);
    kt.push_back(p->A722_kt[n]);
    kr.push_back(p->A722_kr[n]);
    zc.push_back(p->A722_zc[n]);
    }
}

// the segments a-b and c-d cross (proper intersection or touching)
static bool cross(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy)
{
    auto orient = [](double px, double py, double qx, double qy, double rx, double ry)
    {
        const double v = (qx-px)*(ry-py) - (qy-py)*(rx-px);
        return (v>0.0) - (v<0.0);
    };

    const int o1 = orient(ax,ay,bx,by,cx,cy), o2 = orient(ax,ay,bx,by,dx,dy);
    const int o3 = orient(cx,cy,dx,dy,ax,ay), o4 = orient(cx,cy,dx,dy,bx,by);

    return o1*o2<=0 && o3*o4<=0 && !(o1==0 && o2==0);
}

void seastate_obstacle::build(lexer *p)
{
    fe.assign(size_t(ni)*nj,-1);
    fn.assign(size_t(ni)*nj,-1);
    faces.clear();
    fi.clear();
    fj.clear();
    fdir.clear();

    if(nobs==0)
    return;

    for(int i=imin; i<imin+ni; ++i)
    for(int j=jmin; j<jmin+nj; ++j)
    for(int s=0; s<2; ++s)
    {
    const int i2 = (s==0) ? i+1 : i, j2 = (s==0) ? j : j+1;

        if(i2>=imin+ni || j2>=jmin+nj)
        continue;

    const double x1 = p->XP[i+p->margin], y1 = p->YP[j+p->margin];
    const double x2 = p->XP[i2+p->margin], y2 = p->YP[j2+p->margin];

        for(int n=0; n<nobs; ++n)
        if(cross(x1,y1,x2,y2,xs[n],ys[n],xe[n],ye[n]))
        {
        face f;
        f.obs = n;
        f.alpha = float(std::atan2(ye[n]-ys[n],xe[n]-xs[n]));
        f.kt2 = (kt[n]>=0.0) ? float(kt[n]*kt[n]) : 1.0f;
        f.kr2 = float(std::min(kr[n]*kr[n],1.0-double(f.kt2)));

        const size_t c = size_t(i-imin)*nj + (j-jmin);
        (s==0 ? fe : fn)[c] = int(faces.size());
        faces.push_back(f);
        fi.push_back(i);
        fj.push_back(j);
        fdir.push_back(s);
        break;
        }
    }
}

void seastate_obstacle::update(lexer *p, fdm_seastate *e)
{
    const double pi = 3.14159265358979323846;
    const double ga = 2.6, gb = 0.15;
    seastate_param sp;

    auto hs = [&](int i, int j)
    {
        if(e->wet(i,j)!=1 || e->N->spec(i,j)==nullptr)
        return 0.0;
        sp.compute(*e->grid,e->N->spec(i,j));
        return sp.Hs;
    };

    for(size_t f=0; f<faces.size(); ++f)
    {
    const int n = faces[f].obs;

        if(kt[n]>=0.0)
        continue;

    const int i = fi[f], j = fj[f];
    const int i2 = (fdir[f]==0) ? i+1 : i, j2 = (fdir[f]==0) ? j : j+1;
    const double H = std::max(hs(i,j),hs(i2,j2));

    double t = 1.0;

        if(H>1.0e-6)
        {
        const double wl = p->wd + 0.5*(e->eta(i,j) + e->eta(i2,j2));
        const double r = (zc[n]-wl)/H;

        if(r>=ga-gb)
        t = 0.0;
        else if(r>-gb-ga)
        t = 0.5*(1.0 - std::sin(pi/(2.0*ga)*(r+gb)));
        }

    faces[f].kt2 = float(t*t);
    faces[f].kr2 = float(std::min(kr[n]*kr[n],1.0-t*t));
    }
}

const seastate_obstacle::face *seastate_obstacle::east(int i, int j) const
{
    if(nobs==0 || i<imin || i>=imin+ni || j<jmin || j>=jmin+nj)
    return nullptr;

    const int k = fe[size_t(i-imin)*nj + (j-jmin)];
    return k<0 ? nullptr : &faces[k];
}

const seastate_obstacle::face *seastate_obstacle::north(int i, int j) const
{
    if(nobs==0 || i<imin || i>=imin+ni || j<jmin || j>=jmin+nj)
    return nullptr;

    const int k = fn[size_t(i-imin)*nj + (j-jmin)];
    return k<0 ? nullptr : &faces[k];
}
