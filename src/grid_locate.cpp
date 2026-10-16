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

#include"grid.h"
#include<cmath>
#include<algorithm>

// Physical points on the curvilinear grid (G 50): XP/YP and posc_i/posf_i are the logical box of DIVEMesh
// there, the physical cell centres and nodes are XC2D/YC2D and XN2D/YN2D. A bucket table over the box of
// the rank holds three kinds of items: the quadrilaterals of four cell centres (interpolation at a point),
// the cells as quadrilaterals of their nodes (the cell that holds a point) and the cell centres (the cells
// near a triangle). The table is built on the first call and freed with the geometry.

struct geo_locate_cache
{
    double x0=0.0, y0=0.0, hb=1.0;     // origin and size of the buckets
    int nbx=1, nby=1;
    double bx0=0.0, bx1=0.0, by0=0.0, by1=0.0;   // box of the nodes of the rank's cells
    std::vector<std::vector<int>> quad, cell, cent;

    int bucket(double x, double y) const
    {
        int a = int(std::floor((x-x0)/hb)), b = int(std::floor((y-y0)/hb));
        a = std::max(0, std::min(nbx-1, a));
        b = std::max(0, std::min(nby-1, b));
        return a*nby + b;
    }

    template<class F> void range(double xs, double xe, double ys, double ye, F f) const
    {
        const int as = std::max(0, int(std::floor((xs-x0)/hb))), ae = std::min(nbx-1, int(std::floor((xe-x0)/hb)));
        const int bs = std::max(0, int(std::floor((ys-y0)/hb))), be = std::min(nby-1, int(std::floor((ye-y0)/hb)));
        for(int a=as; a<=ae; ++a)
        for(int b=bs; b<=be; ++b)
        f(a*nby + b);
    }
};

namespace
{
// inverse of the bilinear map of the quadrilateral 00-10-11-01 (coordinates relative to node 00)
bool bilinear_inverse(double px, double py, double ax, double ay, double bx, double by, double cx, double cy,
                      double &s, double &t)
{
    // P(s,t) = s A + t C + s t (B - A - C);  A = 10-00, B = 11-00, C = 01-00
    const double dx = bx-ax-cx, dy = by-ay-cy;
    s = 0.5; t = 0.5;

    for(int it=0; it<30; ++it)
    {
        const double fx = s*ax + t*cx + s*t*dx - px;
        const double fy = s*ay + t*cy + s*t*dy - py;
        const double j11 = ax + t*dx, j12 = cx + s*dx;
        const double j21 = ay + t*dy, j22 = cy + s*dy;
        const double det = j11*j22 - j12*j21;

        if(std::fabs(det) < 1.0e-300)
        return false;

        const double ds = ( j22*fx - j12*fy)/det;
        const double dt = (-j21*fx + j11*fy)/det;
        s -= ds;
        t -= dt;

        if(std::fabs(ds)+std::fabs(dt) < 1.0e-15)
        break;
    }
    return std::isfinite(s) && std::isfinite(t);
}
}

bool grid::geo_bilinear_inverse(double px, double py, double ax, double ay, double bx, double by,
                                double cx, double cy, double &s, double &t)
{
    return bilinear_inverse(px,py,ax,ay,bx,by,cx,cy,s,t);
}

static geo_locate_cache *geo_build(grid *g)
{
    geo_locate_cache *c = new geo_locate_cache;
    const int ni = g->knox, nj = g->knoy;

    auto xn = [&](int a, int b) { return g->XN2D[g->nij(a,b)]; };
    auto yn = [&](int a, int b) { return g->YN2D[g->nij(a,b)]; };
    auto xc = [&](int a, int b) { return g->XC2D[(a-g->imin)*g->jmax + (b-g->jmin)]; };
    auto yc = [&](int a, int b) { return g->YC2D[(a-g->imin)*g->jmax + (b-g->jmin)]; };

    // the box of the rank: nodes of the cells, and the centres of the quads (one ghost layer)
    double x0=1.0e300, x1=-1.0e300, y0=1.0e300, y1=-1.0e300;
    for(int a=0; a<=ni; ++a)
    for(int b=0; b<=nj; ++b)
    {
        x0 = std::min(x0,xn(a,b)); x1 = std::max(x1,xn(a,b));
        y0 = std::min(y0,yn(a,b)); y1 = std::max(y1,yn(a,b));
    }
    c->bx0 = x0; c->bx1 = x1; c->by0 = y0; c->by1 = y1;

    double ex0=x0, ex1=x1, ey0=y0, ey1=y1;
    for(int a=-1; a<=ni; ++a)
    for(int b=-1; b<=nj; ++b)
    {
        ex0 = std::min(ex0,xc(a,b)); ex1 = std::max(ex1,xc(a,b));
        ey0 = std::min(ey0,yc(a,b)); ey1 = std::max(ey1,yc(a,b));
    }

    const double area = std::max((ex1-ex0)*(ey1-ey0), 1.0e-300);
    c->hb = std::max(std::sqrt(area/double(std::max(1,(ni+2)*(nj+2))))*2.0, 1.0e-12);
    c->x0 = ex0;
    c->y0 = ey0;
    c->nbx = std::max(1, int((ex1-ex0)/c->hb)+1);
    c->nby = std::max(1, int((ey1-ey0)/c->hb)+1);
    c->quad.assign(c->nbx*c->nby, {});
    c->cell.assign(c->nbx*c->nby, {});
    c->cent.assign(c->nbx*c->nby, {});

    auto insert = [&](std::vector<std::vector<int>> &T, int key, double qx0, double qx1, double qy0, double qy1)
    {
        c->range(qx0, qx1, qy0, qy1, [&](int bk) { T[bk].push_back(key); });
    };

    const int W = nj+2;     // key = (a+1)*W + (b+1)

    for(int a=-1; a<ni; ++a)
    for(int b=-1; b<nj; ++b)
    {
        const double X[4] = {xc(a,b), xc(a+1,b), xc(a+1,b+1), xc(a,b+1)};
        const double Y[4] = {yc(a,b), yc(a+1,b), yc(a+1,b+1), yc(a,b+1)};
        insert(c->quad, (a+1)*W + (b+1), *std::min_element(X,X+4), *std::max_element(X,X+4),
               *std::min_element(Y,Y+4), *std::max_element(Y,Y+4));
    }

    for(int a=0; a<ni; ++a)
    for(int b=0; b<nj; ++b)
    {
        const double X[4] = {xn(a,b), xn(a+1,b), xn(a+1,b+1), xn(a,b+1)};
        const double Y[4] = {yn(a,b), yn(a+1,b), yn(a+1,b+1), yn(a,b+1)};
        insert(c->cell, (a+1)*W + (b+1), *std::min_element(X,X+4), *std::max_element(X,X+4),
               *std::min_element(Y,Y+4), *std::max_element(Y,Y+4));
        c->cent[c->bucket(xc(a,b),yc(a,b))].push_back((a+1)*W + (b+1));
    }

    return c;
}

void grid::geo_locate_free()
{
    delete geoloc;
    geoloc = nullptr;
}

bool grid::geo_quad(double x, double y, int &ii, int &jj, double &fx, double &fy)
{
    if(geoloc==nullptr)
    geoloc = geo_build(this);

    const geo_locate_cache &c = *geoloc;
    const int W = knoy+2;

    auto xc = [&](int a, int b) { return XC2D[(a-imin)*jmax + (b-jmin)]; };
    auto yc = [&](int a, int b) { return YC2D[(a-imin)*jmax + (b-jmin)]; };

    auto test = [&](int key, double &s, double &t)
    {
        const int a = key/W - 1, b = key%W - 1;
        const double ox = xc(a,b), oy = yc(a,b);
        return bilinear_inverse(x-ox, y-oy, xc(a+1,b)-ox, yc(a+1,b)-oy, xc(a+1,b+1)-ox, yc(a+1,b+1)-oy,
                                xc(a,b+1)-ox, yc(a,b+1)-oy, s, t);
    };

    const double tol = 1.0e-12;
    int best=-1;
    double bs=0.0, bt=0.0, bd=1.0e300;

    auto scan = [&](const std::vector<int> &L)
    {
        for(int key : L)
        {
            double s,t;
            if(!test(key,s,t))
            continue;

            // distance of (s,t) from the unit square: 0 inside
            const double d = std::max(0.0, std::max(-s, s-1.0)) + std::max(0.0, std::max(-t, t-1.0));
            if(d < bd || (d==bd && key < best))
            {
                bd = d; best = key; bs = s; bt = t;
            }
        }
    };

    scan(c.quad[c.bucket(x,y)]);

    // outside the quads of the bucket: all quads (rare: points outside the rank)
    if(best<0 || bd>tol)
    for(const auto &L : c.quad)
    scan(L);

    if(best<0)
    {
        ii = 0; jj = 0; fx = fy = 0.0;
        return false;
    }

    ii = best/W - 1;
    jj = best%W - 1;
    fx = std::max(0.0, std::min(1.0, bs));
    fy = std::max(0.0, std::min(1.0, bt));
    return bd<=tol;
}

bool grid::geo_cell(double x, double y, int &ii, int &jj)
{
    if(geoloc==nullptr)
    geoloc = geo_build(this);

    const geo_locate_cache &c = *geoloc;
    const int W = knoy+2;

    if(x<c.bx0 || x>c.bx1 || y<c.by0 || y>c.by1)
    return false;

    auto xn = [&](int a, int b) { return XN2D[nij(a,b)]; };
    auto yn = [&](int a, int b) { return YN2D[nij(a,b)]; };

    for(int key : c.cell[c.bucket(x,y)])
    {
        const int a = key/W - 1, b = key%W - 1;
        const double ox = xn(a,b), oy = yn(a,b);
        double s,t;
        if(!bilinear_inverse(x-ox, y-oy, xn(a+1,b)-ox, yn(a+1,b)-oy, xn(a+1,b+1)-ox, yn(a+1,b+1)-oy,
                             xn(a,b+1)-ox, yn(a,b+1)-oy, s, t))
        continue;

        // half-open: the face s=1 belongs to the next cell; the outer faces of the rank are closed so that
        // a point on the last face of the domain still has a cell
        const bool sin = s>=0.0 && (s<1.0 || (s==1.0 && a+origin_i==gknox-1));
        const bool tin = t>=0.0 && (t<1.0 || (t==1.0 && b+origin_j==gknoy-1));

        if(sin && tin)
        {
            ii = a;
            jj = b;
            return true;
        }
    }
    return false;
}

void grid::geo_cells(double x0, double x1, double y0, double y1, std::vector<int> &ij)
{
    if(geoloc==nullptr)
    geoloc = geo_build(this);

    const geo_locate_cache &c = *geoloc;
    const int W = knoy+2;
    ij.clear();

    c.range(x0, x1, y0, y1, [&](int bk)
    {
        for(int key : c.cent[bk])
        {
            const int a = key/W - 1, b = key%W - 1;
            const int q = (a-imin)*jmax + (b-jmin);
            if(XC2D[q]>=x0 && XC2D[q]<=x1 && YC2D[q]>=y0 && YC2D[q]<=y1)
            {
                ij.push_back(a);
                ij.push_back(b);
            }
        }
    });
}

void grid::geo_box(double &x0, double &x1, double &y0, double &y1)
{
    if(geoloc==nullptr)
    geoloc = geo_build(this);

    x0 = geoloc->bx0; x1 = geoloc->bx1;
    y0 = geoloc->by0; y1 = geoloc->by1;
}

double grid::geo_width(int ii, int jj) const
{
    const int q = (ii-imin)*jmax + (jj-jmin);
    const double A = AREA2D[q];
    const double l1 = 0.5*(LEN1[q] + LEN1[q-jmax]);
    const double l2 = 0.5*(LEN2[q] + LEN2[q-1]);
    return 0.5*(A/l1 + A/l2);
}

void grid::geo_grad(int ii, int jj, double fE, double fW, double fN, double fS, double &fx, double &fy) const
{
    auto qq = [&](int a, int b) { return (a-imin)*jmax + (b-jmin); };
    const double xe = XC2D[qq(ii+1,jj)] - XC2D[qq(ii-1,jj)], ye = YC2D[qq(ii+1,jj)] - YC2D[qq(ii-1,jj)];
    const double xn = XC2D[qq(ii,jj+1)] - XC2D[qq(ii,jj-1)], yn = YC2D[qq(ii,jj+1)] - YC2D[qq(ii,jj-1)];
    const double fe = fE - fW, fn = fN - fS;
    const double det = xe*yn - ye*xn;

    if(std::fabs(det)>0.0)
    {
    fx = (fe*yn - ye*fn)/det;
    fy = (xe*fn - fe*xn)/det;
    }
    else
    fx = fy = 0.0;
}
