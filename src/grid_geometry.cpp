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

// Horizontal geometry layer: 2D nodes, cell centroids and areas, face normal vectors.
// Nothing here is read by a solver yet; on Cartesian grids the metrics are the exact products of
// the 1D spacings (the values the kernels use today), so a later switch of the kernels to the
// geometry layer keeps the Cartesian results bit for bit.

void geo2d::centroid(double x00, double y00, double x10, double y10, double x01, double y01, double x11, double y11,
                     double &xc, double &yc)
{
    // two triangles 00-10-11 and 00-11-01, area weighted, relative to node 00 (large coordinates)
    const double ax=x10-x00, ay=y10-y00, bx=x11-x00, by=y11-y00, cx=x01-x00, cy=y01-y00;
    const double a1 = 0.5*(ax*by - ay*bx);
    const double a2 = 0.5*(bx*cy - by*cx);
    const double a = a1+a2;

    if(std::fabs(a)>0.0)
    {
        xc = x00 + (a1*(ax+bx) + a2*(bx+cx))/(3.0*a);
        yc = y00 + (a1*(ay+by) + a2*(by+cy))/(3.0*a);
    }
    else
    {
        xc = 0.25*(x00+x10+x01+x11);
        yc = 0.25*(y00+y10+y01+y11);
    }
}

void grid::geometry_alloc()
{
    const int nn = (imax+1)*(jmax+1);
    const int ns = imax*jmax;

    if(XN2D!=nullptr && nn==geo_nnode && ns==geo_nslice)
    return;

    geometry_free();

    geo_nnode = nn;
    geo_nslice = ns;

    XN2D = new double[nn]();
    YN2D = new double[nn]();
    XC2D = new double[ns]();
    YC2D = new double[ns]();
    AREA2D = new double[ns]();
    SX1 = new double[ns]();
    SY1 = new double[ns]();
    SX2 = new double[ns]();
    SY2 = new double[ns]();
}

void grid::geometry_free()
{
    delete [] XN2D; delete [] YN2D;
    delete [] XC2D; delete [] YC2D; delete [] AREA2D;
    delete [] SX1; delete [] SY1; delete [] SX2; delete [] SY2;

    XN2D=YN2D=XC2D=YC2D=AREA2D=SX1=SY1=SX2=SY2=nullptr;
    geo_nnode = geo_nslice = 0;
}

void grid::geometry_cartesian_nodes()
{
    for(int ii=imin; ii<=imin+imax; ++ii)
    for(int jj=jmin; jj<=jmin+jmax; ++jj)
    {
        XN2D[nij(ii,jj)] = XN[ii+marge];
        YN2D[nij(ii,jj)] = YN[jj+marge];
    }
}

void grid::geometry_metrics()
{
    for(int ii=imin; ii<imin+imax; ++ii)
    for(int jj=jmin; jj<jmin+jmax; ++jj)
    {
        const int q = (ii-imin)*jmax + (jj-jmin);

        if(geo_curv==0)
        {
            // exact Cartesian values: the 1D spacings and centres of gridspacing()
            const int a = ii+marge, b = jj+marge;

            XC2D[q] = XP[a];
            YC2D[q] = YP[b];
            AREA2D[q] = DXN[a]*DYN[b];
            SX1[q] = DYN[b];
            SY1[q] = 0.0;
            SX2[q] = 0.0;
            SY2[q] = DXN[a];
        }
        else
        {
            const double x00=XN2D[nij(ii,jj)],   y00=YN2D[nij(ii,jj)];
            const double x10=XN2D[nij(ii+1,jj)], y10=YN2D[nij(ii+1,jj)];
            const double x01=XN2D[nij(ii,jj+1)], y01=YN2D[nij(ii,jj+1)];
            const double x11=XN2D[nij(ii+1,jj+1)], y11=YN2D[nij(ii+1,jj+1)];

            geo2d::centroid(x00,y00,x10,y10,x01,y01,x11,y11,XC2D[q],YC2D[q]);
            AREA2D[q] = geo2d::area(x00,y00,x10,y10,x01,y01,x11,y11);

            // face i+1/2: edge 10 -> 11, normal to the right of it (+i side)
            SX1[q] =  (y11-y10);
            SY1[q] = -(x11-x10);

            // face j+1/2: edge 01 -> 11, normal to the left of it (+j side)
            SX2[q] = -(y11-y01);
            SY2[q] =  (x11-x01);
        }
    }
}

int grid::geometry_check(double &err, double &gcl) const
{
    int mismatch = 0;
    err = 0.0;
    gcl = 0.0;

    if(XN2D==nullptr)
    return -1;

    if(geo_curv==0)
    for(int ii=imin; ii<=imin+imax; ++ii)
    for(int jj=jmin; jj<=jmin+jmax; ++jj)
    if(XN2D[nij(ii,jj)]!=XN[ii+marge] || YN2D[nij(ii,jj)]!=YN[jj+marge])
    ++mismatch;

    auto sl = [&](int ii, int jj) { return (ii-imin)*jmax + (jj-jmin); };

    for(int ii=imin; ii<imin+imax; ++ii)
    for(int jj=jmin; jj<jmin+jmax; ++jj)
    {
        const int q = sl(ii,jj);

        const double x00=XN2D[nij(ii,jj)],   y00=YN2D[nij(ii,jj)];
        const double x10=XN2D[nij(ii+1,jj)], y10=YN2D[nij(ii+1,jj)];
        const double x01=XN2D[nij(ii,jj+1)], y01=YN2D[nij(ii,jj+1)];
        const double x11=XN2D[nij(ii+1,jj+1)], y11=YN2D[nij(ii+1,jj+1)];

        double xc,yc;
        geo2d::centroid(x00,y00,x10,y10,x01,y01,x11,y11,xc,yc);
        const double A = geo2d::area(x00,y00,x10,y10,x01,y01,x11,y11);
        const double L = std::sqrt(std::fabs(A));    // length scale of the cell

        if(L>0.0)
        {
            err = std::max(err, std::fabs(AREA2D[q]-A)/std::fabs(A));
            err = std::max(err, std::max(std::fabs(XC2D[q]-xc),std::fabs(YC2D[q]-yc))/L);
            err = std::max(err, std::max(std::fabs(SX1[q]-(y11-y10)),std::fabs(SY1[q]+(x11-x10)))/L);
            err = std::max(err, std::max(std::fabs(SX2[q]+(y11-y01)),std::fabs(SY2[q]-(x11-x01)))/L);
        }
        else
        err = std::max(err,1.0);

        // closure (geometric conservation): the outward face vectors of a cell sum to zero
        if(ii>imin && jj>jmin)
        {
            const int qw = sl(ii-1,jj), qs = sl(ii,jj-1);
            const double rx = SX1[q]-SX1[qw] + SX2[q]-SX2[qs];
            const double ry = SY1[q]-SY1[qw] + SY2[q]-SY2[qs];
            if(L>0.0)
            gcl = std::max(gcl, std::max(std::fabs(rx),std::fabs(ry))/L);
        }
    }

    return mismatch;
}
