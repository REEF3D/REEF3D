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

#include"lexer.h"
#include"ghostcell.h"
#include"geo_mesh.h"
#include<cmath>
#include<algorithm>

// horizontal geometry layer of the level-0 grid (grid_geometry.cpp), with its consistency check,
// and the check of a curvilinear river corridor grid (CURV) of the geometry file
void lexer::geometry_ini(ghostcell *pgc)
{
    geo_curv = 0;

    geometry_alloc();
    geometry_cartesian_nodes();
    geometry_metrics();

    double err=0.0, gcl=0.0;
    int mismatch = geometry_check(err,gcl);

    mismatch = pgc->globalisum(mismatch);
    err = pgc->globalmax(err);
    gcl = pgc->globalmax(gcl);

    if(mpirank==0)
    {
    cout<<"horizontal geometry: "<<(geo_curv==0?"Cartesian":"curvilinear")
        <<", metric check "<<err<<", closure "<<gcl<<endl;

    if(mismatch!=0)
    cout<<"!!! horizontal geometry: "<<mismatch<<" nodes differ from XN/YN !!!"<<endl;
    }

    if(gridgeo!=nullptr && gridgeo->curv_ni>0 && mpirank==0)
    curv_report();
}

// CURV: quadrilateral metrics of the corridor grid (folded cells, areas, orthogonality, stretching)
void lexer::curv_report()
{
    const geo_mesh &G = *gridgeo;
    const int ni = G.curv_ni, nj = G.curv_nj;
    auto id = [&](int ii, int jj) { return ii + (ni+1)*jj; };
    const vector<double> &X = G.curv_x, &Y = G.curv_y;

    int folded = 0;
    double amin = 1.0e20, amax = 0.0, nonorth = 0.0, ratio = 1.0;

    for(int jj=0; jj<nj; ++jj)
    for(int ii=0; ii<ni; ++ii)
    {
        const double x00=X[id(ii,jj)],   y00=Y[id(ii,jj)];
        const double x10=X[id(ii+1,jj)], y10=Y[id(ii+1,jj)];
        const double x01=X[id(ii,jj+1)], y01=Y[id(ii,jj+1)];
        const double x11=X[id(ii+1,jj+1)], y11=Y[id(ii+1,jj+1)];

        const double A = geo2d::area(x00,y00,x10,y10,x01,y01,x11,y11);

        if(A<=0.0)
        ++folded;

        amin = std::min(amin,A);
        amax = std::max(amax,A);

        // angle between the cell's i and j directions (midpoints of opposite edges)
        const double ax = 0.5*(x10+x11-x00-x01), ay = 0.5*(y10+y11-y00-y01);
        const double bx = 0.5*(x01+x11-x00-x10), by = 0.5*(y01+y11-y00-y10);
        const double la = std::sqrt(ax*ax+ay*ay), lb = std::sqrt(bx*bx+by*by);

        if(la>0.0 && lb>0.0)
        {
            const double c = std::fabs(ax*bx+ay*by)/(la*lb);
            nonorth = std::max(nonorth, std::asin(std::min(c,1.0))*180.0/std::acos(-1.0));
        }

        // stretching: area of the next cell along i and across j
        auto areaij = [&](int a, int b)
        {
            return geo2d::area(X[id(a,b)],Y[id(a,b)],X[id(a+1,b)],Y[id(a+1,b)],
                               X[id(a,b+1)],Y[id(a,b+1)],X[id(a+1,b+1)],Y[id(a+1,b+1)]);
        };

        if(A>0.0)
        {
            if(ii+1<ni) { const double An=areaij(ii+1,jj); if(An>0.0) ratio = std::max(ratio,std::sqrt(std::max(A/An,An/A))); }
            if(jj+1<nj) { const double An=areaij(ii,jj+1); if(An>0.0) ratio = std::max(ratio,std::sqrt(std::max(A/An,An/A))); }
        }
    }

    cout<<"CURV corridor grid: "<<ni<<" x "<<nj<<" cells, area "<<amin<<" .. "<<amax<<" m2, "
        <<folded<<" folded, non-orthogonality max "<<nonorth<<" deg, neighbour size ratio max "<<ratio
        <<" (not used by the solvers yet)"<<endl;

    if(folded>0)
    cout<<"!!! CURV corridor grid: "<<folded<<" folded cells !!!"<<endl;
}
