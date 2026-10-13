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
// G 50 1/2: curvilinear grid, the 2D nodes of the global logical box (knox x knoy cells of the
// DIVEMesh grid) from the CURV section or from curvgrid.dat; XN/YN stay the logical box
void lexer::geometry_ini(ghostcell *pgc)
{
    geo_curv = 0;

    geometry_alloc();
    geometry_cartesian_nodes();

    if(G50>0)
    curv_nodes(pgc);

    // multi-block grid: block 0 with its ring is this grid, the other blocks are level-0 grids of
    // the mesh refinement (reefamr_block.cpp)
    geo_mb = 0;
    if(G50==1 && gridgeo!=nullptr && gridgeo->blk.size()>1)
    geo_mb = int(gridgeo->blk.size());

    // the blocks have their own index directions: walls act the same along i and j (G 52 1)
    if(geo_mb>0)
    G52 = 1;

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

    if(geo_curv==1)
    {
        // interior cells must not be folded; multi-block: the ring of block 0 holds the cells of its
        // neighbours, at a singular vertex (three or five blocks around a node) the corner of the ring
        // has no cells: they are degenerate and become solid
        int folded = 0, ringdeg = 0;
        double amin = 1.0e20, amax = 0.0;
        const geo_mesh::blk_block *B0 = (geo_mb>0) ? &gridgeo->blk[0] : nullptr;

        for(int ii=0; ii<knox; ++ii)
        for(int jj=0; jj<knoy; ++jj)
        {
            const int q = (ii-imin)*jmax + (jj-jmin);
            const double A = AREA2D[q];
            if(B0!=nullptr)
            {
                const int I = ii+origin_i, J = jj+origin_j;
                const bool ring = I<B0->ring[0] || I>=B0->ring[0]+B0->ni || J<B0->ring[2] || J>=B0->ring[2]+B0->nj;
                if(ring)
                {
                    if(A<=0.0)
                    {
                        ++ringdeg;
                        flagslice4[q] = -1;
                    }
                    continue;
                }
            }
            if(A<=0.0)
            ++folded;
            amin = std::min(amin,A);
            amax = std::max(amax,A);
        }

        folded = pgc->globalisum(folded);
        ringdeg = pgc->globalisum(ringdeg);
        amin = pgc->globalmin(amin);
        amax = pgc->globalmax(amax);

        if(ringdeg>0 && mpirank==0)
        cout<<"multi-block grid: "<<ringdeg<<" degenerate cells in the ring of block 0 (singular vertices), solid"<<endl;

        if(mpirank==0)
        cout<<"curvilinear grid (G 50 "<<G50<<"): "<<gknox<<" x "<<gknoy<<" cells, area "<<amin<<" .. "<<amax<<" m2"<<endl;

        if(folded>0)
        {
            if(mpirank==0)
            cout<<"!!! curvilinear grid: "<<folded<<" folded cells !!!"<<endl;
            exit(1);
        }
    }

    // (multi-block: the CURV section is block 0 with its ring, checked above)
    if(gridgeo!=nullptr && gridgeo->curv_ni>0 && mpirank==0 && geo_mb==0)
    curv_report();
}

// global 2D nodes (and bed levels at the nodes) of the curvilinear grid; the rank's nodes with
// three ghost layers, ghost nodes outside the domain by point reflection at the boundary nodes
void lexer::curv_nodes(ghostcell *pgc)
{
    int ni=0, nj=0;
    vector<double> gx,gy,gz;

    if(G50==1)
    {
        if(gridgeo==nullptr || gridgeo->curv_ni<=0)
        {
            if(mpirank==0)
            cout<<"!!! G 50 1: no CURV section in DIVEMesh_Grid/grid-geometry.dat (DIVEMesh R 1 1) !!!"<<endl;
            exit(1);
        }

        ni = gridgeo->curv_ni;
        nj = gridgeo->curv_nj;
        gx = gridgeo->curv_x;
        gy = gridgeo->curv_y;
        gz = gridgeo->curv_zb;
    }

    if(G50==2)
    {
        // curvgrid.dat: ni nj, then x y zb of node (i,j) for j=0..nj, i=0..ni (i fastest)
        ifstream in("curvgrid.dat");

        if(!(in>>ni>>nj) || ni<1 || nj<1)
        {
            if(mpirank==0)
            cout<<"!!! G 50 2: curvgrid.dat missing or without \"ni nj\" !!!"<<endl;
            exit(1);
        }

        const size_t nn = size_t(ni+1)*size_t(nj+1);
        gx.resize(nn); gy.resize(nn); gz.resize(nn);

        for(size_t q=0; q<nn; ++q)
        if(!(in>>gx[q]>>gy[q]>>gz[q]))
        {
            if(mpirank==0)
            cout<<"!!! G 50 2: curvgrid.dat ends at node "<<q<<" of "<<nn<<" !!!"<<endl;
            exit(1);
        }
    }

    if(ni!=gknox || nj!=gknoy)
    {
        if(mpirank==0)
        cout<<"!!! G 50: curvilinear grid "<<ni<<" x "<<nj<<" cells, the DIVEMesh grid "<<gknox<<" x "<<gknoy
            <<": the DIVEMesh grid must be the logical box of the curvilinear grid !!!"<<endl;
        exit(1);
    }

    if(ni<margin || nj<margin)
    {
        if(mpirank==0)
        cout<<"!!! G 50: at least "<<margin<<" cells in i and j !!!"<<endl;
        exit(1);
    }

    // the global nodes stay for the AMR patches (placed patches: nodes of other ranks)
    curv_ni = ni;
    curv_nj = nj;
    curv_gx.swap(gx);
    curv_gy.swap(gy);
    curv_gz.swap(gz);

    geo_zbn.assign(geo_nnode,0.0);

    for(int ii=imin; ii<=imin+imax; ++ii)
    for(int jj=jmin; jj<=jmin+jmax; ++jj)
    curv_node(origin_i+ii, origin_j+jj, XN2D[nij(ii,jj)], YN2D[nij(ii,jj)], geo_zbn[nij(ii,jj)]);

    geo_curv = 1;
}

// G 50: global node (K,L), outside the domain by point reflection at the boundary nodes (in i at the
// boundary row, then in j); bed level of the nearest boundary node
void lexer::curv_node(int K, int L, double &x, double &y, double &z) const
{
    const int ni=curv_ni, nj=curv_nj;
    const vector<double> &gx=curv_gx, &gy=curv_gy, &gz=curv_gz;
    auto gid = [&](int a, int b) { return size_t(a) + size_t(ni+1)*size_t(b); };

    int Kr=K, Lr=L;
    int Kb=K, Lb=L;

    if(K<0)  { Kb=0;  Kr=-K; }
    if(K>ni) { Kb=ni; Kr=2*ni-K; }
    if(L<0)  { Lb=0;  Lr=-L; }
    if(L>nj) { Lb=nj; Lr=2*nj-L; }

    Kr = std::min(std::max(Kr,0),ni);
    Lr = std::min(std::max(Lr,0),nj);

    double xr=gx[gid(Kr,Lr)], yr=gy[gid(Kr,Lr)];

    if(Kb!=K)
    {
        xr = 2.0*gx[gid(Kb,Lr)] - xr;
        yr = 2.0*gy[gid(Kb,Lr)] - yr;
    }

    if(Lb!=L)
    {
        const int Kc = std::min(std::max(K,0),ni);
        double xb=gx[gid(Kc,Lb)], yb=gy[gid(Kc,Lb)];

        if(Kb!=K)
        {
            xb = 2.0*gx[gid(Kb,Lb)] - gx[gid(Kr,Lb)];
            yb = 2.0*gy[gid(Kb,Lb)] - gy[gid(Kr,Lb)];
        }

        xr = 2.0*xb - xr;
        yr = 2.0*yb - yr;
    }

    x = xr;
    y = yr;
    z = gz[gid(std::min(std::max(K,0),ni),std::min(std::max(L,0),nj))];
}

// G 50: bed level of the cells from the bed levels of the curvilinear nodes (after grid_solids)
void lexer::geometry_bed()
{
    if(geo_curv!=1 || geo_zbn.empty())
    return;

    for(int ii=0; ii<knox; ++ii)
    for(int jj=0; jj<knoy; ++jj)
    {
        const int q = (ii-imin)*jmax + (jj-jmin);
        const double zb = 0.25*(geo_zbn[nij(ii,jj)] + geo_zbn[nij(ii+1,jj)] + geo_zbn[nij(ii,jj+1)] + geo_zbn[nij(ii+1,jj+1)]);

        topobed[q] = zb;
        bed[q] = MAX(global_zmin,zb);
    }
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
        <<(G50==1 ? " (the grid of the run, G 50 1)" : " (not the grid of the run; G 50 1 uses it)")<<endl;

    if(folded>0)
    cout<<"!!! CURV corridor grid: "<<folded<<" folded cells !!!"<<endl;
}
