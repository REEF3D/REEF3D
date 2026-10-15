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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#ifndef GRID_H_
#define GRID_H_

#include "increment.h"
#include <vector>

class ghostcell;

class grid : virtual public increment
{
public:
    grid() = default;
    virtual ~grid() = default;

    void assign_margin();
    void sigma_coord_ini();

    void gridspacing(ghostcell *pgc);

    // Non-Uniform Mesh
    double *XN,*YN,*ZN; // Nodal coordinates
    double *XP,*YP,*ZP; // Cell center coordinates
    double *DXN,*DYN,*DZN; // Nodal grid spacing
    double *DXP,*DYP,*DZP; // Cell center grid spacing
    double *ZSN,*ZSP;
    double DXM,DYD,DXD;
    double DYM,DZM;
    
    // origin and rotation
    double global_orig_x,global_orig_y;
    double alpha_grid;

    // boundary conditions
    int *IO,*IOSL;
    int *DF,*DF1,*DF2,*DF3,*DFF;
    int *DFBED=nullptr;   // sediment cell flag (2D), allocated in flagini() or, for SFLOW, by the sediment module

    bool i_dir,j_dir,k_dir;
    double x_dir,y_dir,z_dir;

    int **gcin, **gcout;
    int gcin_count, gcout_count;

    // maxcoor
    double xcoormax,xcoormin,ycoormax,ycoormin,zcoormax,zcoormin;
    double maxlength;

    int knox,knoy,knoz;
    const int margin = 3;

    double originx,originy,originz;
    double endx,endy,endz;
    double global_xmin,global_ymin,global_zmin;
    double global_xmax,global_ymax,global_zmax;
    int origin_i, origin_j, origin_k;
    int gknox,gknoy,gknoz;

    int imin,imax,jmin,jmax,kmin,kmax,kmaxF;

    double dx,dy,dz;

    // ---- horizontal geometry layer (grid_geometry.cpp)
    // 2D nodes and the cell and face metrics of the slice, the basis of a curvilinear grid in the
    // horizontal plane.  On the Cartesian grids of now (geo_curv=0) they are filled from XN/YN and
    // no solver reads them; the AMR patches subdivide the level-0 nodes bilinearly.
    //   node (i,j), i=imin..imin+imax, j=jmin..jmin+jmax : XN2D/YN2D[nij(i,j)]
    //   cell (i,j), slice index IJ                       : XC2D/YC2D centroid, AREA2D area
    //   face i+1/2 of cell (i,j) (as u, flagslice1)      : SX1/SY1 = normal (+i side) * face length
    //   face j+1/2 of cell (i,j) (as v, flagslice2)      : SX2/SY2 = normal (+j side) * face length
    //                                                      LEN1/LEN2 face lengths
    int geo_curv = 0;                   // 0: Cartesian product grid, 1: curvilinear (G 50, SFLOW)
    int geo_mb = 0;                     // G 50 1 with blocks (DIVEMesh R 20 1): number of blocks, 0 one grid
    double *XN2D=nullptr, *YN2D=nullptr;
    double *XC2D=nullptr, *YC2D=nullptr, *AREA2D=nullptr;
    double *SX1=nullptr, *SY1=nullptr, *SX2=nullptr, *SY2=nullptr;
    double *LEN1=nullptr, *LEN2=nullptr;
    int geo_nnode=0, geo_nslice=0;      // allocated sizes

    int nij(int ii, int jj) const { return (ii-imin)*(jmax+1) + (jj-jmin); }

    // ---- physical points on the curvilinear grid (grid_locate.cpp; built on the first call, valid as long as
    //      the nodes; usable on Cartesian grids as well)
    //   geo_quad:  the quadrilateral of the cell centres (i,j),(i+1,j),(i+1,j+1),(i,j+1) that holds (x,y), i in
    //              [-1,knox-1], and the bilinear coordinates fx, fy in [0,1] (the Cartesian posf_i/j and
    //              1-wa, 1-wb of the interpolations); outside the centres of the rank: the nearest quad,
    //              fx, fy clamped, false
    //   geo_cell:  the cell (i,j) of the rank (0..knox-1, 0..knoy-1) that holds (x,y), half-open in both
    //              index directions (a point on a face belongs to one cell); false: none
    //   geo_cells: the cells of the rank whose centre lies in the box [x0,x1]x[y0,y1] (pairs i,j)
    //   geo_box:   the box of the nodes of the rank's cells
    struct geo_locate_cache *geoloc=nullptr;
    bool geo_quad(double x, double y, int &i, int &j, double &fx, double &fy);
    bool geo_cell(double x, double y, int &i, int &j);
    void geo_cells(double x0, double x1, double y0, double y1, std::vector<int> &ij);
    void geo_box(double &x0, double &x1, double &y0, double &y1);
    //   geo_width: the mean width of cell (i,j), 0.5 (area/mean length of the i-faces + area/mean length of the j-faces)
    //   geo_grad:  the physical gradient at cell (i,j) from the values of its four neighbours (chain rule on the
    //              cell centres, as sflow_curv::grad5)
    double geo_width(int i, int j) const;
    //   cellarea:  the horizontal area of cell (i,j): AREA2D on the curvilinear grid, DXN DYN else (bit for bit)
    double cellarea(int ii, int jj) const
    { return geo_curv==1 ? AREA2D[(ii-imin)*jmax + (jj-jmin)] : DXN[ii+marge]*DYN[jj+marge]; }
    void geo_grad(int i, int j, double fE, double fW, double fN, double fS, double &fx, double &fy) const;
    void geo_locate_free();
    //   geo_bilinear_inverse: (s,t) of the point p in the quadrilateral 00-10-11-01 (coordinates relative to node 00)
    static bool geo_bilinear_inverse(double px, double py, double ax, double ay, double bx, double by,
                                     double cx, double cy, double &s, double &t);

    void geometry_alloc();
    void geometry_free();
    void geometry_cartesian_nodes();    // XN2D/YN2D = (XN,YN)
    void geometry_metrics();            // cell and face metrics from the 2D nodes
    int geometry_check(double &err, double &gcl) const;   // mismatches of nodes vs XN/YN (Cartesian),
                                                          // max deviation of the stored metrics from the
                                                          // quadrilateral formulas, max closure residual
};

// quadrilateral metrics of a cell from its corner nodes 00,10,01,11 (i,j order), shared by the
// geometry layer and the checks of external curvilinear grids (CURV)
namespace geo2d
{
    inline double area(double x00, double y00, double x10, double y10, double x01, double y01, double x11, double y11)
    {
        return 0.5*((x11-x00)*(y01-y10) - (y11-y00)*(x01-x10));
    }

    void centroid(double x00, double y00, double x10, double y10, double x01, double y01, double x11, double y11,
                  double &xc, double &yc);
}

#endif
