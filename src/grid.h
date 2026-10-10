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
    int geo_curv = 0;                   // 0: Cartesian product grid, 1: curvilinear (no solver yet)
    double *XN2D=nullptr, *YN2D=nullptr;
    double *XC2D=nullptr, *YC2D=nullptr, *AREA2D=nullptr;
    double *SX1=nullptr, *SY1=nullptr, *SX2=nullptr, *SY2=nullptr;
    int geo_nnode=0, geo_nslice=0;      // allocated sizes

    int nij(int ii, int jj) const { return (ii-imin)*(jmax+1) + (jj-jmin); }

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
