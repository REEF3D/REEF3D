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

    double *RN,*SN,*TN; // Temporary arrays
    double *RP,*SP,*TP; // Temporary arrays
    double DRM,DSM,DTM;
    double *DRDXN,*DSDYN,*DTDZN;
    double *DRDXP,*DSDYP,*DTDZP;
    
    // origin and rotation
    double global_orig_x,global_orig_y;
    double alpha_grid;

    // boundary conditions
    int *IO,*IOSL;
    int *DF,*DF1,*DF2,*DF3;
    int *DFBED;

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
};

#endif
