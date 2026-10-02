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

#ifndef GEO_RAYCAST_H_
#define GEO_RAYCAST_H_

#include"increment.h"
#include<vector>

class lexer;

using namespace std;

// Geometry core: ray-cast kernels for triangulated surfaces.
//
// The kernels were taken from the 6DOF floating body geometry and are shared by
// the 6DOF objects (CFD, NHFLOW, FNPF), the NHFLOW solid forcing and the solids
// and topography of the grid file (all hydrodynamic modules).
//
// Triangles are passed as tri_x/tri_y/tri_z[n][3], n in [ts,te).
// Cell fields use the IJK layout of the fdm fields (field.V, p->flag4).
//
// Cartesian grid (CFD, PTF, grid solids)
//   cart_io    inside/outside by ray parity along x (dir 0), y (1) and z (2);
//              each call marks IO=-1 where the entity encloses the cell centre
//   cart_dist  unsigned distance along the grid lines to the crossings,
//              LS = min(LS, |R - x|); dir 2 optionally raises a bed level
//   cart_vertexdist  distance to the triangle vertices (X 188 2)
//
// Sigma grid (NHFLOW)
//   sigma_io       inside/outside by vertical ray parity on the sigma nodes ZSP
//   sigma_band     exact point-triangle distance in a band around the surface
//
// Columns (FNPF)
//   column_hits    sorted, unique z-crossings of the vertical ray of every column

struct geo_cart
{
    // cells of the computation [is,ie) x [js,je) x [ks,ke)
    int is,ie,js,je,ks,ke;

    // true:  the cells are the subdomain interior, cell search with p->posc_*,
    //        only active cells (flag4>0), IO valid in one ghost layer (6DOF)
    // false: the cells may include ghost cells, cell search on the node arrays
    //        incl. ghost nodes, all cells (grid solids)
    bool interior;

    // skip triangles without a vertex inside the global domain (6DOF)
    bool clip;

    // 1: inside where both crossing counts are odd, 2: where both are even
    int raymode;

    // ray placement against rays through edges and vertices
    // 0: ray tilted by psi in the transverse directions (6DOF)
    // 1: ray shifted by psi in the transverse directions (DIVEMesh, grid solids):
    //    unlike the tilt it also leaves the diagonals of axis-aligned faces
    int raystyle;

    double psi;     // ray offset against rays through edges and vertices
    double ext;     // ray extension beyond the global domain
};

struct geo_ray
{
    double Px,Py,Pz,Qx,Qy,Qz;
};

class geo_raycast : public increment
{
public:
    geo_raycast(lexer*);
    virtual ~geo_raycast();

    // 6DOF set-up: interior cells, clipped triangles
    static geo_cart interior(lexer*);

    // grid solids: interior plus the ghost cells of the fields
    static geo_cart extended(lexer*, double dxm);

    // Cartesian
    void cart_io(lexer*, int dir, double **tri_x, double **tri_y, double **tri_z, int ts, int te,
                 const geo_cart&, int *IO, int *CL, int *CR);

    void cart_dist(lexer*, int dir, double **tri_x, double **tri_y, double **tri_z, int ts, int te,
                   const geo_cart&, const int *IO, double *LS, double *bed=nullptr);

    void cart_vertexdist(lexer*, double **tri_x, double **tri_y, double **tri_z, int ts, int te, double *LS);

    // sigma
    void sigma_io(lexer*, double **tri_x, double **tri_y, double **tri_z, int ts, int te,
                  double DSM, int raymode, int *IO, int *CL, int *CR);

    void sigma_band(lexer*, double **tri_x, double **tri_y, double **tri_z, int ts, int te,
                    double band, double *LS);

    static double dist2_tri(double Px, double Py, double Pz,
                            double Ax, double Ay, double Az,
                            double Bx, double By, double Bz,
                            double Cx, double Cy, double Cz);

    // columns: hits[i*knoy+j], i in [0,knox), j in [0,knoy)
    void column_hits(lexer*, double **tri_x, double **tri_y, double **tri_z, int ts, int te,
                     double dxm, vector<vector<double> > &hits);

private:
    // cell index search, Cartesian
    int cell(lexer*, int dir, double x, const geo_cart&);
    static int cell_ext(const double *N, int lo, int hi, double x);
    void bracket(lexer*, int d, double s, double e, const geo_cart&, int &is, int &ie);
    void setray(lexer*, int dir, int i, int j, int k, const geo_cart&, geo_ray&);

    bool checkin(lexer*, int dir, double Ax, double Ay, double Az,
                 double Bx, double By, double Bz, double Cx, double Cy, double Cz);

    const double epsi;
};

#endif
