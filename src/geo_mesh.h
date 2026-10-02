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

#ifndef GEO_MESH_H_
#define GEO_MESH_H_

#include<vector>
#include<string>

class lexer;

using namespace std;

// Geometry core: triangulated objects.
//
// The solids (S) and topography (T) of DIVEMesh arrive as surface triangles in
// DIVEMesh_Grid/grid-geometry.dat (grid format v2) and are kept here for the whole run.
// Triangles are stored as tri_x/tri_y/tri_z[n][3], the layout of the 6DOF objects,
// so all geometry kernels (geo_raycast) work on both.

struct geo_object
{
    int role;       // 1 solid, 2 topo, 3 zero-thickness plate (wall surfaces only)
    int keyword;    // DIVEMesh keyword: 1 STL, 10 box, 11 box array, 32 cylinder_y, ...
    int index;      // entity number of the keyword
    int raymode;    // 1: inside where the crossing counts are odd, 2: even
    int invert;     // inside/outside of the field swapped after this object (S 9 2, T 9 2)
    int ts, te;     // triangle range
    vector<double> param;
};

class geo_mesh
{
public:
    geo_mesh();
    virtual ~geo_mesh();

    static const int role_solid = 1;
    static const int role_topo = 2;
    static const int role_plate = 3;

    // DIVEMesh_Grid/grid-geometry.dat
    void read(lexer*, const char *name);

    int count(int role) const;
    int tricount(int role) const;

    vector<geo_object> obj;

    double **tri_x, **tri_y, **tri_z;
    int ntri;

    // GHDR
    int solidread, toporead;
    int geodat;     // geodat bed level (GEOB) belongs to: 0 none, 1 solid, 2 topo
    double dxm;     // DIVEMesh mean cell size

private:
    void allocate(int);
    void release();
};

#endif
