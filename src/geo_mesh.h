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

    // CURV (optional): river corridor grid of DIVEMesh (R 1), read and checked, not used by a solver yet
    //   nodes (i,j), i=0..curv_ni along the river from the inflow, j=0..curv_nj from the right bank,
    //   index i + (curv_ni+1)*j; curv_par = {R 2 wl0, R 3 wl, R 3 margin, R 5 open, raster spacing}
    int curv_ni, curv_nj;
    vector<double> curv_par, curv_x, curv_y, curv_zb;

    // BLKS (optional): multi-block grid of DIVEMesh (R 20 1); block 0 is the grid of the run (CURV,
    // with a ring of blk_ring cells on its interface faces, blk[0].ring), blocks 1.. carry their
    // nodes; interface segments: face, lo, hi (cells along it), neighbour block and face, cell
    // transform i' = sign0*c[perm0]+off0, j' = sign1*c[perm1]+off1, c = (i,j)
    struct blk_seg { int face, lo, hi, nbr, nface, perm[2], sign[2], off[2]; };
    struct blk_block
    {
        int ni=0, nj=0;
        int ring[4] = {0,0,0,0};
        vector<blk_seg> seg;
        vector<double> x, y, zb;
        int orig=0, oi=0, oj=0;     // G 55: the block of DIVEMesh this piece is cut from, its first cell there
    };
    vector<blk_block> blk;
    int blk_ring = 0;

    // G 55 m (geo_mesh_split.cpp): the blocks 1.. of DIVEMesh cut into pieces along i (j for a block wider in j
    // without open faces) of at most (cells of the blocks 1..) / (m ranks) cells each, at least 16 cells deep;
    // a piece is a block of its own (the first piece keeps the number of its block, the others follow the
    // blocks of DIVEMesh), the interfaces between the pieces and to the neighbours rebuilt cell by cell
    void split_blocks(lexer*);
    int nblk_orig = 0;              // the blocks of DIVEMesh (blk.size() before the cut)
    vector<int> orig_ni, orig_nj;   // their sizes
    // the piece that holds face f (0 -i, 1 +i, 2 -j, 3 +j) of block b of DIVEMesh whole; -1: cut across it
    int face_piece(int b, int f) const;

    // G 56 a dz (geo_mesh_split.cpp): the load of a cell relative to a wet cell, a when its bed lies above the
    // still water level + dz (F 60, else F 62 with F 65); without G 56 (a = 1) or a still water level every cell 1
    static double load_level(const lexer*);
    static bool load_weighted(const lexer*);
    static double cell_load(const lexer*, double zb);
    double blk_cell_load(const lexer*, int b, int ci, int cj) const;   // blocks 1..: the bed is the mean of the cell's nodes
    double blk_load(const lexer*, int b) const;

private:
    void allocate(int);
    void release();
};

#endif
