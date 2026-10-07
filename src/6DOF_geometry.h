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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#ifndef SIXDOF_GEOMETRY_H_
#define SIXDOF_GEOMETRY_H_

#include"increment.h"
#include<vector>
#include<Eigen/Dense>

class lexer;
class ghostcell;

using namespace std;

//  Surface geometry of a 6DOF body: the triangulated hull and its pose transformation.
//
//  Solver independent (lexer for the input keys and the grid spacing used to size the
//  triangles, ghostcell for the MPI reductions; no fdm). Used by sixdof_obj for all
//  hydrodynamic models and by the ship module.
//
//  Triangles: tri_x/tri_y/tri_z[n][3] in the inertial frame (current pose),
//             tri_x0/tri_y0/tri_z0[n][3] relative to the centre of gravity in the body frame.
//  Entities:  one per input object (X 110 box, X 131-133 cylinders, X 153 wedge_sym, X 163 wedge,
//             X 164 hexahedron, X 170/172 wavemakers, X 180 STL), triangles tstart[e] .. tend[e]-1.
//
//  Set-up:  allocate() -> primitives() / read_stl() (and the wavemaker objects of sixdof_obj)
//           -> orient() -> refine()   ... then store_body_frame(c) once the CoG is known
//  Motion:  transform(R, c)

class sixdof_geometry : public increment
{
public:

    sixdof_geometry(int);
    virtual ~sixdof_geometry();

    // ---- set-up
    void allocate(lexer*);                      // triangle arrays for all input objects
    void primitives(lexer*, ghostcell*);        // X 110, 131-133, 153, 163, 164
    void read_stl(lexer*, ghostcell*);          // X 180: floating.stl / floating-<id>.stl
    void orient(lexer*, ghostcell*);            // consistent outward orientation (ray test)
    void refine(lexer*, ghostcell*);            // X 185: split (1-3) or remesh (4) to the grid spacing

    // ---- pose
    void store_body_frame(const Eigen::Vector3d&);                  // tri0 = tri - c
    void rotate(double, double, double, const Eigen::Vector3d&);    // rotate tri about c (Euler angles)
    void transform(const Eigen::Matrix3d&, const Eigen::Vector3d&); // tri = R tri0 + c

    // rotation of a point about (x0,y0,z0) with the Euler angles phi, theta, psi
    static void rotation_tri(double, double, double, double&, double&, double&,
                             const double&, const double&, const double&);

    // ---- properties
    double volume() const;                      // enclosed volume (divergence theorem)

    // ---- data
    const int id;                               // body number
    double **tri_x, **tri_y, **tri_z, **tri_x0, **tri_y0, **tri_z0;
    int tricount;
    int entity_count, entity_sum;
    int *tstart, *tend;

    // X 185 4: the remesher closed holes of the surface (volume and mass properties change)
    bool holes_closed;
    
    // horizontal spacing factor of the triangle sizing (finest level of a mesh refinement zone)
    double amr_hfac;

private:
    // primitives
    void box(lexer*, int);
	void cylinder_x(lexer*, int);
	void cylinder_y(lexer*, int);
	void cylinder_z(lexer*, int);
	void wedge_sym(lexer*, int);
    void wedge(lexer*, int);
    void hexahedron(lexer*, int);

    // refinement
    void geometry_refinement(lexer*, ghostcell*);
	void geometry_remesh(lexer*, ghostcell*);
	void create_triangle(double&,double&,double&,double&,double&,double&,double&,double&,double&,const double&,const double&,const double&);
    vector<vector<double> > tri_x_r, tri_y_r, tri_z_r;

    // orientation
    int *tri_switch, *tri_switch_local;
    int tricount_local, *tricount_local_list, *tricount_local_displ;
    int tricount_switch_total;

    int n, q;
};

#endif
