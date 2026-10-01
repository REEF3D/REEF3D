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

#ifndef SIXDOF_REMESH_H_
#define SIXDOF_REMESH_H_

// Adaptive isotropic re-triangulation of a closed (or open) oriented surface
// triangulation, used for the 6DOF body surface (X 185 4).
//
// The input triangle soup is welded into an indexed mesh, T-junctions are repaired
// and sharp features (dihedral angle > feature angle, open boundaries, non-manifold
// edges) are detected and chained into feature curves between corner vertices.
// The surface is then remeshed towards a target edge length h(x) given by a sizing
// function (the local fluid grid spacing):
//
//   split    edges longer  than 4/3 h
//   collapse edges shorter than 4/5 h (link condition, no normal flips, no new long edges)
//   flip     edges to drive the vertex valence to its ideal value (6 inside, 4 on an open
//            boundary, angle-sum/60deg at corners)
//   relax    vertices tangentially towards the h^-4 weighted centroid of their 1-ring
//            (crease vertices along their feature curve) and project them back onto the
//            original surface (resp. its feature curves)
//
// followed by a Delaunay (max-min-angle) flip pass and final relaxation sweeps.
// Sharp edges and corners are preserved exactly, triangle orientation is preserved and
// all vertices remain on the original surface. The class is self-contained (std only)
// and deterministic.
//
// Reference: M. Botsch, L. Kobbelt, "A remeshing approach to multiresolution modeling",
// Proc. Eurographics Symposium on Geometry Processing, 2004.

#include<array>
#include<vector>
#include<functional>
#include<ostream>

class sixdof_remesh
{
public:
    using vec3 = std::array<double,3>;
    using sizing_func = std::function<double(double,double,double)>;

    struct params
    {
        double feature_angle = 30.0;   // [deg] dihedral angle above which an edge is a sharp feature
        int iterations = 10;           // split/collapse/flip/relax cycles
        int smooth_final = 3;          // final Delaunay flip + relaxation sweeps
        double merge_tol = 1.0e-7;     // vertex welding tolerance relative to the bounding box diagonal
        double tjunction_tol = 1.0e-6; // T-junction detection tolerance relative to the bounding box diagonal
        long max_tri = 2000000;        // abort (and keep the input) if the estimated triangle count exceeds this
    };

    struct stats
    {
        int ntri_in=0, ntri_out=0, nvert_out=0;
        int n_feature_edges=0, n_corners=0, n_chains=0;
        int n_boundary_edges=0, n_nonmanifold_edges=0, n_inconsistent_edges=0, n_tjunctions=0;
        long ntri_estimate=0;
        double area_in=0.0, area_out=0.0, vol_in=0.0, vol_out=0.0;
        double minangle_in=0.0, minangle_out=0.0;           // [deg] smallest angle in the mesh
        double q_mean_in=0.0, q_min_in=0.0;                  // q = 4 sqrt(3) A / sum(l^2), 1 = equilateral
        double q_mean_out=0.0, q_min_out=0.0;
        double frac_q05_in=0.0, frac_q05_out=0.0;            // fraction of triangles with q < 0.5
        double Lh_mean=0.0, Lh_min=0.0, Lh_max=0.0;          // edge length / target length
        bool ok=false;
    };

    // in : triangle soup, 3 consecutive vertices per triangle, consistently oriented
    // out: remeshed triangle soup with the same orientation
    // h  : target edge length at (x,y,z), > 0
    // returns false (out = in) if the mesh could not be remeshed
    bool remesh(const std::vector<vec3> &in, std::vector<vec3> &out, const sizing_func &h, const params &prm, stats &st);

    static void print_stats(std::ostream &os, const stats &st);
};

#endif
