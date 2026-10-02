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

#ifndef GEO_PRIMITIVE_H_
#define GEO_PRIMITIVE_H_

// Geometry core: triangulated primitives (6DOF objects, NHFLOW solids).
// Every function appends its surface triangles to tri_x/tri_y/tri_z[tricount][3].

class geo_primitive
{
public:
    // number of circumferential segments: circumference/ds, at least nmin
    static int segments(double r, double ds, int nmin);

    static void box(double**, double**, double**, int &tricount,
                    double xs, double xe, double ys, double ye, double zs, double ze);

    static void cylinder_x(double**, double**, double**, int &tricount,
                           double ym, double zm, double x1, double x2, double r, int snum);
    static void cylinder_y(double**, double**, double**, int &tricount,
                           double xm, double zm, double y1, double y2, double r, int snum);
    static void cylinder_z(double**, double**, double**, int &tricount,
                           double xm, double ym, double z1, double z2, double r, int snum);

    static void wedge_sym(double**, double**, double**, int &tricount,
                          double xs, double xe, double ys, double ye, double zs, double ze);

    // 6 vertices x1 y1 z1 ... x6 y6 z6
    static void wedge(double**, double**, double**, int &tricount, const double *v);

    // 8 vertices x1 y1 z1 ... x8 y8 z8
    static void hexahedron(double**, double**, double**, int &tricount, const double *v);

    static void sphere(double**, double**, double**, int &tricount,
                       double xm, double ym, double zm, double r, int snum);

    static void wedge_x(double**, double**, double**, int &tricount,
                        double xs, double xe, double ys, double ye, double zs, double ze);
    static void wedge_y(double**, double**, double**, int &tricount,
                        double xs, double xe, double ys, double ye, double zs, double ze);
    static void wedge_z(double**, double**, double**, int &tricount,
                        double xs, double xe, double ys, double ye, double zs, double ze);
};

#endif
