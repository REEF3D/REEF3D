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

#ifndef FNPF_ICE_FLOE_H_
#define FNPF_ICE_FLOE_H_

#include<vector>

using namespace std;

// Rigid ice floe: prism of thickness h over a convex planform.
// Body frame: origin at the centre of mass (planform centroid, mid-thickness),
// planform vertices bx,by counter-clockwise, bottom face at z_body = -h/2.
// Obstacles (type 1) are fixed convex polygons that take part in contact only.

struct fnpf_ice_floe
{
    int id = 0;
    int type = 0;                   // 0: floe, 1: fixed obstacle (contact only, load logged)

    // geometry, body frame
    vector<double> bx,by;
    double h = 1.0;
    double rho = 917.0;
    double area = 0.0;
    double mass = 0.0;
    double Ib[3][3] = {{0.0}};      // inertia tensor about the centre of mass, body frame
    double rbound = 0.0;            // planform bounding radius
    double width2D = 1.0;           // 2D flume: width carried by one cell row (area / planform length)
    double wmin = 0.0;              // minimum caliper width of the planform

    // state
    double x[3] = {0.0};            // centre of mass
    double q[4] = {1.0,0.0,0.0,0.0};// orientation quaternion (w,x,y,z), body -> world
    double v[3] = {0.0};            // velocity of the centre of mass
    double w[3] = {0.0};            // angular velocity, world frame

    // lid parameters
    double klid = 0.0;              // lid stiffness [N/m^3]
    double clid = 0.0;              // lid damping   [Ns/m^3]
    double pw = 0.0;                // weight per unit area rho_i*g*h [Pa]

    // per-step quantities, filled by the footprint and the force integration
    double R[3][3] = {{0.0}};       // rotation matrix, body -> world
    double nb[3] = {0.0,0.0,1.0};   // upward normal of the bottom face, world
    double pb[3] = {0.0};           // centre of the bottom face, world
    vector<double> wx,wy;           // bottom face vertices projected on the horizontal plane
    double bbox[4] = {0.0};         // xmin xmax ymin ymax of wx,wy
    double tap = 0.0;               // footprint edge taper half width of this floe
    double omega = 0.0;             // highest lid-spring frequency (heave, roll/pitch) of the footprint

    double F[3] = {0.0};            // hydrodynamic force  (lid + drag), current RK stage
    double M[3] = {0.0};            // hydrodynamic moment about the centre of mass, world
    double Awet = 0.0;              // loaded footprint area

    double acc[3] = {0.0};          // acceleration of the centre of mass, last RK stage
    double alp[3] = {0.0};          // angular acceleration, world, last RK stage

    double Fc[2] = {0.0};           // contact force (obstacles: ice load), world, planar
    int ncontact = 0;
};

#endif
